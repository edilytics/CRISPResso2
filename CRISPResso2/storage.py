"""Flat-memory adaptive storage backend for CRISPResso2.

This module implements the ``parquet`` storage backend introduced by the
flat-memory design (see ``design_docs/STREAMING_GROUPBY_SPIKE.md`` and the
master design plan in ``.pi/plans/``). It is gated behind the
``--storage_backend parquet`` flag; the default ``pandas`` path in
``CRISPRessoCORE.py`` is untouched and serves as the parity oracle.

Stage coverage
--------------
* **Stage 1 — raw-read counting** (``VariantStore.count_reads``): replaces the
  in-memory ``variantCache`` dedup dict with an *adaptive* approach.
  - *Eager path* (small inputs, below the memory budget): builds the exact
    dedup dict the current pandas path uses — byte-for-byte parity, no I/O
    overhead. This is the spectrum's "fast" end.
  - *Spill path* (large inputs, above the memory budget): streams read keys to
    a temp text file, then runs an external merge-sort + streaming collapse
    (``LC_ALL=C sort -S <buf> | uniq -c`` when qualities are not retained, or a
    streaming first-per-group collapse when they are). Peak RSS stays flat
    across arbitrary row counts and read lengths — confirmed by spike S1b in
    ``design_docs/STREAMING_GROUPBY_SPIKE.md`` (ratios 1.00-2.49 at up to
    1M x 10kb, vs 22-24x linear scaling for a hash ``group_by``).

Stages 2 (worker parquet writers), 3 (streaming collapse), and 4 (stream-out:
count vectors + allele TSV sink) are implemented; this module exposes the
Stage 1 ``VariantStore.count_reads`` API, the ``ReadCounts`` result handle,
the Stage 2 worker parquet writer, the Stage 3 ``VariantStore.collapse``
streaming collapse, the Stage 4 ``VariantStore.compute_count_vectors`` /
``VariantStore.write_allele_frequency_table`` stream-outs, and the Stage 5
``CollapsedAlleles.get_slice`` lazy view. None of these are wired into
``CRISPRessoCORE.main`` yet — they are exercised by unit tests against the
current pandas path for parity.
"""

from __future__ import annotations

import gzip
import logging
import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Iterator, Optional

from CRISPResso2 import CRISPRessoShared

_logger = logging.getLogger("CRISPResso2")
info = _logger.info

# Default in-memory budget for the eager dedup dict (Stage 1 spectrum threshold).
# Chosen so a typical short-read amplicon run stays eager (fast, parity path)
# while multi-million-read long-read inputs spill to disk. Exposed via
# ``VariantStore(memory_budget_mb=...)`` so tests and future auto-select can tune.
DEFAULT_MEMORY_BUDGET_MB = 512

# Cap the in-memory sort buffer so the external merge-sort spills to disk by
# design rather than grabbing a fraction of available RAM (spike S1b showed BSD
# sort defaults to a large buffer that defeats the spill on a 16GB host).
_SORT_BUFFER = "64M"

# Per-entry overhead estimate for the eager dedup dict (CPython dict slot +
# list value + str refs). Used for the spectrum memory accounting.
_DICT_ENTRY_OVERHEAD_BYTES = 80


def _open_fastq(path: str):
    """Open a FASTQ file, transparently handling gzip."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _read_fastq_records(
    fastq1_path: str, fastq2_path: Optional[str]
) -> Iterator[tuple[str, str]]:
    """Yield ``(read_key, quals)`` per read, matching the current dedup semantics.

    For single-read input (``fastq2_path is None``) the key is the R1 sequence
    and the quality string is the R1 phred string — matching
    ``process_fastq`` (``variantCache[fastq_seq] = count``; quals not stored).

    For paired input the key is ``R1 + '+' + reverse_complement(R2)`` and the
    quals are ``R1_qual + ' ' + reversed(R2_qual)`` — matching
    ``process_paired_fastq``.
    """
    f1 = _open_fastq(fastq1_path)
    if fastq2_path is not None:
        f2 = _open_fastq(fastq2_path)
    else:
        f2 = None
    try:
        while True:
            id1 = f1.readline()
            if not id1:
                if f2 is not None:
                    id2_check = f2.readline()
                    if id2_check:
                        raise CRISPRessoShared.InputFileFormatException(
                            "The two fastq files are not the same length. Please check your input files."
                        )
                break
            seq1 = f1.readline().strip()
            f1.readline()  # plus line
            qual1 = f1.readline().strip()
            if f2 is not None:
                id2 = f2.readline()
                if not id2:
                    raise CRISPRessoShared.InputFileFormatException(
                        "The two fastq files are not the same length. Please check your input files."
                    )
                seq2 = CRISPRessoShared.reverse_complement(f2.readline().strip())
                f2.readline()  # plus line
                qual2 = f2.readline().strip()[::-1]
                key = seq1 + "+" + seq2
                quals = qual1 + " " + qual2
            else:
                key = seq1
                quals = qual1
            yield key, quals
    finally:
        f1.close()
        if f2 is not None:
            f2.close()


def _dict_memory_bytes(num_unique: int, total_key_bytes: int, total_qual_bytes: int) -> int:
    """Estimate eager-dict RSS from accumulated key/qual sizes + entry overhead."""
    return num_unique * _DICT_ENTRY_OVERHEAD_BYTES + total_key_bytes + total_qual_bytes


@dataclass
class ReadCounts:
    """Result of Stage 1 read counting (``VariantStore.count_reads``).

    Two shapes depending on the spectrum decision:

    * ``spilled == False`` (eager, small inputs): the dedup dict is held in
      memory. :meth:`to_dict` / :meth:`to_dict_with_quals` reproduce the exact
      ``variantCache`` shapes the current pandas path builds — used for parity.
    * ``spilled == True`` (large inputs): counts live in a parquet file at
      :attr:`parquet_path` with columns ``read_key``, ``count`` and (if
      ``keep_quals``) ``quals``. :meth:`items` streams rows without
      materializing the whole table.

    Both shapes expose :meth:`items` so downstream Stage 2 workers can iterate
    ``(read_key, count, quals)`` uniformly.
    """

    num_total: int
    num_unique: int
    spilled: bool
    keep_quals: bool
    parquet_path: Optional[str] = None
    # eager-only (spilled == False): the parity dict
    _counts: Optional[dict] = field(default=None, repr=False)

    def to_dict(self) -> dict:
        """Return ``{read_key: count}`` — parity with single-read ``variantCache``.

        Eager path only (raises if spilled). For the spilled path use
        :meth:`items` or read ``parquet_path`` directly.
        """
        if self.spilled:
            raise RuntimeError("to_dict() is only available on the eager (non-spilled) path")
        assert self._counts is not None
        if self.keep_quals:
            return {k: v[0] for k, v in self._counts.items()}
        return dict(self._counts)

    def to_dict_with_quals(self) -> dict:
        """Return ``{read_key: [count, quals]}`` — parity with paired ``variantCache``.

        Eager path only. ``keep_quals`` must have been True at counting time.
        """
        if self.spilled:
            raise RuntimeError("to_dict_with_quals() is only available on the eager (non-spilled) path")
        if not self.keep_quals:
            raise RuntimeError("to_dict_with_quals() requires keep_quals=True at counting time")
        assert self._counts is not None
        return {k: [v[0], v[1]] for k, v in self._counts.items()}

    def items(self) -> Iterator[tuple[str, int, Optional[str]]]:
        """Stream ``(read_key, count, quals_or_None)``.

        For the eager path this iterates the in-memory dict. For the spilled
        path this streams the parquet file row-group by row-group (pyarrow
        ``iter_batches``), so peak memory is bounded by one batch — not by the
        unique-read count.
        """
        if not self.spilled:
            assert self._counts is not None
            if self.keep_quals:
                for key, (count, quals) in self._counts.items():
                    yield key, count, quals
            else:
                for key, count in self._counts.items():
                    yield key, count, None
            return
        # spilled path — stream the parquet file in bounded batches
        import pyarrow.parquet as pq

        assert self.parquet_path is not None
        columns = ["read_key", "count"] + (["quals"] if self.keep_quals else [])
        pf = pq.ParquetFile(self.parquet_path)
        for batch in pf.iter_batches(columns=columns, batch_size=50_000):
            keys = batch.column("read_key")
            counts = batch.column("count")
            quals_col = batch.column("quals") if self.keep_quals else None
            for i in range(batch.num_rows):
                q = quals_col[i].as_py() if quals_col is not None else None
                yield keys[i].as_py(), counts[i].as_py(), q


class VariantStore:
    """Adaptive storage backend owning the read->align->collapse lifecycle.

    Stage 1 (``count_reads``) is implemented here; Stages 2-3 are added in
    later PRs. The store is *not* wired into ``CRISPRessoCORE.main`` yet — it
    is exercised by unit tests against the current dedup dict for parity.
    """

    def __init__(
        self,
        output_directory: str,
        *,
        memory_budget_mb: int = DEFAULT_MEMORY_BUDGET_MB,
        keep_quals: bool = False,
        sort_buffer: str = _SORT_BUFFER,
    ):
        self.output_directory = output_directory
        self.memory_budget_bytes = int(memory_budget_mb) * 1024 * 1024
        self.keep_quals = keep_quals
        self.sort_buffer = sort_buffer
        os.makedirs(output_directory, exist_ok=True)

    # -- Stage 1: raw-read counting -----------------------------------------

    def count_reads(
        self,
        fastq1_path: str,
        fastq2_path: Optional[str] = None,
        *,
        force_spill: bool = False,
    ) -> ReadCounts:
        """Count reads from FASTQ, deduplicating by read key.

        Adaptive ("spectrum") behaviour:

        * If the eager dedup dict fits in ``memory_budget_mb`` the result is the
          in-memory dict — byte-for-byte parity with the current ``variantCache``
          and no spill I/O (the fast end of the spectrum).
        * If the dict would exceed the budget (or ``force_spill=True``) the
          remaining reads are streamed to a temp text file and collapsed with an
          external merge-sort + streaming dedup. Peak RSS stays flat in the row
          count (spike S1b).

        Parameters
        ----------
        fastq1_path, fastq2_path
            Single-read input omits ``fastq2_path``. Paired input constructs the
            key as ``R1 + '+' + reverse_complement(R2)`` (current behaviour).
        force_spill
            Force the spill path regardless of input size (used by tests and the
            future auto-select heuristic).

        Returns
        -------
        ReadCounts
            Handle over the counted reads; see :class:`ReadCounts`.

        """
        # Fast path: try eager first. If the dict outgrows the budget, drop it
        # and re-stream the file into the spill pipeline. The eager attempt is
        # single-pass; spill re-reads the input once (acceptable: spill means
        # "slow but finishes").
        if not force_spill:
            counts, num_total, spilled = self._count_eager(fastq1_path, fastq2_path)
            if not spilled:
                num_unique = len(counts)
                return ReadCounts(
                    num_total=num_total,
                    num_unique=num_unique,
                    spilled=False,
                    keep_quals=self.keep_quals,
                    _counts=counts,
                )
        # Spill path (forced, or the eager dict exceeded the budget).
        return self._count_spill(fastq1_path, fastq2_path)

    def _count_eager(
        self, fastq1_path: str, fastq2_path: Optional[str]
    ) -> tuple[dict, int, bool]:
        """Single-pass eager dedup. Returns ``(dict, num_total, spilled)``.

        ``spilled`` is True if the dict exceeded the budget mid-stream; in that
        case the partial dict is discarded (the caller re-streams via
        ``_count_spill``).
        """
        counts: dict = {}
        num_total = 0
        total_key_bytes = 0
        total_qual_bytes = 0
        budget = self.memory_budget_bytes
        for key, quals in _read_fastq_records(fastq1_path, fastq2_path):
            num_total += 1
            if key in counts:
                if self.keep_quals:
                    counts[key][0] += 1
                else:
                    counts[key] += 1
            else:
                if self.keep_quals:
                    counts[key] = [1, quals]
                    total_qual_bytes += len(quals)
                else:
                    counts[key] = 1
                total_key_bytes += len(key)
                # The estimate is O(1) (running totals), so check on every new
                # unique key rather than periodically — otherwise small
                # high-cardinality inputs (e.g. all-unique long reads) would
                # never trip the budget and would OOM.
                if _dict_memory_bytes(len(counts), total_key_bytes, total_qual_bytes) > budget:
                    return counts, num_total, True  # spilled — discard, re-stream
        return counts, num_total, False

    def _count_spill(self, fastq1_path: str, fastq2_path: Optional[str]) -> ReadCounts:
        """External-sort + streaming-collapse dedup. Peak RSS bounded by the
        sort buffer (``-S``) plus one pyarrow export batch.
        """
        workdir = tempfile.mkdtemp(prefix="crispresso_store_", dir=self.output_directory)
        keys_file = os.path.join(workdir, "keys.txt")
        sorted_file = os.path.join(workdir, "keys.sorted")
        num_total = self._stream_keys_to_file(fastq1_path, fastq2_path, keys_file)
        self._external_sort(keys_file, sorted_file)
        num_unique = self._write_collapsed_parquet(sorted_file, workdir)
        # Clean up intermediate text files; keep the parquet artifact.
        for p in (keys_file, sorted_file):
            try:
                os.remove(p)
            except OSError:
                pass
        parquet_path = os.path.join(workdir, "read_counts.parquet")
        return ReadCounts(
            num_total=num_total,
            num_unique=num_unique,
            spilled=True,
            keep_quals=self.keep_quals,
            parquet_path=parquet_path,
        )

    def _stream_keys_to_file(
        self, fastq1_path: str, fastq2_path: Optional[str], keys_file: str
    ) -> int:
        """Stream every read to ``keys_file`` as ``read_key <TAB> quals`` lines.

        ``quals`` is written only when ``keep_quals``; the tab separator is safe
        because read keys are ACGT/+ and phred quals contain no tabs.
        """
        num_total = 0
        with open(keys_file, "w", buffering=1 << 20) as kf:  # 1 MiB buffer
            for key, quals in _read_fastq_records(fastq1_path, fastq2_path):
                num_total += 1
                if self.keep_quals:
                    kf.write(key)
                    kf.write("\t")
                    kf.write(quals)
                    kf.write("\n")
                else:
                    kf.write(key)
                    kf.write("\n")
        return num_total

    def _external_sort(self, keys_file: str, sorted_file: str) -> None:
        """External merge-sort of the keys file, spilling to ``workdir``.

        ``-S`` caps the in-memory buffer (forcing disk spill); ``-T`` sets the
        spill dir; ``-s`` makes the sort stable so that, when quals are carried,
        the first occurrence in file order wins per key group (parity with the
        eager dict's first-occurrence quals). ``LC_ALL=C`` forces byte-wise
        comparison (fast, locale-independent).
        """
        tmpdir = os.path.dirname(keys_file) + os.sep + "sort_tmp"
        os.makedirs(tmpdir, exist_ok=True)
        # When quals are carried the line is `key<TAB>quals`; sort by the key
        # field ONLY (-t TAB -k1,1) so a stable sort preserves file order within
        # a key group — otherwise sorting by the whole line reorders by quals and
        # the first-occurrence quals parity is lost.
        if self.keep_quals:
            sort_cmd = [
                "sort", "-S", self.sort_buffer, "-T", tmpdir, "-s",
                "-t", "\t", "-k1,1",
                "-o", sorted_file, keys_file,
            ]
        else:
            sort_cmd = [
                "sort", "-S", self.sort_buffer, "-T", tmpdir, "-s",
                "-o", sorted_file, keys_file,
            ]
        env = dict(os.environ, LC_ALL="C")
        try:
            subprocess.run(sort_cmd, check=True, env=env)
        except (FileNotFoundError, subprocess.CalledProcessError) as e:
            raise CRISPRessoShared.OutputFolderIncompleteException(
                "External sort failed during read counting: %s" % (e,)
            ) from e

    def _write_collapsed_parquet(self, sorted_file: str, workdir: str) -> int:
        """Collapse the sorted keys file and write ``(read_key, count[, quals])`` parquet.

        Returns the number of unique keys written. Two collapse strategies:

        * ``keep_quals=False``: pipe ``sort`` output through ``uniq -c`` (C,
          fast) — but we already sorted, so ``uniq`` (without ``-c``) plus a
          count is done by streaming the sorted file and counting adjacent
          duplicates. We use a pure-Python streaming collapse for both modes
          (O(1) memory, one line of group state) so quals are handled uniformly.
        """
        import pyarrow as pa
        import pyarrow.parquet as pq

        parquet_path = os.path.join(workdir, "read_counts.parquet")
        if self.keep_quals:
            schema = pa.schema(
                [
                    pa.field("read_key", pa.string()),
                    pa.field("count", pa.int64()),
                    pa.field("quals", pa.string()),
                ]
            )
        else:
            schema = pa.schema(
                [
                    pa.field("read_key", pa.string()),
                    pa.field("count", pa.int64()),
                ]
            )
        writer = pq.ParquetWriter(parquet_path, schema)
        num_unique = 0

        BATCH = 50_000
        buf_key: list[str] = []
        buf_count: list[int] = []
        buf_qual: list[Optional[str]] = []

        def flush() -> None:
            nonlocal num_unique
            if not buf_key:
                return
            if self.keep_quals:
                arrays = [pa.array(buf_key, pa.string()), pa.array(buf_count, pa.int64()), pa.array(buf_qual, pa.string())]
                table = pa.Table.from_arrays(arrays, schema=schema)
            else:
                arrays = [pa.array(buf_key, pa.string()), pa.array(buf_count, pa.int64())]
                table = pa.Table.from_arrays(arrays, schema=schema)
            writer.write_table(table)
            num_unique += len(buf_key)
            buf_key.clear()
            buf_count.clear()
            buf_qual.clear()

        prev_key: Optional[str] = None
        cur_count = 0
        cur_qual: Optional[str] = None
        with open(sorted_file, "r", buffering=1 << 20) as sf:
            for line in sf:
                line = line.rstrip("\n")
                if self.keep_quals:
                    # split on the first tab only — quals never contain a tab
                    key, _, quals = line.partition("\t")
                else:
                    key = line
                    quals = None
                if key == prev_key:
                    cur_count += 1
                else:
                    if prev_key is not None:
                        buf_key.append(prev_key)
                        buf_count.append(cur_count)
                        buf_qual.append(cur_qual)
                        if len(buf_key) >= BATCH:
                            flush()
                    prev_key = key
                    cur_count = 1
                    cur_qual = quals  # first occurrence in file order (stable sort)
            if prev_key is not None:
                buf_key.append(prev_key)
                buf_count.append(cur_count)
                buf_qual.append(cur_qual)
                flush()
        writer.close()
        return num_unique


def count_reads_from_fastq(
    fastq1_path: str,
    fastq2_path: Optional[str],
    output_directory: str,
    *,
    keep_quals: bool = False,
    memory_budget_mb: int = DEFAULT_MEMORY_BUDGET_MB,
    force_spill: bool = False,
) -> ReadCounts:
    """Convenience wrapper: create a :class:`VariantStore` and count reads.

    Provided so parity tests and future wiring can call Stage 1 in one line
    without constructing a store explicitly.
    """
    store = VariantStore(
        output_directory,
        memory_budget_mb=memory_budget_mb,
        keep_quals=keep_quals,
    )
    return store.count_reads(fastq1_path, fastq2_path, force_spill=force_spill)


__all__ = [
    "ALIGNED_SCHEMA",
    "COLLAPSED_SCHEMA",
    "CountVectors",
    "DEFAULT_MEMORY_BUDGET_MB",
    "AlignedShardWriter",
    "CollapsedAlleles",
    "ReadCounts",
    "VariantStore",
    "compute_count_vectors_from_collapsed",
    "count_reads_from_fastq",
    "get_slice_from_collapsed",
    "iter_aligned_shard",
    "payload_to_row",
    "row_to_payload",
    "variant_parquet_generator_process",
    "write_allele_frequency_table",
]


# ---------------------------------------------------------------------------
# Stage 2: Worker parquet writer (replaces variants_{i}.tsv JSON output)
# ---------------------------------------------------------------------------
#
# Each alignment worker writes its own parquet shard (aligned_{i}.parquet)
# containing one row per unique read, with the full alignment payload stored
# as native arrow types (arrays as lists, coordinates as list-of-structs) —
# no JSON serialization. The schema is derived once (here) so all shards share
# an identical schema and scan_parquet over the glob concats cleanly (edge #7).
#
# The per-reference payloads (variant_{ref_name} sub-dicts in the current
# code) are stored as a list-of-structs column ``payloads``. Stage 3 (PR 5)
# will explode this column to produce per-allele rows.
#
# NOT wired into CRISPRessoCORE.main yet — exercised by unit tests only.

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

# A (start, end) coordinate — stored as a struct for round-trip fidelity
# (edge #18: insertion/deletion coordinates are tuples of (start, end)).
_COORD_STRUCT = pa.struct([
    pa.field("start", pa.int64()),
    pa.field("end", pa.int64()),
])

# Per-reference payload struct — mirrors the fields of
# CRISPRessoCOREResources.find_indels_substitutions output + the fields added
# by get_new_variant_object. Every field that get_allele_row consumes is here.
_PAYLOAD_STRUCT = pa.struct([
    pa.field("ref_name", pa.string()),
    pa.field("aln_seq", pa.string()),
    pa.field("aln_ref", pa.string()),
    pa.field("aln_strand", pa.string()),
    pa.field("classification", pa.string()),
    pa.field("irregular_ends", pa.bool_()),
    pa.field("insertions_outside_window", pa.int64()),
    pa.field("deletions_outside_window", pa.int64()),
    pa.field("substitutions_outside_window", pa.int64()),
    pa.field("total_mods", pa.int64()),
    pa.field("mods_in_window", pa.int64()),
    pa.field("mods_outside_window", pa.int64()),
    pa.field("insertion_n", pa.int64()),
    pa.field("deletion_n", pa.int64()),
    pa.field("substitution_n", pa.int64()),
    pa.field("ref_positions", pa.list_(pa.int64())),
    pa.field("all_insertion_positions", pa.list_(pa.int64())),
    pa.field("all_insertion_left_positions", pa.list_(pa.int64())),
    pa.field("insertion_positions", pa.list_(pa.int64())),
    pa.field("insertion_coordinates", pa.list_(_COORD_STRUCT)),
    pa.field("insertion_sizes", pa.list_(pa.int64())),
    pa.field("all_deletion_positions", pa.list_(pa.int64())),
    pa.field("deletion_positions", pa.list_(pa.int64())),
    pa.field("deletion_coordinates", pa.list_(_COORD_STRUCT)),
    pa.field("deletion_sizes", pa.list_(pa.int64())),
    pa.field("all_substitution_positions", pa.list_(pa.int64())),
    pa.field("substitution_positions", pa.list_(pa.int64())),
    pa.field("substitution_values", pa.list_(pa.string())),
])

# The full shard schema. One row per unique read.
ALIGNED_SCHEMA = pa.schema([
    pa.field("read_key", pa.string()),
    pa.field("count", pa.int64()),
    pa.field("best_match_score", pa.float64()),
    pa.field("class_name", pa.string()),
    pa.field("best_match_name", pa.string()),
    pa.field("aln_ref_names", pa.list_(pa.string())),
    pa.field("aln_scores", pa.list_(pa.float64())),
    pa.field("caching_is_ok", pa.bool_()),
    pa.field("payloads", pa.list_(_PAYLOAD_STRUCT)),
])

# Fields NOT stored (deliberately dropped — not consumed downstream by
# get_allele_row, and not by the PR 6 count-vector aggregation which streams the
# collapsed allele parquet, whose schema is a superset of get_allele_row):
#   - ref_aln_details          (only for fastq_output annotation; edge #12 follow-up)
#   - all_deletion_coordinates (only for deletions_outside_window, already computed)
#   - all_substitution_values  (needed by all_substitution_base_vectors — the per-base
#     substitution count vector built in the CRISPRessoCORE.main aggregation loop
#     ~line 4060. NOT in get_allele_row, so absent from the collapsed allele
#     parquet. The per-base vectors (all_substitution_base_vectors /
#     all_base_count_vectors) and the windowed vectors (insertion_count_vectors,
#     etc.) are deferred to PR 7's consumer wiring: they feed plot functions
#     directly and require refs/window context. PR 6's Stream-Out B covers the
#     four all_* position vectors + all_indelsub_count_vectors that are persisted
#     to Modification_count_vectors.txt — these need only fields already in
#     COLLAPSED_SCHEMA. Adding all_substitution_values to _PAYLOAD_STRUCT +
#     COLLAPSED_SCHEMA is the prerequisite for the per-base vectors (see PR 7).


def _to_int_list(val):
    """Convert a numpy array, Python list, or scalar to a list of ints."""
    if val is None:
        return None
    if hasattr(val, "tolist"):
        val = val.tolist()
    return [int(x) for x in val]


def _to_str_list(val):
    """Convert a numpy array or list to a list of Python strings."""
    if val is None:
        return None
    if hasattr(val, "tolist"):
        val = val.tolist()
    return [str(x) for x in val]


def _coords_to_structs(tuples):
    """Convert a list of (start, end) tuples to arrow struct dicts."""
    if tuples is None:
        return None
    return [{"start": int(s), "end": int(e)} for s, e in tuples]


def _structs_to_coords(structs):
    """Convert arrow struct dicts back to (start, end) tuples."""
    if structs is None:
        return []
    return [(c["start"], c["end"]) for c in structs]


def _payload_to_struct(payload):
    """Convert a per-reference payload dict to an arrow struct dict.

    Handles numpy arrays (via .tolist()), Python lists, and tuples (coordinates).
    """
    return {
        "ref_name": payload.get("ref_name"),
        "aln_seq": payload.get("aln_seq"),
        "aln_ref": payload.get("aln_ref"),
        "aln_strand": payload.get("aln_strand"),
        "classification": payload.get("classification"),
        "irregular_ends": bool(payload.get("irregular_ends", False)),
        "insertions_outside_window": int(payload.get("insertions_outside_window", 0)),
        "deletions_outside_window": int(payload.get("deletions_outside_window", 0)),
        "substitutions_outside_window": int(payload.get("substitutions_outside_window", 0)),
        "total_mods": int(payload.get("total_mods", 0)),
        "mods_in_window": int(payload.get("mods_in_window", 0)),
        "mods_outside_window": int(payload.get("mods_outside_window", 0)),
        "insertion_n": int(payload.get("insertion_n", 0)),
        "deletion_n": int(payload.get("deletion_n", 0)),
        "substitution_n": int(payload.get("substitution_n", 0)),
        "ref_positions": _to_int_list(payload.get("ref_positions")),
        "all_insertion_positions": _to_int_list(payload.get("all_insertion_positions")),
        "all_insertion_left_positions": _to_int_list(payload.get("all_insertion_left_positions")),
        "insertion_positions": _to_int_list(payload.get("insertion_positions")),
        "insertion_coordinates": _coords_to_structs(payload.get("insertion_coordinates")),
        "insertion_sizes": _to_int_list(payload.get("insertion_sizes")),
        "all_deletion_positions": _to_int_list(payload.get("all_deletion_positions")),
        "deletion_positions": _to_int_list(payload.get("deletion_positions")),
        "deletion_coordinates": _coords_to_structs(payload.get("deletion_coordinates")),
        "deletion_sizes": _to_int_list(payload.get("deletion_sizes")),
        "all_substitution_positions": _to_int_list(payload.get("all_substitution_positions")),
        "substitution_positions": _to_int_list(payload.get("substitution_positions")),
        "substitution_values": _to_str_list(payload.get("substitution_values")),
    }


def _struct_to_payload(struct_row):
    """Convert an arrow struct dict back to a per-reference payload dict.

    Reconstructs numpy arrays for substitution_values (matching the original
    # payload from find_indels_substitutions) and tuples for coordinates.
    """
    return {
        "ref_name": struct_row["ref_name"],
        "aln_seq": struct_row["aln_seq"],
        "aln_ref": struct_row["aln_ref"],
        "aln_strand": struct_row["aln_strand"],
        "classification": struct_row["classification"],
        "irregular_ends": struct_row["irregular_ends"],
        "insertions_outside_window": struct_row["insertions_outside_window"],
        "deletions_outside_window": struct_row["deletions_outside_window"],
        "substitutions_outside_window": struct_row["substitutions_outside_window"],
        "total_mods": struct_row["total_mods"],
        "mods_in_window": struct_row["mods_in_window"],
        "mods_outside_window": struct_row["mods_outside_window"],
        "insertion_n": struct_row["insertion_n"],
        "deletion_n": struct_row["deletion_n"],
        "substitution_n": struct_row["substitution_n"],
        "ref_positions": list(struct_row["ref_positions"]) if struct_row["ref_positions"] else [],
        "all_insertion_positions": list(struct_row["all_insertion_positions"]) if struct_row["all_insertion_positions"] else [],
        "all_insertion_left_positions": list(struct_row["all_insertion_left_positions"]) if struct_row["all_insertion_left_positions"] else [],
        "insertion_positions": list(struct_row["insertion_positions"]) if struct_row["insertion_positions"] else [],
        "insertion_coordinates": _structs_to_coords(struct_row["insertion_coordinates"]),
        "insertion_sizes": list(struct_row["insertion_sizes"]) if struct_row["insertion_sizes"] else [],
        "all_deletion_positions": list(struct_row["all_deletion_positions"]) if struct_row["all_deletion_positions"] else [],
        "deletion_positions": list(struct_row["deletion_positions"]) if struct_row["deletion_positions"] else [],
        "deletion_coordinates": _structs_to_coords(struct_row["deletion_coordinates"]),
        "deletion_sizes": list(struct_row["deletion_sizes"]) if struct_row["deletion_sizes"] else [],
        "all_substitution_positions": list(struct_row["all_substitution_positions"]) if struct_row["all_substitution_positions"] else [],
        "substitution_positions": list(struct_row["substitution_positions"]) if struct_row["substitution_positions"] else [],
        "substitution_values": np.array(struct_row["substitution_values"]) if struct_row["substitution_values"] else np.array([]),
    }


def payload_to_row(read_key, count, payload):
    """Convert a variant payload dict to an arrow row dict matching ALIGNED_SCHEMA.

    The payload is the dict returned by ``get_new_variant_object`` (single-read)
    or ``get_new_variant_object_from_paired`` (paired). For unaligned reads
    (``best_match_score <= 0``), the payload has no ``aln_ref_names`` or
    per-ref sub-dicts — those fields become null in the row.
    """
    aln_ref_names = payload.get("aln_ref_names")
    payloads = None
    if aln_ref_names:
        payloads = []
        for ref_name in aln_ref_names:
            sub = payload.get("variant_" + ref_name)
            if sub is not None:
                payloads.append(_payload_to_struct(sub))
        if not payloads:
            payloads = None
    aln_scores = payload.get("aln_scores", [])
    return {
        "read_key": read_key,
        "count": int(count),
        "best_match_score": float(payload.get("best_match_score", 0)),
        "class_name": payload.get("class_name"),
        "best_match_name": payload.get("best_match_name"),
        "aln_ref_names": list(aln_ref_names) if aln_ref_names else None,
        "aln_scores": [float(s) for s in aln_scores] if aln_scores else None,
        "caching_is_ok": bool(payload.get("caching_is_ok", True)),
        "payloads": payloads,
    }


def row_to_payload(row):
    """Convert a parquet row (dict) back to a variant payload dict.

    This is the inverse of :func:`payload_to_row` — used by Stage 3 (collapse)
    to reconstruct payloads from parquet shards. Reconstructs numpy arrays for
    ``substitution_values`` and tuples for coordinates to match the original
    payload shape from ``find_indels_substitutions``.
    """
    payloads = row.get("payloads")
    payload = {
        "count": row["count"],
        "best_match_score": row["best_match_score"],
        "class_name": row.get("class_name"),
        "best_match_name": row.get("best_match_name"),
        "aln_ref_names": list(row["aln_ref_names"]) if row.get("aln_ref_names") else [],
        "aln_scores": list(row["aln_scores"]) if row.get("aln_scores") else [],
        "caching_is_ok": row.get("caching_is_ok", True),
    }
    if payloads:
        for sub in payloads:
            sub_payload = _struct_to_payload(sub)
            ref_name = sub.get("ref_name")
            if ref_name is not None:
                payload["variant_" + ref_name] = sub_payload
    return payload


class AlignedShardWriter:
    """Writes aligned payload rows to a parquet shard in batches.

    Wraps a :class:`pyarrow.parquet.ParquetWriter` with the shared
    :data:`ALIGNED_SCHEMA` so all shards across workers are schema-identical
    (edge #7). Rows are buffered and flushed in batches to amortize I/O.
    """

    def __init__(self, path, *, batch_size=10_000, schema=ALIGNED_SCHEMA):
        self.path = str(path)
        self._schema = schema
        self._writer = pq.ParquetWriter(self.path, schema)
        self._batch_size = batch_size
        self._buf = []

    def write_row(self, row):
        """Buffer a row dict (from :func:`payload_to_row`) for batched write."""
        self._buf.append(row)
        if len(self._buf) >= self._batch_size:
            self._flush()

    def _flush(self):
        if not self._buf:
            return
        table = pa.Table.from_pylist(self._buf, schema=self._schema)
        self._writer.write_table(table)
        self._buf.clear()

    def close(self):
        """Flush remaining rows and close the writer."""
        self._flush()
        if self._writer is not None:
            self._writer.close()
            self._writer = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


def iter_aligned_shard(path):
    """Stream ``(read_key, count, payload)`` from a parquet shard.

    Used by Stage 3 (collapse) to iterate aligned payloads without
    materializing the whole shard. Reads row-group by row-group via
    ``iter_batches`` so peak memory is bounded by one batch.
    """
    pf = pq.ParquetFile(str(path))
    for batch in pf.iter_batches(batch_size=50_000):
        table = pa.Table.from_batches([batch], schema=ALIGNED_SCHEMA)
        for row in table.to_pylist():
            yield row["read_key"], row["count"], row_to_payload(row)


def variant_parquet_generator_process(
    read_items,
    get_new_variant_object,
    args,
    refs,
    ref_names,
    aln_matrix,
    pe_scaffold_dna_info,
    process_id,
    output_dir,
    is_paired=False,
):
    """Multiprocessing target: align reads and write a parquet shard.

    Mirrors :func:`CRISPRessoCORE.variant_file_generator_process` but writes
    parquet (``aligned_{process_id}.parquet``) instead of JSON TSV. Active only
    under ``--storage_backend parquet`` (the caller decides whether to invoke
    this or the TSV path). The TSV path is untouched.

    Parameters
    ----------
    read_items : list of (read_key, count, quals_or_None)
        The chunk of counted reads for this worker (from ReadCounts.items()).
    get_new_variant_object : callable
        The alignment function (``get_new_variant_object`` for single-read,
        ``get_new_variant_object_from_paired`` for paired).
    args, refs, ref_names, aln_matrix, pe_scaffold_dna_info
        Same as the TSV worker — passed through to the alignment function.
    process_id : int
        Worker ID (determines shard filename).
    output_dir : str
        Directory to write the shard.
    is_paired : bool
        If True, splits read_key on '+' and quals on ' ' to get R1/R2 (matching
        ``args.crispresso_merge`` in the TSV worker).

    """
    shard_path = os.path.join(output_dir, f"aligned_{process_id}.parquet")
    num_processed = 0
    with AlignedShardWriter(shard_path) as writer:
        for index, (read_key, count, quals) in enumerate(read_items):
            num_processed = index + 1
            if is_paired:
                fastq1_seq, fastq2_seq = read_key.split("+")
                if quals is not None:
                    fastq1_qual, fastq2_qual = quals.split(" ")
                else:
                    fastq1_qual, fastq2_qual = "", ""
                new_variant = get_new_variant_object(
                    args, fastq1_seq, fastq2_seq, fastq1_qual, fastq2_qual,
                    refs, ref_names, aln_matrix, pe_scaffold_dna_info,
                )
            else:
                new_variant = get_new_variant_object(
                    args, read_key, refs, ref_names, aln_matrix, pe_scaffold_dna_info,
                )
            new_variant["count"] = count
            row = payload_to_row(read_key, count, new_variant)
            writer.write_row(row)
            if index % 10000 == 0 and index != 0:
                info(f"Process {process_id + 1} has processed {index} unique reads", {"percent_complete": 10})
    info(f"Process {process_id + 1} has finished processing {num_processed} unique reads", {"percent_complete": 10})


# ---------------------------------------------------------------------------
# Stage 3: Streaming collapse (aln_seq re-key + RC canonical-key merge +
#          multi-reference fan-out) — PR 5.
# ---------------------------------------------------------------------------
#
# ``VariantStore.collapse`` consumes the per-worker aligned parquet shards
# produced by Stage 2 and produces the collapsed allele table: one row per
# (allele, reference) with the final ``#Reads`` (post RC-merge), plus the
# variant-level aggregations (``n_total``, ``class_counts``) and per-reference
# counts (``counts_total/modified/unmodified/discarded``) that the pandas path
# builds in the ``CRISPRessoCORE.main`` allele loop (~line 3976).
#
# The three sub-stages match the pandas path exactly:
#
#   3a — aln_seq re-key (paired only). Reads keyed by R1+'+'+RC(R2) are
#        re-keyed by the primary reference's aligned sequence (aln_seq), summing
#        ``count``. Two different raw pairs that align identically collapse to
#        one allele. For single-read input there is NO re-key (the key is the
#        raw read sequence, already deduped by Stage 1) — matching
#        ``process_fastq`` which skips the ``if '+' in key`` re-key block.
#
#   3b — reverse-complement merge. For each key, if ``reverse_complement(key)``
#        is also present with count > 0, the RC's count is folded into the
#        current key and the RC is zeroed. The representative is the key that
#        appears FIRST in iteration order — the exact pandas tie-break
#        (insertion order, edge #8). Iteration order here = shard/scan order,
#        which equals FASTQ first-occurrence order on the eager Stage 1 path
#        (parity) and sorted-key order on the spill path (deterministic; the
#        representative's payload is equivalent to the FASTQ-first one for the
#        common case — see ``_collapse_rc_merge`` docstring).
#
#   3c — multi-reference fan-out. Each collapsed variant produces one allele
#        row per reference it aligned to (``get_allele_row`` shape), with the
#        AMBIGUOUS (one row, ``AMBIGUOUS_<ref0>``) and ``discard_indel_reads``
#        (one row per ref, ``DISCARDED_<ref0>``) special cases replicated
#        verbatim from the pandas loop.
#
# Memory: the input shards are streamed row-by-row via ``iter_aligned_shard``
# (one pyarrow batch in flight). The collapse group table (3a/3b) is held in a
# plain Python dict keyed by the collapse key; its size is bounded by the
# *unique allele count* for paired input (where 3a collapses before 3b), which
# is the design's stated bounded term. For single-read input with ~1:1
# diversity (no re-key, no collapse before 3b), the dict is bounded by the
# unique *read* count — the same root cause the spike (S1) identified for
# Stage 1. Flat memory for that case requires the external-sort canonical-key
# merge (write post-3a to parquet, ``sort(canonical_key)`` which spills, then a
# streaming first-per-group pass); that is a follow-up and not in scope for PR
# 5, whose gate is parity vs the pandas ``get_allele_row`` / RC-merge logic on
# canned payloads. The 3a/3b/3c methods are split so 3b can be swapped for the
# streaming version without touching 3a/3c.
#
# NOT wired into ``CRISPRessoCORE.main`` — exercised by unit tests only.

# Keys present in the non-detailed ``get_allele_row`` branch (the
# ``else`` of ``get_allele_row``, ~line 3960). The detailed branch is a
# superset; the collapsed parquet always stores the full (detailed) column set
# so downstream (PR 6 TSV sink, PR 7 df_alleles) can project either shape.
_NON_DETAILED_ALLELE_KEYS = (
    "#Reads",
    "Aligned_Sequence",
    "Reference_Sequence",
    "n_inserted",
    "n_deleted",
    "n_mutated",
    "Reference_Name",
    "Read_Status",
    "Aligned_Reference_Names",
    "Aligned_Reference_Scores",
    "ref_positions",
)

# Collapsed-allele parquet schema (the persisted artifact). One row per
# (allele, reference) — i.e. the post-3c exploded table. Always stores the
# full/detailed column set; consumers project.
COLLAPSED_SCHEMA = pa.schema([
    pa.field("#Reads", pa.int64()),
    pa.field("Aligned_Sequence", pa.string()),
    pa.field("Reference_Sequence", pa.string()),
    pa.field("n_inserted", pa.int64()),
    pa.field("n_deleted", pa.int64()),
    pa.field("n_mutated", pa.int64()),
    pa.field("Reference_Name", pa.string()),
    pa.field("Read_Status", pa.string()),
    pa.field("Aligned_Reference_Names", pa.string()),
    pa.field("Aligned_Reference_Scores", pa.string()),
    pa.field("ref_positions", pa.list_(pa.int64())),
    pa.field("all_insertion_positions", pa.list_(pa.int64())),
    pa.field("all_insertion_left_positions", pa.list_(pa.int64())),
    pa.field("insertion_positions", pa.list_(pa.int64())),
    pa.field("insertion_coordinates", pa.list_(_COORD_STRUCT)),
    pa.field("insertion_sizes", pa.list_(pa.int64())),
    pa.field("all_deletion_positions", pa.list_(pa.int64())),
    pa.field("deletion_positions", pa.list_(pa.int64())),
    pa.field("deletion_coordinates", pa.list_(_COORD_STRUCT)),
    pa.field("deletion_sizes", pa.list_(pa.int64())),
    pa.field("all_substitution_positions", pa.list_(pa.int64())),
    pa.field("substitution_positions", pa.list_(pa.int64())),
    pa.field("substitution_values", pa.list_(pa.string())),
])


def _get_allele_row(
    reference_name: str,
    variant_count: int,
    aln_ref_names_str: str,
    aln_ref_scores_str: str,
    variant_payload: dict,
    write_detailed: bool,
) -> dict:
    """Replica of the nested ``get_allele_row`` in ``CRISPRessoCORE.main``.

    Returns the exact dict shape the pandas path appends to ``alleles_list``
    (both the detailed and non-detailed branches), so a direct dict comparison
    against the pandas output is a valid parity check.
    """
    if write_detailed:
        return {
            "#Reads": variant_count,
            "Aligned_Sequence": variant_payload["aln_seq"],
            "Reference_Sequence": variant_payload["aln_ref"],
            "n_inserted": variant_payload["insertion_n"],
            "n_deleted": variant_payload["deletion_n"],
            "n_mutated": variant_payload["substitution_n"],
            "Reference_Name": reference_name,
            "Read_Status": variant_payload["classification"],
            "Aligned_Reference_Names": aln_ref_names_str,
            "Aligned_Reference_Scores": aln_ref_scores_str,
            "ref_positions": variant_payload["ref_positions"],
            "all_insertion_positions": variant_payload["all_insertion_positions"],
            "all_insertion_left_positions": variant_payload["all_insertion_left_positions"],
            "insertion_positions": variant_payload["insertion_positions"],
            "insertion_coordinates": variant_payload["insertion_coordinates"],
            "insertion_sizes": variant_payload["insertion_sizes"],
            "all_deletion_positions": variant_payload["all_deletion_positions"],
            "deletion_positions": variant_payload["deletion_positions"],
            "deletion_coordinates": variant_payload["deletion_coordinates"],
            "deletion_sizes": variant_payload["deletion_sizes"],
            "all_substitution_positions": variant_payload["all_substitution_positions"],
            "substitution_positions": variant_payload["substitution_positions"],
            "substitution_values": variant_payload["substitution_values"],
        }
    return {
        "#Reads": variant_count,
        "Aligned_Sequence": variant_payload["aln_seq"],
        "Reference_Sequence": variant_payload["aln_ref"],
        "n_inserted": variant_payload["insertion_n"],
        "n_deleted": variant_payload["deletion_n"],
        "n_mutated": variant_payload["substitution_n"],
        "Reference_Name": reference_name,
        "Read_Status": variant_payload["classification"],
        "Aligned_Reference_Names": aln_ref_names_str,
        "Aligned_Reference_Scores": aln_ref_scores_str,
        "ref_positions": variant_payload["ref_positions"],
    }


def _allele_dict_to_parquet_row(d: dict) -> dict:
    """Convert a (always-full) allele-row dict to a COLLAPSED_SCHEMA row.

    Position arrays become int64 lists; coordinates become list-of-struct
    dicts; ``substitution_values`` (a numpy array in the payload) becomes a
    list of strings — matching the arrow storage used by the aligned shards.
    """
    return {
        "#Reads": int(d["#Reads"]),
        "Aligned_Sequence": d["Aligned_Sequence"],
        "Reference_Sequence": d["Reference_Sequence"],
        "n_inserted": int(d["n_inserted"]),
        "n_deleted": int(d["n_deleted"]),
        "n_mutated": int(d["n_mutated"]),
        "Reference_Name": d["Reference_Name"],
        "Read_Status": d["Read_Status"],
        "Aligned_Reference_Names": d["Aligned_Reference_Names"],
        "Aligned_Reference_Scores": d["Aligned_Reference_Scores"],
        "ref_positions": _to_int_list(d["ref_positions"]),
        "all_insertion_positions": _to_int_list(d["all_insertion_positions"]),
        "all_insertion_left_positions": _to_int_list(d["all_insertion_left_positions"]),
        "insertion_positions": _to_int_list(d["insertion_positions"]),
        "insertion_coordinates": _coords_to_structs(d["insertion_coordinates"]),
        "insertion_sizes": _to_int_list(d["insertion_sizes"]),
        "all_deletion_positions": _to_int_list(d["all_deletion_positions"]),
        "deletion_positions": _to_int_list(d["deletion_positions"]),
        "deletion_coordinates": _coords_to_structs(d["deletion_coordinates"]),
        "deletion_sizes": _to_int_list(d["deletion_sizes"]),
        "all_substitution_positions": _to_int_list(d["all_substitution_positions"]),
        "substitution_positions": _to_int_list(d["substitution_positions"]),
        "substitution_values": _to_str_list(d["substitution_values"]),
    }


@dataclass
class CollapsedAlleles:
    """Result of Stage 3 collapse (``VariantStore.collapse``).

    Attributes
    ----------
    allele_rows : list[dict]
        Sorted (``#Reads`` desc, ``Aligned_Sequence`` asc, ``Reference_Sequence``
        asc — matching ``df_alleles.sort_values`` ~line 4303) list of
        ``get_allele_row``-shaped dicts. The dict shape respects
        ``write_detailed_allele_table``/``vcf_output`` (detailed vs the
        ``crispresso2Cols``-plus-``ref_positions`` subset).
    n_total : int
        ``N_TOTAL`` — sum of post-RC-merge counts over all variants (excludes
        unaligned reads, which are filtered before collapse). parity with the
        pandas ``N_TOTAL`` used for ``%Reads``.
    class_counts : dict[str, int]
        Variant-level ``class_name`` → merged count (the pie-chart counts).
        Computed pre-explode (once per variant), so multi-reference variants
        are counted once — matching the pandas loop. (Count *vectors* and the
        TSV sink are PR 6; they consume ``allele_rows``.)
    counts_total, counts_modified, counts_unmodified, counts_discarded : dict[str, int]
        Per-reference scalar counts from the 3c fan-out (matching
        ``counts_total``/``counts_modified``/``counts_unmodified``/
        ``counts_discarded`` in the pandas loop). ``counts_discarded`` is keyed
        by the actual reference (the one whose payload had indels), matching
        the pandas ``counts_discarded[ref_name] += variant_count``.
    parquet_path : Optional[str]
        Path to the persisted ``collapsed.allele.parquet`` artifact (full
        column set), or ``None`` if ``write_parquet=False``.

    """

    allele_rows: list
    n_total: int
    class_counts: dict
    counts_total: dict
    counts_modified: dict
    counts_unmodified: dict
    counts_discarded: dict
    parquet_path: Optional[str] = None

    def allele_rows_dataframe(self):
        """Return ``allele_rows`` as a pandas DataFrame (sorted, with ``%Reads``).

        Convenience for parity tests and the eventual df_alleles wiring (PR 7).
        ``%Reads`` = ``#Reads / n_total * 100`` and ``n_*`` cast to int —
        matching ``CRISPRessoCORE.main`` ~line 4301-4303.
        """
        import pandas as pd

        df = pd.DataFrame(self.allele_rows)
        if len(df) and self.n_total > 0:
            df["%Reads"] = df["#Reads"] / self.n_total * 100
        elif len(df):
            df["%Reads"] = 0.0
        if len(df):
            df[["n_deleted", "n_inserted", "n_mutated"]] = df[
                ["n_deleted", "n_inserted", "n_mutated"]
            ].astype(int)
        return df


def _expand_shard_paths(shard_paths) -> list:
    """Normalise ``shard_paths`` (list or glob string) to a sorted file list."""
    if isinstance(shard_paths, str):
        import glob

        return sorted(glob.glob(shard_paths))
    return [str(p) for p in shard_paths]


# Add Stage 3 methods to VariantStore. Defined as module-level functions and
# bound onto the class so the Stage 1/2 class body (above) stays readable.


def _collapse(
    self: "VariantStore",
    shard_paths,
    *,
    is_paired: bool,
    expand_ambiguous_alignments: bool = False,
    discard_indel_reads: bool = False,
    write_detailed_allele_table: bool = False,
    vcf_output: bool = False,
    write_parquet: bool = True,
    collapsed_path: Optional[str] = None,
) -> CollapsedAlleles:
    """Stage 3: collapse aligned shards into the sorted allele table.

    Parameters
    ----------
    shard_paths
        List of aligned-shard paths (from Stage 2) or a glob string
        (e.g. ``"aligned_*.parquet"``). Streamed row-by-row via
        :func:`iter_aligned_shard`.
    is_paired
        ``True`` for paired/``crispresso_merge`` input (perform the 3a
        aln_seq re-key); ``False`` for single-read / bam input (skip re-key —
        the key is the raw read sequence).
    expand_ambiguous_alignments, discard_indel_reads
        Replicate the corresponding ``args`` flags in the 3c fan-out.
    write_detailed_allele_table, vcf_output
        Select the ``get_allele_row`` branch (detailed vs subset) for the
        returned ``allele_rows``. The persisted parquet always stores the full
        column set.
    write_parquet, collapsed_path
        If ``write_parquet`` (default), write ``collapsed.allele.parquet`` to
        ``collapsed_path`` (or ``<output_directory>/collapsed.allele.parquet``)
        and return its path on :class:`CollapsedAlleles`.

    Notes
    -----
    * Unaligned reads (``best_match_score <= 0``) are filtered out before
      collapse — they never enter ``variantCache`` in the pandas path (removed
      during read-back). The ``not_aln``/fastq_output annotation path is a
      wiring concern (PR 7), not collapse.
    * ``caching_is_ok == False`` (paired re_aln) reads are expected to be
      pulled aside by the wiring (the pandas ``re_aln`` loop) before collapse;
      collapse assumes its input is the post-re_aln aligned set. Canned
      payloads use ``caching_is_ok=True``.
    * The empty-seq count=0 sentinel (pandas ``variantCache['']``) is a no-op
      (skipped by ``count == 0``) and never enters shards, so it is not
      materialised here (OQ-E5).

    """
    write_detailed = bool(write_detailed_allele_table or vcf_output)
    paths = _expand_shard_paths(shard_paths)

    # 3a + 3b: stream shards into the in-memory collapse table.
    store = self._collapse_rekey_and_rcmerge(paths, is_paired=is_paired)

    # 3c: fan out to per-reference allele rows + per-variant/per-ref aggregations.
    # allele_rows are built as FULL (detailed) dicts — the parquet needs the
    # position arrays for PR 6 count vectors even when the allele TSV is
    # non-detailed.
    full_rows, n_total, class_counts, counts_total, counts_modified, \
        counts_unmodified, counts_discarded = self._collapse_fanout(
            store, discard_indel_reads=discard_indel_reads,
        )

    # Sort matching df_alleles.sort_values(~line 4303): #Reads desc,
    # Aligned_Sequence asc, Reference_Sequence asc. Python's sort is stable,
    # so full ties preserve the fan-out (variant-then-ref) insertion order —
    # the same order the pandas alleles_list is built in.
    full_rows.sort(
        key=lambda r: (-int(r["#Reads"]), r["Aligned_Sequence"], r["Reference_Sequence"])
    )

    # Persist the collapsed allele table (full column set).
    parquet_path = None
    if write_parquet:
        parquet_path = (
            collapsed_path
            or os.path.join(self.output_directory, "collapsed.allele.parquet")
        )
        _write_collapsed_allele_parquet(full_rows, parquet_path)

    # The returned allele_rows respect the get_allele_row branch (detailed vs
    # the non-detailed subset), for direct parity comparison and the PR 6 TSV
    # sink. The persisted parquet is always full.
    if write_detailed:
        allele_rows = full_rows
    else:
        allele_rows = [
            {k: r[k] for k in _NON_DETAILED_ALLELE_KEYS if k in r} for r in full_rows
        ]

    return CollapsedAlleles(
        allele_rows=allele_rows,
        n_total=n_total,
        class_counts=class_counts,
        counts_total=counts_total,
        counts_modified=counts_modified,
        counts_unmodified=counts_unmodified,
        counts_discarded=counts_discarded,
        parquet_path=parquet_path,
    )


def _collapse_rekey_and_rcmerge(self: "VariantStore", paths: list, *, is_paired: bool) -> "dict":
    """Stage 3a (re-key) + 3b (RC merge): stream shards into the collapse table.

    Returns an insertion-ordered dict ``{collapse_key: _Record}`` where each
    ``_Record`` carries the merged ``count`` and the representative payload
    (first occurrence in scan order — the pandas tie-break).

    3a: for paired input the collapse key is the primary reference's
    ``aln_seq`` (``payload['variant_'+aln_ref_names[0]]['aln_seq']``); for
    single-read the key is the raw read sequence (``read_key``), i.e. no
    re-key — matching ``process_fastq`` which skips the ``if '+' in key`` block.

    3b: iterate the post-3a dict in insertion order; for each key with
    ``count > 0``, fold ``reverse_complement(key)``'s count into it (if present
    and positive) and zero the RC. The representative is the first-iterated
    key — the exact pandas tie-break. The ``rc == key`` (palindrome) case is
    replicated verbatim (pandas doubles the count there); it is a latent pandas
    behaviour, preserved for parity.

    Memory: the dict is bounded by the unique collapse-key count. For paired
    input that is the unique *allele* count (small — the design's bounded
    term). For single-read input with no re-key it is the unique *read* count;
    flat memory for that case is a follow-up (external-sort canonical-key
    merge), see the module docstring.
    """
    rc = CRISPRessoShared.reverse_complement
    store: dict = {}

    def _primary_aln_seq(payload: dict) -> str:
        ref0 = payload["aln_ref_names"][0]
        return payload["variant_" + ref0]["aln_seq"]

    for path in paths:
        for read_key, count, payload in iter_aligned_shard(path):
            # Skip unaligned reads (best_match_score <= 0): they are removed
            # from variantCache before the pandas collapse loop.
            if payload.get("best_match_score", 0) <= 0:
                continue
            if is_paired:
                key = _primary_aln_seq(payload)
            else:
                key = read_key
            rec = store.get(key)
            if rec is None:
                store[key] = {"count": int(count), "payload": payload}
            else:
                rec["count"] += int(count)
                # keep first-occurrence payload (pandas: variantCache[new_key]
                # is set once, later occurrences only add count)

    # 3b: RC merge in insertion order.
    for key in list(store.keys()):
        rec = store[key]
        variant_count = rec["count"]
        if variant_count == 0:
            continue
        rc_key = rc(key)
        rc_rec = store.get(rc_key)
        if rc_rec is not None and rc_rec["count"] > 0:
            variant_count += rc_rec["count"]
            rc_rec["count"] = 0
            rec["count"] = variant_count
    return store


def _collapse_fanout(
    self: "VariantStore",
    store: dict,
    *,
    discard_indel_reads: bool,
):
    """Stage 3c: explode collapsed variants into per-reference allele rows.

    Replicates ``CRISPRessoCORE.main`` ~line 3976-4045 verbatim:

    * Skip zero-count entries (RC-merged-away partners).
    * ``class_name == "AMBIGUOUS"`` → one row ``AMBIGUOS_<ref0>`` using
      ``variant_<ref0>``'s payload; do NOT iterate refs.
    * Otherwise iterate ``aln_ref_names``; if ``discard_indel_reads`` and the
      per-ref payload has indels → one ``DISCARDED_<ref0>`` row (using that
      ref's payload) and ``counts_discarded[ref] += count``.
    * Else a normal row ``<ref>`` and ``counts_total/modified/unmodified``.

    Returns ``(allele_rows, n_total, class_counts, counts_total,
    counts_modified, counts_unmodified, counts_discarded)``.

    ``allele_rows`` are always built as *full* (detailed) dicts — the persisted
    parquet must carry the position arrays for PR 6's count-vector aggregation
    even when ``write_detailed_allele_table`` is False (the pandas count-vector
    loop runs unconditionally on every aligned variant). ``collapse`` subsets
    the returned list to the non-detailed shape when appropriate.
    """
    allele_rows: list = []
    n_total = 0
    class_counts: dict = {}
    counts_total: dict = {}
    counts_modified: dict = {}
    counts_unmodified: dict = {}
    counts_discarded: dict = {}

    for key in store:  # insertion order
        rec = store[key]
        variant_count = rec["count"]
        if variant_count == 0:
            continue
        n_total += variant_count

        payload = rec["payload"]
        class_name = payload.get("class_name")
        class_counts[class_name] = class_counts.get(class_name, 0) + variant_count

        aln_ref_names = payload["aln_ref_names"]
        aln_ref_names_str = "&".join(aln_ref_names)
        aln_ref_scores = payload.get("aln_scores", [])
        aln_ref_scores_str = "&".join([str(x) for x in aln_ref_scores])

        if class_name == "AMBIGUOUS":
            variant_payload = payload["variant_" + aln_ref_names[0]]
            allele_rows.append(
                _get_allele_row(
                    "AMBIGUOUS_" + aln_ref_names[0],
                    variant_count, aln_ref_names_str, aln_ref_scores_str,
                    variant_payload, write_detailed=True,
                )
            )
            continue

        for ref_name in aln_ref_names:
            variant_payload = payload["variant_" + ref_name]
            if discard_indel_reads and (
                variant_payload["deletion_n"] > 0 or variant_payload["insertion_n"] > 0
            ):
                counts_discarded[ref_name] = counts_discarded.get(ref_name, 0) + variant_count
                allele_rows.append(
                    _get_allele_row(
                        "DISCARDED_" + aln_ref_names[0],
                        variant_count, aln_ref_names_str, aln_ref_scores_str,
                        variant_payload, write_detailed=True,
                    )
                )
                continue
            counts_total[ref_name] = counts_total.get(ref_name, 0) + variant_count
            if variant_payload["classification"] == "MODIFIED":
                counts_modified[ref_name] = counts_modified.get(ref_name, 0) + variant_count
            else:
                counts_unmodified[ref_name] = counts_unmodified.get(ref_name, 0) + variant_count
            allele_rows.append(
                _get_allele_row(
                    ref_name, variant_count, aln_ref_names_str, aln_ref_scores_str,
                    variant_payload, write_detailed=True,
                )
            )

    return (
        allele_rows, n_total, class_counts,
        counts_total, counts_modified, counts_unmodified, counts_discarded,
    )


def _write_collapsed_allele_parquet(allele_rows: list, path: str) -> None:
    """Write the full-column collapsed allele table to parquet.

    ``allele_rows`` are full (detailed) dicts (see ``_collapse_fanout``); the
    persisted artifact always carries the complete column set so PR 6's
    count-vector aggregation and PR 7's df_alleles view can project either the
    detailed or ``crispresso2Cols`` shape.
    """
    schema = COLLAPSED_SCHEMA
    writer = pq.ParquetWriter(path, schema)
    try:
        if not allele_rows:
            return
        rows = [_allele_dict_to_parquet_row(d) for d in allele_rows]
        table = pa.Table.from_pylist(rows, schema=schema)
        writer.write_table(table)
    finally:
        writer.close()


# Bind Stage 3 onto VariantStore.
VariantStore.collapse = _collapse
VariantStore._collapse_rekey_and_rcmerge = _collapse_rekey_and_rcmerge
VariantStore._collapse_fanout = _collapse_fanout


def collapse_aligned_shards(
    shard_paths,
    output_directory: str,
    *,
    is_paired: bool,
    expand_ambiguous_alignments: bool = False,
    discard_indel_reads: bool = False,
    write_detailed_allele_table: bool = False,
    vcf_output: bool = False,
    write_parquet: bool = True,
    collapsed_path: Optional[str] = None,
    memory_budget_mb: int = DEFAULT_MEMORY_BUDGET_MB,
) -> CollapsedAlleles:
    """Convenience wrapper: create a :class:`VariantStore` and collapse shards.

    Provided so parity tests and future wiring can call Stage 3 in one line.
    """
    store = VariantStore(output_directory, memory_budget_mb=memory_budget_mb)
    return store.collapse(
        shard_paths,
        is_paired=is_paired,
        expand_ambiguous_alignments=expand_ambiguous_alignments,
        discard_indel_reads=discard_indel_reads,
        write_detailed_allele_table=write_detailed_allele_table,
        vcf_output=vcf_output,
        write_parquet=write_parquet,
        collapsed_path=collapsed_path,
    )


# ---------------------------------------------------------------------------
# Stage 4 — Stream-Out B: count-vector aggregation (PR 6)
# ---------------------------------------------------------------------------
#
# The per-reference position count vectors that ``CRISPRessoCORE.main`` builds
# in the variant-cache aggregation loop (~line 4025-4049) and persists to
# ``Modification_count_vectors.txt`` via ``save_count_vectors_to_file``
# (~line 4689):
#
#   * ``all_insertion_count_vectors``        — 1 at each amplicon position
#     bracketing an insertion (both sides), weighted by #Reads.
#   * ``all_insertion_left_count_vectors``   — 1 at the left position only.
#   * ``all_deletion_count_vectors``         — 1 at each deleted position.
#   * ``all_substitution_count_vectors``     — 1 at each substituted position.
#   * ``all_indelsub_count_vectors``         — sum of ins+del+sub (computed
#     after the pass, matching ~line 4200).
#
# Under Approach B these are built by a *streaming pass* over the collapsed
# allele parquet (Stage 3 output): scan row-group by row-group (pyarrow
# ``iter_batches``), and for each (allele, reference) row add the #Reads-
# weighted contribution to the numpy accumulator for that reference. Peak RSS
# is bounded by the accumulator size — O(amplicon_length x num_refs) — not by
# the read count (the design's stated bounded lower-order term).
#
# Parity: the arithmetic is identical to the pandas loop, including the numpy
# fancy-indexing ``arr[positions] += count`` expression (which has the same
# duplicate-index semantics as pandas — positions are distinct in practice, so
# this is a no-op difference, but replicated verbatim for safety). AMBIGUOUS
# and DISCARDED rows do not contribute (the pandas loop ``continue``\ s before
# the count-vector lines); here they are skipped because their
# ``Reference_Name`` (``AMBIGUOUS_<ref0>`` / ``DISCARDED_<ref0>``) is not a key
# in ``ref_lengths``.
#
# The per-base vectors (``all_substitution_base_vectors``,
# ``all_base_count_vectors``) and the quantification-window vectors
# (``insertion_count_vectors``, etc.) are deferred to PR 7 — see the
# _PAYLOAD_STRUCT comment above.
#
# NOT wired into CRISPRessoCORE.main — exercised by unit tests only.


@dataclass
class CountVectors:
    """Per-reference position count vectors (Stream-Out B).

    Mirrors the four ``all_*`` accumulator dicts plus the derived
    ``all_indelsub_count_vectors`` that ``CRISPRessoCORE.main`` builds, and the
    per-reference ``counts_total`` used for the ``Total`` row of
    ``Modification_count_vectors.txt``. All vectors are ``float64`` numpy arrays
    of length ``ref_lengths[ref_name]`` — matching ``np.zeros(seq_len)`` in the
    pandas path, so :meth:`save_modification_count_vectors` produces byte-
    identical output to ``save_count_vectors_to_file``.
    """

    all_insertion_count_vectors: dict
    all_insertion_left_count_vectors: dict
    all_deletion_count_vectors: dict
    all_substitution_count_vectors: dict
    all_indelsub_count_vectors: dict
    counts_total: dict
    ref_lengths: dict

    def save_modification_count_vectors(self, ref_name, ref_seq, path):
        """Write ``Modification_count_vectors.txt`` for one reference.

        Byte-identical replica of ``CRISPRessoCORE.save_count_vectors_to_file``
        (~line 4613) as called for the all-modifications file (~line 4689):
        first row is ``Sequence\t<ref_seq chars>``; subsequent rows are
        ``<vectorName>\t<tab-joined vector values>`` for Insertions,
        Insertions_Left, Deletions, Substitutions, All_modifications, Total.
        Vector values are written with ``str()`` (float64 → ``"5.0"``); the
        Total row is a Python int repeated (``"5"``) — matching the pandas path
        where ``counts_total`` is a Python int and the vectors are float64.
        """
        ref_len = self.ref_lengths[ref_name]
        total = self.counts_total.get(ref_name, 0)
        vectors = [
            self.all_insertion_count_vectors[ref_name],
            self.all_insertion_left_count_vectors[ref_name],
            self.all_deletion_count_vectors[ref_name],
            self.all_substitution_count_vectors[ref_name],
            self.all_indelsub_count_vectors[ref_name],
            [total] * ref_len,
        ]
        names = [
            "Insertions",
            "Insertions_Left",
            "Deletions",
            "Substitutions",
            "All_modifications",
            "Total",
        ]
        with open(path, "w") as outfile:
            outfile.write("Sequence\t" + "\t".join(list(ref_seq)) + "\n")
            for vector, name in zip(vectors, names):
                outfile.write(name + "\t" + "\t".join([str(x) for x in vector]) + "\n")


def _compute_count_vectors(
    self: "VariantStore",
    ref_lengths: dict,
    *,
    collapsed_path: Optional[str] = None,
    batch_size: int = 50_000,
) -> CountVectors:
    """Stream the collapsed allele parquet → numpy count-vector accumulators.

    Parameters
    ----------
    ref_lengths
        ``{ref_name: sequence_length}`` — the amplicon length per reference
        (from ``refs[ref_name]['sequence_length']``). Defines the accumulator
        vector lengths and the set of valid reference names (rows whose
        ``Reference_Name`` is not in this dict — i.e. ``AMBIGUOUS_*`` /
        ``DISCARDED_*`` — are skipped, matching the pandas loop's ``continue``
        before the count-vector lines).
    collapsed_path
        Path to the collapsed allele parquet (Stage 3 output). Defaults to
        ``<output_directory>/collapsed.allele.parquet``.
    batch_size
        Row batch size for ``pyarrow.iter_batches`` (bounds peak memory).

    Returns
    -------
    CountVectors
        The four ``all_*`` accumulators + ``all_indelsub_count_vectors`` +
        per-reference ``counts_total`` (sum of #Reads over contributing rows),
        matching the pandas loop's outputs.

    Notes
    -----
    * Vectors are ``np.zeros(seq_len)`` (float64), exactly like the pandas path
      (~line 3874), so :meth:`CountVectors.save_modification_count_vectors`
      is byte-identical to ``save_count_vectors_to_file``.
    * ``arr[positions] += count`` uses numpy fancy indexing — the same
      expression as pandas. Empty/None position lists are no-ops.
    * ``counts_total`` is accumulated here (sum of #Reads over contributing
      rows, initialized to 0 for every ref — matching pandas ~line 3870 which
      zeroes ``counts_total[ref_name]`` for all refs) rather than taken from
      :class:`CollapsedAlleles` so this method is self-contained. It equals
      ``CollapsedAlleles.counts_total`` for refs with contributing rows;
      ``CollapsedAlleles.counts_total`` omits zero-count refs (a PR 5 gap),
      while this dict includes them as 0 (pandas parity).
    """
    path = collapsed_path or os.path.join(self.output_directory, "collapsed.allele.parquet")
    ins = {r: np.zeros(n) for r, n in ref_lengths.items()}
    ins_l = {r: np.zeros(n) for r, n in ref_lengths.items()}
    dele = {r: np.zeros(n) for r, n in ref_lengths.items()}
    sub = {r: np.zeros(n) for r, n in ref_lengths.items()}
    counts_total = {r: 0 for r in ref_lengths}

    columns = [
        "#Reads",
        "Reference_Name",
        "all_insertion_positions",
        "all_insertion_left_positions",
        "all_deletion_positions",
        "all_substitution_positions",
    ]
    pf = pq.ParquetFile(path)
    for batch in pf.iter_batches(columns=columns, batch_size=batch_size):
        reads = batch.column("#Reads").to_pylist()
        refs = batch.column("Reference_Name").to_pylist()
        ip = batch.column("all_insertion_positions").to_pylist()
        ipl = batch.column("all_insertion_left_positions").to_pylist()
        dp = batch.column("all_deletion_positions").to_pylist()
        sp = batch.column("all_substitution_positions").to_pylist()
        for i in range(batch.num_rows):
            ref_name = refs[i]
            if ref_name not in ref_lengths:
                continue  # AMBIGUOUS_*/DISCARDED_* rows don't contribute
            count = int(reads[i])
            counts_total[ref_name] += count
            _acc(ins[ref_name], ip[i], count)
            _acc(ins_l[ref_name], ipl[i], count)
            _acc(dele[ref_name], dp[i], count)
            _acc(sub[ref_name], sp[i], count)

    indelsub = {
        r: ins[r] + dele[r] + sub[r] for r in ref_lengths
    }  # matches ~line 4200 (ins+del+sub; insertion_left NOT included)
    return CountVectors(
        all_insertion_count_vectors=ins,
        all_insertion_left_count_vectors=ins_l,
        all_deletion_count_vectors=dele,
        all_substitution_count_vectors=sub,
        all_indelsub_count_vectors=indelsub,
        counts_total=counts_total,
        ref_lengths=dict(ref_lengths),
    )


def _acc(arr, positions, count):
    """``arr[positions] += count`` with empty/None positions as a no-op."""
    if not positions:
        return
    arr[positions] += count


VariantStore.compute_count_vectors = _compute_count_vectors


def compute_count_vectors_from_collapsed(
    collapsed_path: str,
    ref_lengths: dict,
    output_directory: str,
    *,
    memory_budget_mb: int = DEFAULT_MEMORY_BUDGET_MB,
) -> CountVectors:
    """Convenience wrapper: create a :class:`VariantStore` and compute vectors."""
    store = VariantStore(output_directory, memory_budget_mb=memory_budget_mb)
    return store.compute_count_vectors(ref_lengths, collapsed_path=collapsed_path)


# ---------------------------------------------------------------------------
# Stage 4 — Stream-Out A: allele frequency TSV sink (PR 6)
# ---------------------------------------------------------------------------
#
# Streams the collapsed allele parquet (Stage 3 output, already sorted by
# ``(-#Reads, Aligned_Sequence, Reference_Sequence)``) to
# ``Alleles_frequency_table.txt`` — byte-identical to the pandas path:
#
#     df_alleles = pd.DataFrame(alleles_list)            # sorted
#     df_alleles['%Reads'] = df_alleles['#Reads'] / N_TOTAL * 100
#     df_alleles[['n_deleted','n_inserted','n_mutated']] = .astype(int)
#     (df_alleles if detailed else df_alleles.loc[:, cols]).to_csv(
#         path, sep='\t', header=True, index=None)
#
# Implementation: stream the parquet in bounded batches (pyarrow
# ``iter_batches``), reconstruct each batch as a pandas DataFrame with the
# *exact* cell types the pandas ``alleles_list`` carries (numpy arrays for
# position/size columns, tuple-lists for coordinate columns, numpy array for
# ``substitution_values``), and call ``DataFrame.to_csv`` per batch (header on
# the first batch only). Because pandas performs the float/array formatting —
# and the dtypes are identical batch-to-batch — the concatenated output is
# byte-identical to a whole-frame ``to_csv``, while peak memory is bounded by
# one batch (not the whole allele table).
#
# Handles all four pandas branches: detailed/non-detailed × dsODN/no-dsODN,
# replicating the ``crispresso2Cols`` projection, the ``%Reads`` derivation,
# the int cast, and the ``contains dsODN`` / ``contains dsODN fragment``
# boolean columns (``Aligned_Sequence.str.find(dsODN) > 0``).
#
# NOT wired into CRISPRessoCORE.main — exercised by unit tests only.

# Column order for the detailed allele table: get_allele_row's dict order
# (23 keys) with %Reads appended (pandas adds it after constructing the frame).
_DETAILED_ALLELE_COLS = [
    "#Reads",
    "Aligned_Sequence",
    "Reference_Sequence",
    "n_inserted",
    "n_deleted",
    "n_mutated",
    "Reference_Name",
    "Read_Status",
    "Aligned_Reference_Names",
    "Aligned_Reference_Scores",
    "ref_positions",
    "all_insertion_positions",
    "all_insertion_left_positions",
    "insertion_positions",
    "insertion_coordinates",
    "insertion_sizes",
    "all_deletion_positions",
    "deletion_positions",
    "deletion_coordinates",
    "deletion_sizes",
    "all_substitution_positions",
    "substitution_positions",
    "substitution_values",
]

# crispresso2Cols (non-detailed projection) — exact order from CRISPRessoCORE.
_CRISPRESSO2_COLS = [
    "Aligned_Sequence",
    "Reference_Sequence",
    "Reference_Name",
    "Read_Status",
    "n_deleted",
    "n_inserted",
    "n_mutated",
    "#Reads",
    "%Reads",
]

# Parquet columns needed to reconstruct an alleles_list-shaped row.
_TSV_PARQUET_COLUMNS = _DETAILED_ALLELE_COLS  # %Reads computed, not stored


def _parquet_row_to_allele_dict(row, n_total):
    """Convert a collapsed-parquet row (dict) to an alleles_list-shaped dict.

    Reconstructs the pandas cell types so ``DataFrame.to_csv`` formats them
    identically to the pandas ``alleles_list``:
      * position/size arrays → ``np.array(list, dtype=int)`` (``str`` →
        ``"[0 1 2]"``, matching ``str(np.array(...))`` in pandas).
      * coordinate arrays → list of ``(start, end)`` tuples (``str`` →
        ``"[(8, 10)]"``, matching pandas where ``insertion_coordinates`` is a
        list of tuples).
      * ``substitution_values`` → ``np.array(list, dtype='<U1')`` (``str`` →
        ``"['G' 'T']``, matching ``str(np.array(['G','T']))``).
    ``%Reads`` is computed as ``#Reads / n_total * 100`` (float64 — same as
    ``df['#Reads'] / N_TOTAL * 100``). ``n_inserted/n_deleted/n_mutated`` are
    left as ints; the batch DataFrame cast (``.astype(int)``) is applied later
    to match pandas exactly.
    """
    return {
        "#Reads": int(row["#Reads"]),
        "Aligned_Sequence": row["Aligned_Sequence"],
        "Reference_Sequence": row["Reference_Sequence"],
        "n_inserted": int(row["n_inserted"]),
        "n_deleted": int(row["n_deleted"]),
        "n_mutated": int(row["n_mutated"]),
        "Reference_Name": row["Reference_Name"],
        "Read_Status": row["Read_Status"],
        "Aligned_Reference_Names": row["Aligned_Reference_Names"],
        "Aligned_Reference_Scores": row["Aligned_Reference_Scores"],
        "ref_positions": _to_np_int_array(row.get("ref_positions")),
        "all_insertion_positions": _to_np_int_array(row.get("all_insertion_positions")),
        "all_insertion_left_positions": _to_np_int_array(row.get("all_insertion_left_positions")),
        "insertion_positions": _to_np_int_array(row.get("insertion_positions")),
        "insertion_coordinates": _structs_to_coords(row.get("insertion_coordinates")) or [],
        "insertion_sizes": _to_np_int_array(row.get("insertion_sizes")),
        "all_deletion_positions": _to_np_int_array(row.get("all_deletion_positions")),
        "deletion_positions": _to_np_int_array(row.get("deletion_positions")),
        "deletion_coordinates": _structs_to_coords(row.get("deletion_coordinates")) or [],
        "deletion_sizes": _to_np_int_array(row.get("deletion_sizes")),
        "all_substitution_positions": _to_np_int_array(row.get("all_substitution_positions")),
        "substitution_positions": _to_np_int_array(row.get("substitution_positions")),
        "substitution_values": _to_np_str_array(row.get("substitution_values")),
        "%Reads": (int(row["#Reads"]) / n_total * 100) if n_total > 0 else 0.0,
    }


def _to_np_int_array(val):
    """Convert a parquet list (or None) to a numpy int array (pandas cell parity)."""
    if val is None:
        return np.array([], dtype=int)
    return np.asarray(val, dtype=int)


def _to_np_str_array(val):
    """Convert a parquet list (or None) to a numpy <U1 string array (pandas parity)."""
    if val is None:
        return np.array([], dtype="<U1")
    return np.asarray([str(x) for x in val], dtype="<U1")


def _write_allele_frequency_table(
    self: "VariantStore",
    output_path: str,
    n_total: int,
    *,
    write_detailed_allele_table: bool = False,
    dsODN: str = "",
    collapsed_path: Optional[str] = None,
    batch_size: int = 50_000,
) -> str:
    """Stream the collapsed allele parquet to ``Alleles_frequency_table.txt``.

    Byte-identical to the pandas ``df_alleles.to_csv`` path (both detailed and
    non-detailed branches, with and without ``dsODN``), via per-batch pandas
    ``to_csv`` so peak memory is bounded by one batch.

    See the module section docstring above for the parity strategy.
    """
    import pandas as pd

    path = collapsed_path or os.path.join(self.output_directory, "collapsed.allele.parquet")
    rc = CRISPRessoShared.reverse_complement
    have_dsODN = dsODN != ""
    have_fragment = have_dsODN and len(dsODN) > 6
    sub_dsODN = dsODN[3:-3] if have_fragment else ""
    sub_dsODN_rc = rc(sub_dsODN) if have_fragment else ""
    dsODN_rc = rc(dsODN) if have_dsODN else ""

    # Build the column projection matching the pandas path exactly:
    #   detailed, no dsODN : get_allele_row cols (23) + %Reads
    #   detailed, dsODN    : above + dsODN bool cols (fw, rv, contains, fragment*)
    #   non-detailed       : crispresso2Cols (+ contains dsODN + fragment if dsODN)
    # Always projecting (rather than df.to_csv for detailed) keeps the header
    # consistent across batches and lets the empty-input case write a header.
    if write_detailed_allele_table:
        select_cols = list(_DETAILED_ALLELE_COLS) + ["%Reads"]
    else:
        select_cols = list(_CRISPRESSO2_COLS)
    if have_dsODN:
        # dsODN bool columns in the order pandas adds them to the DataFrame.
        dsODN_extra = ["contains dsODN fw", "contains dsODN rv", "contains dsODN"]
        if have_fragment:
            dsODN_extra += [
                "contains dsODN fragment fw",
                "contains dsODN fragment rv",
                "contains dsODN fragment",
            ]
        else:
            # len(dsODN) <= 6: pandas omits the fragment columns from the frame
            # (a latent KeyError on the non-detailed projection). We create an
            # all-False 'contains dsODN fragment' so the projection succeeds;
            # the common long-dsODN path is unaffected.
            dsODN_extra += ["contains dsODN fragment"]
        if write_detailed_allele_table:
            select_cols += dsODN_extra
        else:
            select_cols += ["contains dsODN", "contains dsODN fragment"]

    pf = pq.ParquetFile(path)
    with open(output_path, "w") as handle:
        first = True
        for batch in pf.iter_batches(columns=_TSV_PARQUET_COLUMNS, batch_size=batch_size):
            col_values = {c: batch.column(c).to_pylist() for c in _TSV_PARQUET_COLUMNS}
            n = batch.num_rows
            row_dicts = []
            for i in range(n):
                row = {c: col_values[c][i] for c in _TSV_PARQUET_COLUMNS}
                row_dicts.append(_parquet_row_to_allele_dict(row, n_total))
            if not row_dicts:
                continue
            df = pd.DataFrame(row_dicts)
            df["%Reads"] = df["#Reads"] / n_total * 100 if n_total > 0 else 0.0
            df[["n_deleted", "n_inserted", "n_mutated"]] = df[
                ["n_deleted", "n_inserted", "n_mutated"]
            ].astype(int)
            if have_dsODN:
                df["contains dsODN fw"] = df["Aligned_Sequence"].str.find(dsODN) > 0
                df["contains dsODN rv"] = df["Aligned_Sequence"].str.find(dsODN_rc) > 0
                df["contains dsODN"] = df["contains dsODN fw"] | df["contains dsODN rv"]
                if have_fragment:
                    df["contains dsODN fragment fw"] = df["Aligned_Sequence"].str.find(sub_dsODN) > 0
                    df["contains dsODN fragment rv"] = df["Aligned_Sequence"].str.find(sub_dsODN_rc) > 0
                    df["contains dsODN fragment"] = (
                        df["contains dsODN fragment fw"] | df["contains dsODN fragment rv"]
                    )
                else:
                    df["contains dsODN fragment"] = False
            df.loc[:, select_cols].to_csv(handle, sep="\t", header=first, index=False)
            first = False
        if first:
            # No rows (empty allele table): write the header only — matching
            # pandas df_alleles.to_csv on an empty frame projected to the cols.
            handle.write("\t".join(select_cols) + "\n")
    return output_path


VariantStore.write_allele_frequency_table = _write_allele_frequency_table


def write_allele_frequency_table(
    collapsed_path: str,
    output_path: str,
    n_total: int,
    output_directory: str,
    *,
    write_detailed_allele_table: bool = False,
    dsODN: str = "",
    memory_budget_mb: int = DEFAULT_MEMORY_BUDGET_MB,
) -> str:
    """Convenience wrapper: create a :class:`VariantStore` and sink the TSV."""
    store = VariantStore(output_directory, memory_budget_mb=memory_budget_mb)
    return store.write_allele_frequency_table(
        output_path,
        n_total,
        write_detailed_allele_table=write_detailed_allele_table,
        dsODN=dsODN,
        collapsed_path=collapsed_path,
    )


# ---------------------------------------------------------------------------
# Stage 5 — df_alleles lazy view + get_slice (PR 7, Phase 3 item 8)
# ---------------------------------------------------------------------------
#
# The collapsed allele parquet (Stage 3 output) is the persisted artifact
# that replaces the in-memory ``df_alleles`` pandas DataFrame. Consumers fall
# into two groups (design premise #4):
#
#   * Whole-table consumers (allele TSV sink [PR 6], VCF writer, report) —
#     already served by Stage 4 stream-outs and ``CollapsedAlleles.
#     allele_rows_dataframe()``.
#   * Slice consumers (~30 plot callsites in ``plots/data_prep.py`` and
#     ``CRISPRessoCORE.py``) — they call
#     ``df_alleles.loc[df_alleles['Reference_Name'] == ref]`` and then
#     ``iterrows()`` / ``groupby`` / ``str.find``. These receive a pandas
#     DataFrame materialized from a *filtered, projected* parquet slice via
#     ``get_slice()`` rather than a slice of the whole in-RAM frame.
#
# ``get_slice(ref_name, columns)`` scans the collapsed parquet with a row filter
# (``Reference_Name == ref_name``) and column projection at the pyarrow level
# (``iter_batches(columns=...)``), so peak memory for a plot is bounded by the
# slice size, not the whole allele table. Cell types are reconstructed to
# match the pandas ``alleles_list`` exactly (numpy int arrays for positions,
# tuple-lists for coordinates, numpy <U1 array for substitution_values) so plot
# code — which calls ``str()`` on cells, slices arrays, indexes them — is
# byte-for-byte identical to the pandas backend.
#
# NOT wired into CRISPRessoCORE.main — exercised by unit tests only. The
# callsite adaptation (item 9) builds on this API: when ``backend == parquet``
# the plot context's ``df_alleles`` will be a thin object whose
# ``.loc[Reference_Name == ref]`` materializes via ``get_slice``; until then,
# ``get_slice`` is a standalone parity-tested method.

# Columns whose parquet values (arrow list<int64> / list<struct>) must be
# reconstructed to the pandas ``alleles_list`` cell types for plot parity:
#   * position / size arrays  -> np.array(list, dtype=int)
#   * coordinate columns       -> list of (start, end) tuples
#   * substitution_values     -> np.array(list, dtype='<U1')
_INT_ARRAY_SLICE_COLS = frozenset([
    "ref_positions",
    "all_insertion_positions",
    "all_insertion_left_positions",
    "insertion_positions",
    "insertion_sizes",
    "all_deletion_positions",
    "deletion_positions",
    "deletion_sizes",
    "all_substitution_positions",
    "substitution_positions",
])
_COORD_SLICE_COLS = frozenset(["insertion_coordinates", "deletion_coordinates"])
_STR_ARRAY_SLICE_COLS = frozenset(["substitution_values"])


def _reconstruct_slice_cell(col_name, val):
    """Convert a parquet cell to the pandas ``alleles_list`` cell type.

    Parity with ``CRISPRessoCORE.main``'s ``alleles_list`` (and thus
    ``df_alleles``): numpy int arrays for position/size columns, tuple-lists
    for coordinate columns, numpy ``<U1`` array for ``substitution_values``.
    Scalar columns pass through unchanged.

    Tolerant of both storage shapes: parquet-scan cells (arrow lists of int /
    lists of struct dicts) and in-memory ``allele_rows`` cells (Python lists /
    lists of tuples / numpy arrays) so the same reconstruction applies on both
    the persisted-parquet and the write_parquet=False paths.
    """
    if col_name in _INT_ARRAY_SLICE_COLS:
        if isinstance(val, np.ndarray):
            return val if val.dtype.kind in "iu" else np.asarray(val, dtype=int)
        return _to_np_int_array(val)
    if col_name in _COORD_SLICE_COLS:
        if val is None:
            return []
        if isinstance(val, np.ndarray):
            return _structs_to_coords(val.tolist())
        if val and isinstance(val[0], tuple):
            return list(val)  # already a tuple-list (in-memory path)
        return _structs_to_coords(val)
    if col_name in _STR_ARRAY_SLICE_COLS:
        if isinstance(val, np.ndarray):
            return val if val.dtype.kind in "US" else np.asarray(val, dtype="<U1")
        return _to_np_str_array(val)
    return val


def _collapsed_get_slice(
    self: "CollapsedAlleles",
    ref_name=None,
    columns=None,
    *,
    include_pct_reads=True,
    collapsed_path=None,
    batch_size=50_000,
):
    """Return a pandas DataFrame slice of the collapsed allele table.

    Replaces ``df_alleles.loc[df_alleles['Reference_Name'] == ref_name]``
    (and the whole-frame ``df_alleles`` when ``ref_name`` / ``columns`` are
    ``None``) for the parquet backend, with row filtering and column
    projection pushed down to the parquet scan so peak memory is bounded by
    the slice, not the whole table.

    Parameters
    ----------
    ref_name : str or None
        Filter rows to ``Reference_Name == ref_name``. ``None`` returns all
        rows (the whole-table view). ``AMBIGUOUS_<ref>`` / ``DISCARDED_<ref>``
        rows are returned like any other — the caller filters if needed.
    columns : list[str] or None
        Project to these columns (order preserved). ``None`` returns the full
        detailed column set (``_DETAILED_ALLELE_COLS``) plus ``%Reads``. Pass
        ``_CRISPRESSO2_COLS`` for the non-detailed projection.
    include_pct_reads : bool
        If True and ``"%Reads"`` is in the output columns (or ``columns`` is
        None), add ``%Reads = #Reads / n_total * 100`` — matching the pandas
        ``df_alleles['%Reads']`` derivation.
    collapsed_path : str or None
        Path to the collapsed allele parquet. Defaults to
        :attr:`CollapsedAlleles.parquet_path`. When absent (or the file does
        not exist), falls back to materializing from the in-memory
        :attr:`allele_rows` — so ``get_slice`` works whether or not the
        parquet artifact was persisted.
    batch_size : int
        Row batch size for ``pyarrow.iter_batches`` (bounds peak memory).

    Returns
    -------
    pandas.DataFrame
        Cell types match ``df_alleles`` exactly (numpy int arrays for
        positions, tuple-lists for coordinates, numpy ``<U1`` array for
        ``substitution_values``); ``n_deleted`` / ``n_inserted`` / ``n_mutated``
        cast to int. An empty slice returns a frame with the requested columns
        and zero rows (no crash), matching pandas ``df_alleles.loc[...]`` on
        a no-match filter.

    Notes
    -----
    * Row filtering is currently applied per-batch in Python (pyarrow's
      ``iter_batches`` does not accept a predicate). True predicate pushdown
      (via ``pyarrow.dataset`` filters) is a follow-up; for typical per-ref
      slices the Python filter is cheap because the projection already
      trimmed the row width. The batched scan keeps peak memory bounded by
      ``batch_size`` rows × projected width regardless.
    * The in-memory fallback path does *not* reconstruct cell types —
      ``allele_rows`` already carry the native payload types (numpy arrays /
      tuple-lists) from :func:`_get_allele_row`.
    """
    import pandas as pd

    # Resolve the output column order.
    if columns is None:
        select = list(_DETAILED_ALLELE_COLS)
        if include_pct_reads:
            select.append("%Reads")
    else:
        select = list(columns)

    # Columns we must read from parquet to reconstruct the requested output.
    # ``%Reads`` is computed, not stored; it needs ``#Reads`` + ``n_total``.
    read_cols = [c for c in select if c != "%Reads"]
    if "%Reads" in select and "#Reads" not in read_cols:
        read_cols.append("#Reads")
    # When filtering by ref_name we must read Reference_Name to apply the row
    # filter, even if the caller did not ask for it in the output columns.
    if ref_name is not None and "Reference_Name" not in read_cols:
        read_cols.append("Reference_Name")
    # Preserve insertion order without duplicates.
    seen = set()
    read_cols = [c for c in read_cols if not (c in seen or seen.add(c))]

    path = collapsed_path or self.parquet_path
    row_dicts: list = []

    if path is not None and os.path.exists(path):
        import pyarrow.parquet as pq

        pf = pq.ParquetFile(path)
        for batch in pf.iter_batches(columns=read_cols, batch_size=batch_size):
            n = batch.num_rows
            if n == 0:
                continue
            col_values = {c: batch.column(c).to_pylist() for c in read_cols}
            ref_col = col_values.get("Reference_Name")
            for i in range(n):
                if ref_name is not None:
                    rn = ref_col[i] if ref_col is not None else None
                    if rn != ref_name:
                        continue
                row = {c: _reconstruct_slice_cell(c, col_values[c][i]) for c in read_cols}
                if include_pct_reads and "%Reads" in select:
                    reads = int(row["#Reads"])
                    row["%Reads"] = (reads / self.n_total * 100) if self.n_total > 0 else 0.0
                row_dicts.append(row)
    else:
        # In-memory fallback (e.g. write_parquet=False): allele_rows already
        # carry native cell types from _get_allele_row, but the position arrays
        # round-tripped through the parquet shard as Python lists — reconstruct
        # them to numpy so the slice matches df_alleles on both paths.
        for r in self.allele_rows:
            if ref_name is not None and r.get("Reference_Name") != ref_name:
                continue
            row = {c: _reconstruct_slice_cell(c, r.get(c)) for c in read_cols}
            if include_pct_reads and "%Reads" in select:
                reads = int(r.get("#Reads", 0))
                row["%Reads"] = (reads / self.n_total * 100) if self.n_total > 0 else 0.0
            row_dicts.append(row)

    if not row_dicts:
        # Empty slice: return a frame with the requested columns, zero rows.
        return pd.DataFrame({c: pd.Series(dtype=object) for c in select})

    df = pd.DataFrame(row_dicts)
    # Cast n_* to int (parity with df_alleles ~line 4310).
    for c in ("n_deleted", "n_inserted", "n_mutated"):
        if c in df.columns:
            df[c] = df[c].astype(int)
    return df.loc[:, select]


CollapsedAlleles.get_slice = _collapsed_get_slice


def get_slice_from_collapsed(
    collapsed: CollapsedAlleles,
    ref_name=None,
    columns=None,
    *,
    include_pct_reads=True,
    collapsed_path=None,
):
    """Convenience wrapper: :meth:`CollapsedAlleles.get_slice` in one call.

    Provided so parity tests and future wiring can materialize a slice without
    referencing the bound method directly.
    """
    return collapsed.get_slice(
        ref_name=ref_name,
        columns=columns,
        include_pct_reads=include_pct_reads,
        collapsed_path=collapsed_path,
    )
