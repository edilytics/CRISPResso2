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

Stages 2 (worker parquet writers) and 3 (streaming collapse) live in later PRs;
this module currently exposes only the Stage 1 ``VariantStore.count_reads`` API
and the ``ReadCounts`` result handle. Nothing here is wired into
``CRISPRessoCORE.main`` yet.
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
    "DEFAULT_MEMORY_BUDGET_MB",
    "AlignedShardWriter",
    "ReadCounts",
    "VariantStore",
    "count_reads_from_fastq",
    "iter_aligned_shard",
    "payload_to_row",
    "row_to_payload",
    "variant_parquet_generator_process",
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
# get_allele_row or the count-vector aggregation):
#   - ref_aln_details      (only for fastq_output annotation; edge #12 follow-up)
#   - all_deletion_coordinates (only for deletions_outside_window, already computed)
#   - all_substitution_values  (not in get_allele_row; substitution_values is)


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
