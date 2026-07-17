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
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Iterator, Optional

from CRISPResso2 import CRISPRessoShared

_logger = logging.getLogger("CRISPResso2")
info = _logger.info

# Default in-memory budget for the eager dedup dict (Stage 1 spectrum threshold).
# The eager path holds the full dedup dict in RAM — inherently O(unique_reads),
# so for the flat-memory success criterion (SC #1) the budget must be small
# enough that only genuinely-small inputs stay eager. 128 MB keeps ~50k short
# reads (or ~10k × 2 kb long reads) eager — the "fast end of the spectrum"
# where the eager path has no measurable slowdown vs. the current pandas path
# (SC #4) — while spilling the multi-million-read long-read regime to the
# external-sort path (peak RSS bounded by the sort buffer, not by read count).
# lowered from 512 MB after the PR 7 end-to-end RSS benchmark
# (scripts/bench_pipeline_memory.py) showed 1 M × 200 bp high-diversity reads
# (~280 MB dict) stayed eager and scaled linearly. Exposed via
# ``VariantStore(memory_budget_mb=...)`` so tests and future auto-select can tune.
DEFAULT_MEMORY_BUDGET_MB = 128

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
        for batch in pf.iter_batches(columns=columns, batch_size=_adaptive_batch_size(pf)):
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
    "AlleleAggregates",
    "CountVectors",
    "DEFAULT_MEMORY_BUDGET_MB",
    "AlignedShardWriter",
    "CollapsedAlleles",
    "ReadCounts",
    "VariantStore",
    "aggregate_alleles_from_collapsed",
    "compute_aln_stats_and_homology_from_shards",
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
    # all_substitution_values: per-position substituted base across the WHOLE
    # amplicon (length = ref_len; '-' where not substituted). NOT consumed by
    # get_allele_row (so absent from the detailed allele TSV) but required by
    # the per-base substitution count vectors (all_substitution_base_vectors)
    # built in the CRISPRessoCORE.main aggregation loop. Added here + to
    # COLLAPSED_SCHEMA so PR 7's streaming aggregation can build those vectors
    # from the collapsed parquet.
    pa.field("all_substitution_values", pa.list_(pa.string())),
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
# get_allele_row or the PR 7 streaming aggregation):
#   - ref_aln_details          (only for the HDR/prime-editing re-alignment
#     block and fastq_output annotation; edge #12 follow-up. The parquet
#     backend raises NotImplementedError on those paths until wired.)
#   - all_deletion_coordinates (only for deletions_outside_window, already
#     computed as a scalar in the payload)
# all_substitution_values IS stored (in _PAYLOAD_STRUCT + COLLAPSED_SCHEMA)
# as of PR 7 — it is needed by the per-base substitution count vectors
# (all_substitution_base_vectors) built in the CRISPRessoCORE.main aggregation
# loop. It is NOT in get_allele_row / _DETAILED_ALLELE_COLS, so it does not
# appear in the detailed allele TSV or the get_slice default projection; it is
# consumed only by the streaming aggregation (aggregate_alleles).


def _to_int_list(val):
    """Convert a numpy array, Python list, or scalar to an int64 numpy array.

    Returns a ``numpy.ndarray`` (NOT a Python list) so pyarrow's
    ``Table.from_pylist`` builds the list<int64> column directly from the
    contiguous buffer instead of materializing a Python list of boxed ints
    first (a ~3x transient reduction per position-array column — see
    ``design_docs/PARQUET_MEMORY_PROFILE.md`` item 3). Callers that consume
    the value only via pyarrow are unaffected; the parquet round-trip reader
    (``_struct_to_payload``) does ``list(...)`` on the way back out, so
    read-back values remain Python lists.
    """
    if val is None:
        return None
    import numpy as _np
    if hasattr(val, "tolist"):  # numpy array
        return _np.asarray(val, dtype=_np.int64)
    return _np.asarray([int(x) for x in val], dtype=_np.int64)


def _to_str_list(val):
    """Convert a numpy array or list to a numpy array of Python strings.

    Returns a ``numpy.ndarray`` (object dtype) for the same reason as
    :func:`_to_int_list` — pyarrow builds the list<string> column from the
    array without a boxed-Python-str list intermediate.
    """
    if val is None:
        return None
    import numpy as _np
    if hasattr(val, "tolist"):
        val = val.tolist()
    return _np.asarray([str(x) for x in val], dtype=object)


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
    The per-ref sub-payload emitted by ``find_indels_substitutions`` is a Cython
    ``ResultsSlotsDict`` (not a plain ``dict``); coerce it via ``vars()`` so the
    ``.get()`` calls below work uniformly on both the Cython and plain-dict
    shapes (the latter is what the unit-test canned payloads use).
    """
    if not isinstance(payload, dict):
        payload = vars(payload)
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
        "all_substitution_values": _to_str_list(payload.get("all_substitution_values")),
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
        "all_substitution_values": np.array(struct_row["all_substitution_values"]) if struct_row["all_substitution_values"] else np.array([]),
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


def iter_aligned_shard(path, batch_size=None, *, columns=None):
    """Stream ``(read_key, count, payload)`` from a parquet shard.

    Used by Stage 3 (collapse) and the homology/stat stream-outs to iterate
    aligned payloads without materializing the whole shard. Reads in batches
    so peak memory is bounded by one batch.

    ``batch_size`` defaults to an *adaptive* value: the parquet metadata's
    uncompressed row-group size is used to estimate bytes/row, and the batch
    is sized to target ~50 MB of arrow data per batch. This keeps peak RSS
    flat for long-read amplicons (a fixed 50k-row batch would hold ~8 GB of
    2 kb-read payloads). Pass an explicit ``batch_size`` to override.

    ``columns`` selects which top-level columns to read off disk (parquet is
    columnar, so unnamed columns are never read). Pass a minimal set for
    passes that only need scalars — e.g. the canonical-key projection needs
    only ``["read_key", "count", "best_match_score"]`` and must NOT pay to
    read/decompress the amplicon-length ``payloads`` struct. ``None`` (default)
    reads all columns (full payload) for backward compatibility. NOTE: pyarrow
    cannot project nested struct leaves (``payloads.aln_seq`` is silently
    ignored), so ``payloads`` is all-or-nothing; passes that need *any* nested
    field must read the whole ``payloads`` column.
    """
    pf = pq.ParquetFile(str(path))
    if batch_size is None:
        batch_size = _adaptive_batch_size(pf)
    for batch in pf.iter_batches(batch_size=batch_size, columns=columns):
        # Iterate row-by-row via column-cell access (``.as_py()``) rather than
        # ``batch.to_pylist()``: the latter materializes the WHOLE batch as
        # Python objects at once, which blows up ~10-50x for string-array
        # columns (all_substitution_values is amplicon-length single-char
        # strings/row). Row-by-row keeps only one row's Python objects in
        # flight — peak RSS is the (compact) arrow batch + one payload.
        names = batch.schema.names
        cols = [batch.column(n) for n in names]
        for i in range(batch.num_rows):
            row = {names[j]: cols[j][i].as_py() for j in range(len(names))}
            yield row["read_key"], row["count"], row_to_payload(row)


def _adaptive_batch_size(pf, target_bytes=50_000_000):
    """Pick a row batch size that targets ~``target_bytes`` of arrow data.

    Estimates bytes/row from the first row group's *uncompressed* column sizes
    (the in-memory arrow footprint, not the smaller compressed on-disk size).
    Falls back to 50k rows when metadata is unavailable. Used by every
    streaming parquet scan in this module (shard iteration, count vectors,
    allele aggregation, TSV sink, get_slice) so peak RSS stays flat for
    long-read amplicons regardless of which stage is running.
    """
    try:
        meta = pf.metadata
        if meta is None or meta.num_row_groups == 0:
            return 50_000
        rg = meta.row_group(0)
        uncompressed = sum(
            rg.column(j).total_uncompressed_size for j in range(rg.num_columns)
        )
        bytes_per_row = max(1, uncompressed // max(1, rg.num_rows))
        return max(500, min(50_000, target_bytes // bytes_per_row))
    except Exception:
        return 50_000


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
    # all_substitution_values (see _PAYLOAD_STRUCT): per-position substituted
    # base across the whole amplicon; needed by the streaming per-base
    # substitution count vectors. NOT part of the detailed allele TSV
    # (_DETAILED_ALLELE_COLS) — the TSV sink and get_slice default projection
    pa.field("all_substitution_values", pa.list_(pa.string())),
])

# COLLAPSED_SCHEMA + a synthetic ``row_idx`` column. Used by the partitioned
# gather (``design_docs/LARGE_ALLELE_ORDER_PERF.md``) to write the unsorted
# full-row parquet in Step A and seek rows back by index in Step D.
_COLLAPSED_SCHEMA_WITH_ROW_IDX = COLLAPSED_SCHEMA.append(
    pa.field("row_idx", pa.int64())
)


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
        "all_substitution_values": _to_str_list(d.get("all_substitution_values")),
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
    aln_stats, aln_homology, not_aln_homology : Optional
        Alignment statistics + homology dicts (the outputs of
        ``compute_aln_stats_and_homology_from_shards``), fused into the
        single-read Pass 2 scan so homology free-rides on the payloads read
        that fanout already needs (no separate shard pass). Populated by the
        single-read branch; ``None`` for the paired/empty paths (which don't
        compute homology here).

    """

    allele_rows: list
    n_total: int
    class_counts: dict
    counts_total: dict
    counts_modified: dict
    counts_unmodified: dict
    counts_discarded: dict
    parquet_path: Optional[str] = None
    aln_stats: Optional[dict] = None
    aln_homology: Optional[dict] = None
    not_aln_homology: Optional[dict] = None

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
    keep_allele_rows: bool = True,
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

    # Dispatch: single-read input routes through the streaming external-sort
    # canonical-key merge (flat peak RSS regardless of unique-read count —
    # see ``_collapse_streaming_single_read``). Paired input stays on the
    # in-memory path (bounded by the unique allele count after 3a re-key).
    if not is_paired:
        return self._collapse_streaming_single_read(
            paths,
            discard_indel_reads=discard_indel_reads,
            write_detailed=write_detailed,
            write_parquet=write_parquet,
            collapsed_path=collapsed_path,
            expand_ambiguous_alignments=expand_ambiguous_alignments,
            keep_allele_rows=keep_allele_rows,
        )

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
    # sink. The persisted parquet is always full. When ``keep_allele_rows`` is
    # False (the wired bench/streaming-TSV path — df_alleles skipped) the
    # in-memory copy is pure dead weight: ``get_slice`` reads the parquet
    # artifact directly, so return ``[]`` and drop ``full_rows`` before the
    # aggregation stage allocates (design_docs/PARQUET_MEMORY_PROFILE.md
    # item 2).
    if not keep_allele_rows:
        allele_rows = []
        del full_rows
    elif write_detailed:
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
            allele_row = _get_allele_row(
                "AMBIGUOUS_" + aln_ref_names[0],
                variant_count, aln_ref_names_str, aln_ref_scores_str,
                variant_payload, write_detailed=True,
            )
            allele_row["all_substitution_values"] = variant_payload.get("all_substitution_values")
            allele_rows.append(allele_row)
            continue

        for ref_name in aln_ref_names:
            variant_payload = payload["variant_" + ref_name]
            if discard_indel_reads and (
                variant_payload["deletion_n"] > 0 or variant_payload["insertion_n"] > 0
            ):
                counts_discarded[ref_name] = counts_discarded.get(ref_name, 0) + variant_count
                allele_row = _get_allele_row(
                    "DISCARDED_" + aln_ref_names[0],
                    variant_count, aln_ref_names_str, aln_ref_scores_str,
                    variant_payload, write_detailed=True,
                )
                allele_row["all_substitution_values"] = variant_payload.get("all_substitution_values")
                allele_rows.append(allele_row)
                continue
            counts_total[ref_name] = counts_total.get(ref_name, 0) + variant_count
            if variant_payload["classification"] == "MODIFIED":
                counts_modified[ref_name] = counts_modified.get(ref_name, 0) + variant_count
            else:
                counts_unmodified[ref_name] = counts_unmodified.get(ref_name, 0) + variant_count
            allele_row = _get_allele_row(
                ref_name, variant_count, aln_ref_names_str, aln_ref_scores_str,
                variant_payload, write_detailed=True,
            )
            allele_row["all_substitution_values"] = variant_payload.get("all_substitution_values")
            allele_rows.append(allele_row)

    return (
        allele_rows, n_total, class_counts,
        counts_total, counts_modified, counts_unmodified, counts_discarded,
    )


# ---------------------------------------------------------------------------
# Stage 3 (single-read): streaming external-sort canonical-key merge.
#
# Replaces the in-memory ``store`` dict (``_collapse_rekey_and_rcmerge``) for
# the single-read path. The in-memory store is bounded by the unique *allele*
# count for paired input (3a re-key collapses first), but for single-read input
# there is no re-key, so the store holds one entry per unique *read* — the
# dominant memory term that OOMs long-read / high-diversity data (see
# ``design_docs/STREAMING_SINGLE_READ_COLLAPSE.md``: ~3.5 GB at 1M×200 bp,
# ~40× the Stage 1 dedup dict). This streaming path keeps peak RSS bounded by
# one payload + the (tiny) group index, regardless of unique-read count.
#
# Algorithm (mirrors Stage 1's S1b-proven external merge-sort pattern):
#   1. Pass 1 — stream aligned shards, compute ``canonical_key =
#      min(read_key, rc(read_key))`` per row, write a text projection
#      ``canonical_key \t seq_no \t count`` (seq_no = monotonic scan order,
#      used for the tie-break). Bounded by one shard batch in flight.
#   2. External merge-sort the text projection by ``canonical_key`` (system
#      ``sort -s -S <buf>``; spike ``bench_polars_sort_spill.py`` confirmed
#      polars ``sort`` does NOT spill — ratio 33 — while system ``sort`` stays
#      flat — ratio 1.24 from 100k→1M).
#   3. Streaming first-per-group collapse — stream the sorted text; adjacent
#      rows share a canonical_key. For each group: sum ``count``, keep the
#      representative = the row with the smallest ``seq_no`` (first in scan
#      order — the exact pandas insertion-order tie-break, since stable sort
#      preserves ascending seq_no within a group). Palindrome doubling
#      (``canonical_key == rc(canonical_key)``) is replicated verbatim for
#      byte-for-byte parity with the pandas path's latent double-count.
#   4. Sort the groups by ``rep_seq_no`` → scan order (pandas insertion order).
#   5. Pass 2 — re-stream the aligned shards (same deterministic order); for
#      each representative (seq_no match), reconstruct the payload and run
#      ``_collapse_fanout`` on a single-entry store. This produces the same
#      allele rows + aggregations as the in-memory path, one variant at a time.
#   6. Final ordering — the allele rows must be sorted by
#      ``(-#Reads, Aligned_Sequence, Reference_Sequence)`` (matching
#      ``df_alleles.sort_values`` ~line 4303) for TSV/get_slice parity. If the
#      allele set fits the memory budget (small inputs: tests, basic test),
#      sort in memory and populate ``allele_rows``. Otherwise use the
#      partitioned external gather from ``LARGE_ALLELE_ORDER_PERF.md``: write
#      full rows to an unsorted parquet + a small sort-key projection,
#      external-sort ONLY the key projection, then gather full rows into
#      in-memory buckets (one bucket at a time) and emit each in order. The
#      full row data never passes through a global text sort — peak RSS is
#      bounded by one bucket + accumulators, and the 1M x 2 kb run completes
#      (the JSON/TSV-carry global sort that preceded it timed out).
#
# Parity holds by construction: ``seq_no`` encodes scan order, so the
# representative (min seq_no) is the first-occurrence-in-scan-order key — the
# exact pandas tie-break — and sorting groups by rep_seq_no reproduces the
# pandas store insertion order, so the final stable sort's tie-break is
# identical. Paired input stays on the in-memory path (genuinely bounded by
# allele count); only single-read routes here.


def _ensure_nofile(needed: int) -> None:
    """Raise the soft ``RLIMIT_NOFILE`` so the partitioned gather can keep
    ``needed`` bucket spill files open simultaneously.

    No-ops if the soft limit already suffices or the platform lacks
    ``resource``. Never raises — a failure here is non-fatal (the gather will
    surface a clear ``OSError: Too many open files`` if it truly cannot
    proceed, and ``n_buckets`` is sized to stay within the hard limit on the
    platforms CRISPResso2 supports: macOS hard ≥ 1024, Linux ≥ 4096).
    """
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        target = int(min(hard, needed))
        if target > soft:
            resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))
    except (ImportError, OSError, ValueError):
        pass


def _run_external_sort(input_file: str, output_file: str, workdir: str,
                       key_args: list, *, sort_buffer: str = _SORT_BUFFER) -> None:
    """External merge-sort a text file, spilling to ``workdir``.

    ``key_args`` is the list of sort key specifiers (e.g. ``["-k1,1"]`` or
    ``["-k1,1", "-k2,2", "-k3,3"]``). ``-s`` stable so equal-key rows preserve
    file order (parity: the first occurrence in scan order wins per group).
    ``-S`` caps the in-memory buffer (forces disk spill); ``-T`` sets the spill
    dir; ``LC_ALL=C`` forces byte-wise comparison (fast, locale-independent,
    and matches Python's code-point sort for ACGT/-/+ strings).
    """
    tmpdir = os.path.join(workdir, "sort_tmp")
    os.makedirs(tmpdir, exist_ok=True)
    sort_cmd = ["sort", "-S", sort_buffer, "-T", tmpdir, "-s", "-t", "\t"]
    sort_cmd += key_args
    sort_cmd += ["-o", output_file, input_file]
    env = dict(os.environ, LC_ALL="C")
    try:
        subprocess.run(sort_cmd, check=True, env=env)
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        raise CRISPRessoShared.OutputFolderIncompleteException(
            "External sort failed during single-read collapse: %s" % (e,)
        ) from e


def _join_ints(val) -> str:
    """Encode a list/numpy array of ints as a space-joined string for TSV carry.

    Empty/None → empty string. Negatives are fine (``-8 -8 8 9``) — the space
    separator is unambiguous for ints.
    """
    if val is None:
        return ""
    if hasattr(val, "tolist"):
        val = val.tolist()
    return " ".join(str(int(x)) for x in val)


def _split_ints(s: str) -> list:
    """Inverse of :func:`_join_ints`. Empty string → empty list."""
    if not s:
        return []
    return [int(x) for x in s.split()]


def _join_coords(val) -> str:
    """Encode a list of (start, end) coordinate tuples as ``start,end`` ";"-joined.

    Comma separates start/end within a coord; semicolon separates coords. Handles
    negative positions unambiguously (``-8,10;12,14``). Empty/None → empty.
    Accepts either (start, end) tuples/lists or arrow struct dicts
    (``{"start": .., "end": ..}``) so it can encode rows read back from a
    collapsed parquet (Step D of the partitioned gather) without a separate
    struct→tuple conversion.
    """
    if val is None:
        return ""
    parts = []
    for c in val:
        if isinstance(c, dict):
            parts.append(f"{int(c['start'])},{int(c['end'])}")
        else:
            parts.append(f"{int(c[0])},{int(c[1])}")
    return ";".join(parts)


def _split_coords(s: str) -> list:
    """Inverse of :func:`_join_coords`. Empty string → empty list of tuples."""
    if not s:
        return []
    return [tuple(int(x) for x in c.split(",")) for c in s.split(";")]


def _join_strs(val) -> str:
    """Encode a list/numpy array of strings as a space-joined string for TSV carry.

    Used for ``substitution_values`` / ``all_substitution_values`` — single-char
    ACGT/'-' values, so the space separator is unambiguous. Empty/None → empty.
    """
    if val is None:
        return ""
    if hasattr(val, "tolist"):
        val = val.tolist()
    return " ".join(str(x) for x in val)


def _split_strs(s: str) -> list:
    """Inverse of :func:`_join_strs`. Empty string → empty list."""
    if not s:
        return []
    return s.split()


# TSV field layout for the external-sort row carry (one line per allele row).
# Fields 1-4 are the sort keys (ascending sort = the final allele order);
# fields 5+ are the row data. No JSON — compact tab/space/semicolon encoding.
#   1: {10**15 - #Reads:015d}   (ascending sort = descending #Reads)
#   2: Aligned_Sequence         (also row data — no duplication)
#   3: Reference_Sequence        (also row data)
#   4: {row_idx:015d}            (tiebreaker = emission order = pandas insertion order)
#   5: #Reads, 6: n_inserted, 7: n_deleted, 8: n_mutated,
#   9: Reference_Name, 10: Read_Status,
#  11: Aligned_Reference_Names, 12: Aligned_Reference_Scores,
#  13-16: ref_positions, all_insertion_positions, all_insertion_left_positions, insertion_positions,
#  17: insertion_coordinates, 18: insertion_sizes,
#  19-20: all_deletion_positions, deletion_positions,
#  21: deletion_coordinates, 22: deletion_sizes,
#  23-24: all_substitution_positions, substitution_positions,
#  25-26: substitution_values, all_substitution_values
_TSV_ROW_FIELDS = (
    "#Reads", "n_inserted", "n_deleted", "n_mutated",
    "Reference_Name", "Read_Status",
    "Aligned_Reference_Names", "Aligned_Reference_Scores",
    "ref_positions", "all_insertion_positions", "all_insertion_left_positions",
    "insertion_positions", "insertion_coordinates", "insertion_sizes",
    "all_deletion_positions", "deletion_positions",
    "deletion_coordinates", "deletion_sizes",
    "all_substitution_positions", "substitution_positions",
    "substitution_values", "all_substitution_values",
)


def _allele_row_to_tsv(row: dict, row_idx: int) -> str:
    """Serialize a full allele-row dict to one TSV line for the external sort.

    The first four fields are the sort keys (neg_reads, Aligned_Sequence,
    Reference_Sequence, row_idx); the rest is the row data. Handles numpy arrays
    (via ``.tolist()``) and coordinate tuples. The string columns (aln_seq,
    ref_seq, Reference_Name, Read_Status, Aligned_Reference_Names/Scores) never
    contain tabs/spaces/semicolons (ACGT/'-'/'&'/digits only), so the encoding
    is unambiguous. ~2-5x smaller and ~2-3x faster to round-trip than JSON
    (no recursive tokenizer; arrays are ``" ".join(map(str,...))``).
    """
    reads = int(row["#Reads"])
    parts = [
        f"{10**15 - reads:015d}",
        row["Aligned_Sequence"],
        row["Reference_Sequence"],
        f"{row_idx:015d}",
        str(reads),
        str(int(row["n_inserted"])),
        str(int(row["n_deleted"])),
        str(int(row["n_mutated"])),
        row["Reference_Name"],
        row["Read_Status"],
        row["Aligned_Reference_Names"],
        row["Aligned_Reference_Scores"],
        _join_ints(row["ref_positions"]),
        _join_ints(row["all_insertion_positions"]),
        _join_ints(row["all_insertion_left_positions"]),
        _join_ints(row["insertion_positions"]),
        _join_coords(row["insertion_coordinates"]),
        _join_ints(row["insertion_sizes"]),
        _join_ints(row["all_deletion_positions"]),
        _join_ints(row["deletion_positions"]),
        _join_coords(row["deletion_coordinates"]),
        _join_ints(row["deletion_sizes"]),
        _join_ints(row["all_substitution_positions"]),
        _join_ints(row["substitution_positions"]),
        _join_strs(row["substitution_values"]),
        _join_strs(row["all_substitution_values"]),
    ]
    return "\t".join(parts)


def _tsv_line_to_allele_row(line: str) -> dict:
    """Inverse of :func:`_allele_row_to_tsv` → a row dict for ``_allele_dict_to_parquet_row``.

    Skips the 4 sort-key fields; reconstructs plain Python lists/tuples/strings
    for the row-data fields. ``_allele_dict_to_parquet_row`` then converts to
    the arrow-native COLLAPSED_SCHEMA row (int64 lists, coord structs, etc.).
    """
    f = line.split("\t")
    return {
        "Aligned_Sequence": f[1],
        "Reference_Sequence": f[2],
        "#Reads": int(f[4]),
        "n_inserted": int(f[5]),
        "n_deleted": int(f[6]),
        "n_mutated": int(f[7]),
        "Reference_Name": f[8],
        "Read_Status": f[9],
        "Aligned_Reference_Names": f[10],
        "Aligned_Reference_Scores": f[11],
        "ref_positions": _split_ints(f[12]),
        "all_insertion_positions": _split_ints(f[13]),
        "all_insertion_left_positions": _split_ints(f[14]),
        "insertion_positions": _split_ints(f[15]),
        "insertion_coordinates": _split_coords(f[16]),
        "insertion_sizes": _split_ints(f[17]),
        "all_deletion_positions": _split_ints(f[18]),
        "deletion_positions": _split_ints(f[19]),
        "deletion_coordinates": _split_coords(f[20]),
        "deletion_sizes": _split_ints(f[21]),
        "all_substitution_positions": _split_ints(f[22]),
        "substitution_positions": _split_ints(f[23]),
        "substitution_values": _split_strs(f[24]),
        "all_substitution_values": _split_strs(f[25]),
    }


def _streaming_collapse_groups(sorted_keys_file: str):
    """Stream a canonical-key-sorted text file → one (key, total_count, rep_seq_no) per group.

    Input lines: ``canonical_key \t seq_no \t count``. Stable sort guarantees
    ascending seq_no within a group, so the first row is the representative
    (min seq_no = first in scan order — the pandas tie-break). Palindrome
    doubling (``key == rc(key)``) is applied for parity with the pandas path.
    Returns a list in sorted-canonical-key order.
    """
    rc = CRISPRessoShared.reverse_complement
    groups = []
    prev_key = None
    total = 0
    rep_seq = None
    with open(sorted_keys_file, "r", buffering=1 << 20) as sf:
        for line in sf:
            line = line.rstrip("\n")
            if not line:
                continue
            key, _, rest = line.partition("\t")
            seq_s, _, count_s = rest.partition("\t")
            seq = int(seq_s)
            cnt = int(count_s)
            if key != prev_key:
                if prev_key is not None:
                    if prev_key == rc(prev_key):
                        total *= 2
                    groups.append((prev_key, total, rep_seq))
                prev_key = key
                total = cnt
                rep_seq = seq
            else:
                total += cnt
        if prev_key is not None:
            if prev_key == rc(prev_key):
                total *= 2
            groups.append((prev_key, total, rep_seq))
    return groups


def _accumulate_agg(dst_class_counts, dst_counts_total, dst_counts_modified,
                    dst_counts_unmodified, dst_counts_discarded,
                    class_counts, counts_total, counts_modified,
                    counts_unmodified, counts_discarded):
    """Merge one ``_collapse_fanout`` call's aggregation deltas into the running totals."""
    for k, v in class_counts.items():
        dst_class_counts[k] = dst_class_counts.get(k, 0) + v
    for k, v in counts_total.items():
        dst_counts_total[k] = dst_counts_total.get(k, 0) + v
    for k, v in counts_modified.items():
        dst_counts_modified[k] = dst_counts_modified.get(k, 0) + v
    for k, v in counts_unmodified.items():
        dst_counts_unmodified[k] = dst_counts_unmodified.get(k, 0) + v
    for k, v in counts_discarded.items():
        dst_counts_discarded[k] = dst_counts_discarded.get(k, 0) + v


def _write_collapsed_allele_parquet_from_tsv(sorted_tsv: str, path: str,
                                              batch_size: int = 10_000) -> int:
    """Stream a sort-key-prefixed TSV file into the collapsed allele parquet.

    Input lines (already sorted by fields 1-4):
    ``{neg_reads}\t{aligned}\t{ref}\t{row_idx}\t{row data...}``. Parses each line
    via :func:`_tsv_line_to_allele_row`, converts to a COLLAPSED_SCHEMA row,
    and writes in batches. Returns the row count.

    The batch size is adaptive: the first line's byte length estimates
    bytes/row, and the batch targets ~32 MB of Python row dicts so peak RSS
    stays flat for long-read amplicons. TSV parse (``split`` + int joins) is
    ~2-3x cheaper than JSON ``loads`` — no recursive tokenizer.
    """
    writer = pq.ParquetWriter(path, COLLAPSED_SCHEMA)
    n = 0
    buf = []
    try:
        with open(sorted_tsv, "r", buffering=1 << 20) as fh:
            # Estimate bytes/row from the first line to size batches by bytes.
            first = fh.readline()
            if not first:
                return 0
            bytes_per_row = max(1, len(first))
            flush_every = max(100, min(batch_size, 32_000_000 // bytes_per_row))
            def _flush():
                nonlocal n
                if not buf:
                    return
                table = pa.Table.from_pylist(buf, schema=COLLAPSED_SCHEMA)
                writer.write_table(table)
                n += len(buf)
                buf.clear()
            line = first.rstrip("\n")
            if line:
                buf.append(_allele_dict_to_parquet_row(_tsv_line_to_allele_row(line)))
                if len(buf) >= flush_every:
                    _flush()
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                buf.append(_allele_dict_to_parquet_row(_tsv_line_to_allele_row(line)))
                if len(buf) >= flush_every:
                    _flush()
            _flush()
    finally:
        writer.close()
    return n


def _partitioned_gather_write_unsorted(
    self: "VariantStore",
    paths: list,
    rep_seq_sorted,
    rep_map: dict,
    *,
    discard_indel_reads: bool,
    est_row_size: int,
    unsorted_parquet: str,
    gather_keys_file: str,
    homology_state: Optional["_AlnStatsHomologyState"] = None,
) -> tuple:
    """Step A of the partitioned gather (``LARGE_ALLELE_ORDER_PERF.md``).

    One streaming pass over the aligned shards: for each representative
    (``seq_no`` match) run ``_collapse_fanout`` and (a) write the FULL row to
    ``unsorted_parquet`` (native arrow, ``_COLLAPSED_SCHEMA_WITH_ROW_IDX``) and
    (b) append the small sort-key projection to ``gather_keys_file``
    (``neg_reads \t Aligned_Sequence \t Reference_Sequence \t row_idx``).
    Returns ``(n_allele_rows, n_total, class_counts, counts_total,
    counts_modified, counts_unmodified, counts_discarded)``.
    """
    n_total = 0
    class_counts: dict = {}
    counts_total: dict = {}
    counts_modified: dict = {}
    counts_unmodified: dict = {}
    counts_discarded: dict = {}

    flush_every = max(100, min(8192, 32_000_000 // max(1, est_row_size)))
    writer = pq.ParquetWriter(unsorted_parquet, _COLLAPSED_SCHEMA_WITH_ROW_IDX)
    try:
        buf: list = []
        row_idx = 0
        seq_no = 0
        rep_iter = iter(rep_seq_sorted)
        next_rep = next(rep_iter, None)
        with open(gather_keys_file, "w", buffering=1 << 20) as kf:
            for path in paths:
                for read_key, count, payload in iter_aligned_shard(path):
                    if homology_state is not None:
                        homology_state.accumulate(read_key, int(count), payload)
                    if next_rep is not None and seq_no == next_rep:
                        ck, tc = rep_map[seq_no]
                        store = {ck: {"count": tc, "payload": payload}}
                        rows, nt, cc, ct, cm, cu, cd = self._collapse_fanout(
                            store, discard_indel_reads=discard_indel_reads)
                        n_total += nt
                        _accumulate_agg(class_counts, counts_total, counts_modified,
                                        counts_unmodified, counts_discarded,
                                        cc, ct, cm, cu, cd)
                        for r in rows:
                            d = _allele_dict_to_parquet_row(r)
                            d["row_idx"] = row_idx
                            buf.append(d)
                            kf.write(
                                f"{10**15 - int(r['#Reads']):015d}\t"
                                f"{r['Aligned_Sequence']}\t"
                                f"{r['Reference_Sequence']}\t"
                                f"{row_idx:015d}\n")
                            row_idx += 1
                            if len(buf) >= flush_every:
                                writer.write_table(pa.Table.from_pylist(
                                    buf, schema=_COLLAPSED_SCHEMA_WITH_ROW_IDX))
                                buf.clear()
                        next_rep = next(rep_iter, None)
                    seq_no += 1
                    # No early break — homology accumulation must visit every row.
        if buf:
            writer.write_table(pa.Table.from_pylist(
                buf, schema=_COLLAPSED_SCHEMA_WITH_ROW_IDX))
            buf.clear()
    finally:
        writer.close()

    return (row_idx, n_total, class_counts, counts_total, counts_modified,
            counts_unmodified, counts_discarded)


def _partitioned_gather_emit(
    self: "VariantStore",
    *,
    n_allele_rows: int,
    est_row_size: int,
    unsorted_parquet: str,
    gather_keys_sorted: str,
    bucket_dir: str,
    parquet_path: str,
) -> None:
    """Steps C–E of the partitioned gather.

    *C* — read ``gather_keys_sorted`` (output order) and build
    ``out_pos[row_idx] = output index``.
    *D* — stream ``unsorted_parquet`` once; scatter each row to a per-bucket TSV
    spill file (``out_pos`` + the row's TSV encoding). Peak RSS = K small write
    buffers + one arrow batch.
    *E* — per bucket (in output order): read, sort in memory by ``out_pos``, and
    write to ``parquet_path`` (batched arrow writes). Peak RSS = one bucket.
    """
    import numpy as _np

    # -- Step C: out_pos[row_idx] = output index (from the sorted key projection).
    out_pos = _np.empty(n_allele_rows, dtype=_np.int64)
    pos = 0
    with open(gather_keys_sorted, "r", buffering=1 << 20) as sf:
        for line in sf:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            out_pos[int(f[3])] = pos
            pos += 1

    # Bucket size: one bucket materialised as Python dicts in Step E should
    # stay within the memory budget (dict overhead ≈ 3× the arrow row).
    chunk_size = max(1, int(self.memory_budget_bytes // max(1, est_row_size * 3)))
    n_buckets = (n_allele_rows + chunk_size - 1) // chunk_size

    os.makedirs(bucket_dir, exist_ok=True)
    bucket_paths = [
        os.path.join(bucket_dir, f"bucket_{b:06d}.tsv") for b in range(n_buckets)
    ]
    # Step D keeps all n_buckets spill files open simultaneously; raise the FD
    # soft limit so this does not hit the (often 256) default.
    _ensure_nofile(n_buckets + 128)

    # -- Step D: partitioned gather (one streaming pass over the unsorted parquet).
    bucket_fhs = [open(p, "w", buffering=1 << 20) for p in bucket_paths]
    try:
        pf = pq.ParquetFile(unsorted_parquet)
        read_batch = _adaptive_batch_size(pf)
        for batch in pf.iter_batches(batch_size=read_batch):
            names = batch.schema.names
            # Column access + per-row .as_py() (mirrors ``iter_aligned_shard``):
            # materialises ONE row's Python objects at a time so peak RSS is the
            # (compact) arrow batch + one payload, not the whole batch as pylist.
            cols = {n: batch.column(n) for n in names}
            row_idx_col = cols["row_idx"]
            data_names = [n for n in names if n != "row_idx"]
            for i in range(batch.num_rows):
                ri = row_idx_col[i].as_py()
                op = int(out_pos[ri])
                b = op // chunk_size
                if b >= n_buckets:  # guard against rounding edge
                    b = n_buckets - 1
                row = {n: cols[n][i].as_py() for n in data_names}
                # ``_allele_row_to_tsv`` encodes coords via ``_join_coords``,
                # which accepts the struct dicts read back from parquet.
                bucket_fhs[b].write(f"{op:015d}\t{_allele_row_to_tsv(row, ri)}\n")
    finally:
        for fh in bucket_fhs:
            try:
                fh.close()
            except OSError:
                pass

    # -- Step E: emit (per bucket, in output order).
    final_writer = pq.ParquetWriter(parquet_path, COLLAPSED_SCHEMA)
    try:
        flush_every = max(100, min(8192, 32_000_000 // max(1, est_row_size)))
        for b in range(n_buckets):
            entries = []
            with open(bucket_paths[b], "r", buffering=1 << 20) as bf:
                for line in bf:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    op_s, _, rest = line.partition("\t")
                    entries.append((int(op_s), rest))
            entries.sort(key=lambda e: e[0])
            buf = []
            for _, rest in entries:
                buf.append(_allele_dict_to_parquet_row(
                    _tsv_line_to_allele_row(rest)))
                if len(buf) >= flush_every:
                    final_writer.write_table(
                        pa.Table.from_pylist(buf, schema=COLLAPSED_SCHEMA))
                    buf.clear()
            if buf:
                final_writer.write_table(
                    pa.Table.from_pylist(buf, schema=COLLAPSED_SCHEMA))
                buf.clear()
            del entries
    finally:
        final_writer.close()


def _collapse_streaming_single_read(
    self: "VariantStore",
    paths: list,
    *,
    discard_indel_reads: bool = False,
    write_detailed: bool = False,
    write_parquet: bool = True,
    collapsed_path: Optional[str] = None,
    expand_ambiguous_alignments: bool = False,
    keep_allele_rows: bool = True,
) -> CollapsedAlleles:
    """Stage 3 for single-read input: streaming external-sort canonical-key merge.

    See the section docstring above for the algorithm and parity argument.
    Produces the same ``CollapsedAlleles`` as the in-memory path, but with peak
    RSS bounded by one payload + the group index (O(unique_reads) scalars) +
    accumulators, instead of one full payload per unique read.

    Homology (``aln_stats`` / ``aln_homology`` / ``not_aln_homology``) is fused
    into the Pass 2 fanout scan — that scan already reads ``payloads`` for
    fanout, so the per-row homology accumulation free-rides on it (no separate
    shard pass). Populated on the returned ``CollapsedAlleles``.
    """
    rc = CRISPRessoShared.reverse_complement
    workdir = tempfile.mkdtemp(prefix="crispresso_collapse_", dir=self.output_directory)
    keys_file = os.path.join(workdir, "keys.txt")
    sorted_keys_file = os.path.join(workdir, "keys.sorted")

    parquet_path = None
    if write_parquet:
        parquet_path = (
            collapsed_path
            or os.path.join(self.output_directory, "collapsed.allele.parquet")
        )

    try:
        # -- Pass 1: stream shards → text projection (canonical_key, seq_no, count).
        # Column-projected to the 3 top-level scalars this pass needs — it must
        # NOT read/decompress the amplicon-length ``payloads`` struct (the
        # dominant shard cost). ``best_match_score`` is a top-level column, so
        # the unaligned-skip needs no payload access. ``est_amplicon_len``
        # (which only sizes the in-memory-vs-external threshold) is taken from
        # ``len(read_key)`` — for single-read input read_key IS the read
        # sequence, ≈ amplicon length — avoiding a dive into ``payloads``.
        est_amplicon_len = 0
        seq_no = 0
        any_aligned = False
        _pass1_columns = ["read_key", "count", "best_match_score"]
        with open(keys_file, "w", buffering=1 << 20) as kf:
            for path in paths:
                for read_key, count, payload in iter_aligned_shard(
                        path, columns=_pass1_columns):
                    if payload.get("best_match_score", 0) <= 0:
                        seq_no += 1  # consume seq_no for seek-back parity
                        continue
                    any_aligned = True
                    rc_read = rc(read_key)
                    canonical_key = read_key if read_key <= rc_read else rc_read
                    if est_amplicon_len == 0:
                        est_amplicon_len = len(read_key)
                    kf.write(canonical_key)
                    kf.write("\t")
                    kf.write(str(seq_no))
                    kf.write("\t")
                    kf.write(str(count))
                    kf.write("\n")
                    seq_no += 1

        # Empty input / all-unaligned: emit an empty parquet (parity with the
        # in-memory path's empty case).
        if not any_aligned:
            if parquet_path is not None:
                _write_collapsed_allele_parquet([], parquet_path)
            # No aligned reads → Pass 2 won't run, so accumulate homology for
            # the (all-unaligned) rows here via the shared state. Rare edge;
            # the normal path fuses this into Pass 2.
            _st = _AlnStatsHomologyState(
                expand_ambiguous_alignments=expand_ambiguous_alignments)
            _unaln_cols = ["read_key", "count", "best_match_score", "aln_scores"]
            for path in paths:
                for read_key, count, payload in iter_aligned_shard(
                        path, columns=_unaln_cols):
                    _st.accumulate(read_key, int(count), payload)
            _as, _ah, _nah = _st.finalize()
            return CollapsedAlleles(
                allele_rows=[], n_total=0, class_counts={},
                counts_total={}, counts_modified={}, counts_unmodified={},
                counts_discarded={}, parquet_path=parquet_path,
                aln_stats=_as, aln_homology=_ah, not_aln_homology=_nah,
            )

        # -- External sort by canonical_key (stable → ascending seq_no within a group).
        _run_external_sort(keys_file, sorted_keys_file, workdir, ["-k1,1"],
                           sort_buffer=self.sort_buffer)

        # -- Streaming first-per-group collapse.
        groups = _streaming_collapse_groups(sorted_keys_file)
        num_groups = len(groups)

        # Sort groups by rep_seq_no → scan order (pandas insertion order), so
        # the final stable sort's tie-break is identical to the in-memory path.
        groups.sort(key=lambda g: g[2])
        rep_map = {g[2]: (g[0], g[1]) for g in groups}
        rep_seq_sorted = sorted(rep_map.keys())

        # Decide in-memory vs external final sort. est_row_size ≈ 10 int64
        # arrays of amplicon length + the two amplicon-length strings.
        est_row_size = max(1, est_amplicon_len) * 10 * 8 + 4096
        use_in_memory = num_groups * est_row_size < self.memory_budget_bytes

        n_total = 0
        class_counts: dict = {}
        counts_total: dict = {}
        counts_modified: dict = {}
        counts_unmodified: dict = {}
        counts_discarded: dict = {}
        # Homology/stats accumulators — fused into the Pass 2 fanout scan
        # (both branches) so homology free-rides on the payloads read that
        # fanout already needs. ``compute_aln_stats_and_homology_from_shards``
        # is NOT called separately on the wired path (fix #3).
        _homology_st = _AlnStatsHomologyState(
            expand_ambiguous_alignments=expand_ambiguous_alignments)

        if use_in_memory:
            # Small input (tests, basic test): hold allele rows in memory,
            # sort, write parquet, populate allele_rows — same shape as the
            # in-memory path. One variant's payload in flight at a time.
            full_rows: list = []
            seq_no = 0
            rep_iter = iter(rep_seq_sorted)
            next_rep = next(rep_iter, None)
            for path in paths:
                for read_key, count, payload in iter_aligned_shard(path):
                    _homology_st.accumulate(read_key, int(count), payload)
                    if next_rep is not None and seq_no == next_rep:
                        ck, tc = rep_map[seq_no]
                        store = {ck: {"count": tc, "payload": payload}}
                        rows, nt, cc, ct, cm, cu, cd = self._collapse_fanout(
                            store, discard_indel_reads=discard_indel_reads)
                        full_rows.extend(rows)
                        n_total += nt
                        _accumulate_agg(class_counts, counts_total, counts_modified,
                                        counts_unmodified, counts_discarded,
                                        cc, ct, cm, cu, cd)
                        next_rep = next(rep_iter, None)
                    seq_no += 1
                    # NOTE: no early break once reps are consumed — homology
                    # accumulation (fused into this scan) must visit EVERY row,
                    # including unaligned reads trailing the last representative.

            full_rows.sort(key=lambda r: (
                -int(r["#Reads"]), r["Aligned_Sequence"], r["Reference_Sequence"]))
            if parquet_path is not None:
                _write_collapsed_allele_parquet(full_rows, parquet_path)
            # When ``keep_allele_rows`` is False (the wired bench/streaming-
            # TSV path) ``allele_rows`` is never consumed — ``get_slice`` reads
            # the persisted parquet directly. Drop ``full_rows`` now so the
            # ~num_groups x amplicon-length transient is freed before Stage 4
            # allocates (design_docs/PARQUET_MEMORY_PROFILE.md item 2).
            if not keep_allele_rows:
                allele_rows = []
                del full_rows
            elif write_detailed:
                allele_rows = full_rows
            else:
                allele_rows = [
                    {k: r[k] for k in _NON_DETAILED_ALLELE_KEYS if k in r}
                    for r in full_rows
                ]
        else:
            # Large input (bench): partitioned external gather — key-sort +
            # bucketed gather. The full row data never passes through a global
            # text sort (the bottleneck that timed out 1M x 2000 bp under the
            # old JSON/TSV-carry global sort). See
            # ``design_docs/LARGE_ALLELE_ORDER_PERF.md``.
            #
            #   A. fanout → write FULL rows to an unsorted parquet (native
            #      arrow, ``_COLLAPSED_SCHEMA_WITH_ROW_IDX``) + a small
            #      ``gather.keys.txt`` (sort-key prefix + row_idx only).
            #   B. external-sort the small key projection (stable, by output
            #      keys — byte-identical order to the old global TSV sort).
            #   C. build ``out_pos[row_idx]`` (output index) from keys.sorted.
            #   D. stream the unsorted parquet once; scatter each row to a
            #      per-bucket TSV spill file (``out_pos`` + row TSV).
            #   E. per bucket (in output order): read, sort by ``out_pos`` in
            #      memory, write to the final parquet. Peak RSS = one bucket.
            #
            # ``allele_rows`` is left empty — the pipeline consumes the parquet
            # via ``get_slice`` (``write_parquet`` is required on this branch).
            unsorted_parquet = os.path.join(workdir, "collapsed.unsorted.parquet")
            gather_keys_file = os.path.join(workdir, "gather.keys.txt")
            gather_keys_sorted = os.path.join(workdir, "gather.keys.sorted")
            bucket_dir = os.path.join(workdir, "buckets")

            (n_allele_rows, n_total, class_counts, counts_total,
             counts_modified, counts_unmodified, counts_discarded) = \
                _partitioned_gather_write_unsorted(
                    self, paths, rep_seq_sorted, rep_map,
                    discard_indel_reads=discard_indel_reads,
                    est_row_size=est_row_size,
                    unsorted_parquet=unsorted_parquet,
                    gather_keys_file=gather_keys_file,
                    homology_state=_homology_st)

            if parquet_path is not None and n_allele_rows > 0:
                # Step B — sort the small key projection (sort keys + row_idx
                # only; ~4 kB/row at 2 kb amplicon vs ~26 kB/row for full TSV).
                _run_external_sort(gather_keys_file, gather_keys_sorted, workdir,
                                   ["-k1,1", "-k2,2", "-k3,3", "-k4,4"],
                                   sort_buffer=self.sort_buffer)
                # Steps C–E — bucket assignment + gather + emit.
                _partitioned_gather_emit(
                    self,
                    n_allele_rows=n_allele_rows,
                    est_row_size=est_row_size,
                    unsorted_parquet=unsorted_parquet,
                    gather_keys_sorted=gather_keys_sorted,
                    bucket_dir=bucket_dir,
                    parquet_path=parquet_path)
            elif parquet_path is not None:
                # Fan-out produced zero rows (all representatives discarded).
                _write_collapsed_allele_parquet([], parquet_path)
            allele_rows = []

        _as, _ah, _nah = _homology_st.finalize()
        return CollapsedAlleles(
            allele_rows=allele_rows,
            n_total=n_total,
            class_counts=class_counts,
            counts_total=counts_total,
            counts_modified=counts_modified,
            counts_unmodified=counts_unmodified,
            counts_discarded=counts_discarded,
            parquet_path=parquet_path,
            aln_stats=_as, aln_homology=_ah, not_aln_homology=_nah,
        )
    finally:
        # All temp artifacts live under ``workdir`` (keys files, the unsorted
        # parquet, and the per-bucket spill dir). ``shutil.rmtree`` removes
        # the bucket subdirectory and every spill file in one go (edges #15/#21).
        shutil.rmtree(workdir, ignore_errors=True)


def _write_collapsed_allele_parquet(allele_rows: list, path: str) -> None:
    """Write the full-column collapsed allele table to parquet.

    ``allele_rows`` are full (detailed) dicts (see ``_collapse_fanout``); the
    persisted artifact always carries the complete column set so PR 6's
    count-vector aggregation and PR 7's df_alleles view can project either the
    detailed or ``crispresso2Cols`` shape.

    Writes in *batched* fashion (convert + ``write_table`` per batch, then
    clear) rather than materializing the whole ``rows`` list + a single
    ``pa.Table.from_pylist(rows)`` of the entire table at once. At 2000 bp
    that full-table arrow materialization is the dominant transient
    (~one full copy of every position array alive simultaneously);
    batching bounds peak RSS by one batch regardless of allele count
    (``design_docs/PARQUET_MEMORY_PROFILE.md`` item 1).
    """
    schema = COLLAPSED_SCHEMA
    writer = pq.ParquetWriter(path, schema)
    try:
        if not allele_rows:
            return
        # Batch sized on a byte budget so long-read amplicons (large
        # per-row position arrays) still bound peak RSS. ~16 MB of row dicts
        # per batch is a handful of thousand rows at 2 kb.
        _BATCH_BYTES = 16_000_000
        batch: list = []
        batch_bytes = 0
        for d in allele_rows:
            row = _allele_dict_to_parquet_row(d)
            batch.append(row)
            # estimate this row's dict bytes from the position-array columns
            # (numpy arrays expose ``.nbytes``; Python lists/scalars fall
            # back to a coarse 8 B/element estimate, scalars → 8 B).
            for v in row.values():
                if v is None:
                    continue
                nb = getattr(v, "nbytes", None)
                if nb is not None:
                    batch_bytes += nb
                elif isinstance(v, (list, tuple)):
                    batch_bytes += len(v) * 8
                else:
                    batch_bytes += 8
            if batch_bytes >= _BATCH_BYTES:
                writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                batch.clear()
                batch_bytes = 0
        if batch:
            writer.write_table(pa.Table.from_pylist(batch, schema=schema))
            batch.clear()
    finally:
        writer.close()


# Bind Stage 3 onto VariantStore.
VariantStore.collapse = _collapse
VariantStore._collapse_rekey_and_rcmerge = _collapse_rekey_and_rcmerge
VariantStore._collapse_fanout = _collapse_fanout
VariantStore._collapse_streaming_single_read = _collapse_streaming_single_read


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
    keep_allele_rows: bool = True,
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
        keep_allele_rows=keep_allele_rows,
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
    batch_size: Optional[int] = None,
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
    if batch_size is None:
        batch_size = _adaptive_batch_size(pf)
    for batch in pf.iter_batches(columns=columns, batch_size=batch_size):
        # Row-by-row cell access (not batch.column(...).to_pylist()) so the
        # amplicon-length position-array columns don't materialize as whole-
        # batch Python lists at once (~50x peak-RSS blowup for string/list
        # columns per design_docs/PARQUET_MEMORY_PROFILE.md item 6).
        reads_col = batch.column("#Reads")
        refs_col = batch.column("Reference_Name")
        ip_col = batch.column("all_insertion_positions")
        ipl_col = batch.column("all_insertion_left_positions")
        dp_col = batch.column("all_deletion_positions")
        sp_col = batch.column("all_substitution_positions")
        for i in range(batch.num_rows):
            ref_name = refs_col[i].as_py()
            if ref_name not in ref_lengths:
                continue  # AMBIGUOUS_*/DISCARDED_* rows don't contribute
            count = int(reads_col[i].as_py())
            counts_total[ref_name] += count
            _acc(ins[ref_name], ip_col[i].as_py(), count)
            _acc(ins_l[ref_name], ipl_col[i].as_py(), count)
            _acc(dele[ref_name], dp_col[i].as_py(), count)
            _acc(sub[ref_name], sp_col[i].as_py(), count)

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
# Stage 4b — Streaming per-allele aggregation (PR 7, Phase 3 item 9)
# ---------------------------------------------------------------------------
#
# ``CRISPRessoCORE.main``'s allele loop (~line 3973-4180) does more than the
# four ``all_*`` position vectors + ``class_counts`` covered by Stage 4 / the
# ``CollapsedAlleles`` pre-fanout outputs. In a single pass over
# ``variantCache`` it also builds, per reference:
#
#   * the windowed count vectors (``insertion_count_vectors``, etc.),
#   * the per-base vectors (``all_substitution_base_vectors``,
#     ``all_base_count_vectors``),
#   * the length vectors (``insertion_length_vectors``, ``deletion_length_vectors``),
#   * the per-allele-size Counters (``inserted_n_dicts``, etc.),
#   * the composite modification counts (``counts_only_*``, ``counts_*_and_*``),
#   * the exon / splicing / frameshift counts + ``hists_inframe`` / ``hists_frameshift``.
#
# Under the parquet backend these are built by a *streaming pass* over the
# collapsed allele parquet (Stage 3 output): scan row-group by row-group, and
# for each (allele, reference) row run the per-ref aggregation body verbatim.
# The pre-fanout parts (``class_counts``, ``N_TOTAL``,
# ``counts_total/modified/unmodified/discarded``) are already computed by
# ``CollapsedAlleles`` (collapse) and are NOT recomputed here; AMBIGUOUS and
# DISCARDED rows are skipped (the pandas loop ``continue``\ s before the per-ref
# body for them).
#
# Parity note: every accumulator is integer arithmetic (exact in float64 up to
# 2^53) or a Counter (commutative), so the result is independent of row order —
# the collapsed parquet's ``(-#Reads, Aligned, Ref)`` sort order yields identical
# totals to the pandas loop's ``variantCache`` insertion order. This is the key
# fact that lets the streaming pass match the pandas path byte-for-byte.
#
# Scope: replicates the FULL per-ref body including the exon/splicing/frameshift
# logic (gated by ``refs[ref_name]['contains_coding_seq']``), so it is correct
# for coding-seq amplicons too. The HDR / prime-editing re-alignment block
# (~line 4230) is NOT replicated — it needs ``ref_aln_details`` (raw alignment
# strings) which are deliberately not stored (see _PAYLOAD_STRUCT comment); the
# parquet backend raises ``NotImplementedError`` on that path until wired.


@dataclass
class AlleleAggregates:
    """Per-reference aggregation outputs from the streaming allele pass.

    Mirrors the dict/numpy/Counter structures ``CRISPRessoCORE.main`` builds in
    the allele loop, so the downstream plot/quantification code consumes them
    unchanged under the parquet backend. All fields are keyed by ``ref_name``
    (or ``ref_name + "_" + nuc`` for the per-base vectors), matching the pandas
    path's initialization (~line 3850).
    """

    all_insertion_count_vectors: dict
    all_insertion_left_count_vectors: dict
    all_deletion_count_vectors: dict
    all_substitution_count_vectors: dict
    all_indelsub_count_vectors: dict
    all_substitution_base_vectors: dict
    all_base_count_vectors: dict
    insertion_count_vectors: dict
    deletion_count_vectors: dict
    substitution_count_vectors: dict
    indelsub_count_vectors: dict
    insertion_count_vectors_noncoding: dict
    deletion_count_vectors_noncoding: dict
    substitution_count_vectors_noncoding: dict
    insertion_length_vectors: dict
    deletion_length_vectors: dict
    inserted_n_dicts: dict
    deleted_n_dicts: dict
    substituted_n_dicts: dict
    effective_len_dicts: dict
    hists_inframe: dict
    hists_frameshift: dict
    counts_insertion: dict
    counts_deletion: dict
    counts_substitution: dict
    counts_only_insertion: dict
    counts_only_deletion: dict
    counts_only_substitution: dict
    counts_insertion_and_deletion: dict
    counts_insertion_and_substitution: dict
    counts_deletion_and_substitution: dict
    counts_insertion_and_deletion_and_substitution: dict
    counts_modified_frameshift: dict
    counts_modified_non_frameshift: dict
    counts_non_modified_non_frameshift: dict
    counts_splicing_sites_modified: dict
    substitution_base_vectors: dict  # windowed subset, populated in finalize()


@dataclass
class _AggState:
    """Mutable accumulator state for the streaming per-allele pass.

    Held in a single dataclass so :func:`_aggregate_one_row` can mutate it
    without a long argument list, and so :func:`_finalize_aggregates` can derive
    the post-loop vectors (``all_indelsub``, ``indelsub``, ``substitution_base_vectors``)
    in one place — matching the pandas post-loop block (~line 4184-4190).
    """

    refs: dict
    ref_names: list
    ignore_insertions: bool
    ignore_deletions: bool
    ignore_substitutions: bool
    all_insertion_count_vectors: dict = field(default_factory=dict)
    all_insertion_left_count_vectors: dict = field(default_factory=dict)
    all_deletion_count_vectors: dict = field(default_factory=dict)
    all_substitution_count_vectors: dict = field(default_factory=dict)
    insertion_count_vectors: dict = field(default_factory=dict)
    deletion_count_vectors: dict = field(default_factory=dict)
    substitution_count_vectors: dict = field(default_factory=dict)
    insertion_count_vectors_noncoding: dict = field(default_factory=dict)
    deletion_count_vectors_noncoding: dict = field(default_factory=dict)
    substitution_count_vectors_noncoding: dict = field(default_factory=dict)
    insertion_length_vectors: dict = field(default_factory=dict)
    deletion_length_vectors: dict = field(default_factory=dict)
    all_substitution_base_vectors: dict = field(default_factory=dict)
    all_base_count_vectors: dict = field(default_factory=dict)
    inserted_n_dicts: dict = field(default_factory=dict)
    deleted_n_dicts: dict = field(default_factory=dict)
    substituted_n_dicts: dict = field(default_factory=dict)
    effective_len_dicts: dict = field(default_factory=dict)
    hists_inframe: dict = field(default_factory=dict)
    hists_frameshift: dict = field(default_factory=dict)
    counts_insertion: dict = field(default_factory=dict)
    counts_deletion: dict = field(default_factory=dict)
    counts_substitution: dict = field(default_factory=dict)
    counts_only_insertion: dict = field(default_factory=dict)
    counts_only_deletion: dict = field(default_factory=dict)
    counts_only_substitution: dict = field(default_factory=dict)
    counts_insertion_and_deletion: dict = field(default_factory=dict)
    counts_insertion_and_substitution: dict = field(default_factory=dict)
    counts_deletion_and_substitution: dict = field(default_factory=dict)
    counts_insertion_and_deletion_and_substitution: dict = field(default_factory=dict)
    counts_modified_frameshift: dict = field(default_factory=dict)
    counts_modified_non_frameshift: dict = field(default_factory=dict)
    counts_non_modified_non_frameshift: dict = field(default_factory=dict)
    counts_splicing_sites_modified: dict = field(default_factory=dict)


def _init_agg_state(refs: dict, ref_names: list, args) -> _AggState:
    """Initialize the per-ref accumulator structures.

    Verbatim replica of ``CRISPRessoCORE.main``'s init block (~line 3850-3855):
    one zeroed numpy vector per ref for each position vector, one Counter per
    ref for each size/effective-length/hist dict, and zero-initialized scalars
    for every count bucket.
    """
    from collections import Counter

    st = _AggState(
        refs=refs,
        ref_names=ref_names,
        ignore_insertions=getattr(args, "ignore_insertions", False),
        ignore_deletions=getattr(args, "ignore_deletions", False),
        ignore_substitutions=getattr(args, "ignore_substitutions", False),
    )
    for ref_name in ref_names:
        n = refs[ref_name]["sequence_length"]
        st.all_insertion_count_vectors[ref_name] = np.zeros(n)
        st.all_insertion_left_count_vectors[ref_name] = np.zeros(n)
        st.all_deletion_count_vectors[ref_name] = np.zeros(n)
        st.all_substitution_count_vectors[ref_name] = np.zeros(n)
        st.insertion_count_vectors[ref_name] = np.zeros(n)
        st.deletion_count_vectors[ref_name] = np.zeros(n)
        st.substitution_count_vectors[ref_name] = np.zeros(n)
        st.insertion_count_vectors_noncoding[ref_name] = np.zeros(n)
        st.deletion_count_vectors_noncoding[ref_name] = np.zeros(n)
        st.substitution_count_vectors_noncoding[ref_name] = np.zeros(n)
        st.insertion_length_vectors[ref_name] = np.zeros(n)
        st.deletion_length_vectors[ref_name] = np.zeros(n)
        for nuc in ("A", "C", "G", "T", "N"):
            st.all_substitution_base_vectors[ref_name + "_" + nuc] = np.zeros(n)
        for nuc in ("A", "C", "G", "T", "N", "-"):
            st.all_base_count_vectors[ref_name + "_" + nuc] = np.zeros(n)
        st.inserted_n_dicts[ref_name] = Counter()
        st.deleted_n_dicts[ref_name] = Counter()
        st.substituted_n_dicts[ref_name] = Counter()
        st.effective_len_dicts[ref_name] = Counter()
        st.hists_inframe[ref_name] = Counter()
        st.hists_inframe[ref_name][0] = 0
        st.hists_frameshift[ref_name] = Counter()
        st.hists_frameshift[ref_name][0] = 0
        for attr in (
            "counts_insertion", "counts_deletion", "counts_substitution",
            "counts_only_insertion", "counts_only_deletion", "counts_only_substitution",
            "counts_insertion_and_deletion", "counts_insertion_and_substitution",
            "counts_deletion_and_substitution", "counts_insertion_and_deletion_and_substitution",
            "counts_modified_frameshift", "counts_modified_non_frameshift",
            "counts_non_modified_non_frameshift", "counts_splicing_sites_modified",
        ):
            getattr(st, attr)[ref_name] = 0
    return st


def _agg_acc(arr, positions, count):
    """``arr[positions] += count`` with empty/None positions as a no-op."""
    if positions is None or len(positions) == 0:
        return
    arr[positions] += count


def _aggregate_one_row(st: _AggState, ref_name: str, variant_count: int, p: dict) -> None:
    """Run the per-ref aggregation body for one collapsed allele row.

    ``p`` is a payload-shaped dict (keys matching ``variant_payload`` in the
    pandas loop: ``aln_seq``, ``ref_positions``, ``all_insertion_positions``,
    ``insertion_n``, ``all_substitution_values``, ``insertion_coordinates``,
    etc.) built by :func:`_collapsed_row_to_payload`. This is a verbatim
    replica of ``CRISPRessoCORE.main``'s per-ref body (~line 4006-4180) minus
    the ``counts_total/modified/unmodified`` lines (handled by collapse) and
    the ``alleles_list.append`` (df_alleles comes from collapse / get_slice).
    """
    refs = st.refs
    this_effective_len = refs[ref_name]["sequence_length"]

    this_has_insertions = False
    _agg_acc(st.all_insertion_count_vectors[ref_name], p["all_insertion_positions"], variant_count)
    _agg_acc(st.all_insertion_left_count_vectors[ref_name], p["all_insertion_left_positions"], variant_count)

    if not st.ignore_insertions:
        st.inserted_n_dicts[ref_name][p["insertion_n"]] += variant_count
        _agg_acc(st.insertion_count_vectors[ref_name], p["insertion_positions"], variant_count)
        this_effective_len = this_effective_len + p["insertion_n"]
        if p["insertion_n"] > 0:
            st.counts_insertion[ref_name] += variant_count
            this_has_insertions = True

    this_has_deletions = False
    _agg_acc(st.all_deletion_count_vectors[ref_name], p["all_deletion_positions"], variant_count)
    if not st.ignore_deletions:
        st.deleted_n_dicts[ref_name][p["deletion_n"]] += variant_count
        _agg_acc(st.deletion_count_vectors[ref_name], p["deletion_positions"], variant_count)
        this_effective_len = this_effective_len - p["deletion_n"]
        if p["deletion_n"] > 0:
            st.counts_deletion[ref_name] += variant_count
            this_has_deletions = True

    st.effective_len_dicts[ref_name][this_effective_len] += variant_count

    this_has_substitutions = False
    _agg_acc(st.all_substitution_count_vectors[ref_name], p["all_substitution_positions"], variant_count)

    if not st.ignore_substitutions:
        st.substituted_n_dicts[ref_name][p["substitution_n"]] += variant_count
        _agg_acc(st.substitution_count_vectors[ref_name], p["substitution_positions"], variant_count)
        if p["substitution_n"] > 0:
            st.counts_substitution[ref_name] += variant_count
            this_has_substitutions = True

        nucs = ["A", "T", "C", "G", "N"]
        all_sub_values = p["all_substitution_values"]
        all_sub_positions = p["all_substitution_positions"]
        for nuc in nucs:
            isNuc = [n == nuc for n in all_sub_values]
            if np.sum(isNuc) > 0:
                locs = np.array(all_sub_positions)[isNuc]
                st.all_substitution_base_vectors[ref_name + "_" + nuc][locs] += variant_count

    if this_has_deletions:
        if this_has_insertions:
            if this_has_substitutions:
                st.counts_insertion_and_deletion_and_substitution[ref_name] += variant_count
            else:
                st.counts_insertion_and_deletion[ref_name] += variant_count
        elif this_has_substitutions:
            st.counts_deletion_and_substitution[ref_name] += variant_count
        else:
            st.counts_only_deletion[ref_name] += variant_count
    elif this_has_insertions:
        if this_has_substitutions:
            st.counts_insertion_and_substitution[ref_name] += variant_count
        else:
            st.counts_only_insertion[ref_name] += variant_count
    elif this_has_substitutions:
        st.counts_only_substitution[ref_name] += variant_count

    # set all_base_count_vectors
    aln_seq = p["aln_seq"]
    ref_pos = p["ref_positions"]
    for i in range(len(aln_seq)):
        if ref_pos[i] < 0:
            continue
        nuc = aln_seq[i]
        st.all_base_count_vectors[ref_name + "_" + nuc][ref_pos[i]] += variant_count

    exon_len_mods = refs[ref_name]["exon_len_mods"]
    tot_exon_len_mod = sum(exon_len_mods)
    if this_has_insertions or this_has_deletions or this_has_substitutions or tot_exon_len_mod != 0:
        exon_positions = refs[ref_name]["exon_positions"]
        splicing_positions = refs[ref_name]["splicing_positions"]
        insertion_coordinates = p["insertion_coordinates"]
        insertion_sizes = p["insertion_sizes"]
        insertion_positions = p["insertion_positions"]
        deletion_coordinates = p["deletion_coordinates"]
        deletion_sizes = p["deletion_sizes"]
        deletion_positions = p["deletion_positions"]
        substitution_positions = p["substitution_positions"]

        length_modified_positions_exons = []
        current_read_exons_modified = False
        current_read_spliced_modified = False

        for idx_ins, (ins_start, ins_end) in enumerate(insertion_coordinates):
            st.insertion_length_vectors[ref_name][ins_start] += (insertion_sizes[idx_ins] * variant_count)
            st.insertion_length_vectors[ref_name][ins_end] += (insertion_sizes[idx_ins] * variant_count)

            if refs[ref_name]["contains_coding_seq"]:
                if set(exon_positions).intersection((ins_start, ins_end)):
                    current_read_exons_modified = True
                    length_modified_positions_exons.append((insertion_sizes[idx_ins]))

        for idx_del, (del_start, del_end) in enumerate(deletion_coordinates):
            st.deletion_length_vectors[ref_name][list(range(del_start, del_end))] += (deletion_sizes[idx_del] * variant_count)

        if refs[ref_name]["contains_coding_seq"]:
            del_positions_to_append = sorted(set(exon_positions).intersection(set(deletion_positions)))
            if del_positions_to_append:
                current_read_exons_modified = True
                length_modified_positions_exons.append(-len(del_positions_to_append))

            if set(exon_positions).intersection(substitution_positions):
                current_read_exons_modified = True

            if set(splicing_positions).intersection(deletion_positions):
                current_read_spliced_modified = True
            if set(splicing_positions).intersection(insertion_positions):
                current_read_spliced_modified = True
            if set(splicing_positions).intersection(substitution_positions):
                current_read_spliced_modified = True

            if current_read_spliced_modified:
                st.counts_splicing_sites_modified[ref_name] += variant_count

            if tot_exon_len_mod != 0:
                effective_length = sum(length_modified_positions_exons) + tot_exon_len_mod
                if (effective_length % 3) == 0:
                    st.counts_modified_non_frameshift[ref_name] += variant_count
                    st.hists_inframe[ref_name][effective_length] += variant_count
                else:
                    st.counts_modified_frameshift[ref_name] += variant_count
                    st.hists_frameshift[ref_name][effective_length] += variant_count
            elif current_read_exons_modified:
                if not length_modified_positions_exons:
                    st.counts_modified_non_frameshift[ref_name] += variant_count
                    st.hists_inframe[ref_name][0] += variant_count
                else:
                    effective_length = sum(length_modified_positions_exons)
                    if (effective_length % 3) == 0:
                        st.counts_modified_non_frameshift[ref_name] += variant_count
                        st.hists_inframe[ref_name][effective_length] += variant_count
                    else:
                        st.counts_modified_frameshift[ref_name] += variant_count
                        st.hists_frameshift[ref_name][effective_length] += variant_count
            else:
                st.counts_non_modified_non_frameshift[ref_name] += variant_count
                _agg_acc(st.insertion_count_vectors_noncoding[ref_name], insertion_positions, variant_count)
                _agg_acc(st.deletion_count_vectors_noncoding[ref_name], deletion_positions, variant_count)
                _agg_acc(st.substitution_count_vectors_noncoding[ref_name], substitution_positions, variant_count)
                st.hists_inframe[ref_name][0] += variant_count
    elif tot_exon_len_mod != 0:
        if (tot_exon_len_mod % 3) == 0:
            st.counts_modified_non_frameshift[ref_name] += variant_count
            st.hists_inframe[ref_name][tot_exon_len_mod] += variant_count
        else:
            st.counts_modified_frameshift[ref_name] += variant_count
            st.hists_frameshift[ref_name][tot_exon_len_mod] += variant_count


def _finalize_aggregates(st: _AggState) -> AlleleAggregates:
    """Derive the post-loop vectors and pack into :class:`AlleleAggregates`.

    Verbatim replica of ``CRISPRessoCORE.main``'s post-loop block (~line
    4184-4190): ``all_indelsub = all_ins + all_del + all_sub``; ``indelsub =
    ins + del + sub`` (windowed); ``substitution_base_vectors`` = the windowed
    subset of ``all_substitution_base_vectors`` via ``refs[ref]['include_idxs']``.
    """
    all_indelsub = {}
    indelsub = {}
    substitution_base_vectors = {}
    for ref_name in st.ref_names:
        all_indelsub[ref_name] = (
            st.all_insertion_count_vectors[ref_name]
            + st.all_deletion_count_vectors[ref_name]
            + st.all_substitution_count_vectors[ref_name]
        )
        indelsub[ref_name] = (
            st.insertion_count_vectors[ref_name]
            + st.deletion_count_vectors[ref_name]
            + st.substitution_count_vectors[ref_name]
        )
        this_include_idx = st.refs[ref_name]["include_idxs"]
        for nuc in ("A", "C", "G", "T", "N"):
            substitution_base_vectors[ref_name + "_" + nuc] = [
                st.all_substitution_base_vectors[ref_name + "_" + nuc][x] for x in this_include_idx
            ]
    return AlleleAggregates(
        all_insertion_count_vectors=st.all_insertion_count_vectors,
        all_insertion_left_count_vectors=st.all_insertion_left_count_vectors,
        all_deletion_count_vectors=st.all_deletion_count_vectors,
        all_substitution_count_vectors=st.all_substitution_count_vectors,
        all_indelsub_count_vectors=all_indelsub,
        all_substitution_base_vectors=st.all_substitution_base_vectors,
        all_base_count_vectors=st.all_base_count_vectors,
        insertion_count_vectors=st.insertion_count_vectors,
        deletion_count_vectors=st.deletion_count_vectors,
        substitution_count_vectors=st.substitution_count_vectors,
        indelsub_count_vectors=indelsub,
        insertion_count_vectors_noncoding=st.insertion_count_vectors_noncoding,
        deletion_count_vectors_noncoding=st.deletion_count_vectors_noncoding,
        substitution_count_vectors_noncoding=st.substitution_count_vectors_noncoding,
        insertion_length_vectors=st.insertion_length_vectors,
        deletion_length_vectors=st.deletion_length_vectors,
        inserted_n_dicts=st.inserted_n_dicts,
        deleted_n_dicts=st.deleted_n_dicts,
        substituted_n_dicts=st.substituted_n_dicts,
        effective_len_dicts=st.effective_len_dicts,
        hists_inframe=st.hists_inframe,
        hists_frameshift=st.hists_frameshift,
        counts_insertion=st.counts_insertion,
        counts_deletion=st.counts_deletion,
        counts_substitution=st.counts_substitution,
        counts_only_insertion=st.counts_only_insertion,
        counts_only_deletion=st.counts_only_deletion,
        counts_only_substitution=st.counts_only_substitution,
        counts_insertion_and_deletion=st.counts_insertion_and_deletion,
        counts_insertion_and_substitution=st.counts_insertion_and_substitution,
        counts_deletion_and_substitution=st.counts_deletion_and_substitution,
        counts_insertion_and_deletion_and_substitution=st.counts_insertion_and_deletion_and_substitution,
        counts_modified_frameshift=st.counts_modified_frameshift,
        counts_modified_non_frameshift=st.counts_modified_non_frameshift,
        counts_non_modified_non_frameshift=st.counts_non_modified_non_frameshift,
        counts_splicing_sites_modified=st.counts_splicing_sites_modified,
        substitution_base_vectors=substitution_base_vectors,
    )


def _collapsed_row_to_payload(row: dict) -> dict:
    """Convert a collapsed-parquet row (dict) to a payload-shaped dict.

    Maps the collapsed parquet's column names to the ``variant_payload`` keys
    that :func:`_aggregate_one_row` (a verbatim replica of the pandas per-ref
    body) expects: ``Aligned_Sequence`` → ``aln_seq``, ``Read_Status`` →
    ``classification``, ``n_inserted`` → ``insertion_n``, etc. Position arrays
    are reconstructed to numpy int arrays and coordinate columns to tuple-lists
    so the fancy-indexing / iteration in the per-ref body behaves identically to
    the pandas path (where ``variant_payload`` comes from
    ``find_indels_substitutions`` emitting numpy arrays).
    """
    return {
        "aln_seq": row["Aligned_Sequence"],
        "aln_ref": row["Reference_Sequence"],
        "classification": row["Read_Status"],
        "insertion_n": int(row["n_inserted"]),
        "deletion_n": int(row["n_deleted"]),
        "substitution_n": int(row["n_mutated"]),
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
        "all_substitution_values": _to_np_str_array(row.get("all_substitution_values")),
    }


def _aggregate_alleles(
    self: "VariantStore",
    refs: dict,
    args,
    ref_names: list,
    *,
    collapsed_path: Optional[str] = None,
    batch_size: Optional[int] = None,
) -> AlleleAggregates:
    """Stream the collapsed allele parquet → per-reference aggregation outputs.

    Replaces the per-ref body of ``CRISPRessoCORE.main``'s allele loop for the
    parquet backend. Scans the collapsed parquet row-group by row-group
    (``pyarrow.iter_batches``), reconstructs each row as a payload-shaped dict,
    and runs :func:`_aggregate_one_row` for every normal (non-AMBIGUOUS,
    non-DISCARDED) row. Peak memory is bounded by one batch + the accumulator
    vectors (O(amplicon_length × num_refs × num_base_types) — the design's
    stated bounded lower-order term), not by the allele count.

    Parameters
    ----------
    refs
        ``{ref_name: ref_dict}`` — must carry ``sequence_length``,
        ``include_idxs``, ``exon_len_mods``, ``exon_positions``,
        ``splicing_positions``, ``contains_coding_seq`` (same structures the
        pandas loop reads from ``refs[ref_name]``).
    args
        The CRISPResso args namespace (reads ``ignore_insertions`` /
        ``ignore_deletions`` / ``ignore_substitutions``).
    ref_names
        Ordered list of reference names (defines accumulator initialization).
    collapsed_path
        Path to the collapsed allele parquet (Stage 3 output). Defaults to
        ``<output_directory>/collapsed.allele.parquet``.
    batch_size
        Row batch size for ``pyarrow.iter_batches`` (bounds peak memory).

    Returns
    -------
    AlleleAggregates
        All per-ref accumulator dicts/vectors/Counters, with the post-loop
        derived vectors (``all_indelsub``, ``indelsub``,
        ``substitution_base_vectors``) finalized — matching the pandas loop's
        outputs for direct consumption by the downstream plot / quantification
        code.
    """
    path = collapsed_path or os.path.join(self.output_directory, "collapsed.allele.parquet")
    st = _init_agg_state(refs, ref_names, args)

    columns = [
        "#Reads",
        "Reference_Name",
        "Aligned_Sequence",
        "Reference_Sequence",
        "Read_Status",
        "n_inserted",
        "n_deleted",
        "n_mutated",
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
        "all_substitution_values",
    ]
    pf = pq.ParquetFile(path)
    if batch_size is None:
        # Cap the adaptive batch so the 21 in-flight ChunkedArrays (several of
        # them amplicon-length list columns) don't land the whole collapsed
        # table in one batch at small scale. The row-by-row ``.as_py()`` access
        # is batch-size-independent for correctness; a 1000-row cap bounds the
        # arrow decompression buffer (design_docs/PARQUET_MEMORY_PROFILE.md
        # item 5).
        batch_size = min(_adaptive_batch_size(pf), 1000)
    for batch in pf.iter_batches(columns=columns, batch_size=batch_size):
        n = batch.num_rows
        if n == 0:
            continue
        # Row-by-row cell access (not batch.to_pylist()) so string-array
        # columns (all_substitution_values) don't blow up peak RSS ~50x.
        col_arrays = {c: batch.column(c) for c in columns}
        ref_col = col_arrays["Reference_Name"]
        reads_col = col_arrays["#Reads"]
        for i in range(n):
            ref_name = ref_col[i].as_py()
            # AMBIGUOUS_* / DISCARDED_* rows do not contribute to the per-ref
            # aggregations — the pandas loop ``continue``\ s before the per-ref
            # body for them (class_counts / counts_discarded are handled by
            # collapse). Their Reference_Name is not in ref_names.
            if ref_name not in st.refs:
                continue
            row = {c: col_arrays[c][i].as_py() for c in columns}
            p = _collapsed_row_to_payload(row)
            _aggregate_one_row(st, ref_name, int(reads_col[i].as_py()), p)

    return _finalize_aggregates(st)


VariantStore.aggregate_alleles = _aggregate_alleles


@dataclass
class _AlnStatsHomologyState:
    """Mutable accumulators for ``aln_stats`` + homology dicts, built per-row.

    Extracted so the SAME per-row logic runs either as a standalone shard
    scan (``compute_aln_stats_and_homology_from_shards``) OR fused into
    collapse's Pass 2 (the fanout scan that already reads ``payloads``),
    avoiding a second full payloads read for homology.
    """
    expand_ambiguous_alignments: bool = False
    N_TOT_READS: int = 0
    N_CACHED_ALN: int = 0
    N_CACHED_NOTALN: int = 0
    N_COMPUTED_ALN: int = 0
    N_COMPUTED_NOTALN: int = 0
    N_GLOBAL_SUBS: int = 0
    N_SUBS_OUTSIDE_WINDOW: int = 0
    N_MODS_IN_WINDOW: int = 0
    N_MODS_OUTSIDE_WINDOW: int = 0
    N_READS_IRREGULAR_ENDS: int = 0
    READ_LENGTH: int = 0
    aln_homology: dict = field(default_factory=dict)
    not_aln_homology: dict = field(default_factory=dict)

    def accumulate(self, read_key, count: int, payload):
        """Fold one shard row (pre-merge) into the stats + homology dicts."""
        self.N_TOT_READS += count
        best = payload.get("best_match_score", 0)
        if best <= 0:
            self.N_COMPUTED_NOTALN += 1
            self.N_CACHED_NOTALN += (count - 1)
            self.not_aln_homology[read_key] = {
                "aln_scores": list(payload.get("aln_scores", []) or []),
                "count": count,
            }
            return
        self.N_COMPUTED_ALN += 1
        self.N_CACHED_ALN += (count - 1)
        self.aln_homology[read_key] = {
            "aln_scores": list(payload.get("aln_scores", []) or []),
            "count": count,
        }
        aln_ref_names = payload.get("aln_ref_names", []) or []
        if len(aln_ref_names) == 1 or self.expand_ambiguous_alignments:
            for name in aln_ref_names:
                sub = payload.get("variant_" + name)
                if sub is None:
                    continue
                if self.READ_LENGTH == 0:
                    self.READ_LENGTH = len(sub["aln_seq"])
                self.N_GLOBAL_SUBS += (sub["substitution_n"] + sub["substitutions_outside_window"]) * count
                self.N_SUBS_OUTSIDE_WINDOW += sub["substitutions_outside_window"] * count
                self.N_MODS_IN_WINDOW += sub["mods_in_window"] * count
                self.N_MODS_OUTSIDE_WINDOW += sub["mods_outside_window"] * count
                if sub["irregular_ends"]:
                    self.N_READS_IRREGULAR_ENDS += count

    def finalize(self):
        """Return ``(aln_stats, aln_homology, not_aln_homology)``."""
        aln_stats = {
            "N_TOT_READS": self.N_TOT_READS,
            "N_CACHED_ALN": self.N_CACHED_ALN,
            "N_CACHED_NOTALN": self.N_CACHED_NOTALN,
            "N_COMPUTED_ALN": self.N_COMPUTED_ALN,
            "N_COMPUTED_NOTALN": self.N_COMPUTED_NOTALN,
            "N_GLOBAL_SUBS": self.N_GLOBAL_SUBS,
            "N_SUBS_OUTSIDE_WINDOW": self.N_SUBS_OUTSIDE_WINDOW,
            "N_MODS_IN_WINDOW": self.N_MODS_IN_WINDOW,
            "N_MODS_OUTSIDE_WINDOW": self.N_MODS_OUTSIDE_WINDOW,
            "N_READS_IRREGULAR_ENDS": self.N_READS_IRREGULAR_ENDS,
            "READ_LENGTH": self.READ_LENGTH,
        }
        return aln_stats, self.aln_homology, self.not_aln_homology


def compute_aln_stats_and_homology_from_shards(
    shard_paths,
    num_unique_reads: int,
    *,
    expand_ambiguous_alignments: bool = False,
):
    """Stream aligned shards → (aln_stats, aln_homology, not_aln_homology).

    Replaces the read-back half of ``CRISPRessoCORE.process_fastq`` (~line
    1910-1985) for the parquet backend. One streaming pass over the per-worker
    shards (``iter_aligned_shard``) computes:

    * ``aln_stats`` — the ``N_TOT_READS`` / ``N_COMPUTED_ALN`` / ``N_CACHED_ALN``
      / ``N_GLOBAL_SUBS`` / ... dict that ``CRISPResso_mapping_statistics.txt``
      and the guardrails consume. Uses the PRE-merge per-unique-read counts
      (the shards are pre-collapse), matching the pandas read-back which counts
      each unique read once before the 3b RC-merge folds partner counts to 0.
    * ``aln_homology`` — ``{read_key: {'aln_scores': list, 'count': int}}`` for
      aligned reads (``best_match_score > 0``), for ``get_and_save_homology_scores``.
    * ``not_aln_homology`` — same shape for unaligned reads
      (``best_match_score <= 0``).

    Notes
    -----
    * ``READ_LENGTH`` is taken from the first aligned read's primary ``aln_seq``
      in scan order — matching the pandas read-back (which sets it from the
      first aligned read in worker-file order). For typical amplicon data every
      read is the same length, so this is order-independent.
    * Homology counts here are PRE-RC-merge (both RC partners carry their
      original counts). The pandas ``variantCache`` at homology time is
      POST-merge (folded partners have ``count = 0``). These differ only when
      the input contains reverse-complement read pairs; the wiring handles
      that case by building homology from the post-merge variant store when
      needed (see the ``CollapsedAlleles`` variants parquet — TODO if the basic
      test reveals RC pairs).
    * In the wired parquet path this standalone scan is NOT used — homology is
      fused into collapse's Pass 2 (``_collapse_streaming_single_read``) to
      avoid a second ``payloads`` read. This function is kept for the exported
      API and standalone use; it shares the per-row logic via
      :class:`_AlnStatsHomologyState`.
    """
    st = _AlnStatsHomologyState(expand_ambiguous_alignments=expand_ambiguous_alignments)
    for path in _expand_shard_paths(shard_paths):
        for read_key, count, payload in iter_aligned_shard(path):
            st.accumulate(read_key, int(count), payload)
    return st.finalize()


def aggregate_alleles_from_collapsed(
    collapsed_path: str,
    refs: dict,
    args,
    ref_names: list,
    output_directory: str,
    *,
    memory_budget_mb: int = DEFAULT_MEMORY_BUDGET_MB,
) -> AlleleAggregates:
    """Convenience wrapper: create a :class:`VariantStore` and aggregate."""
    store = VariantStore(output_directory, memory_budget_mb=memory_budget_mb)
    return store.aggregate_alleles(
        refs, args, ref_names, collapsed_path=collapsed_path,
    )


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
    batch_size: Optional[int] = None,
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
    if batch_size is None:
        batch_size = _adaptive_batch_size(pf)
    # Read only the columns the output projection actually needs (plus #Reads
    # for %Reads). The non-detailed TSV needs only scalar columns — reading
    # the amplicon-length position/substitution arrays would blow up peak RSS
    # ~50x via to_pylist() for no benefit (they're never emitted).
    _stored_cols = set(COLLAPSED_SCHEMA.names)
    read_cols = [c for c in select_cols if c in _stored_cols]
    if "#Reads" not in read_cols:
        read_cols.append("#Reads")
    with open(output_path, "w") as handle:
        first = True
        for batch in pf.iter_batches(columns=read_cols, batch_size=batch_size):
            col_arrays = {c: batch.column(c) for c in read_cols}
            n = batch.num_rows
            row_dicts = []
            for i in range(n):
                row = {c: col_arrays[c][i].as_py() for c in read_cols}
                if write_detailed_allele_table:
                    row_dicts.append(_parquet_row_to_allele_dict(row, n_total))
                else:
                    # Non-detailed TSV: only the scalar crispresso2Cols are
                    # emitted — skip the amplicon-length array conversion
                    # (and the array column read) entirely.
                    row_dicts.append({
                        "#Reads": int(row["#Reads"]),
                        "Aligned_Sequence": row["Aligned_Sequence"],
                        "Reference_Sequence": row["Reference_Sequence"],
                        "n_inserted": int(row["n_inserted"]),
                        "n_deleted": int(row["n_deleted"]),
                        "n_mutated": int(row["n_mutated"]),
                        "Reference_Name": row["Reference_Name"],
                        "Read_Status": row["Read_Status"],
                    })
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
    ``df_alleles``): ``find_indels_substitutions`` (``CRISPRessoCOREResources.pyx``)
    emits the position/size arrays as **Python lists** (``ref_positions=[]`` +
    ``.append``), the coordinate arrays as lists of tuples, and ONLY
    ``substitution_values`` / ``all_substitution_values`` as numpy arrays
    (``np.array(...)``). Plot code calls ``.index()`` on ``ref_positions``
    (~line 790, 1653), which requires a list — so position/size columns must
    round-trip as lists, not numpy. Scalar columns pass through unchanged.

    Tolerant of both storage shapes: parquet-scan cells (arrow lists of int /
    lists of struct dicts) and in-memory ``allele_rows`` cells (Python lists /
    lists of tuples / numpy arrays) so the same reconstruction applies on both
    the persisted-parquet and the write_parquet=False paths.
    """
    if col_name in _INT_ARRAY_SLICE_COLS:
        # find_indels emits Python lists; preserve that (plot code .index()es them).
        if val is None:
            return []
        if isinstance(val, np.ndarray):
            return val.tolist()
        return [int(x) for x in val]
    if col_name in _COORD_SLICE_COLS:
        if val is None:
            return []
        if isinstance(val, np.ndarray):
            return _structs_to_coords(val.tolist())
        if val and isinstance(val[0], tuple):
            return list(val)  # already a tuple-list (in-memory path)
        return _structs_to_coords(val)
    if col_name in _STR_ARRAY_SLICE_COLS:
        # substitution_values / all_substitution_values are numpy arrays in the
        # real payload (find_indels emits np.array(...)).
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
    batch_size=None,
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
        if batch_size is None:
            batch_size = _adaptive_batch_size(pf)
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
