# Fix #3 — Streaming single-read collapse (external-sort canonical-key merge)

**Status:** IMPLEMENTED (2026-07-07) — streaming single-read collapse live;
SC #1 verified at 2000 bp (RSS(100k)/RSS(10k) = 2.20 < 3).
**Blocks:** Success Criterion #1 (flat peak RSS).
**Supersedes:** the master-plan Stage 3b spec (polars `group_by(canonical_key)`),
which the spike (`STREAMING_GROUPBY_SPIKE.md`) falsified for high-cardinality
keys and which PR 5's in-memory `_collapse_rekey_and_rcmerge` only implements
for the bounded-cardinality case.

## The problem (confirmed by `scripts/bench_pipeline_memory.py`, 2026-07-07)

`_collapse_rekey_and_rcmerge` (CRISPResso2/storage.py) builds an in-memory
`store: dict[collapse_key, {count, payload}]`. For **paired** input the collapse
key is the primary `aln_seq` (3a re-key) → unique *allele* count ≪ read count →
the dict is bounded (the design's stated bounded term). For **single-read**
input there is **no 3a re-key** (`key = read_key`), and 3b only collapses
reverse-complement pairs — so for high-diversity single-read data (PacBio HiFi /
Nanopore, ~1:1 diversity, few RC pairs) the store holds **one entry per unique
read**, each carrying the full payload (~10 Python lists of ~amplicon-length
ints ≈ 12–16 KB/read).

End-to-end benchmark (200 bp amplicon, ~1:1 diversity, n_processes=2):

| reads (input → unique) | parquet peak RSS |
|------------------------:|-----------------:|
| 10k → 5.2k              | 436 MB           |
| 100k → 44k             | 1308 MB          |
| 1M → 290k              | 3519 MB          |

The collapse store is the dominant term (~12 KB × 290k ≈ 3.5 GB at 1M — ~40×
the Stage 1 dedup dict). SC #1 fails: RSS(1M)/RSS(10k) = 8.08 (target < 3).
Fixes #1 (stream worker chunks) and #2 (Stage 1 budget 512→128 MB) are kept but
barely move the number — they address the two *smaller* in-memory copies, not
this one.

## The fix

Replace the in-memory `store` with a streaming **external-sort canonical-key
merge** for the single-read path. Paired stays in-memory (it is genuinely
bounded by allele count). The approach mirrors Stage 1's S1b-proven pattern
(`STREAMING_GROUPBY_SPIKE.md`: external merge-sort → streaming first-per-group
collapse → flat RSS regardless of unique count).

### Algorithm (single-read path only)

Input: the aligned shard parquets (`aligned_{i}.parquet`, one row per unique
read, schema `ALIGNED_SCHEMA`), produced by Stage 2.

1. **Compute the canonical key.** For each shard row, `canonical_key =
   min(read_key, reverse_complement(read_key))`. (For single-read, 3a is a
   no-op — the collapse key IS the read_key — and 3b folds RC pairs via the
   canonical key. `read_key` is the shard's `read_key` column.)
2. **Spill to a sortable temp parquet.** Stream the shards via
   `iter_aligned_shard` (bounded batches) and write one row per read to
   `collapsed.sort.{tmp}.parquet` with schema = `ALIGNED_SCHEMA` + a
   `canonical_key: string` column + a monotonic `seq_no: int64` (the original
   scan order, for the tie-break). The payload columns are copied through
   unchanged (native arrow lists/structs — no JSON). This file is large (≈
   unique-read count × payload size) but lives on **disk**, not RAM.
3. **External-sort by `canonical_key`.** Sort `collapsed.sort.{tmp}.parquet` by
   `canonical_key` into `collapsed.sorted.{tmp}.parquet`. The sort MUST spill
   to disk (external merge-sort) so peak RSS is bounded by the sort buffer, not
   the row count. **Mechanism to confirm with a 15-min spike** (see below);
   fallback is proven.
4. **Streaming first-per-group collapse.** Stream `collapsed.sorted.{tmp}.parquet`
   row-group by row-group (`pq.iter_batches`). Adjacent rows share a
   `canonical_key` (post-sort). For each group: sum `count`, keep the
   representative payload = the row with the **smallest `seq_no`** in the group
   (see tie-break below). Emit one collapsed (count, payload) record per group.
   Peak RSS = one group + the accumulators (the same bounded terms
   `aggregate_alleles` already uses).
5. **Fan-out (3c).** Feed each collapsed record to the existing
   `_collapse_fanout` per-variant logic (AMBIGUOUS / DISCARDED / per-ref
   explode → `_get_allele_row`), writing directly to `collapsed.allele.parquet`.
   Unchanged from PR 5.

Output is identical to the in-memory path: the same `collapsed.allele.parquet`,
`CollapsedAlleles` (allele_rows, n_total, class_counts, counts_*), just produced
without holding the whole store.

### Tie-break / parity

The pandas path's representative for an RC pair is the **first-occurrence in
scan order** (insertion order). PR 5's in-memory store replicates this exactly
(first occurrence wins). The streaming version picks the **min `seq_no`** per
canonical-key group — which IS first-occurrence-in-scan-order, so parity holds
**by construction** (`seq_no` encodes scan order). Note: on the Stage 1 *spill*
path, scan order is already sorted-key order (the dedup parquet is sorted), so
the streaming collapse's "sorted-canonical-key-first" is consistent with the
existing spill-path behavior — no NEW parity divergence. The eager-Stage-1 path
keeps FASTQ-first scan order via `seq_no`. Edge #8 (palindrome: `rc(seq)==seq`)
is handled naturally — the group has one member, count unchanged.

### Sort mechanism (the one open implementation question)

The spike (`STREAMING_GROUPBY_SPIKE.md` S1b) proved **system `sort`** (external
merge-sort, `-S` buffer cap) stays flat on multi-kb string keys. polars
streaming `group_by` does NOT (Premise #3 falsified). Whether **polars streaming
`sort`** spills is unverified.

- **Preferred (verify first):** `pl.scan_parquet(tmp).sort("canonical_key").sink_parquet(sorted)`.
  If a 15-min spike (1M × 2 kb `canonical_key`, measure peak RSS vs row count)
  shows RSS flat, use it — cleanest (one polars call, native parquet round-trip).
- **Fallback (proven):** shell out to system `sort` on a text projection
  `canonical_key \t seq_no \t shard_row_offset`, `-S 64M -T <workdir> -s` (stable),
  then a streaming pass that reads the sorted text and pulls each representative
  payload from the temp parquet by `seq_no`/offset. This is exactly Stage 1's
  pattern (`_count_spill`) so the plumbing exists. Slower (text + seek-back) but
  guaranteed flat.

Pick the mechanism by the spike; the surrounding algorithm is mechanism-agnostic.

### Integration points

- **Dispatch:** in `_collapse`, branch on `is_paired`. Paired → existing
  `_collapse_rekey_and_rcmerge` (in-memory, bounded). Single-read → new
  `_collapse_streaming_single_read` (this design).
- **`CollapsedAlleles` outputs:** identical — the streaming path produces the
  same `allele_rows` / `n_total` / `class_counts` / `counts_*` / `parquet_path`.
- **`files_to_remove`:** add the temp sort parquets (`collapsed.sort.*`,
  `collapsed.sorted.*`); cleanup in a `finally` (edge #15/#21).
- **`re_aln` reads** (`caching_is_ok=False`, paired only — edge #13): unaffected;
  single-read has no `re_aln` loop.
- **Empty input / all-unaligned:** the streaming sort over 0 rows must emit an
  empty `collapsed.allele.parquet` (parity with the in-memory path's empty case).

### Test plan

- **P0 parity:** canned single-read shards through BOTH `_collapse` paths
  (in-memory oracle + streaming) → identical `CollapsedAlleles` (allele_rows,
  n_total, class_counts, counts_*). Cover: no RC pairs; one RC pair (tie-break
  = min `seq_no`); palindrome; AMBIGUOUS; DISCARDED; multi-ref fan-out.
- **P0 SC #1 re-run:** `scripts/bench_pipeline_memory.py` at 10k/100k/1M (200 bp)
  AND a 2–5 kb amplicon variant (to exercise Stage 1 spill alongside #3) →
  assert RSS(1M)/RSS(10k) < 3.
- **P1:** empty input; single-read; all-unaligned.
- **P1:** `make basic` (FANC.Cas9) still byte-identical (the basic test is
  single-read — it now exercises the streaming collapse).

### Out of scope for this fix

- Paired-input collapse (already bounded; no change).
- The homology RC-pair count semantics (latent parity risk #2 in the master
  plan) — that needs the *post-merge* variant store persisted, which is a
  separate artifact (OQ-E3). The streaming collapse does not by itself fix
  homology parity; homology still reads the pre-merge shards. (If the homology
  fix and this fix both need a post-merge variants parquet, they should share
  one artifact — note for the implementer.)
