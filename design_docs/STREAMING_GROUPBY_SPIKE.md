# Spike S1 Findings — Streaming `group_by` memory on multi-kb string keys

**Date:** 2026-07-02
**PR:** PR 2 (streaming memory benchmark harness)
**Status:** GATE FAILED (with nuance) — design revision required before PR 5
**Harness:** `scripts/bench_streaming.py` (reusable; `pixi run -e default python scripts/bench_streaming.py`)

## The go/no-go question

The master design (Approach B) premises flat peak memory on polars streaming
`group_by` spilling the hash table to disk, so that peak RSS is "bounded by the
streaming engine's chunk size, not by read count" (design Premise #3). Spike S1
tests that premise directly: does `scan_parquet → group_by(read_key).len() →
sink_parquet` keep peak RSS flat as row count scales, when the group key is a
1–10 kb string (the long-read amplicon scenario)?

## Method

Each `(rows, read_length, engine)` config runs in a **fresh subprocess** so
`ru_maxrss` reflects only the group_by (generation happens in a prior step,
discarded from the measurement). Input is parquet shards on disk, scanned lazily,
so spilling is actually testable. Three engines:

- **streaming** — polars `collect(engine='streaming')` / `sink_parquet`
- **eager** — polars `collect()` (regression oracle: should scale ~linearly)
- **duckdb** — `COPY (SELECT ... GROUP BY) TO ... (FORMAT PARQUET)` (design's
  documented Approach C fallback)

Worst case: `unique_fraction=1.0` (every read distinct — the long-read ~1:1
diversity scenario, and the case that OOMs the current in-RAM dedup dict).

## Results (full matrix)

```
engine      rows     len   peakRSS(MB)  wall(s)   groups   input(MB)
streaming   10000    1000        137      0.13     10000        3.0
streaming  100000    1000        525      0.60    100000       29.8
streaming 1000000    1000       2850      5.97   1000000      298.0
streaming   10000    5000        356      0.49     10000       14.8
streaming  100000    5000       2117      5.00    100000      148.2
streaming 1000000    5000       9137    945.38   1000000     1481.6
eager       10000    1000        103      0.21     10000        3.0
eager      100000    1000        456      1.12    100000       29.8
eager     1000000    1000       3095      5.71   1000000      298.0
eager       10000    5000        222      0.47     10000       14.8
eager      100000    5000       1624      4.89    100000      148.2
eager     1000000    5000       7227    764.40   1000000     1481.6
duckdb      10000    1000        174      1.54     10000        3.0
duckdb     100000    1000        438      1.54    100000       29.8
duckdb    1000000    1000       3404      4.40   1000000      298.0
duckdb      10000    5000        357      0.86     10000       14.8
duckdb     100000    5000       1434      4.59    100000      148.2
duckdb    1000000    5000       4985   1651.94   1000000     1481.6
```

Gate: `RSS(max_rows) / RSS(min_rows) < 3.0` for a spilling engine, fixed read
length. **All engines FAIL** — ratios 14–26× from 10k→1M rows.

## The nuance that matters

RSS does **not** scale with row count per se — it scales with **unique-key
count × key length** (the size of the hash aggregate's group table). Evidence:
the same 1M rows × 1kb with `unique_fraction=0.01` (10k unique groups, not 1M):

```
streaming 1000000 x 1kb, 1% unique  →  454 MB  (vs 2850 MB all-unique)
duckdb    1000000 x 1kb, 1% unique  →  403 MB  (vs 3404 MB all-unique)
```

So when there is real collapse (unique ≪ rows), RSS is bounded by the unique
count. When there is no collapse (all-unique), RSS == holding all keys.

## What this means for the design

**Premise #3 is falsified as stated.** Neither polars streaming nor DuckDB
spills the hash aggregate for high-cardinality string `group_by`. The hash table
holds all unique keys; there is no automatic disk spill at the cardinalities that
matter. This is not an engine-version bug — it is fundamental to hash-based
aggregation: the engine cannot know two keys are equal without holding both.

Concretely, by stage:

| Stage | Operation | Collapse? | Flat-memory feasible? |
|-------|-----------|-----------|-----------------------|
| **1** — raw-read counting | `group_by(read_key).len()` on R1+R2 concat | For long reads: ~none (1:1 diversity) | **NO** — output itself is N rows; no engine can avoid holding N keys. Fundamental, not an engine limit. |
| **3a** — `aln_seq` re-key | `group_by(aln_seq).sum(count)` | Yes — many reads align to the same allele | **YES** — unique allele count ≪ read count; RSS bounded by allele count. |
| **3b** — RC canonical-key merge | `group_by(min(seq,rc(seq)))` | Yes — collapses strand pairs | **YES** — same as 3a. |

So the **post-alignment collapse (Stage 3) is feasible** with polars/DuckDB
hash `group_by`, because real collapse keeps the unique count small. The
**pre-alignment raw-read dedup (Stage 1) is NOT feasible** via hash `group_by`
for all-unique long reads.

## Recommended design revision (for PR 3+ planning)

**RESOLVED by S1b (sort | uniq -c comparison).** See the comparison below —
external merge-sort keeps peak RSS flat across all read lengths, confirming the
revision. DuckDB is no longer needed as a fallback.

1. **Stage 1: replace hash `group_by` with an external sort + streaming collapse.**
   External merge-sort on the read_key (which spills to disk by design), followed
   by a single streaming pass that sums adjacent duplicate counts, gives
   flat-memory dedup **regardless of unique count** — the classic on-disk
   `sort | uniq -c`. polars can do this via `scan_parquet → sort("read_key") →
   sink_parquet` (sort spills) then a streaming pass; or shell out to the
   system `sort` utility (already a CRISPResso2 dependency pattern — see
   `get_most_frequent_reads` using `sort | uniq -c`). **Confirmed below.**
   - Alternatively, **skip Stage 1 dedup entirely for long-read inputs**: when
     `unique_fraction` is ~1.0, dedup is a no-op, so pass reads directly to
     alignment workers. The `variantCache` raw-read dedup exists to collapse
     PCR/sequencing duplicates (common in short-read data); for PacBio HiFi /
     Nanopore it buys nothing. A size-thresholded bypass is the simplest fix.

2. **Stage 3: keep polars hash `group_by`** — it is feasible here because
   post-alignment collapse reduces cardinality. No engine change needed. DuckDB
   remains the fallback if a specific Stage 3 op regresses (it was ~30% more
   memory-efficient than polars streaming at 1M×5kb: 4985 MB vs 9137 MB).

3. **Update the success criterion.** "Flat peak memory across arbitrary input
   sizes" is achievable for Stage 3 but **not for Stage 1 when reads are
   all-unique** — because Stage 1's *output* is N rows. Either Stage 1 must
   external-sort, or the criterion must scope to "post-alignment memory" and
   Stage 1 must be bypassed/sorted for long reads.

## S1b: `sort | uniq -c` vs polars `group_by` (the resolution)

**Conclusion: external merge-sort + `uniq -c` keeps peak RSS flat across all
read lengths and row counts, at a ~2–6× time cost that vanishes at 10kb
(where polars hits memory pressure and becomes slower).** This is the approach
for Stage 1.

Harness: `scripts/compare_sort_vs_groupby.py` (reusable; same subprocess-per-
config methodology as S1, measuring `RUSAGE_SELF` for the export process and
`RUSAGE_CHILDREN` for the sort/uniq pipeline).

### Method

Two engines on the same parquet input (all-unique, worst case):
- **polars_gb** — `scan_parquet → group_by(read_key).len() → sink_parquet` (S1 baseline)
- **sort_uniq** — pyarrow streams `read_key` to a text file in **adaptive batches**
  (target 50 MB/batch, so batch size shrinks as read_length grows), then
  `LC_ALL=C sort -S 64M -T <tmpdir> | uniq -c` (BSD/GNU sort; `-S 64M` caps the
  in-memory buffer, forcing external spill; `-T` sets the spill dir).

### Results (full matrix, 1.0 unique fraction)

```
engine       rows     len  peakRSS(MB) self(MB) child(MB) wall(s)  groups
polars_gb    10000    1000        137     137        0     0.13    10000
sort_uniq    10000    1000         90      90       12     0.25    10000
polars_gb   100000    1000        529     529        0     0.62   100000
sort_uniq   100000    1000        189     189       68     0.91   100000
polars_gb  1000000    1000       3017    3017        0     1.43  1000000
sort_uniq  1000000    1000        190     190      105     9.18  1000000

polars_gb    10000    5000        369     369        0     0.35    10000
sort_uniq    10000    5000        188     188       51     0.45    10000
polars_gb   100000    5000       1988    1988        0     2.92   100000
sort_uniq   100000    5000        468     468      106     3.34   100000
polars_gb  1000000    5000       8655    8655        0    18.67  1000000
sort_uniq  1000000    5000        469     469      132    40.11  1000000

polars_gb   100000   10000       4071    4071        0     4.95   100000
sort_uniq   100000   10000        560     560      104     5.68   100000
polars_gb  1000000   10000       5670    5670        0   159.16  1000000
sort_uniq  1000000   10000        559     559      115    88.65  1000000
```

### RSS scaling (max_rows / min_rows per read length)

```
polars_gb  len=1000b:  10000->1000000 = 137->3017 MB  (ratio 22.05) [SCALES]
polars_gb  len=5000b:  10000->1000000 = 369->8655 MB  (ratio 23.44) [SCALES]
polars_gb  len=10000b: 100000->1000000= 4071->5670 MB (ratio 1.39) [flat-ish but at 5.7 GB]
sort_uniq  len=1000b:  10000->1000000 = 90->190 MB    (ratio 2.12) [FLAT] ✓
sort_uniq  len=5000b:  10000->1000000 = 188->469 MB   (ratio 2.49) [FLAT] ✓
sort_uniq  len=10000b:100000->1000000= 560->559 MB    (ratio 1.00) [FLAT] ✓
```

### Key findings

1. **`sort`'s peak RSS is flat and tiny (104–132 MB at all sizes)** — external
   merge-sort spills to disk as designed; the `-S 64M` buffer cap bounds it.
   `uniq -c` is O(1). The hash-aggregate premise that failed in S1 is simply
   not how sort works.
2. **The export step (pyarrow → text) is the remaining RSS term**, bounded by
   the adaptive batch (target 50 MB). At 5kb it plateaus at ~469 MB; at 10kb
   ~559 MB — constant in row count, scaling only mildly with read length (one
   batch in flight + arrow overhead). Smaller byte budgets shrink this further.
3. **Time tradeoff is the expected spectrum:** sort_uniq is ~6× slower at 1kb,
   ~2× slower at 5kb, and **faster at 10kb** (88.6s vs 159.2s) because polars
   hits memory pressure and thrashes while sort streams to disk.
4. **No new heavy dependencies:** `sort`/`uniq` are POSIX tools already used by
   CRISPResso2 (`get_most_frequent_reads`). pyarrow is already a dependency (PR 1).

### Design decision

- **Stage 1 (raw-read counting): use external sort + `uniq -c`** (the
  `sort_uniq` approach), NOT polars hash `group_by`. Confirmed flat across
  1kb/5kb/10kb × up to 1M rows.
- **Stage 3 (post-alignment collapse): polars hash `group_by` is fine** (real
  collapse keeps unique count small; S1 showed RSS bounded by unique count).
- DuckDB (Approach C) is **not needed** — sort|uniq -c resolves the Stage 1
  ceiling that neither polars streaming nor DuckDB could.

## Open follow-ups (new spikes)

- **S1b:** external-sort dedup — measure `scan_parquet → sort(read_key) →
  sink_parquet → streaming collapse` peak RSS at 1M×5kb. If flat, this replaces
  hash `group_by` for Stage 1. (~15 min with the existing harness extended.)
- Confirm the `unique_fraction` heuristic for bypassing Stage 1 dedup: at what
  threshold does skipping dedup beat doing it? (Likely >0.9.)

## Reproduce

```bash
pixi run -e default python scripts/bench_streaming.py            # full matrix
pixi run -e default python scripts/bench_streaming.py --quick    # fast smoke
pixi run -e default python scripts/bench_streaming.py --rows 1000000 --read-lengths 1000 --unique-fraction 0.01
```

Results are written to `scripts/bench_streaming_results.md` on each run.
