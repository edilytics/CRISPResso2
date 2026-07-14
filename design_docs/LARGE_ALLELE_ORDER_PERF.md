# Performance: large-allele final ordering (external sort without JSON carry)

**Status:** DESIGN (2026-07-07) — not yet implemented.
**Blocks:** Success Criterion #2 (long-read completion) at scale — RSS is flat
(`STREAMING_SINGLE_READ_COLLAPSE.md`), but the 1M × 2000 bp run times out
because the final allele ordering carries ~12.8 GB of JSON through a global
external sort. This doc removes that bottleneck.
**Supersedes:** the external-sort branch of `_collapse_streaming_single_read`
step 5 (the JSON-carry `output.jsonl` → `sort` → `parquet` path), keeping the
in-memory branch for small inputs (tests, `make basic test`).

## The problem (confirmed by `scripts/bench_pipeline_memory.py`, 2026-07-07)

After Fix #3 (`STREAMING_SINGLE_READ_COLLAPSE.md`), peak RSS is flat for
single-read input at 2000 bp: SC #1 PASSES (RSS(100k)/RSS(10k) = 2.20 < 3).
But the 1M × 2000 bp run **times out** (5400 s wall, target < 5400):

| reads (×2000 bp) | unique allele rows | parquet peak RSS | wall   | status |
|----------------:|-------------------:|----------------:|-------:|:-------|
| 10k             | ~7.7k              | 837 MB          | 112 s  | OK     |
| 100k            | ~44k               | 1838 MB         | 759 s  | OK     |
| 1M              | ~290k              | ~1.8 GB (flat)  | >5400 s| TIMEOUT |

RSS is flat (~1.8 GB, same order as 100k) — the streaming collapse (Fix #3)
succeeded. The timeout is **purely the final allele-ordering step** in
`_collapse_streaming_single_read`'s external branch: it writes every collapsed
allele row as JSON to `output.jsonl` (~290k × ~44 kB ≈ 12.8 GB — the full payload
including amplicon-length position arrays, JSON-encoded), runs
`system sort` over the 12.8 GB, then parses the sorted JSON back into parquet.
The sort does ~5–9 merge passes over 12.8 GB (64 MB buffer) = ~64–128 GB of
I/O, plus ~25 GB of JSON serialize/parse CPU. At 200 MB/s disk + JSON CPU that
is ~20 min — over the timeout.

The root cause is **carrying the full row data through a global text sort**.
The sort only needs the key `(-#Reads, Aligned_Sequence, Reference_Sequence)`
plus a stable tiebreaker; the amplicon-length position arrays are dead weight
through every merge pass and through the JSON round-trip.

## Spikes (2026-07-07) — what does NOT spill

The obvious "use a spill-capable parquet sort instead of JSON-carry" was
spiked and falsified, mirroring the earlier canonical-key-sort spike
(`bench_polars_sort_spill.py`):

| engine (sort large parquet rows by a key) | 1M × 2 kb peak RSS | spills? |
|-------------------------------------------|-------------------:|:--------|
| polars `scan_parquet().sort().sink_parquet()` | 4367 MB (ratio 33 vs 10k) | NO (eager) |
| DuckDB `COPY (... ORDER BY ...) TO parquet` (`memory_limit=512MB`, `preserve_insertion_order=false`, `threads=2`) | 5867 MB (OOM) | NO (eager) |
| system `sort` (external merge-sort, `-S` buffer cap) | ~410 MB (ratio 1.24) | **YES** |

Neither polars `sort` nor DuckDB `ORDER BY` spills for large-row sorts — both
materialize the table eagerly and OOM. **System `sort` remains the only
spill-capable mechanism.** So the fix must keep system `sort` but stop feeding
it the full rows.

## The fix — key-sort + partitioned gather (bucket-by-output-position)

Replace the global JSON-carry sort with a classic **partitioned external
gather**: sort the *small* key projection (sort keys + a row index) with system
`sort`, then gather the full rows from an unsorted parquet into per-bucket
spill files and emit each bucket in order. The full row data never passes
through a global sort — it is read once, scattered to buckets, and each bucket
is sorted in memory (it fits).

### Algorithm (replaces the external branch of step 5)

Inputs: the collapsed groups `(canonical_key, total_count, rep_seq_no)`,
sorted by `rep_seq_no` (scan order), plus the aligned shards (Pass 2 input).

**Step A — fanout + key capture (one streaming pass over the shards).**
Re-stream the shards; for each representative (`rep_seq_no` match) run
`_collapse_fanout` → 1+ allele rows. Assign each emitted row a monotonic
`row_idx`. Write the FULL row to `collapsed.unsorted.parquet` (arrow,
`COLLAPSED_SCHEMA` + a `row_idx: int64` column) — compact native arrow, no
JSON. Simultaneously append to `keys.txt`:
`{neg_reads:015d}\t{Aligned_Sequence}\t{Reference_Sequence}\t{row_idx}` — the
sort key + tiebreaker only (~4 kB/row at 2 kb amplicon, vs ~44 kB/row JSON).
The aln_seq/ref_seq are read from the row being emitted (no separate seek-back).

**Step B — key sort (small).** `system sort -s -S <buf> -k1,1 -k2,2 -k3,3
keys.txt` → `keys.sorted`. Stable (`-s`) so equal keys preserve `row_idx`
(fanout/scan order = the pandas tie-break). At 1M × 2 kb this is ~1.2 GB
through the sort (vs 12.8 GB for JSON-carry) — ~10× less I/O per merge pass.

**Step C — bucket assignment.** Read `keys.sorted` (output order). Build
`out_pos[row_idx]` = output index (numpy int64 array, N elements ≈ 2.3 MB at
290k). `chunk_size = max(1, memory_budget_bytes // est_row_size)`;
`bucket = out_pos[row_idx] // chunk_size`. Each bucket holds `chunk_size` rows
and fits in memory.

**Step D — partitioned gather (one streaming pass over the unsorted parquet).**
Stream `collapsed.unsorted.parquet` (row_idx 0..N-1, adaptive batches). For
each row, look up its `bucket` and append `(out_pos, full_row)` to that
bucket's **spill file** (`bucket_{b}.txt`, compact TSV: `out_pos` + a
single-line encoding of the row). K spill files are open simultaneously (K ≈
N/chunk_size; manageable file-handle count, each with a modest write buffer).
Peak RSS = K write buffers (small); temp disk = N rows total (structured,
~1.5× arrow, NOT through any global sort).

**Step E — emit (per bucket, in order).** For `b = 0..K-1`: read
`bucket_{b}.txt` (chunk_size rows), sort in memory by `out_pos`, strip
`out_pos`, write to the final `collapsed.allele.parquet` (batched pyarrow
writes). Peak RSS = one bucket (chunk_size rows). Buckets are emitted in
order → the final parquet is globally sorted by
`(-#Reads, Aligned_Sequence, Reference_Sequence, row_idx)` — byte-identical
to the pandas path's `df_alleles.sort_values` (the `row_idx` tiebreaker
encodes fanout/scan order, matching the stable `-s` sort the JSON-carry used).

**Step F — cleanup.** Remove `collapsed.unsorted.parquet`, `keys.txt`,
`keys.sorted`, and `bucket_*.txt` in a `finally` (edges #15/#21).

### I/O / time estimate (1M × 2000 bp, ~290k allele rows, ~30 kB arrow/row)

| step                  | data moved                 | I/O (×merge passes)   | est. wall |
|-----------------------|----------------------------|-----------------------|----------:|
| A — fanout + write    | 8.7 GB parquet + 1.2 GB keys| 1× each               | ~1.5 min  |
| B — key sort          | 1.2 GB keys                | ~2 passes (256 MB buf)| ~0.5 min  |
| D — gather (read+spill)| 8.7 GB parq → 13 GB spill  | 1× read + 1× write    | ~3 min    |
| E — emit (read+write) | 13 GB spill → 8.7 GB parq  | 1× read + 1× write    | ~3 min    |
| **total**             |                            |                       | **~8 min**|

vs the JSON-carry path's ~20 min (12.8 GB × ~6 passes global sort + ~25 GB
JSON serialize/parse). ~2.5–3× faster — brings 1M × 2000 bp comfortably under
the timeout, while keeping RSS flat (bounded by one bucket + accumulators).

### Quick partial win (independent, ship first if desired)

Bumping the system-sort buffer from `-S 64M` to `-S 256M`–`1G` (still constant,
so still "flat" per SC #1) cuts the JSON-carry merge passes ~3× and is a
one-line change. It does NOT remove the JSON serialize/parse CPU cost, so it
is a partial improvement (~30–40% faster, not the 3× above). Ship as a
standalone quick win or fold into the partitioned-gather PR; the
partitioned-gather design supersedes it either way.

**STATUS (2026-07-07): the JSON→TSV row-carry swap is IMPLEMENTED** (commit
pending). The external branch now writes compact TSV lines (sort-key prefix +
row data, no JSON — arrays are space/semicolon-joined; ~2–5× smaller per row
and ~2–3× faster to round-trip than JSON). Measured at 100k × 2000 bp:
RSS(100k)/RSS(10k) improved 2.20 → 1.84 (smaller temps → lower peak RSS);
byte-identical parity retained (the 4 forced-external tests pass; `make basic
test` passes). The 1M × 2000 bp run still times out — TSV shrinks the carried
bytes ~1.7× but the full payload (incl. amplicon-length position arrays) still
passes through the global sort, so 1M completion still needs the
partitioned-gather design below. TSV is also the format the gather uses for
its key projection + bucket spills, so this swap is a prerequisite win that
the gather builds on.

### Parity

The `row_idx` tiebreaker makes the order deterministic AND identical to the
current JSON-carry stable sort: `row_idx` = fanout emission order = scan order
= the pandas insertion order that `sort -s` preserved. The output is
byte-identical to the in-memory path (and thus the pandas oracle) on canned
payloads — the existing forced-external parity tests
(`test_streaming_external_*`) cover this once the external branch is swapped.

### Integration points

- **Dispatch:** unchanged — `_collapse` routes single-read to
  `_collapse_streaming_single_read`; the `use_in_memory` threshold still
  picks the in-memory branch for small inputs (tests, `make basic test`).
  Only the `else` (external) branch changes: replace the `output.jsonl` →
  `sort` → `_write_collapsed_allele_parquet_from_jsonl` sequence with
  steps A–E above.
- **`CollapsedAlleles` outputs:** identical (`allele_rows` stays `[]` on the
  external branch — the pipeline consumes the parquet via `get_slice` / the
  streaming TSV sink; `n_total` / `class_counts` / `counts_*` built in step A).
- **Temp lifecycle:** `files_to_remove` / `finally` extended with
  `collapsed.unsorted.parquet`, `keys.*`, `bucket_*.txt` (edges #15/#21).
- **`_allele_dict_to_parquet_row` / `COLLAPSED_SCHEMA`:** reused as-is for the
  unsorted and final parquet writes; the compact TSV spill encoding needs a
  new (small) row ↔ text pair (symmetric to `_allele_row_to_json` but TSV).
- **No new dependencies:** system `sort` (already used) + pyarrow (already a
  dep). DuckDB is NOT needed (spike falsified its ORDER BY for large rows).

### Test plan

- **P0 parity:** swap the external branch; the existing forced-external tests
  (`test_streaming_external_matches_in_memory_rc_pair`,
  `test_streaming_external_palindrome_doubles`,
  `test_streaming_external_fuzz_matches_in_memory`,
  `test_streaming_external_empty_input`) must still pass byte-identical.
- **P0 SC #2 (completion):** `scripts/bench_pipeline_memory.py` at
  1M × 2000 bp completes (wall < 5400 s, ideally < 1800 s) with flat RSS.
- **P1:** `make basic test` still byte-identical (the in-memory branch serves it).
- **P1:** a 100k × 5 kb run (exercises the partitioned gather with ~5 kB keys
  and many buckets) completes with flat RSS.

### Out of scope

- The in-memory branch (small inputs) — unchanged.
- Paired-input collapse — unchanged (bounded by allele count).
- The canonical-key sort (step 2 of the collapse) — already the small
  `keys.txt` projection; this doc only changes the FINAL ordering (step 5).
- The Stage 5 `get_slice` whole-frame materialization for real (plots-on) runs
  — a separate follow-up (logged in `bench_pipeline_memory_results.md`); this
  doc is the long-read *completion* perf fix, not the plots-on RSS fix.
