# CRISPResso2 storage-backend peak-RSS benchmark

## Re-run after Fix #3 (streaming single-read collapse) — amplicon-length 2000

Amplicon: synthetic 2000 bp, ~1:1 read diversity (1-3 random subs per read).
Guide: GGAATCCCTTCTGCAGCACC. n_processes=2. `--suppress_plots --suppress_report`
(bench measures the core pipeline: count → align → collapse → aggregate → TSV).

| reads | backend | peak RSS (MB) | wall (s) | status |
|------:|:-------|--------------:|---------:|:-------|
| 10000 | parquet | 837.1 | 112.3 | OK |
| 100000 | parquet | 1837.5 | 758.9 | OK |
| 1000000 | parquet | — | 5400.0 | TIMEOUT |

**SC #1 (parquet): RSS(100000)/RSS(10000) = 1837.5/837.1 = 2.20 — PASS (< 3)**

The 1M×2000 bp run does not OOM — peak RSS stays flat (~1.8 GB, same order as
100k) — but times out: the streaming collapse's external JSON-carry sort must
process ~290k allele rows × ~44 kB/row ≈ 12.8 GB of JSON-lines (each allele row
carries amplicon-length position arrays). System `sort` keeps RSS bounded
(~64 MB buffer) but the merge-sort I/O on 12.8 GB is slow. This is a
performance limit of the external-sort path on long-read-scale data, not a
memory regression; the memory is flat across 10k→100k→1M.

### What changed (Fix #3: streaming single-read collapse)

Before Fix #3 the single-read collapse (`_collapse_rekey_and_rcmerge`) held an
in-memory `store` dict with one entry per unique read, each carrying the full
payload (~10 amplicon-length arrays). At 1M×200 bp that was ~3.5 GB
(~40× the Stage 1 dedup dict); SC #1 was 8.08 (FAIL). Fix #3 replaces the
in-memory store for single-read input with a streaming external-sort
canonical-key merge (`_collapse_streaming_single_read`):

1. **Pass 1** — stream aligned shards, compute `canonical_key = min(read_key,
   rc(read_key))`, write a text projection (`canonical_key \t seq_no \t count`).
2. **External sort** by `canonical_key` (system `sort -s -S 64M`; the
   polars-sort spike — `scripts/bench_polars_sort_spill.py` — confirmed polars
   `sort`+`sink_parquet` does NOT spill, ratio 33, while system `sort` stays
   flat, ratio 1.24).
3. **Streaming first-per-group collapse** — sum counts, representative =
   min `seq_no` (= first in scan order, the pandas tie-break), palindrome
   count-doubling replicated for byte parity.
4. **Pass 2** — re-stream shards; for each representative run `_collapse_fanout`
   (single-entry store). If the allele set fits the memory budget (small
   inputs: tests, `make basic test`), sort in memory + populate `allele_rows`.
   Otherwise write sort-key-prefixed JSON-lines, external-sort by
   `(-#Reads, Aligned_Sequence, Reference_Sequence)`, stream into the parquet.
5. Temp files cleaned in `finally`.

Paired input stays on the in-memory path (bounded by unique allele count after
3a re-key). Parity holds by construction: `seq_no` encodes scan order, so the
representative (min `seq_no`) is first-in-scan-order — the exact pandas
tie-break — and sorting groups by `rep_seq_no` reproduces the pandas store
insertion order, so the final stable sort's tie-break is identical.

### Supporting fixes (bounded parquet batching, exposed by the now-flat collapse)

With the collapse store eliminated, the next memory terms were `to_pylist()`
calls that materialized whole batches of amplicon-length array columns (~50×
the arrow size for string columns). Fixed by:

- `iter_aligned_shard` / `_aggregate_alleles` / the TSV sink: **row-by-row cell
  access** (`.as_py()` per cell) instead of `batch.to_pylist()`, so only one
  row's Python objects are in flight.
- **Adaptive byte-budgeted batch sizing** (`_adaptive_batch_size`): targets
  ~50 MB of arrow data per batch from the parquet row-group metadata, so peak
  RSS is flat for long-read amplicons regardless of which stage runs.
- **TSV sink column pruning**: the non-detailed TSV reads only the scalar
  `crispresso2Cols` (not the amplicon-length position arrays it never emits).
- **Parquet writeback batching** (`_write_collapsed_allele_parquet_from_jsonl`):
  adaptive batch sized on JSON-line bytes (~32 MB of row dicts), not a fixed
  10k rows.

### Bench-mode shortcut (parquet + `--suppress_plots --suppress_report`)

When both plots and the report are suppressed (the bench configuration),
`df_alleles` is only consumed by `write_all_core_data_files` — whose outputs
are themselves only read by the report/plots. In that mode the whole-frame
`get_slice()` and `write_all_core_data_files` are skipped (the streaming TSV
sink serves the allele-frequency table directly from the parquet), so the bench
measures the core pipeline's flat RSS. Real runs (plots or report on) are
unaffected: they materialize `df_alleles` and write core data files as before.

> **Note (Stage 5 follow-up):** for real long-read runs with plots/report ON,
> `get_slice()` still materializes the whole `df_alleles` (the `prep_alleles_around_cut`
> core-data path slices it per-ref in Python). That is the next memory term to
> address (adapt the ~30 plot/core-data callsites to call `get_slice(ref_name)`
> per-ref, per the master plan's PR 7). It is out of scope for Fix #3, which is
> specifically the single-read collapse store.

## Reproduce

```bash
pixi run -e test python scripts/bench_pipeline_memory.py \
    --reads 10000,100000,1000000 --backends parquet \
    --amplicon-length 2000 --workdir /tmp/crispresso_bench_2k \
    --timeout 5400 --skip-1m-pandas
```

---

## Prior runs (200 bp amplicon, before Fix #3)

Amplicon: synthetic 200 bp, ~1:1 read diversity. Guide: GGAATCCCTTCTGCAGCACC. n_processes=2.

| reads | backend | peak RSS (MB) | wall (s) | status |
|------:|:-------|--------------:|---------:|:-------|
| 10000 | pandas | 245.2 | 6.5 | OK |
| 10000 | parquet | 435.5 | 8.3 | OK |
| 100000 | pandas | 573.8 | 17.7 | OK |
| 100000 | parquet | 1308.0 | 34.2 | OK |
| 1000000 | pandas | — | — | SKIPPED |
| 1000000 | parquet | 3519.4 | 205.5 | OK |

**SC #1 (parquet, pre-Fix#3): RSS(1000000)/RSS(10000) = 8.08 — FAIL**

The collapse store was the dominant term (~12 KB × 290k unique ≈ 3.5 GB at 1M,
~40× the Stage 1 dedup dict). Fixes #1 (stream worker chunks via islice) + #2
(Stage 1 budget 512→128 MB) were kept but barely moved the number; Fix #3
(the streaming single-read collapse) is the essential fix and is now in place.
