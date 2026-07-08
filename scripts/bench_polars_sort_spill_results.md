# polars sort spill spike

**Date:** 2026-07-07 15:49:31
**Matrix:** rows=[10000, 100000, 1000000] x key_lengths=[2000]bp x engines=['polars_sort', 'system_sort']
**Key cardinality:** all-unique (worst case — no sort collapse)
**Gate:** peak RSS(max rows) / RSS(min rows) < 3.0 for a given key length

Context: `design_docs/STREAMING_SINGLE_READ_COLLAPSE.md` step 3 — the one
open implementation question. polars `sort`+`sink_parquet` is preferred
(cleanest); system `sort` (external merge-sort, `-S 64M`) is the proven
fallback (spike S1b).

## Results

```
engine              rows    len peakRSS(MB)  self(MB)  child(MB)  wall(s)     n_out   ok
-----------------------------------------------------------------------------------------------
polars_sort        10000   2000       132.5     132.5        0.0     0.24     10000    Y
system_sort        10000   2000       120.2     120.2       41.9     0.47     10000    Y
polars_sort       100000   2000       668.2     668.2        0.0     1.11    100000    Y
system_sort       100000   2000       331.1     331.1       79.9     2.34    100000    Y
polars_sort      1000000   2000      4367.4    4367.4        0.0     3.58   1000000    Y
system_sort      1000000   2000       409.8     409.8      115.7    24.99   1000000    Y

## RSS scaling (max_rows / min_rows per key length)

  polars_sort    len=2000b: 10000->1000000 rows = 132->4367 MB (ratio 32.96) [SCALES]
  system_sort    len=2000b: 10000->1000000 rows = 120->410 MB (ratio 3.41) [SCALES]

## Reproduce

```bash
pixi run -e default python scripts/bench_polars_sort_spill.py
pixi run -e default python scripts/bench_polars_sort_spill.py --quick
```
