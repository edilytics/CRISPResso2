# Streaming Memory Benchmark (PR 2 / Spike S1)

**Date:** 2026-07-03 10:17:49
**Matrix:** rows=[1000, 10000, 100000] x read_lengths=[1000, 2000]bp x engines=['streaming', 'eager', 'duckdb']
**Unique fraction:** 1.0 (1.0 = worst-case all-unique)
**Gate:** peak RSS(max rows) / RSS(min rows) < 3.0 for the streaming engine

## Verdict

**FAIL** — streaming len=1000b: RSS 1000->100000 = 80->526 MB (ratio 6.58) [FAIL] | streaming len=2000b: RSS 1000->100000 = 86->898 MB (ratio 10.43) [FAIL] | duckdb len=1000b: RSS 1000->100000 = 134->446 MB (ratio 3.33) [FAIL] | duckdb len=2000b: RSS 1000->100000 = 138->766 MB (ratio 5.53) [FAIL]

## Results

```
engine          rows   len(b)  peakRSS(MB)  wall(s)    groups  input(MB)   ok
-----------------------------------------------------------------------------
streaming       1000     1000         79.9     0.08      1000        0.3    Y
streaming      10000     1000        137.4     0.12     10000        3.0    Y
streaming     100000     1000        525.8      0.6    100000       29.8    Y
streaming       1000     2000         86.1     0.08      1000        0.6    Y
streaming      10000     2000        196.6     0.19     10000        5.9    Y
streaming     100000     2000        897.7     3.15    100000       59.4    Y
eager           1000     1000         82.4     0.13      1000        0.3    Y
eager          10000     1000        104.0      0.2     10000        3.0    Y
eager         100000     1000        456.5     1.14    100000       29.8    Y
eager           1000     2000         88.1     0.16      1000        0.6    Y
eager          10000     2000        132.7     0.28     10000        5.9    Y
eager         100000     2000        774.4      2.0    100000       59.4    Y
duckdb          1000     1000        134.0     1.02      1000        0.3    Y
duckdb         10000     1000        172.7     0.85     10000        3.0    Y
duckdb        100000     1000        445.7     1.41    100000       29.8    Y
duckdb          1000     2000        138.5     0.82      1000        0.6    Y
duckdb         10000     2000        215.0     0.82     10000        5.9    Y
duckdb        100000     2000        766.1     2.51    100000       59.4    Y
```

## Interpretation

- The **streaming** engine must keep peak RSS roughly flat as row count grows (that is the whole point of the parquet backend — flat memory across arbitrary input sizes).
- The **eager** engine is the regression oracle: its RSS should scale ~linearly with row count, demonstrating that the test actually exercises the dedup at scale.
- If streaming RSS tracks eager RSS (linear growth), polars is falling back to eager on multi-kb string columns. Per the design, the fallback is DuckDB for the RC-dedup stage (Stage 3b) only.

## How to reproduce

```bash
pixi run -e default python scripts/bench_streaming.py            # full matrix
pixi run -e default python scripts/bench_streaming.py --quick    # fast smoke
```
