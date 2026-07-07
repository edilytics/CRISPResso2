# CRISPResso2 storage-backend peak-RSS benchmark

Amplicon: synthetic 200 bp, ~1:1 read diversity (1-3 random subs per read). Guide: GGAATCCCTTCTGCAGCACC. n_processes=2.

| reads | backend | peak RSS (MB) | wall (s) | status |
|------:|:-------|--------------:|---------:|:-------|
| 10000 | pandas | 246.5 | 6.4 | OK |
| 10000 | parquet | 453.5 | 8.6 | OK |
| 100000 | pandas | 567.8 | 17.4 | OK |
| 100000 | parquet | 1360.0 | 33.5 | OK |
| 1000000 | pandas | — | — | SKIPPED |
| 1000000 | parquet | 3629.7 | 206.6 | OK |

**SC #1 (parquet): RSS(1000000)/RSS(10000) = 8.00 — FAIL**
