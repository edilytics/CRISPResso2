# CRISPResso2 storage-backend peak-RSS benchmark

Amplicon: synthetic 200 bp, ~1:1 read diversity (1-3 random subs per read). Guide: GGAATCCCTTCTGCAGCACC. n_processes=2.

| reads | backend | peak RSS (MB) | wall (s) | status |
|------:|:-------|--------------:|---------:|:-------|
| 10000 | pandas | 245.2 | 6.5 | OK |
| 10000 | parquet | 435.5 | 8.3 | OK |
| 100000 | pandas | 573.8 | 17.7 | OK |
| 100000 | parquet | 1308.0 | 34.2 | OK |
| 1000000 | pandas | — | — | SKIPPED |
| 1000000 | parquet | 3519.4 | 205.5 | OK |

**SC #1 (parquet): RSS(1000000)/RSS(10000) = 8.08 — FAIL**

## Re-run after fixes #1 (stream chunks via islice) + #2 (budget 512→128 MB)

| reads | parquet RSS before | parquet RSS after | delta |
|------:|-------------------:|------------------:|------:|
| 10000 | 454 MB | 436 MB | -18 MB |
| 100000 | 1360 MB | 1308 MB | -52 MB |
| 1000000 | 3630 MB | 3519 MB | -111 MB |

SC #1 still FAILS: RSS(1M)/RSS(10k) = 8.08. Fixes #1+#2 are correct
improvements (less parent memory during worker chunking; spill engages for
genuinely large high-cardinality inputs) but barely move the end-to-end
number. Root cause:

- **Stage 1 did NOT spill** on this benchmark. The 200 bp amplicon yields far
  fewer unique reads than input reads (collisions): 10k→5.2k, 100k→44k,
  1M→290k unique. The dedup dict for 290k × 200 bp keys is ~81 MB, under the
  128 MB budget, so it stays eager. Fix #2 engages for the real long-read
  target (e.g. 5 kb × 100k unique ≈ 500 MB > 128 MB → spills) but the 200 bp
  benchmark under-exercises it.
- **The collapse store is the dominant term.** `_collapse_rekey_and_rcmerge`
  builds an in-memory `store` dict with one entry per unique *aligned* read,
  each holding the full payload (~10 Python lists of ~200 ints ≈ 12-16 KB/read).
  At 1M input (290k unique) that is ~3.5 GB — ~40× the Stage 1 dict. This is
  copy #3, and it is the one that matters.

**Conclusion: fix #3 (external-sort canonical-key merge for single-read
collapse, replacing the in-memory `store`) is the essential fix for SC #1.**
Fixes #1+#2 are kept (they are correct and help the long-read spill case) but
are insufficient on their own. A follow-up benchmark with a longer amplicon
(2-5 kb, to exercise the Stage 1 spill) and after fix #3 is needed to
re-validate SC #1.
