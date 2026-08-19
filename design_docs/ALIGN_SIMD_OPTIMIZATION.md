# CRISPResso2Align.global_align — micro-optimization + SIMD investigation

**Date:** 2026-07-23 (benchmarked + implemented)
**Branch / worktree:** `align-simd` → `../CRISPResso2.align-simd`
**Question:** `global_align` is the hottest kernel in the pipeline (every read is
aligned to every reference, forward + reverse, often thousands of times). Can
we speed it up with SIMD? Where should we actually spend effort?

**Status: small, safe win landed. ~5–10% faster, byte-identical output. The
honest finding is that *true* SIMD here is a large rewrite, and the data says
memory traffic — not ALU — is the bigger lever at long-read scale.**

## The kernel

`global_align` is Gotoh affine-gap Needleman–Wunsch. Per call it allocates
**six** `(max_i+1) × (max_j+1)` int32 matrices — M/I/J scores and M/I/J
pointers — and fills them in three near-identical blocks (inner cells, last
column, last row). Each cell computes, for the three matrices:

```
I[i,j] = max(M[i,j-1] + gap_open + gi[i],  I[i,j-1] + gap_extend + gi[i])
J[i,j] = max(M[i-1,j] + gap_open + gi[i-1], J[i-1,j] + gap_extend)
M[i,j] = max(M[i-1,j-1], I[i-1,j-1], J[i-1,j-1]) + matrix[ci,cj]
```

(`gi` = `gap_incentive`, the CRISPR cut-site gap preference.)

## Method

Built a reproducible harness in `scripts/` (all in the worktree):

| script | purpose |
|--------|---------|
| `bench_align.py` | per-call wall time (median/mean/p95/min) over N iterations at representative amplicon/read sizes, with a determinism self-check. Reports **ns/cell** to separate compute cost from alloc/page-fault cost. |
| `align_golden.py --save/--check` | correctness gate: a fixed 15-case battery (identical seqs, indels, mismatches, N bases, gap_incentive on/off, short+long, asymmetric) asserting **byte-identical** output to a baseline snapshot. |
| `compare_bench.py` | A/B diff of bench JSONs (averages per-size medians across multiple runs; reports the within-version spread as the noise floor and flags `real` vs `noise`). |

A/B protocol: 3 bench runs (30 iters each) per version, version swapped by
`git stash`-ing the `.pyx` and force-recompiling. Build: `clang -Ofast` on
Apple arm64. Correctness gate run on every rebuild: `align_golden --check` +
`pytest tests/unit_tests/test_CRISPResso2Align.py` (16 tests).

Reproduce:
```bash
pixi run -e test python setup.py build_ext --inplace
pixi run -e test python scripts/align_golden.py --save scripts/align_golden.json   # baseline
pixi run -e test python scripts/bench_align.py --iters 30 --json bench.json
# ...edit .pyx, rebuild...
pixi run -e test python scripts/align_golden.py --check scripts/align_golden.json
pixi run -e test python scripts/compare_bench.py --before base*.json --after opt*.json
```

## What changed (the small win)

Two zero-risk, readability-preserving edits to `CRISPResso2Align.pyx`, no
algorithm change:

1. **`@cython.wraparound(False)`** on `global_align`. Without it Cython emits a
   negative-index guard branch on every memoryview access; removing it lets the
   compiler drop those branches. (Deliberately *not* adding `cdivision(True)` —
   there is no integer division in the hot path, so it would be a misleading
   no-op.)
2. **Hoist `gap_incentive[i]` / `gap_incentive[i-1]` into scalars `gi_i`/`gi_im1`
   outside the inner `j` loop.** These are loop-invariant in `j` but live in a
   *separate* buffer, so the original inner loop re-read a different cache line
   twice per iteration. The hoist makes the inner loop touch only the score/
   pointer matrices (stride-1, hot). Applied to the inner block and the
   last-row block (both have loop-invariant `gi`); the last-column block varies
   `gi` with `i`, so it's left alone.

## Measured impact (A/B, 3 runs × 30 iters, median of medians)

| amp | read | before (ms) | after (ms) | delta | ns/cell b → a |
|----:|----:|-----------:|----------:|------:|--------------|
| 150 | 150 |   0.157 |   0.144 |  −8.2% | 6.98 → 6.41 |
| 200 | 200 |   0.298 |   0.270 |  −9.3% | 7.44 → 6.75 |
| 500 | 500 |   2.067 |   1.858 | −10.1% | 8.27 → 7.43 |
|1000 |1000 |   9.317 |   8.761 |  −6.0% | 9.32 → 8.76 |
|2000 |2000 |  34.448 |  32.775 |  −4.9% | 8.61 → 8.19 |
|5000 |5000 | 222.460 | 214.518 |  −3.6% | 8.90 → 8.58 |

**~5–10% faster, consistent in direction at every size**, biggest at the
mid sizes. Correctness: `align_golden --check` **PASS** (fingerprint
`43662296acce2dcf`, byte-identical) and **16/16** unit tests pass on the
optimized build.

## Why "SIMD" didn't mean SIMD here

The compiler (`-Ofast`) already tries to auto-vectorize but **can't**, and it's
not a compiler-quality problem — it's the algorithm. The I-matrix recurrence
`I[i,j] = … I[i,j-1] …` is a loop-carried dependency along `j`: every cell needs
its left neighbor. That serial chain is the textbook reason Gotoh SW/NW does not
vectorize row-by-row. The win above is from removing redundant memory reads and
guard branches, not from packed SIMD instructions.

## What the data says about where the real speedup is

`ns/cell` **rises with size**: ~6.4 at 150 bp (working set 0.5 MB, fits in L2)
to ~8.6–9.0 at ≥1000 bp (22–572 MB of matrices, blows the cache). That rise is
the signature of a **memory-bandwidth** bottleneck, not compute. Six int32
matrices = 24 bytes streamed per cell.

Implications for follow-on work:

- **At long-read scale the lever is memory traffic, not ALU.** The three
  *pointer* matrices are written every cell only to be consumed once during
  traceback. Options: pack pointers into a byte array (3 bytes→1), or drop them
  entirely and **recompute argmax during traceback** (3 matrices instead of 6 →
  roughly half the bandwidth). This targets the actual bottleneck and is lower
  risk than restructuring the recurrence.
- **True SIMD (anti-diagonal/wavefront, or Farrar striped)** is a substantial
  rewrite and would mainly help the *compute-bound* short-read regime (≤200 bp),
  where calls are already sub-millisecond. For CRISPR's typical 150–200 bp
  amplicons the absolute time saved per call is small; for long reads the
  bandwidth wall limits the SIMD payoff. The per-position `gap_incentive` also
  complicates the striped technique.
- **Buffer reuse** (the six DP matrices re-`np.empty`'d every call) was profiled
  separately (`scripts/bench_align_matrix_reuse.py` in the main worktree). It
  matters at the page-fault level for cold calls but the hot steady-state cost is
  the DP fill, which is what we improved here.

## Recommendation

1. **Land this change.** It's a free, safe, readability-positive ~5–10%.
2. **Next, if pursuing long-read throughput: reduce DP memory traffic** — start
   with a pointer-free traceback prototype and benchmark it (re-use this harness)
   before committing. Expected to help the ≥1 kb sizes where `ns/cell` is highest.
3. **Defer true SIMD** (anti-diagonal/striped) until/unless the short-read
   per-call cost is shown to dominate a real end-to-end run; it's a large effort
   with a size-dependent, bandwidth-capped payoff.

---

# Part 2: pointer-free traceback (landed)

**Date:** 2026-07-23
**Status: shipped.** Recommendation #2 above panned out far better than
expected: dropping the three pointer matrices and recomputing argmax during
traceback is **~2.5× faster** (not the modest halving predicted), byte-identical
outputs, and it surfaced + fixed a pre-existing crash bug. `global_align` (the
name every caller uses) is now the pointer-free implementation; the original is
retained as `global_align_pointers`.

## The change

`global_align_pointers` (the original) writes **six** int32 matrices per call —
M/I/J scores *and* M/I/J pointers — purely so the O(N+M) traceback can follow
stored pointers. The pointers are written every cell in the O(N*M) fill but
read once. The new `global_align` stores **only the three score matrices**;
traceback recomputes each cell's argmax from the scores on the fly. The
min_score sentinels on the border make the border behave correctly with no
special-casing; the only care is reproducing the fill's exact tie-breaking and
the position-dependent `gap_open`/`gap_extend` (inner cells use `gap_open`, the
last row/col use `gap_extend`).

The fill also got cleaner: with no pointers to store, each cell is a plain `max`
rather than argmax-with-side-effect.

## Correctness: differential property test (hypothesis)

The risk with rewriting traceback is subtle tie-breaking drift. We guard it with
a **property-based differential test** in `tests/unit_tests/test_CRISPResso2Align.py`
(`test_global_align_ptrfree_matches_pointers`): Hypothesis generates thousands
of (seqj, seqi, gap_incentive, gap_open, gap_extend) inputs across three scoring
matrices and asserts `global_align == global_align_pointers` byte-for-byte. It is
part of the normal unit-test suite (`pixi run -e test test`), caches failing
seeds, and `@example`-pins edge cases (length-1, indels, all-different,
formerly-crashing inputs). Hypothesis was added as a test dependency
(`pixi.toml`).

## Pre-existing crash bug found and fixed (min_score sentinel)

The property test immediately surfaced a **pre-existing** bug in the original
`global_align`: on very short sequences and/or large `gap_incentive` at the
border, `min_score = gap_open*max_i*max_j` is too weak a sentinel, so an optimal
path routes through border cells whose pointer matrices are uninitialized
(`np.empty`) — the traceback reads garbage and **raises or segfaults** (e.g.
`global_align("A","C",...)` crashed; ~10% of random length-1 inputs, ~0.8% of
length-2). Real CRISPR data (long reads/amplicons) never hit it because the
sentinel is hugely negative there.

Fix: use a fixed near-`-inf` sentinel (`min_score = -1000000000`) for the
"impossible" border cells. Verified non-regressive: the 15-case golden
(`scripts/align_golden.py`, fingerprint `43662296acce2dcf`) is unchanged and all
**757 unit tests pass**; the formerly-crashing inputs now return sensible
alignments (e.g. `"A"/"C"` mismatch instead of a crash).

## Measured impact (A/B, 3 runs × 30 iters, clang -Ofast, arm64)

`global_align_pointers` (before) vs `global_align` pointer-free (after):

| amp | ns/cell before → after | delta | bufs (MB) |
|----:|------------------------|------:|----------:|
| 150 | 5.43 → 2.71 | −50% | 0.5 vs 0.3 |
| 200 | 6.82 → 2.70 | −60% | 0.9 vs 0.5 |
| 500 | 7.70 → 2.92 | −62% | 5.7 vs 2.9 |
|1000 | 8.66 → 3.35 | −61% | 22.9 vs 11.4 |
|2000 | 8.16 → 3.15 | −61% | 91.6 vs 45.8 |
|5000 | 8.50 → 3.38 | −60% | 572 vs 286 |

**~2.5× faster at every size** (≈60% reduction), well above the noise floor.
The win is bigger than the predicted "half the traffic" because the pointer
writes were pure overhead in the hot fill (the score already determined the
pointer), so removing them cuts both memory traffic *and* instructions, and
shrinks the working set so more stays in cache. Peak RSS per call also roughly
halves (3 matrices instead of 6).

This lands the production speedup with **no changes to `CRISPRessoCORE.py`** —
every existing `CRISPResso2Align.global_align(...)` call site is now ~2.5× faster.

## Reproduce

```bash
pixi run -e test install
pixi run -e test python setup.py build_ext --inplace
pixi run -e test python -m pytest tests/unit_tests/test_CRISPResso2Align.py::test_global_align_ptrfree_matches_pointers -v
pixi run -e test python scripts/bench_align.py --func global_align         --json fast.json
pixi run -e test python scripts/bench_align.py --func global_align_pointers --json ref.json
pixi run -e test python scripts/compare_bench.py --before ref.json --after fast.json
```
