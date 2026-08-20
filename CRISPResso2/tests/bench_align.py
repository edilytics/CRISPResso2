#!/usr/bin/env python3
r"""Benchmark CRISPResso2Align.global_align per-call wall time at representative
amplicon/read sizes.

No Cython recompile required — measures the currently-built extension.

WHY THIS EXISTS
---------------
global_align is the hot inner kernel of the whole pipeline (every read is
aligned to every reference, forward + reverse). Before changing the .pyx we
need a reproducible number to beat and a way to tell compute cost apart from
allocation/page-fault cost. This script reports:

    T/call  : median / mean / p95 / min over N timed iterations
    ns/cell : median / (max_i * max_j)  -> DP cells per second
             (stays ~constant when compute-bound; drops when alloc/fault bound)

It also self-checks that every iteration returns a byte-identical result
(global_align is deterministic) so a flaky timing run is never mistaken for
a correctness change.

Correctness across rebuilds is handled separately by CRISPResso2/tests/align_golden.py.

USAGE
-----
    pixi run -e test python CRISPResso2/tests/bench_align.py \
        --sizes 150,200,500,1000,2000,5000 --iters 20 --json baseline.json
"""

import argparse
import gc
import json
import os
import time

import numpy as np

from CRISPResso2 import CRISPResso2Align

_PKG_DIR = os.path.dirname(CRISPResso2Align.__file__)
_DEFAULT_MATRIX = os.path.join(_PKG_DIR, "EDNAFULL")


def random_dna(length, rng):
    bases = np.array(["A", "C", "G", "T"])
    idx = rng.integers(0, 4, size=length)
    return "".join(bases[idx])


def time_calls(seqj, seqi, matrix, gi, gap_open, gap_extend, n_iters, align_fn):
    """Return (list_of_per_iter_times, all_identical).

    Two untimed warmups prime the page cache and Python dispatch path so the
    measured iterations reflect steady-state per-call cost (the regime that
    matters when the pipeline aligns thousands of reads in a row).
    """
    call = lambda: align_fn(
        seqj, seqi, matrix=matrix, gap_incentive=gi,
        gap_open=gap_open, gap_extend=gap_extend,
    )
    call()  # warmup 1: heat DP pages
    call()  # warmup 2: stabilize dispatch

    # Time each call individually. global_align is deterministic, so every
    # iteration must return a byte-identical (align_j, align_i, score); we flag
    # the run BAD if not (would indicate a real correctness problem, not noise).
    gc.disable()
    t0 = time.perf_counter()
    first = call()
    times = [time.perf_counter() - t0]
    all_identical = True
    for _ in range(n_iters - 1):
        t0 = time.perf_counter()
        res = call()
        times.append(time.perf_counter() - t0)
        if (res[0], res[1], res[2]) != (first[0], first[1], first[2]):
            all_identical = False
    gc.enable()
    return times, all_identical


def fmt_ms(t):
    return f"{t * 1000:.3f} ms"


def run(sizes, rect_ratio, iters, gap_open, gap_extend, matrix_path, seed, json_path,
        func_name="global_align"):
    rng = np.random.default_rng(seed)
    matrix = CRISPResso2Align.read_matrix(matrix_path)
    align_fn = getattr(CRISPResso2Align, func_name)

    print(f"function    : {func_name}")
    print(f"matrix      : {matrix_path}")
    print(f"gap_open={gap_open} gap_extend={gap_extend} iters={iters} seed={seed}")
    if rect_ratio != 1.0:
        print(f"rectangular: read_len = amp_len * {rect_ratio}")
    print()

    hdr = (
        f"{'read_len':>8} {'amp_len':>7} {'cells(M)':>9} {'bufs(MB)':>9} "
        f"{'median':>11} {'mean':>11} {'p95':>11} {'min':>11} {'ns/cell':>9} {'det':>4}"
    )
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for amp_len in sizes:
        read_len = max(1, round(amp_len * rect_ratio))
        max_i, max_j = amp_len, read_len
        seqi = random_dna(max_i, rng)
        seqj = random_dna(max_j, rng)
        gi = np.zeros(max_i + 1, dtype=np.int64)

        cells = max_i * max_j
        # global_align (pointer-free) stores 3 int32 score matrices;
        # global_align_pointers (reference) stores 6 (3 scores + 3 pointers).
        n_matrices = 3 if func_name == "global_align" else 6
        bufs_mb = n_matrices * (max_i + 1) * (max_j + 1) * 4 / (1024 * 1024)

        times, ok = time_calls(
            seqj, seqi, matrix, gi, gap_open, gap_extend, iters, align_fn,
        )
        med = float(np.median(times))
        mean = float(np.mean(times))
        p95 = float(np.percentile(times, 95))
        mn = float(np.min(times))
        ns_per_cell = med * 1e9 / cells if cells else float("nan")

        print(
            f"{read_len:>8} {amp_len:>7} {cells / 1e6:>9.2f} {bufs_mb:>9.1f} "
            f"{fmt_ms(med):>11} {fmt_ms(mean):>11} {fmt_ms(p95):>11} {fmt_ms(mn):>11} "
            f"{ns_per_cell:>9.2f} {'ok' if ok else 'BAD':>4}"
        )
        rows.append({
            "read_len": read_len, "amp_len": amp_len,
            "cells_M": cells / 1e6, "bufs_MB": bufs_mb,
            "median_ms": med * 1000, "mean_ms": mean * 1000,
            "p95_ms": p95 * 1000, "min_ms": mn * 1000,
            "ns_per_cell": ns_per_cell, "deterministic": ok,
        })
        del seqi, seqj, gi
        gc.collect()

    print()
    print("Reading the table:")
    print("  ns/cell ~constant across sizes => compute-bound (DP recurrence dominates)")
    print("  ns/cell rising at small sizes  => alloc/page-fault overhead per call")
    print("  det=BAD => non-deterministic output; results not comparable")

    if json_path:
        out = {
            "function": func_name,
            "matrix": matrix_path, "gap_open": gap_open, "gap_extend": gap_extend,
            "iters": iters, "seed": seed, "rect_ratio": rect_ratio,
            "rows": rows,
        }
        with open(json_path, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nwrote {json_path}")


def parse_sizes(s):
    return [int(p.strip()) for p in s.split(",") if p.strip()]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", default="150,200,500,1000,2000,5000", type=parse_sizes,
                    help="amplicon lengths (bp) [default 150,200,500,1000,2000,5000]")
    ap.add_argument("--rect", dest="rect_ratio", type=float, default=1.0,
                    help="read_len / amp_len [default 1.0 = square]")
    ap.add_argument("--iters", type=int, default=20)
    ap.add_argument("--gap-open", type=int, default=-1)
    ap.add_argument("--gap-extend", type=int, default=-1)
    ap.add_argument("--matrix", default=_DEFAULT_MATRIX)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--func", default="global_align",
                    choices=["global_align", "global_align_pointers"],
                    help="which alignment kernel to benchmark")
    ap.add_argument("--json", default=None, help="dump results as JSON to this path")
    args = ap.parse_args(argv)
    run(args.sizes, args.rect_ratio, args.iters,
        args.gap_open, args.gap_extend, args.matrix, args.seed, args.json, args.func)


if __name__ == "__main__":
    main()
