#!/usr/bin/env python3
"""Correctness gate for CRISPResso2Align.global_align.

Runs a fixed battery of probe inputs (identical seqs, indels, mismatches, N
bases, gap_incentive on/off, short + long, asymmetric lengths) and either:

    --save PATH  : write all (inputs, outputs) as JSON   [capture baseline]
    --check PATH : recompute and assert byte-identical to PATH [verify a rebuild]

This is the gate you run after editing the .pyx and rebuilding: a PASS means
the optimized kernel produces exactly the same alignments as the baseline,
cell for cell. Use it alongside `pytest tests/unit_tests/test_CRISPResso2Align.py`.

Determinism is global_align's invariant, so byte-equality is the right check
(no fuzzy score comparison).

USAGE
-----
    # 1. on baseline build:
    pixi run -e test python scripts/align_golden.py --save scripts/align_golden.json
    # 2. edit .pyx, rebuild, then:
    pixi run -e test python scripts/align_golden.py --check scripts/align_golden.json
"""

import argparse
import hashlib
import json
import os
import sys

import numpy as np

from CRISPResso2 import CRISPResso2Align

_PKG_DIR = os.path.dirname(CRISPResso2Align.__file__)
_DEFAULT_MATRIX = os.path.join(_PKG_DIR, "EDNAFULL")


def random_dna(length, seed):
    rng = np.random.default_rng(seed)
    bases = np.array(["A", "C", "G", "T"])
    return "".join(bases[rng.integers(0, 4, size=length)])


def random_dna_with_n(length, seed, n_frac=0.1):
    rng = np.random.default_rng(seed)
    bases = np.array(["A", "C", "G", "T", "N"])
    return "".join(bases[rng.integers(0, 5, size=length)])


def build_battery(matrix_path):
    """Return a list of (name, kwargs-for-global_align) probe cases.

    kwargs keys: seqj (read), seqi (ref), matrix, gap_incentive, gap_open,
    gap_extend. The reference is seqi, the read is seqj (matches global_align).
    """
    matrix = CRISPResso2Align.read_matrix(matrix_path)
    m = matrix  # alias for brevity
    gi0 = lambda n: np.zeros(n, dtype=np.int64)

    battery = [
        # --- anchored against the known-good unit-test expectations ---
        ("identical_4mer", dict(
            seqj="ATTA", seqi="ATTA", matrix=m, gap_incentive=gi0(5),
            gap_open=-1, gap_extend=-1)),
        ("insertion_incentive", dict(
            seqj="ATTTA", seqi="ATTA", matrix=m,
            gap_incentive=np.array([0, 1, 0, 0, 0], dtype=np.int64),
            gap_open=-1, gap_extend=-1)),
        ("insertion_incentive_pos2", dict(
            seqj="ATTTA", seqi="ATTA", matrix=m,
            gap_incentive=np.array([0, 0, 1, 0, 0], dtype=np.int64),
            gap_open=-1, gap_extend=-1)),
        # --- indels without incentive ---
        ("single_insertion", dict(
            seqj="ACGTACGT", seqi="ACGTACGT", matrix=m, gap_incentive=gi0(9),
            gap_open=-1, gap_extend=-1)),
        ("read_longer_than_ref", dict(
            seqj="ACGTAAAAACGT", seqi="ACGTACGT", matrix=m, gap_incentive=gi0(9),
            gap_open=-1, gap_extend=-1)),
        ("read_shorter_than_ref", dict(
            seqj="ACGTACGT", seqi="ACGTAAAAACGT", matrix=m, gap_incentive=gi0(13),
            gap_open=-1, gap_extend=-1)),
        # --- N base handling ---
        ("ref_has_N", dict(
            seqj="ACGTACGT", seqi="ACGNACGT", matrix=m, gap_incentive=gi0(9),
            gap_open=-1, gap_extend=-1)),
        ("read_has_N", dict(
            seqj="ACGNACGT", seqi="ACGTACGT", matrix=m, gap_incentive=gi0(9),
            gap_open=-1, gap_extend=-1)),
        # --- nonzero gap penalties that matter ---
        ("gap_open_extend_penalty", dict(
            seqj="GGGGACGTCCCC", seqi="GGGGTTTTCCCC", matrix=m,
            gap_incentive=gi0(13), gap_open=-20, gap_extend=-2)),
        # --- representative CRISPR-scale, deterministic ---
        ("square_200", dict(
            seqj=random_dna(200, 11), seqi=random_dna(200, 11), matrix=m,
            gap_incentive=gi0(201), gap_open=-1, gap_extend=-1)),
        ("square_200_incentive_cut", dict(
            seqj=random_dna(200, 12), seqi=random_dna(200, 12), matrix=m,
            gap_incentive=_incentive_around(201, cut=100), gap_open=-1,
            gap_extend=-1)),
        ("square_1000", dict(
            seqj=random_dna(1000, 13), seqi=random_dna(1000, 13), matrix=m,
            gap_incentive=gi0(1001), gap_open=-1, gap_extend=-1)),
        ("rect_read_long_300x500", dict(
            seqj=random_dna(500, 14), seqi=random_dna(300, 14), matrix=m,
            gap_incentive=gi0(301), gap_open=-1, gap_extend=-1)),
        ("rect_read_short_500x300", dict(
            seqj=random_dna(300, 15), seqi=random_dna(500, 15), matrix=m,
            gap_incentive=gi0(501), gap_open=-1, gap_extend=-1)),
        ("with_N_200", dict(
            seqj=random_dna_with_n(200, 16), seqi=random_dna_with_n(200, 17),
            matrix=m, gap_incentive=gi0(201), gap_open=-1, gap_extend=-1)),
    ]
    return battery


def _incentive_around(n, cut, width=5, val=10):
    """gap_incentive that prefers gaps near a cut site (CRISPR cut)."""
    gi = np.zeros(n, dtype=np.int64)
    lo = max(0, cut - width)
    hi = min(n, cut + width + 1)
    gi[lo:hi] = val
    return gi


def run_battery():
    results = []
    for name, kwargs in build_battery(_DEFAULT_MATRIX):
        # seqj (read) and seqi (ref) are passed POSITIONALLY: global_align's
        # params are literally named pystr_seqj/pystr_seqi, and CORE.py calls it
        # the same way. matrix/gap_* go by keyword.
        a, b, score = CRISPResso2Align.global_align(
            kwargs["seqj"], kwargs["seqi"],
            matrix=kwargs["matrix"], gap_incentive=kwargs["gap_incentive"],
            gap_open=kwargs["gap_open"], gap_extend=kwargs["gap_extend"])
        results.append({
            "name": name,
            "inputs": {
                "seqj": kwargs["seqj"], "seqi": kwargs["seqi"],
                "gap_open": kwargs["gap_open"], "gap_extend": kwargs["gap_extend"],
                "gap_incentive": [int(x) for x in kwargs["gap_incentive"]],
            },
            "output": {"align_j": a, "align_i": b, "score": float(score)},
        })
    return results


def fingerprint(results):
    """Stable short hash over all outputs (ignores inputs) for a quick summary."""
    h = hashlib.sha256()
    for r in results:
        o = r["output"]
        h.update(r["name"].encode())
        h.update(b"\x00")
        h.update(o["align_j"].encode())
        h.update(b"\x01")
        h.update(o["align_i"].encode())
        h.update(b"\x02")
        h.update(repr(o["score"]).encode())
        h.update(b"\x03")
    return h.hexdigest()[:16]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--save", metavar="PATH", help="write baseline outputs as JSON")
    g.add_argument("--check", metavar="PATH", help="assert outputs match saved JSON")
    args = ap.parse_args(argv)

    results = run_battery()
    fp = fingerprint(results)

    if args.save:
        with open(args.save, "w") as fh:
            json.dump({"fingerprint": fp, "cases": results}, fh, indent=2)
        print(f"saved {len(results)} cases  fingerprint={fp}  -> {args.save}")
        return 0

    # --check
    with open(args.check) as fh:
        saved = json.load(fh)
    saved_fp = saved["fingerprint"]
    by_name = {c["name"]: c for c in saved["cases"]}

    mismatches = []
    for r in results:
        exp = by_name.get(r["name"])
        if exp is None:
            mismatches.append((r["name"], "case missing from saved file"))
            continue
        eo = exp["output"]
        ro = r["output"]
        if eo["align_j"] != ro["align_j"] or eo["align_i"] != ro["align_i"] \
                or abs(eo["score"] - ro["score"]) > 1e-9:
            mismatches.append((r["name"], "output differs"))

    print(f"checked {len(results)} cases  saved_fp={saved_fp}  now_fp={fp}")
    if mismatches:
        print(f"FAIL: {len(mismatches)} case(s) differ:")
        for name, why in mismatches:
            print(f"  - {name}: {why}")
        return 1
    if saved_fp != fp:
        print("FAIL: all cases match individually but aggregate fingerprint differs "
              "(should not happen).")
        return 1
    print("PASS: byte-identical to baseline on all cases.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
