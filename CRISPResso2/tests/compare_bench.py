#!/usr/bin/env python3
"""Compare bench_align.py JSON outputs across runs/versions.

Given one or more --before and --after JSON files (e.g. 3 bench runs each to
average out OS-scheduling noise), print per-size median-of-medians and the
delta, with the per-run spread so the noise floor is visible.

A delta that is larger than the within-version spread is real; one that is
inside the spread is noise.
"""

import argparse
import json
import statistics
import sys


def load(paths):
    """Return {amp_len: (read_len, [median_ms across provided files])}."""
    by_size = {}
    for p in paths:
        data = json.load(open(p))
        for row in data["rows"]:
            key = row["amp_len"]
            by_size.setdefault(key, {"read_len": row["read_len"], "meds": [],
                                     "ns": []})
            by_size[key]["meds"].append(row["median_ms"])
            by_size[key]["ns"].append(row["ns_per_cell"])
    return by_size


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--before", nargs="+", required=True,
                    help="baseline JSON file(s) (glob expanded by shell)")
    ap.add_argument("--after", nargs="+", required=True,
                    help="optimized JSON file(s)")
    args = ap.parse_args(argv)

    before = load(args.before)
    after = load(args.after)

    hdr = (f"{'amp':>5} {'read':>5} "
           f"{'before(ms)':>11} {'b spread':>9} "
           f"{'after(ms)':>10} {'a spread':>9} "
           f"{'delta':>8} {'ns/cell b':>9} {'ns/cell a':>9} {'verdict':>8}")
    print(hdr)
    print("-" * len(hdr))

    for amp in sorted(before.keys() & after.keys()):
        b = before[amp]
        a = after[amp]
        b_med = statistics.median(b["meds"])
        a_med = statistics.median(a["meds"])
        # spread = max - min of per-run medians (ms), i.e. the noise floor
        b_spread = (max(b["meds"]) - min(b["meds"])) if len(b["meds"]) > 1 else 0.0
        a_spread = (max(a["meds"]) - min(a["meds"])) if len(a["meds"]) > 1 else 0.0
        delta_pct = 100 * (a_med - b_med) / b_med
        b_ns = statistics.median(b["ns"])
        a_ns = statistics.median(a["ns"])
        # "real" if the speedup exceeds the larger within-version spread
        verdict = "noise" if abs(a_med - b_med) <= max(b_spread, a_spread) else "real"

        print(f"{amp:>5} {b['read_len']:>5} "
              f"{b_med:>10.3f} {b_spread:>8.3f} "
              f"{a_med:>9.3f} {a_spread:>8.3f} "
              f"{delta_pct:>+7.1f}% {b_ns:>9.2f} {a_ns:>9.2f} {verdict:>8}")
    print()
    print("delta: negative = faster (after < before). "
          "spread = max-min of per-run medians (noise floor). "
          "verdict = 'noise' if |delta| <= spread.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
