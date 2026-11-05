"""Per-input worker for ab_vs_master.py.

Reads ONE JSON request from stdin, runs global_align from whichever
CRISPResso2 build it was launched with, and prints ONE JSON result. Run in a
fresh subprocess per input so a segfault in the master build (the pre-fix
weak-sentinel crash) is isolated to a single input and cannot abort the run.

A Python-level exception is reported as {"ok": false, "error": ...}; a native
crash (segfault) produces NO stdout and a non-zero/signal exit code, which the
driver treats as "crashed".

Request:  {"seqj","seqi","gap_incentive","gap_open","gap_extend","matrix_path"}
Result:   {"ok": true, "align_j","align_i","score"} | {"ok": false, "error"}
"""
import json
import sys


def main():
    req = json.loads(sys.stdin.read())
    import numpy as np
    from CRISPResso2 import CRISPResso2Align

    matrix = CRISPResso2Align.read_matrix(req["matrix_path"])
    gi = np.array(req["gap_incentive"], dtype=np.int64)
    # The master build prints stray traceback diagnostics to stdout (its
    # `else: print('i: ...')` branch) right before raising. Redirect those to
    # stderr so the ONE result line on stdout is clean JSON the driver can parse.
    real_stdout = sys.stdout
    sys.stdout = sys.stderr
    try:
        try:
            aj, ai, score = CRISPResso2Align.global_align(
                req["seqj"], req["seqi"], matrix=matrix, gap_incentive=gi,
                gap_open=req["gap_open"], gap_extend=req["gap_extend"])
            result = {"ok": True, "align_j": aj, "align_i": ai, "score": score}
        except Exception as e:  # any Python-level error (not a native segfault)
            result = {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        sys.stdout = real_stdout
    print(json.dumps(result))


if __name__ == "__main__":
    main()
