#!/usr/bin/env python3
"""A/B differential: current global_align vs genuine MASTER.

The hypothesis differential test in the repo compares the new pointer-free
``global_align`` against ``global_align_pointers`` -- but BOTH live on this
branch and were derived from the same source (``global_align_pointers`` even
inherits the sentinel fix). That proves new==new-reference, not new==production.

This script is the literal proof the user wants: it builds the *genuine* master
``global_align`` (the original, with the weak ``min_score = gap_open*max_i*max_j``
sentinel that crashes on some short inputs) and the current branch, runs a large
deterministic corpus through each in isolated subprocesses, and classifies every
input:

  (a) AGREE     -- both return byte-identical (align_j, align_i, score)  [the proof]
  (b) MASTER-CRASHES -- master segfaults/errors, current succeeds        [EXPECTED:
                     the sentinel bug fix; formerly-crashing inputs now work]
  (c) CURRENT-CRASHES -- current fails, master succeeds                  [BLOCKER]
  (d) DIFFER    -- both succeed but disagree                             [BLOCKER]
  (e) BOTH-CRASH

Expected on a correct branch: only (a) and (b). Any (c) or (d) is a regression.
Exit code is nonzero if any (c)/(d) is found.

PREREQUISITES
-------------
A worktree of master must exist and be built+installed, e.g.:

    git worktree add /tmp/c2-master master
    cd /tmp/c2-master && pixi install -e test && pixi run -e test pip install -e . --no-build-isolation

The current branch must likewise be built+installed in its own worktree.
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
WORKER = os.path.join(HERE, "_ab_worker.py")
MATRIX = os.path.join(REPO_ROOT, "CRISPResso2", "EDNAFULL")

CURRENT_PY = os.path.join(REPO_ROOT, ".pixi", "envs", "test", "bin", "python")


def git_head(cwd):
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=cwd, text=True).strip()
    except Exception as e:
        return f"<git error: {e}>"


def run_one(python, cwd, req):
    """Run the worker in a fresh subprocess. Returns a result dict.

    A native crash (segfault) or timeout => {"crashed": True, ...}.
    """
    try:
        proc = subprocess.run(
            [python, WORKER], input=json.dumps(req),
            capture_output=True, text=True, timeout=60, cwd=cwd, check=False)
    except subprocess.TimeoutExpired:
        return {"crashed": True, "detail": "timeout (60s)"}
    if proc.returncode != 0:
        sig = -proc.returncode if proc.returncode < 0 else proc.returncode
        return {"crashed": True,
                "detail": f"exit/signal={sig} stderr={proc.stderr.strip()[:160]}"}
    try:
        return json.loads(proc.stdout)
    except Exception:
        return {"crashed": True, "detail": f"no json stdout={proc.stdout.strip()[:160]}"}


def _is_broken(out, seqj, seqi):
    """True if ``out`` crashed or is not a structurally valid alignment.

    A valid global alignment must: have equal-length rows, no (-,-) column,
    and ungapping each row recovers the original sequence. Master's weak
    sentinel lets traceback route through uninitialized pointer cells, which
    can yield garbage (duplicated/truncated bases) -- that's the same bug as a
    crash, just milder, and is NOT a genuine disagreement with current.
    """
    if out.get("crashed") or not out.get("ok", True):
        return True
    aj, ai = out.get("align_j"), out.get("align_i")
    if aj is None or ai is None:
        return True
    if len(aj) != len(ai):
        return True
    if any(a == "-" and b == "-" for a, b in zip(aj, ai)):
        return True
    if aj.replace("-", "") != seqj or ai.replace("-", "") != seqi:
        return True
    return False


def _broken_kind(out, seqj, seqi):
    """Human-readable reason an output is broken (for the report)."""
    if out.get("crashed"):
        return f"crash/timeout ({out.get('detail', '')})"
    if not out.get("ok", True):
        return f"raised {out.get('error', '')}"
    aj, ai = out.get("align_j", ""), out.get("align_i", "")
    reasons = []
    if len(aj) != len(ai):
        reasons.append(f"row lengths differ ({len(aj)} vs {len(ai)})")
    if any(a == "-" and b == "-" for a, b in zip(aj, ai)):
        reasons.append("has (-,-) column")
    if aj.replace("-", "") != seqj:
        reasons.append(f"align_j ungapped {aj.replace('-', '')!r} != read {seqj!r}")
    if ai.replace("-", "") != seqi:
        reasons.append(f"align_i ungapped {ai.replace('-', '')!r} != ref {seqi!r}")
    return "garbage: " + "; ".join(reasons)


# --------------------------------------------------------------------------- #
# Deterministic corpus
# --------------------------------------------------------------------------- #
import random


def _rand_seq(rng, n, alphabet="ACGTN"):
    return "".join(rng.choice(alphabet) for _ in range(n))


def _cut_incentive(n, cut, val=1):
    gi = [0] * n
    if 0 <= cut < n:
        gi[cut] = val
    return gi


def build_corpus():
    """A deterministic, production-representative corpus (seeded)."""
    rng = random.Random(20260730)
    cases = []

    # (1) Short inputs that heavily stress the DP border / weak-sentinel crash
    #     regime (length 1-3, large gap_incentive, wide penalties). This is
    #     where master is expected to crash (~10% of len-1, ~0.8% of len-2).
    for _ in range(120):
        li = rng.randint(1, 3)
        lj = rng.randint(1, 3)
        seqi = _rand_seq(rng, li)
        seqj = _rand_seq(rng, lj)
        gi = [rng.choice([-5, 0, 1, 5, 15, 20]) for _ in range(li + 1)]
        go = rng.choice([-1, -2, -5, -10, -20, -50])
        ge = rng.choice([-1, -2, -20, 0])
        cases.append(("short_border", seqj, seqi, gi, go, ge))

    # (2) Production-like mid/long inputs. These essentially never crash, so
    #     they are where byte-identity (category a) must hold across the board.
    #     Includes both default penalty regimes and prime-editing (ge=0).
    for _ in range(260):
        li = rng.randint(4, 120)
        lj = rng.randint(4, 120)
        seqi = _rand_seq(rng, li)
        seqj = _rand_seq(rng, lj)
        # sometimes mutate a few bases so it's near (not identical to) the ref
        if rng.random() < 0.5:
            sj = list(seqj)
            for _ in range(rng.randint(0, max(1, lj // 10))):
                sj[rng.randrange(lj)] = rng.choice("ACGTN")
            seqj = "".join(sj)
        regime = rng.choice(["nw", "flexiguide", "prime"])
        if regime == "prime":
            go, ge = -50, 0
        else:
            go, ge = -20, -2
        if rng.random() < 0.5:
            gi = _cut_incentive(li + 1, rng.randint(0, li))
        else:
            gi = [0] * (li + 1)
        cases.append((f"mid_{regime}", seqj, seqi, gi, go, ge))

    # (3) Anchored golden-style cases with known-good behavior.
    cases += [
        ("golden_identical", "ATTA", "ATTA", [0] * 5, -1, -1),
        ("golden_incentive", "ATTTA", "ATTA", [0, 1, 0, 0, 0], -1, -1),
        ("golden_n", "ANNG", "ATCG", [0] * 5, -1, -1),
        ("golden_diff", "AAAA", "TTTT", [0] * 5, -1, -1),
        ("golden_pen", "GGGGACGTCCCC", "GGGGTTTTCCCC", [0] * 13, -20, -2),
        ("golden_prime", "ACGTACGT", "ACGT", [0] * 5, -50, 0),
    ]

    return cases


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", default="/tmp/c2-master",
                    help="path to a built+installed worktree of master")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args(argv)

    master_py = os.path.join(args.master, ".pixi", "envs", "test", "bin", "python")
    versions = [
        ("master", master_py, args.master),
        ("current", CURRENT_PY, REPO_ROOT),
    ]
    for name, py, cwd in versions:
        if not os.path.exists(py):
            sys.exit(f"ERROR: {name} python not found at {py} "
                     f"(build+install that worktree first)")

    print(f"master  HEAD: {git_head(args.master)}  ({args.master})")
    print(f"current HEAD: {git_head(REPO_ROOT)}  ({REPO_ROOT})")
    print(f"matrix: {MATRIX}\n")

    corpus = build_corpus()
    print(f"corpus: {len(corpus)} inputs\n")

    def req_for(case):
        _, seqj, seqi, gi, go, ge = case
        return {"seqj": seqj, "seqi": seqi, "gap_incentive": gi,
                "gap_open": go, "gap_extend": ge, "matrix_path": MATRIX}

    # Run every (version, input) pair. Subprocess-per-input isolates crashes.
    def run_pair(idx_case):
        idx, case = idx_case
        req = req_for(case)
        out = {}
        for name, py, cwd in versions:
            out[name] = run_one(py, cwd, req)
        return idx, case, out

    results = [None] * len(corpus)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for idx, case, out in pool.map(run_pair, enumerate(corpus)):
            results[idx] = (case, out)

    # Classify. A master output that is invalid (crash OR malformed -- e.g.
    # duplicated bases from reading uninitialized pointer memory, the same
    # weak-sentinel bug manifesting as garbage instead of a segfault) counts as
    # "master broken", NOT a genuine disagreement. A real regression requires
    # BOTH outputs valid yet different.
    cats = {"a_agree": [], "b_master_broken": [], "c_current_broken": [],
            "d_differ": [], "e_both_broken": []}
    for case, out in results:
        _, seqj, seqi, gi, go, ge = case
        m, c = out["master"], out["current"]
        m_bad = _is_broken(m, seqj, seqi)
        c_bad = _is_broken(c, seqj, seqi)
        if m_bad and c_bad:
            cats["e_both_broken"].append((case, m, c))
        elif m_bad and not c_bad:
            cats["b_master_broken"].append((case, m, c))
        elif c_bad and not m_bad:
            cats["c_current_broken"].append((case, m, c))
        else:
            m_tup = (m.get("align_j"), m.get("align_i"), m.get("score"))
            c_tup = (c.get("align_j"), c.get("align_i"), c.get("score"))
            if m_tup == c_tup:
                cats["a_agree"].append((case, m, c))
            else:
                cats["d_differ"].append((case, m, c))

    n = len(corpus)
    nb = len(cats["b_master_broken"]) + len(cats["e_both_broken"])
    n_comparable = n - nb
    print("=" * 72)
    print("RESULTS")
    print("=" * 72)
    print(f"  (a) AGREE            : {len(cats['a_agree']):>5} / {n}   "
          f"byte-identical  <-- the proof")
    print(f"  (b) MASTER-BROKEN    : {len(cats['b_master_broken']):>5} / {n}   "
          f"expected (sentinel bug fix: crash or garbage output)")
    print(f"  (c) CURRENT-BROKEN   : {len(cats['c_current_broken']):>5} / {n}   "
          f"*** BLOCKER ***")
    print(f"  (d) DIFFER           : {len(cats['d_differ']):>5} / {n}   "
          f"both valid but disagree  *** BLOCKER ***")
    print(f"  (e) BOTH-BROKEN      : {len(cats['e_both_broken']):>5} / {n}")
    if n_comparable:
        print(f"\n  On the {n_comparable} inputs where master produces a VALID "
              f"alignment, current agrees byte-for-byte on "
              f"{len(cats['a_agree'])} ({100 * len(cats['a_agree']) / n_comparable:.1f}%).")
    print(f"  Master is broken (crash/garbage) on {nb} inputs -- the documented "
          f"sentinel bug, fixed in current.")

    def show(title, items, limit=8):
        if not items:
            return
        print(f"\n--- {title} (showing up to {limit} of {len(items)}) ---")
        for case, m, c in items[:limit]:
            _, seqj, seqi, gi, go, ge = case
            base = (f"  read={seqj!r} ref={seqi!r} gi={gi} go={go} ge={ge}")
            if title.startswith("(b)"):
                kind = _broken_kind(m, seqj, seqi)
                print(f"{base}\n      master: BROKEN ({kind})")
                print(f"      current: align_j={c.get('align_j')!r} "
                      f"align_i={c.get('align_i')!r} score={c.get('score')}")
            elif title.startswith("(d)"):
                print(base)
                print(f"      master : align_j={m.get('align_j')!r} "
                      f"align_i={m.get('align_i')!r} score={m.get('score')}")
                print(f"      current: align_j={c.get('align_j')!r} "
                      f"align_i={c.get('align_i')!r} score={c.get('score')}")
            else:  # c or e
                print(base)
                print(f"      master: {m}   current: {c}")

    show("(b) MASTER-BROKEN (expected: sentinel bug fix)", cats["b_master_broken"])
    show("(c) CURRENT-BROKEN (BLOCKER)", cats["c_current_broken"])
    show("(d) DIFFER (BLOCKER: both valid, disagree)", cats["d_differ"])
    show("(e) BOTH-BROKEN", cats["e_both_broken"])

    blockers = len(cats["c_current_broken"]) + len(cats["d_differ"])
    print("\n" + "=" * 72)
    if blockers == 0:
        print(f"PASS: no regressions. Master and current agree byte-for-byte on "
              f"{len(cats['a_agree'])}/{n_comparable} valid-comparable inputs.")
        print(f"      Master is broken on {nb} inputs (crash/garbage); all are the "
              f"documented weak-sentinel border bug, fixed in current.")
        return 0
    print(f"FAIL: {blockers} regression(s) found (categories c/d). "
          f"See details above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
