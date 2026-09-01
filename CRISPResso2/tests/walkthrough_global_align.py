#!/usr/bin/env python3
"""Pedagogical walkthrough of the pointer-free global_align.

Mirrors CRISPResso2Align.pyx global_align exactly (min_score sentinel, the
inner / last-column / last-row fill blocks, and the recomputed-argmax
traceback), but in plain Python and with everything printed: the three DP
matrices and a step-by-step traceback. The final (align_j, align_i, score) is
asserted equal to the compiled global_align so you can trust the trace.

    pixi run -e test python CRISPResso2/tests/walkthrough_global_align.py
"""

import numpy as np

from CRISPResso2 import CRISPResso2Align

NEG = -1_000_000_000  # the near -inf sentinel used for "impossible" border cells


def show(mat, seqi, seqj, name):
    """Pretty-print a DP matrix with seq chars as row/col labels; NEG -> '-inf'."""
    max_i, max_j = len(seqi), len(seqj)
    cols = ["j"] + [f"{c}" for c in ("." + seqj)]
    print(f"\n{name}  (rows i = ref positions, cols j = read positions)")
    print("   " + " ".join(f"{c:>6}" for c in cols))
    for i in range(max_i + 1):
        rowlbl = "." if i == 0 else seqi[i - 1]
        cells = []
        for j in range(max_j + 1):
            v = mat[i][j]
            cells.append("  -inf" if v <= NEG // 2 else f"{v:>6}")
        print(f"{rowlbl}  " + " ".join(cells))


def fill(seqj, seqi, matrix, gi, gap_open, gap_extend):
    """Forward pass: build the M/I/J score matrices (no pointer matrices)."""
    max_i, max_j = len(seqi), len(seqj)
    M = [[NEG] * (max_j + 1) for _ in range(max_i + 1)]
    I = [[NEG] * (max_j + 1) for _ in range(max_i + 1)]
    J = [[NEG] * (max_j + 1) for _ in range(max_i + 1)]

    # border init: M impossible except origin; I/J prefix-gap scores on the
    # edges where they are the only legal continuation.
    M[0][0] = 0
    for j in range(1, max_j + 1):
        I[0][j] = gap_extend * j + gi[0]
    for i in range(1, max_i + 1):
        J[i][0] = gap_extend * i + gi[0]

    def step(i, j, ci, cj, gi_i, gi_im1, gap_pen_i, gap_pen_j):
        # I = best alignment ending in a gap in the REF (read base inserted):
        #   open a gap from M, or extend an existing I-run. Both add gi_i.
        iFromM = gap_pen_i + M[i][j - 1] + gi_i
        iExt = gap_extend + I[i][j - 1] + gi_i
        I[i][j] = iFromM if iFromM > iExt else iExt
        # J = best alignment ending in a gap in the READ (ref base deleted):
        #   open from M (+ gi_im1), or extend an existing J-run (no extra gi).
        jFromM = gap_pen_j + M[i - 1][j] + gi_im1
        jExt = gap_extend + J[i - 1][j]
        J[i][j] = jFromM if jFromM > jExt else jExt
        # M = best alignment ending in a (mis)match: max of the three diagonal
        # predecessors + substitution score. No pointer stored -- just the max.
        best = M[i - 1][j - 1]
        if I[i - 1][j - 1] > best:
            best = I[i - 1][j - 1]
        if J[i - 1][j - 1] > best:
            best = J[i - 1][j - 1]
        M[i][j] = best + matrix[ord(ci), ord(cj)]

    # inner cells: opening a gap costs gap_open ...
    for i in range(1, max_i):
        ci = seqi[i - 1]
        for j in range(1, max_j):
            step(i, j, ci, seqj[j - 1], gi[i], gi[i - 1], gap_open, gap_open)
    # ... last column: gap opening is "free" (use gap_extend) ...
    j = max_j
    for i in range(1, max_i):
        step(i, j, seqi[i - 1], seqj[j - 1], gi[i], gi[i - 1], gap_extend, gap_extend)
    # ... last row: same.
    i = max_i
    for j in range(1, max_j + 1):
        step(i, j, seqi[i - 1], seqj[j - 1], gi[i], gi[i - 1], gap_extend, gap_extend)
    return M, I, J


def traceback(M, I, J, seqj, seqi, gi, gap_open, gap_extend):
    """Pointer-free traceback: recompute each argmax from the score matrices."""
    max_i, max_j = len(seqi), len(seqj)
    i, j = max_i, max_j
    # start: which matrix holds the global optimum at the bottom-right corner?
    if M[i][j] > J[i][j]:
        cur = "M" if M[i][j] > I[i][j] else "I"
    else:
        cur = "J" if J[i][j] > I[i][j] else "I"

    aj, ai = [], []        # built backwards, reversed at the end
    steps = []
    while i > 0 or j > 0:
        if cur == "M":  # diagonal (mis)match move
            mm, ii, jj = M[i - 1][j - 1], I[i - 1][j - 1], J[i - 1][j - 1]
            # recompute which predecessor won (ties -> I, matching the fill)
            if mm > jj:
                nxt = "M" if mm > ii else "I"
            else:
                nxt = "J" if jj > ii else "I"
            aj.append(seqj[j - 1])
            ai.append(seqi[i - 1])
            steps.append((i, j, "M", f"diag; pred M/I/J={mm}/{ii}/{jj} -> {nxt}", seqj[j - 1], seqi[i - 1]))
            i, j = i - 1, j - 1
            cur = nxt
        elif cur == "J":  # up: ref base deleted (gap in read)
            gap_pen = gap_open if (i < max_i and j < max_j) else gap_extend
            jFromM = gap_pen + M[i - 1][j] + gi[i - 1]
            jExt = gap_extend + J[i - 1][j]
            nxt = "M" if jFromM > jExt else "J"
            aj.append("-")
            ai.append(seqi[i - 1])
            steps.append((i, j, "J", f"up;   M-open={jFromM} vs J-ext={jExt} -> {nxt}", "-", seqi[i - 1]))
            i = i - 1
            cur = nxt
        else:             # I: left: read base inserted (gap in ref)
            gap_pen = gap_open if (i < max_i and j < max_j) else gap_extend
            iFromM = gap_pen + M[i][j - 1] + gi[i]
            iExt = gap_extend + I[i][j - 1] + gi[i]
            nxt = "M" if iFromM > iExt else "I"
            aj.append(seqj[j - 1])
            ai.append("-")
            steps.append((i, j, "I", f"left; M-open={iFromM} vs I-ext={iExt} -> {nxt}", seqj[j - 1], "-"))
            j = j - 1
            cur = nxt
    return aj[::-1], ai[::-1], steps


def run(seqj, seqi, gap_open, gap_extend, gap_incentive):
    """Silent core: fill + pointer-free traceback, verified against compiled global_align."""
    M = CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")
    gi = np.array(gap_incentive, dtype=np.int64)
    assert len(gi) == len(seqi) + 1
    ScM, ScI, ScJ = fill(seqj, seqi, M, list(gi), gap_open, gap_extend)
    aj, ai, steps = traceback(ScM, ScI, ScJ, seqj, seqi, list(gi), gap_open, gap_extend)
    align_j, align_i = "".join(aj), "".join(ai)
    matches = sum(1 for a, b in zip(aj, ai) if a == b and a != "-")
    score = round(100 * matches / len(aj), 3)
    exp = CRISPResso2Align.global_align(seqj, seqi, matrix=M, gap_incentive=gi,
                                        gap_open=gap_open, gap_extend=gap_extend)
    assert (align_j, align_i, score) == exp, "MISMATCH vs compiled!"
    return align_j, align_i, score, (ScM, ScI, ScJ), steps


def walkthrough(seqj, seqi, gap_open=-1, gap_extend=-1, gap_incentive=None):
    gi = gap_incentive if gap_incentive is not None else [0] * (len(seqi) + 1)
    print(f"=== align read {seqj!r} to ref {seqi!r}  "
          f"(match=5, mismatch=-4, gap_open={gap_open}, gap_extend={gap_extend}) ===")
    print(f"gap_incentive (len = len(ref)+1): {list(int(x) for x in gi)}")
    align_j, align_i, score, (ScM, ScI, ScJ), steps = \
        run(seqj, seqi, gap_open, gap_extend, gi)
    show(ScM, seqi, seqj, "M (best ending in a match/mismatch)")
    show(ScI, seqi, seqj, "I (best ending in a gap in REF / read insertion)")
    show(ScJ, seqi, seqj, "J (best ending in a gap in READ / ref deletion)")

    print("\nTraceback (start bottom-right, walk to top-left):")
    print(f"  {'i':>2} {'j':>2} {'state':>5}  {'decision':<44} {'read':>4} {'ref':>4}")
    for i, j, st, dec, rj, ri in steps:
        print(f"  {i:>2} {j:>2} {st:>5}  {dec:<44} {rj:>4} {ri:>4}")

    print(f"\nresult: align_j(read)={align_j!r}  align_i(ref)={align_i!r}  score={score}")
    print("(verified byte-identical to compiled global_align)")
    return align_j, align_i, score


def compare(seqj, seqi, configs, gap_open=-1, gap_extend=-1):
    """Show how gap_incentive moves where the gap lands, all verified vs compiled."""
    print(f"\n--- gap_incentive comparison: read {seqj!r} vs ref {seqi!r} ---")
    print(f"  {'gap_incentive':<34}{'read':>10}{'ref':>10}{'score':>7}")
    for label, gi in configs:
        aj, ai, sc = run(seqj, seqi, gap_open, gap_extend, gi)[:3]
        print(f"  {label:<34}{aj:>10}{ai:>10}{sc:>7}")


def cut_column(align_i, cut_after):
    """Screen column (0-based into the aligned strings) right after ref base `cut_after`."""
    ref_seen = 0
    for col, c in enumerate(align_i):
        if c != '-':
            ref_seen += 1
            if ref_seen == cut_after:
                return col + 1
    return len(align_i)


def summarize_indels(align_j, align_i):
    """List (kind, ref_start, ref_end, seq) for each indel; ref coords are 1-indexed."""
    events = []
    ref_pos = 0
    col = 0
    while col < len(align_i):
        if align_j[col] == '-' and align_i[col] != '-':
            start = ref_pos + 1
            seq = ""
            while col < len(align_i) and align_j[col] == '-' and align_i[col] != '-':
                seq += align_i[col]
                ref_pos += 1
                col += 1
            events.append(("deletion", start, ref_pos, seq))
        elif align_i[col] == '-' and align_j[col] != '-':
            seq = ""
            while col < len(align_i) and align_i[col] == '-' and align_j[col] != '-':
                seq += align_j[col]
                col += 1
            events.append(("insertion", ref_pos, ref_pos, seq))
        else:
            ref_pos += 1
            col += 1
    return events


def print_alignment(align_i, align_j, cut_after, title):
    cc = cut_column(align_i, cut_after)
    print(f"\n{title}")
    print(f"  ref : {align_i}")
    print(f"  read: {align_j}")
    print(f"{' ' * (8 + cc)}^  <- cut after ref base {cut_after}")


def crispr_cut_example():
    """A realistic 4bp microhomology-mediated deletion at a CRISPR cut site.

    The amplicon has two identical 'TGCA' copies straddling the cut, so a 4bp
    deletion can be placed on either copy (a score tie). gap_incentive=1 at the
    cut (the production value) anchors the gap to the cut instead of leaving it
    to arbitrary tie-breaking -- exactly what CRISPResso wants so that indels
    across thousands of reads are reported consistently relative to the cut.
    """
    flank5 = "TGCAGTCAGC"         # ref bases 1-10
    mh1 = "TGCA"                  # ref bases 11-14  (microhomology copy 1)
    mh2 = "TGCA"                  # ref bases 15-18  (microhomology copy 2)
    flank3 = "GTACGTACGT"         # ref bases 19-28
    REF = flank5 + mh1 + mh2 + flank3
    CUT = 14                      # cut after ref base 14 (between the two copies)
    READ = flank5 + mh1 + flank3   # copy 2 deleted (the true biological deletion)
    GO, GE = -20, -2              # CRISPResso defaults (proper affine gaps)
    n = len(REF) + 1

    print("=" * 78)
    print("Realistic CRISPR cut-site example: 4bp microhomology-mediated deletion")
    print("=" * 78)
    print(f"  ref  ({len(REF):>2}bp): {REF}")
    print(f"  read ({len(READ):>2}bp): {READ}    (one 'TGCA' copy deleted)")
    print(f"  cut after ref base {CUT} -- between the two identical 'TGCA' copies")
    print(f"  gap_open={GO} gap_extend={GE} (CRISPResso defaults); gap_incentive default = 1")

    # The read deleted copy 2, but copy 1 is identical, so the alignment can
    # legally place the 4bp gap on either copy (a score tie). Show that the
    # cut-site incentive deterministically anchors it to start at the cut.
    def placement(gi):
        aj, ai, sc = run(READ, REF, GO, GE, gi)[:3]
        dels = [(s, e, seq) for kind, s, e, seq in summarize_indels(aj, ai) if kind == "deletion"]
        return "; ".join(f"delete {seq!r} at bases {s}-{e}" for s, e, seq in dels), sc

    print("\nWhere does the 4bp deletion land? (a tie between the two copies)")
    print(f"  {'gap_incentive':<34}{'deletion placement':<30}{'score':>7}")
    for label, k in [("no incentive", None),
                     ("incentive[10]=1  (LEFT of cut)", 10),
                     ("incentive[14]=1  (AT the cut)", 14)]:
        gi = [0] * n if k is None else [1 if i == k else 0 for i in range(n)]
        desc, sc = placement(gi)
        print(f"  {label:<34}{desc:<30}{sc:>7}")
    print("\n  -> With no incentive the gap lands wherever the tie-break happens to put it;")
    print("     the cut-site incentive (production default = 1) anchors it to start at the")
    print("     cut, so indels are reported consistently relative to the cut site across reads.")

    gi_cut = [0] * n
    gi_cut[CUT] = 1
    aj, ai, sc = run(READ, REF, GO, GE, gi_cut)[:3]
    print_alignment(ai, aj, CUT, f"at-cut alignment (gap_incentive[{CUT}]=1), score={sc}")


if __name__ == "__main__":
    # Example 1: read is missing the middle ref base -> a J gap (deletion).
    walkthrough("AG", "ACG")
    print("\n" + "=" * 78 + "\n")
    # Example 2: read has an extra base -> an I gap (insertion / gap in ref).
    walkthrough("ACGT", "ACG")
    print("\n" + "=" * 78 + "\n")
    # Example 3: gap_incentive. A homopolymer makes the single-deletion position
    # ambiguous; gap_incentive at a ref position biases the gap to land there
    # (this is how CRISPResso prefers gaps at the predicted cut site).
    walkthrough("AAA", "AAAA", gap_incentive=[0, 0, 10, 0, 0])
    compare("AAA", "AAAA", [
        ("[0,0,0,0,0]   no incentive", [0, 0, 0, 0, 0]),
        ("[0,10,0,0,0]  incentive at idx 1", [0, 10, 0, 0, 0]),
        ("[0,0,10,0,0]  incentive at idx 2", [0, 0, 10, 0, 0]),
        ("[0,0,0,10,0]  incentive at idx 3", [0, 0, 0, 10, 0]),
        ("[0,0,0,0,10]  incentive at idx 4", [0, 0, 0, 0, 10]),
    ])
    print("\n" + "=" * 78 + "\n")
    crispr_cut_example()
