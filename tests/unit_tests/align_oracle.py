"""Independent optimal-score oracle for CRISPResso2Align.global_align.

This module re-implements the CRISPResso Gotoh affine-gap Needleman-Wunsch
*scoring model* in pure Python, deliberately NOT mirroring the structure of the
production Cython (``CRISPResso2/CRISPResso2Align.pyx``). Its purpose is to be
an independent witness in the differential/property tests:

    rescore(global_align(seqj, seqi, ...)) == dp_optimal_score(seqj, seqi, ...)

i.e. the alignment that ``global_align`` returns is provably *optimal* under
the model. This catches a class of bug that the same-author differential test
(``global_align`` vs ``global_align_pointers``) cannot: if both .pyx
implementations share a recurrence misunderstanding and return the *same*
suboptimal alignment, the differential test passes but the optimal score
computed here would be strictly higher than the returned alignment's score.

Tie-breaking is intentionally NOT modeled here. Optimal alignment *strings*
can legitimately differ when scores tie; only the optimal *score* is unique.
That is what we compare.

THE MODEL (re-derived from CRISPResso2Align.pyx, documented so it can be
audited independently of this code):

States over grid (i, j), i = ref (seqi) index 0..max_i, j = read (seqj) index
0..max_j:
  M[i,j] best alignment ending in a (mis)match consuming ref[i-1], read[j-1]
  I[i,j] best alignment ending in a gap in the REF (read base inserted);
         consumes read[j-1], ref stays. "insertion in ref / deletion in read".
  J[i,j] best alignment ending in a gap in the READ (ref base deleted);
         consumes ref[i-1], read stays. "deletion in ref / insertion in read".

Sentinel NEG = -1_000_000_000 marks "impossible" border cells (no valid
alignment can end there), matching the production ``min_score``.

Border init:
  M[0,0] = 0;  M[0,1:] = NEG;  M[1:,0] = NEG
  I[0,j] = gap_extend * j + gi[0]            for j = 1..max_j   (prefix read overhang)
  I[:,0] = NEG
  J[i,0] = gap_extend * i + gi[0]            for i = 1..max_i   (prefix ref overhang)
  J[0,:] = NEG
  Note the CRISPResso prefix quirk: prefix gaps open with gap_extend (NOT
  gap_open), and gi[0] is added to the whole prefix run (i.e. once, not per
  step) -- see score_alignment for the per-step marginal form.

Inner cells (1 <= i < max_i, 1 <= j < max_j):
  I[i,j] = max(gap_open + M[i,j-1] + gi[i], gap_extend + I[i,j-1] + gi[i])
  J[i,j] = max(gap_open + M[i-1,j] + gi[i-1], gap_extend + J[i-1,j])
  M[i,j] = max(M[i-1,j-1], I[i-1,j-1], J[i-1,j-1]) + matrix[ref,read]

Last column (j == max_j, 1 <= i < max_i) and last row (i == max_i):
  identical to the inner recurrences EXCEPT gap_open is replaced by gap_extend
  in the I-open and J-open terms. (CRISPResso: "for last column and last row,
  ignore gap opening penalty".) The M recurrence is unchanged everywhere.

gap_incentive semantics: gi[i] is added on BOTH opening and extending an I-gap
at ref index i; gi[i-1] is added ONLY when opening a J-gap at cell (i,*) (not
on J-extension). This asymmetry is reproduced exactly.

Optimal score = max(M[max_i,max_j], I[max_i,max_j], J[max_i,max_j]).
"""

NEG = -1_000_000_000


def dp_matrices(seqj, seqi, matrix, gap_incentive, gap_open, gap_extend):
    """Forward DP: return (M, I, J) as lists of lists of ints.

    Pure-python Gotoh fill. ``matrix`` is indexed by ord(char). ``gap_incentive``
    is a sequence of length len(seqi)+1.
    """
    max_i, max_j = len(seqi), len(seqj)
    gi = gap_incentive
    NEG_ = NEG

    M = [[NEG_] * (max_j + 1) for _ in range(max_i + 1)]
    I = [[NEG_] * (max_j + 1) for _ in range(max_i + 1)]
    J = [[NEG_] * (max_j + 1) for _ in range(max_i + 1)]

    M[0][0] = 0
    for j in range(1, max_j + 1):
        I[0][j] = gap_extend * j + gi[0]
    for i in range(1, max_i + 1):
        J[i][0] = gap_extend * i + gi[0]

    def fill_cell(i, j, open_pen):
        ci = ord(seqi[i - 1])
        cj = ord(seqj[j - 1])
        sub = matrix[ci, cj]
        # I = gap in ref (read consumed): open from M or extend from I; both add gi[i]
        i_open = open_pen + M[i][j - 1] + gi[i]
        i_ext = gap_extend + I[i][j - 1] + gi[i]
        I[i][j] = i_open if i_open > i_ext else i_ext
        # J = gap in read (ref consumed): open from M adds gi[i-1]; extend adds nothing
        j_open = open_pen + M[i - 1][j] + gi[i - 1]
        j_ext = gap_extend + J[i - 1][j]
        J[i][j] = j_open if j_open > j_ext else j_ext
        # M = best diagonal predecessor + substitution
        best = M[i - 1][j - 1]
        if I[i - 1][j - 1] > best:
            best = I[i - 1][j - 1]
        if J[i - 1][j - 1] > best:
            best = J[i - 1][j - 1]
        M[i][j] = best + sub

    # inner cells: gap_open for opening
    for i in range(1, max_i):
        for j in range(1, max_j):
            fill_cell(i, j, gap_open)
    # last column: gap_extend for opening
    for i in range(1, max_i):
        fill_cell(i, max_j, gap_extend)
    # last row (includes corner): gap_extend for opening
    for j in range(1, max_j + 1):
        fill_cell(max_i, j, gap_extend)

    return M, I, J


def dp_optimal_score(seqj, seqi, matrix, gap_incentive, gap_open, gap_extend):
    """Optimal Needleman-Wunsch score under the CRISPResso model."""
    if len(seqi) == 0 or len(seqj) == 0:
        # global_align is undefined for empty input (contract: non-empty);
        # degenerate case scored directly: only prefix gaps are possible.
        if len(seqi) == 0 and len(seqj) == 0:
            return 0
        if len(seqi) == 0:
            return gap_extend * len(seqj) + gap_incentive[0]
        return gap_extend * len(seqi) + gap_incentive[0]
    M, I, J = dp_matrices(seqj, seqi, matrix, gap_incentive, gap_open, gap_extend)
    i, j = len(seqi), len(seqj)
    return max(M[i][j], I[i][j], J[i][j])


def score_alignment(align_j, align_i, seqj, seqi, matrix, gap_incentive,
                    gap_open, gap_extend):
    """Re-score a gapped alignment by walking it left-to-right.

    Structurally independent of the DP fill: a single linear pass tracking
    (i, j) consumed and the current gap-run state. ``align_j`` is the read row,
    ``align_i`` is the ref row (matching global_align's return order).

    Reproduces the border rule (gap_extend for opens on the last row/col) and
    the prefix-row/col init quirk (prefix opens use gap_extend; gi[0] applied
    once to the prefix run, i.e. only on its first step).
    """
    assert len(align_j) == len(align_i), "alignment rows must be equal length"
    max_i, max_j = len(seqi), len(seqj)
    gi = gap_incentive
    NEG_ = NEG

    i = 0          # ref chars consumed
    j = 0          # read chars consumed
    state = "START"
    score = 0

    for r, q in zip(align_i, align_j):  # r=ref char/'-', q=read char/'-'
        if r != "-" and q != "-":       # diagonal (mis)match -> M
            score += matrix[ord(seqi[i]), ord(seqj[j])]
            i += 1
            j += 1
            state = "M"
        elif r == "-" and q != "-":     # gap in REF -> I (read consumed)
            # Gotoh constraint: an I-gap may only OPEN from M (or the prefix)
            # or EXTEND an I-run. A direct J->I switch is forbidden -- the
            # production DP's I-matrix only transitions from M or I, never J.
            jn = j + 1                  # arrives at cell (i, jn)
            if i == 0:                  # prefix read overhang: I[0,jn]
                if state == "START":    # first prefix step carries gi[0]
                    score += gap_extend + gi[0]
                elif state == "I":      # extending prefix run
                    score += gap_extend
                else:                   # M/J cannot occur at i==0
                    raise ValueError("illegal transition into prefix I")
            elif state == "I":          # extend an I-run
                score += gap_extend + gi[i]
            elif state == "M":          # open from M
                open_pen = gap_open if (i < max_i and jn < max_j) else gap_extend
                score += open_pen + gi[i]
            else:                       # state == "J": J->I forbidden
                raise ValueError("illegal J->I transition (gaps must be separated by a match)")
            j = jn
            state = "I"
        elif r != "-" and q == "-":     # gap in READ -> J (ref consumed)
            # Gotoh constraint: a J-gap may only OPEN from M (or the prefix)
            # or EXTEND a J-run. A direct I->J switch is forbidden.
            in_ = i + 1                 # arrives at cell (in_, j)
            if j == 0:                  # prefix ref overhang: J[in_,0]
                if state == "START":
                    score += gap_extend + gi[0]
                elif state == "J":      # extending prefix run
                    score += gap_extend
                else:                   # M/I cannot occur at j==0
                    raise ValueError("illegal transition into prefix J")
            elif state == "J":          # extend a J-run (no gi)
                score += gap_extend
            elif state == "M":          # open from M; gi at ref index in_-1
                open_pen = gap_open if (in_ < max_i and j < max_j) else gap_extend
                score += open_pen + gi[in_ - 1]
            else:                       # state == "I": I->J forbidden
                raise ValueError("illegal I->J transition (gaps must be separated by a match)")
            i = in_
            state = "J"
        else:  # r == "-" and q == "-": double-gap column, impossible alignment
            raise ValueError("alignment has a (-,-) column")

    # sanity: the alignment must consume both sequences exactly
    assert i == max_i, f"ref not fully consumed ({i} != {max_i})"
    assert j == max_j, f"read not fully consumed ({j} != {max_j})"
    return score


def brute_optimal_score(seqj, seqi, matrix, gap_incentive, gap_open, gap_extend,
                        cap=7):
    """Brute-force optimal score by enumerating ALL valid alignments.

    A third, structurally independent oracle: explores every gapped alignment
    of the two sequences (choosing match / I-gap / J-gap at each step until
    both are consumed), scores each with score_alignment, returns the max.
    Shares no DP recurrence structure with dp_optimal_score -- it only reuses
    score_alignment's model, so it cross-validates the DP's *transitions*.

    Exponential in the sequence lengths; capped (returns None if too large).
    Use only for tiny inputs (<= cap chars each) to anchor the model.
    """
    if len(seqi) > cap or len(seqj) > cap:
        return None
    max_i, max_j = len(seqi), len(seqj)
    best = [None]

    # alignment built as two lists; scored once complete.
    aj, ai = [], []

    def rec(i, j, state):
        if i == max_i and j == max_j:
            sc = score_alignment("".join(aj), "".join(ai), seqj, seqi,
                                 matrix, gap_incentive, gap_open, gap_extend)
            if best[0] is None or sc > best[0]:
                best[0] = sc
            return
        # match/mismatch (consume both) -- M is reachable from any state
        if i < max_i and j < max_j:
            aj.append(seqj[j]); ai.append(seqi[i])
            rec(i + 1, j + 1, "M")
            aj.pop(); ai.pop()
        # I-gap (gap in ref, consume read) -- legal from START/M/I, NOT from J
        if j < max_j and state != "J":
            aj.append(seqj[j]); ai.append("-")
            rec(i, j + 1, "I")
            aj.pop(); ai.pop()
        # J-gap (gap in read, consume ref) -- legal from START/M/J, NOT from I
        if i < max_i and state != "I":
            aj.append("-"); ai.append(seqi[i])
            rec(i + 1, j, "J")
            aj.pop(); ai.pop()

    rec(0, 0, "START")
    return best[0]
