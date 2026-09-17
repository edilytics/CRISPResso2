"""Unit tests for CRISPResso2Align."""

from CRISPResso2 import CRISPResso2Align
import numpy as np
import pytest
from hypothesis import given, settings, strategies as st, example

# Independent optimal-score oracle (pure-python Gotoh; shares no logic with
# the .pyx). Lives in the test tree only -- never imported by core code.
from align_oracle import dp_optimal_score, score_alignment, brute_optimal_score

ALN_MATRIX = CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")
AA_MATRIX = CRISPResso2Align.read_matrix("./CRISPResso2/BLOSUM62")


# =============================================================================
# Tests for read_matrix
# =============================================================================


def test_read_matrix_ednafull():
    """Test reading EDNAFULL matrix returns correct scores."""
    matrix = CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")
    # Match score for A-A should be 5, mismatch A-T should be -4
    assert matrix[ord('A'), ord('A')] == 5
    assert matrix[ord('A'), ord('T')] == -4


def test_read_matrix_blosum62():
    """Test reading BLOSUM62 matrix returns correct scores."""
    matrix = CRISPResso2Align.read_matrix("./CRISPResso2/BLOSUM62")
    # Self-score for M should be 5 in BLOSUM62
    assert matrix[ord('M'), ord('M')] == 5


# =============================================================================
# Tests for global_align - basic alignment
# =============================================================================


def test_global_align():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "ATTA"
    assert seq2 == "ATTA"
    assert score == 100

def test_global_align_BLOSUM62():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTA', 'ATTA', matrix=AA_MATRIX, gap_incentive=np.array([0,0,0,0,0],dtype=int))
    assert seq1 == 'ATTA'
    assert seq2 == 'ATTA'
    assert score == 100

def test_global_align_gap_incentive_pos1_BLOSUM62():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=AA_MATRIX, gap_incentive=np.array([0,1,0,0,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'A-TTA'
    assert round(score,3) == round(100*4/5.0,3)

def test_global_align_gap_incentive_pos2_BLOSUM62():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=AA_MATRIX, gap_incentive=np.array([0,0,1,0,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'AT-TA'
    assert round(score,3) == round(100*4/5.0,3)

def test_global_align_gap_incentive_pos3_BLOSUM62():
    """General alignment tests."""
    seq1, seq2, score = CRISPResso2Align.global_align('ATTTA', 'ATTA', matrix=AA_MATRIX, gap_incentive=np.array([0,0,0,1,0],dtype=int))
    assert seq1 == 'ATTTA'
    assert seq2 == 'ATT-A'
    assert round(score,3) == round(100*4/5.0,3)


def test_global_align_BLOSUM62():
    """Test amino acid alignment with BLOSUM62 matrix."""
    seq1, seq2, score = CRISPResso2Align.global_align("MRWY", "MRWY", matrix=AA_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "MRWY"
    assert seq2 == "MRWY"
    assert score == 100


def test_global_align_score_perfect_match():
    """Test that perfect match gives score of 100."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "ATCGATCG", "ATCGATCG", matrix=ALN_MATRIX,
        gap_incentive=np.array([0] * 9, dtype=int)
    )
    assert score == 100


def test_global_align_symmetry():
    """Test that alignment is symmetric (swapping inputs gives same alignment)."""
    seq1_a, seq2_a, score_a = CRISPResso2Align.global_align(
        "ATCG", "ATCC", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    seq1_b, seq2_b, score_b = CRISPResso2Align.global_align(
        "ATCC", "ATCG", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    # Scores should be equal
    assert score_a == score_b


# =============================================================================
# Tests for global_align - gap incentives with BLOSUM62
# =============================================================================


def test_global_align_gap_incentive_pos1_BLOSUM62():
    """Test amino acid alignment with gap incentive at position 1."""
    # Use MRRRY vs MRRY - repeated R allows gap at different positions
    seq1, seq2, score = CRISPResso2Align.global_align("MRRRY", "MRRY", matrix=AA_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "MRRRY"
    assert seq2 == "M-RRY"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


def test_global_align_gap_incentive_pos2_BLOSUM62():
    """Test amino acid alignment with gap incentive at position 2."""
    # Use MRRRY vs MRRY - repeated R allows gap at different positions
    seq1, seq2, score = CRISPResso2Align.global_align("MRRRY", "MRRY", matrix=AA_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "MRRRY"
    assert seq2 == "MR-RY"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


def test_global_align_gap_incentive_pos3_BLOSUM62():
    """Test amino acid alignment with gap incentive at position 3."""
    # Use MRRRY vs MRRY - repeated R allows gap at different positions
    seq1, seq2, score = CRISPResso2Align.global_align("MRRRY", "MRRY", matrix=AA_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "MRRRY"
    assert seq2 == "MRR-Y"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


# =============================================================================
# Tests for global_align - gap incentives in sequence 1
# =============================================================================


def test_global_align_gap_incentive_s1():
    """Test the global_align gap incentives for gaps in sequence 1 (the first sequence)."""
    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0, 0], dtype=int))
    #    print('seq1: ' + seq1 + ' seq2: ' + seq2 + ' score ' + str(score))
    assert seq1 == "ATTTA"
    assert seq2 == "ATTTA"
    assert score == 100

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([1, 0, 0, 0, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "ATT-A"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "A-TTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "AT-TA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "ATT-A"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "ATT-A"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTTA", "ATTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 1], dtype=int))
    assert seq1 == "ATTTA"
    assert seq2 == "ATT-A"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "TTTT-"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([1, 0, 0, 0, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "-TTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "T-TTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "TT-TT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "TTT-T"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTTT", "TTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 1], dtype=int))
    assert seq1 == "TTTTT"
    assert seq2 == "TTTT-"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


# =============================================================================
# Tests for global_align - gap incentives in sequence 2
# =============================================================================


def test_global_align_gap_incentive_s2():
    """Test the global_align gap incentives for gaps in sequence 2 (the second sequence)."""
    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([1, 0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "ATT-A"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0, 0], dtype=int))
    assert seq1 == "A-TTA"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "AT-TA"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "ATT-A"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "ATT-A"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("ATTA", "ATTTA", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0, 1], dtype=int))
    assert seq1 == "ATT-A"
    assert seq2 == "ATTTA"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([1, 0, 0, 0, 0, 0], dtype=int))
    assert seq1 == "-TTTT"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 1, 0, 0, 0, 0], dtype=int))
    assert seq1 == "T-TTT"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 1, 0, 0, 0], dtype=int))
    assert seq1 == "TT-TT"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 1, 0, 0], dtype=int))
    assert seq1 == "TTT-T"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 1, 0], dtype=int))
    assert seq1 == "TTTT-"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)

    seq1, seq2, score = CRISPResso2Align.global_align("TTTT", "TTTTT", matrix=ALN_MATRIX, gap_incentive=np.array([0, 0, 0, 0, 0, 1], dtype=int))
    assert seq1 == "TTTT-"
    assert seq2 == "TTTTT"
    assert round(score, 3) == round(100 * 4 / 5.0, 3)


# =============================================================================
# Tests for global_align - edge cases
# =============================================================================


def test_global_align_single_base():
    """Test alignment with single base sequences."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "A", "A", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0], dtype=int)
    )
    assert seq1 == "A"
    assert seq2 == "A"
    assert score == 100


def test_global_align_with_n():
    """Test alignment with N bases."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "ANNG", "ATCG", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    assert seq1 == "A-NNG"
    assert seq2 == "ATC-G"
    assert score == 40.0


def test_global_align_all_n():
    """Test alignment with all N bases."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "NNNN", "NNNN", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    assert seq1 == "NNNN"
    assert seq2 == "NNNN"


def test_global_align_repeated_base():
    """Test alignment with repeated bases (homopolymer)."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "AAAAAAAAAA", "AAAAAAAAAA", matrix=ALN_MATRIX,
        gap_incentive=np.array([0] * 11, dtype=int)
    )
    assert seq1 == "AAAAAAAAAA"
    assert seq2 == "AAAAAAAAAA"
    assert score == 100


def test_global_align_completely_different():
    """Test alignment with completely different sequences."""
    seq1, seq2, score = CRISPResso2Align.global_align(
        "AAAA", "TTTT", matrix=ALN_MATRIX,
        gap_incentive=np.array([0, 0, 0, 0, 0], dtype=int)
    )
    assert seq1 == "---AAAA"
    assert seq2 == "TTTT---"
    assert score == 0.0


# =============================================================================
# Property-based differential test: global_align vs global_align_ptrfree
# =============================================================================
#
# global_align (the production implementation) is pointer-free: it drops the
# three pointer matrices and recomputes each cell's argmax during traceback
# (see CRISPResso2Align.pyx), ~2.5x faster with half the memory. global_align_pointers
# is the original implementation kept as the reference oracle. Their invariant:
# produce the EXACT same (align_j, align_i, score) on every input, including
# identical tie-breaking. Hypothesis searches the input space (sequences, gap
# penalties, per-position gap_incentive, scoring matrix) for a counterexample;
# the @example cases pin known edge cases (identical, indels, all-different,
# gap-incentive, and formerly-crashing short inputs) so they always run.


@st.composite
def _align_inputs(draw):
    """Draw a self-consistent (seqj, seqi, gap_incentive, gap_open, gap_extend).

    gap_incentive must have length len(seqi)+1, so it is derived from the drawn
    seqi rather than drawn independently.

    Parameter ranges span every production code path:
      - needleman_wunsch / flexiguide defaults: gap_open=-20, gap_extend=-2
      - prime editing defaults:               gap_open=-50, gap_extend=0
    so gap_open spans [-50,-1] and gap_extend spans [-20, 0] (0 included: free
    gap extension is the regime where tie-breaking is most stressed).

    Lengths start at 1 (single bases) to exercise the DP border heavily.
    global_align previously crashed on some short / large-gap_incentive inputs
    because its min_score sentinel was too weak (optimal paths reached border
    cells with uninitialized pointer memory -> garbage / segfault). That is now
    fixed with a true -inf sentinel, so the border is well-defined and both
    implementations agree there; the @example cases pin those formerly-buggy
    inputs so a sentinel regression can never silently slip through.
    """
    seqi = draw(st.text("ACGTN", min_size=1, max_size=40))
    seqj = draw(st.text("ACGTN", min_size=1, max_size=40))
    gi_len = len(seqi) + 1
    gi_vals = draw(
        st.lists(st.integers(min_value=-5, max_value=20),
                 min_size=gi_len, max_size=gi_len))
    gap_open = draw(st.integers(min_value=-50, max_value=-1))
    gap_extend = draw(st.integers(min_value=-20, max_value=0))
    return seqj, seqi, np.array(gi_vals, dtype=np.int64), gap_open, gap_extend


@st.composite
def _tiny_inputs(draw):
    """Tiny inputs (len<=6) for exhaustive brute-force cross-checks of the oracle."""
    seqi = draw(st.text("ACGTN", min_size=1, max_size=6))
    seqj = draw(st.text("ACGTN", min_size=1, max_size=6))
    gi_len = len(seqi) + 1
    gi_vals = draw(
        st.lists(st.integers(min_value=-5, max_value=20),
                 min_size=gi_len, max_size=gi_len))
    gap_open = draw(st.integers(min_value=-50, max_value=-1))
    gap_extend = draw(st.integers(min_value=-20, max_value=0))
    return seqj, seqi, np.array(gi_vals, dtype=np.int64), gap_open, gap_extend


# Diverse scoring matrices: the two shipped NCBI matrices plus an extreme
# make_matrix to stress match/mismatch magnitude.
_MATRICES = [
    pytest.param(ALN_MATRIX, id="EDNAFULL"),
    pytest.param(AA_MATRIX, id="BLOSUM62"),
    pytest.param(CRISPResso2Align.make_matrix(10, -10, -4, -1), id="make_matrix_extreme"),
]


@pytest.mark.parametrize("matrix", _MATRICES)
@given(inputs=_align_inputs())
@settings(max_examples=200, deadline=None)
@example(inputs=("A", "C", np.array([0, 0], dtype=np.int64), -1, -1))
@example(inputs=("AA", "C", np.array([0, 0, 0], dtype=np.int64), -1, -1))
@example(inputs=("CAAAA", "AAAAA", np.array([0, 0, 0, 0, 0, 15], dtype=np.int64), -1, -1))
@example(inputs=("ACGTA", "ACGTA", np.array([0, 0, 0, 0, 0, 0], dtype=np.int64), -1, -1))
@example(inputs=("ACGTACGT", "ACGTA", np.array([0, 0, 0, 0, 0, 0], dtype=np.int64), -1, -1))
@example(inputs=("AAAAA", "TTTTT", np.array([0, 0, 0, 0, 0, 0], dtype=np.int64), -1, -1))
def test_global_align_ptrfree_matches_pointers(matrix, inputs):
    """Pointer-free global_align is byte-identical to global_align_pointers."""
    seqj, seqi, gap_incentive, gap_open, gap_extend = inputs
    expected = CRISPResso2Align.global_align_pointers(
        seqj, seqi, matrix=matrix, gap_incentive=gap_incentive,
        gap_open=gap_open, gap_extend=gap_extend)
    actual = CRISPResso2Align.global_align(
        seqj, seqi, matrix=matrix, gap_incentive=gap_incentive,
        gap_open=gap_open, gap_extend=gap_extend)
    assert actual == expected


# =============================================================================
# Property-based tests: output well-formedness, determinism, and optimality
# =============================================================================
#
# These guard against bug classes the differential test above CANNOT see. The
# differential test proves new==old, but both .pyx implementations were written
# by the same author from the same recurrence; a shared misunderstanding (e.g.
# a wrong transition, or a suboptimal traceback) would pass new==old while the
# output is still wrong. The independent pure-python oracle in align_oracle.py
# breaks that coupling:
#
#   - well-formedness: structural invariants on (align_j, align_i) that any
#     valid global alignment must satisfy, regardless of scoring.
#   - optimality: rescore(global_align(...)) == dp_optimal_score(...) -- the
#     returned alignment is provably optimal under an independently derived
#     model. A shared suboptimal-traceback bug would leave rescore < optimum.
#   - determinism: the same input twice yields identical output. The score
#     matrices are np.empty (uninitialized); this catches any path that reads
#     an uninitialized cell (the same class of bug that crashed the old impl).
#
# The oracle itself is anchored by test_oracle_dp_matches_brute_force below
# (exhaustive enumeration on tiny inputs), so we are not trusting a black box.


def _assert_alignment_well_formed(align_j, align_i, seqj, seqi):
    """Structural invariants any valid global alignment must satisfy."""
    assert len(align_j) == len(align_i), (
        "alignment rows differ in length")
    # no column may be a gap in both sequences
    assert not any(a == "-" and b == "-" for a, b in zip(align_j, align_i)), (
        "alignment has a (-,-) column")
    # ungapping each row must recover the original sequence exactly
    assert align_j.replace("-", "") == seqj, (
        f"align_j ungapped != seqj: {align_j!r} vs {seqj!r}")
    assert align_i.replace("-", "") == seqi, (
        f"align_i ungapped != seqi: {align_i!r} vs {seqi!r}")
    # an alignment of (m, n) sequences has between max(m,n) and m+n columns
    lo = max(len(seqi), len(seqj))
    hi = len(seqi) + len(seqj)
    assert lo <= len(align_j) <= hi, (
        f"alignment length {len(align_j)} outside [{lo}, {hi}]")


@pytest.mark.parametrize("matrix", _MATRICES)
@given(inputs=_align_inputs())
@settings(max_examples=300, deadline=None)
@example(inputs=("A", "C", np.array([0, 0], dtype=np.int64), -1, 0))
@example(inputs=("ACGTACGT", "ACGT", np.array([0, 0, 0, 0, 0], dtype=np.int64), -50, 0))
@example(inputs=("AAAA", "TTTT", np.array([0, 0, 0, 0, 0], dtype=np.int64), -50, 0))
def test_global_align_output_well_formed(matrix, inputs):
    """global_align returns a structurally valid, deterministic alignment."""
    seqj, seqi, gap_incentive, gap_open, gap_extend = inputs
    aj, ai, _ = CRISPResso2Align.global_align(
        seqj, seqi, matrix=matrix, gap_incentive=gap_incentive,
        gap_open=gap_open, gap_extend=gap_extend)
    _assert_alignment_well_formed(aj, ai, seqj, seqi)
    # determinism: identical result on a second call (matrices are np.empty)
    aj2, ai2, _ = CRISPResso2Align.global_align(
        seqj, seqi, matrix=matrix, gap_incentive=gap_incentive,
        gap_open=gap_open, gap_extend=gap_extend)
    assert (aj, ai) == (aj2, ai2), "global_align is non-deterministic across calls"


@pytest.mark.parametrize("matrix", _MATRICES)
@given(inputs=_align_inputs())
@settings(max_examples=300, deadline=None)
@example(inputs=("A", "C", np.array([0, 0], dtype=np.int64), -1, 0))
@example(inputs=("ACGTACGT", "ACGT", np.array([0, 0, 0, 0, 0], dtype=np.int64), -50, 0))
@example(inputs=("AAAA", "TTTT", np.array([0, 0, 0, 0, 0], dtype=np.int64), -50, 0))
def test_global_align_is_optimal(matrix, inputs):
    """global_align's returned alignment is provably optimal.

    The returned alignment, re-scored by the independent pure-python oracle,
    must equal the oracle's independently computed optimal NW score. A higher
    oracle optimum would mean global_align returned a suboptimal alignment
    (a shared-recurrence bug that new==old cannot detect).
    """
    seqj, seqi, gap_incentive, gap_open, gap_extend = inputs
    aj, ai, _ = CRISPResso2Align.global_align(
        seqj, seqi, matrix=matrix, gap_incentive=gap_incentive,
        gap_open=gap_open, gap_extend=gap_extend)
    rescored = score_alignment(aj, ai, seqj, seqi, matrix, gap_incentive,
                              gap_open, gap_extend)
    optimum = dp_optimal_score(seqj, seqi, matrix, gap_incentive,
                               gap_open, gap_extend)
    assert rescored == optimum, (
        f"returned alignment scores {rescored} but optimal is {optimum} "
        f"(read={seqj!r} ref={seqi!r} gi={list(gap_incentive)} "
        f"go={gap_open} ge={gap_extend}); global_align returned a suboptimal "
        f"alignment: align_j={aj!r} align_i={ai!r}")


@pytest.mark.parametrize("matrix", _MATRICES)
@given(inputs=_tiny_inputs())
@settings(max_examples=400, deadline=None)
def test_oracle_dp_matches_brute_force(matrix, inputs):
    """Anchor the oracle: its DP optimum == exhaustive brute-force optimum.

    On tiny inputs (len<=6) we enumerate every legal gapped alignment and take
    the max score; it must equal the oracle's DP. This proves the oracle's own
    recurrence is complete/correct independently of global_align, so
    test_global_align_is_optimal is comparing against a trustworthy reference
    rather than a black box.
    """
    seqj, seqi, gap_incentive, gap_open, gap_extend = inputs
    dp = dp_optimal_score(seqj, seqi, matrix, gap_incentive, gap_open, gap_extend)
    bf = brute_optimal_score(seqj, seqi, matrix, gap_incentive, gap_open,
                             gap_extend, cap=6)
    assert bf is not None, "brute force unexpectedly skipped (input too large)"
    assert dp == bf, (
        f"oracle DP={dp} != brute={bf} for read={seqj!r} ref={seqi!r} "
        f"gi={list(gap_incentive)} go={gap_open} ge={gap_extend}: "
        f"the oracle itself is wrong, so optimality checks cannot be trusted")


if __name__ == "__main__":
    # execute only if run as a script
    test_global_align()
    test_global_align_gap_incentive_s1()
    test_global_align_gap_incentive_s2()
