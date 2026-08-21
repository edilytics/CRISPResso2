"""Lightweight formal verification of CRISPResso2's ``global_align`` using Crosshair.

Crosshair (crosshair-tool) symbolically explores the input space of each
property function and reports a concrete counterexample whenever an assertion
can fail. The alignment kernels are Cython-typed (``str`` / ``np.ndarray``
parameters), so they cannot receive symbolic objects directly; each property
uses ``crosshair.core.realize()`` to materialize a concrete witness under the
current path condition, runs the real compiled kernels on it, and asserts
oracle properties on the concrete outputs. Crosshair's solver-driven search
then systematically steers toward corner cases — empty-ish sequences,
degenerate scoring matrices, gap incentives at array boundaries (the exact
class of bug found during the diagonal-kernel port, see _P9_REGRESSIONS).

Ground truth comes from ``tests/unit_tests/align_oracle.py``: an independent
pure-Python Gotoh DP, a linear alignment rescorer, and a brute-force
alignment enumerator — the DP itself is anchored against the brute force in
the unit suite (test_oracle_dp_matches_brute_force), so this module builds on
a validated oracle rather than a second hand-transcription.

Properties (bounded to keep the solver tractable):

1. ``check_matches_reference_impl`` — differential: the production
   diagonal-major kernel vs the retained pointer-based reference.
2. ``check_alignment_structurally_valid`` — well-formed alignment output.
3. ``check_score_is_optimal`` — oracle DP optimum == rescored kernel output.

Run from the repo root:

    .pixi/envs/test/bin/pip install crosshair-tool   # once
    .pixi/envs/test/bin/crosshair check tests/formal/test_global_align_crosshair.py \
        --per_path_timeout=30 --per_condition_timeout=300

The module also runs under plain pytest (concrete regression vectors); the
symbolic properties no-op without crosshair installed.
"""

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ORACLE_DIR = os.path.normpath(os.path.join(_HERE, "..", "unit_tests"))
if _ORACLE_DIR not in sys.path:
    sys.path.insert(0, _ORACLE_DIR)

import align_oracle  # noqa: E402

from CRISPResso2 import CRISPResso2Align  # noqa: E402

try:  # crosshair-only; plain pytest imports work without it
    from crosshair.core import realize as _realize
except ImportError:  # pragma: no cover
    _realize = None

MATRIX = CRISPResso2Align.read_matrix(
    os.path.join(os.path.dirname(CRISPResso2Align.__file__), "EDNAFULL")
)


# ---------------------------------------------------------------------------
# Concrete input builders (symbolic ints -> kernel inputs).
# ---------------------------------------------------------------------------

def _ints_to_dna(v, n):
    """Map a (realized) int's base-4 digits to an ACGT string of length n."""
    v = abs(v)
    s = []
    for _ in range(max(n, 0)):
        s.append("ACGT"[v & 3])
        v >>= 2
    return "".join(s)


def build_matrix(kind, a, b):
    """A concrete (256,256) int64 matrix from a family selector."""
    m = np.full((256, 256), -4, dtype=np.int64)
    if kind == 0:  # simple match/mismatch
        np.fill_diagonal(m, 5)
    elif kind == 1:  # one huge asymmetric entry (lures traceback)
        np.fill_diagonal(m, a)
        m[ord("A")][ord("C")] = 40
        m[ord("C")][ord("A")] = -40
    elif kind == 2:  # small realized values
        np.fill_diagonal(m, a)
        m[ord("G")][ord("T")] = b
    else:  # zero-diagonal world (mismatch-driven)
        np.fill_diagonal(m, 0)
    return m


def build_gi(place, val, n):
    """gap_incentive of length n: flat zeros with one bump.

    place 2 (bump at the array END) is the exact P9 regression site.
    """
    g = [0] * n
    if n:
        if place == 0:
            g[0] = val
        elif place == 1:
            g[n // 2] = val
        elif place == 2:
            g[n - 1] = val
    return np.array(g, dtype=np.int64)


# ---------------------------------------------------------------------------
# Symbolic properties (executed by `crosshair check`).
# ---------------------------------------------------------------------------

def check_matches_reference_impl(w1: int, w2: int, len1: int, len2: int,
                                 kind: int, a: int, b: int,
                                 bump: int, place: int, go: int) -> None:
    """post: _"""
    if _realize is None:
        return
    # lengths 1..6 (the kernel contract requires non-empty inputs)
    len1 = _realize(len1) % 6 + 1
    len2 = _realize(len2) % 6 + 1
    seqj = _ints_to_dna(_realize(w1), len1)
    seqi = _ints_to_dna(_realize(w2), len2)
    kind = _realize(kind) % 4
    a = _realize(a) % 9
    b = _realize(b) % 9
    go = -(_realize(go) % 4)  # 0..-3
    ge = -2
    mat = build_matrix(kind, a, b)
    gi = build_gi(_realize(place) % 3, _realize(bump) % 5, len2 + 1)
    got = CRISPResso2Align.global_align(seqj, seqi, mat, gi, go, ge)
    exp = CRISPResso2Align.global_align_pointers(seqj, seqi, mat, gi, go, ge)
    assert got == exp, (
        f"kernel divergence: {got} != {exp} for seqj={seqj!r} seqi={seqi!r} "
        f"gi={gi.tolist()} go={go} ge={ge} kind={kind}"
    )
    return True


def check_alignment_structurally_valid(w1: int, w2: int, len1: int, len2: int,
                                       kind: int, a: int, bump: int,
                                       place: int, go: int) -> None:
    """post: _"""
    if _realize is None:
        return
    len1 = _realize(len1) % 6 + 1
    len2 = _realize(len2) % 6 + 1
    seqj = _ints_to_dna(_realize(w1), len1)
    seqi = _ints_to_dna(_realize(w2), len2)
    mat = build_matrix(_realize(kind) % 4, _realize(a) % 9, 3)
    gi = build_gi(_realize(place) % 3, _realize(bump) % 5, len2 + 1)
    go = -(_realize(go) % 4)
    aj, ai, score = CRISPResso2Align.global_align(seqj, seqi, mat, gi, go, -2)
    assert len(aj) == len(ai), "alignment rows differ in length"
    for cj, ci in zip(aj, ai):
        assert not (cj == "-" and ci == "-"), "all-gap column"
    assert [c for c in aj if c != "-"] == list(seqj), "read chars out of order"
    assert [c for c in ai if c != "-"] == list(seqi), "ref chars out of order"
    assert score >= 0.0
    return True


def check_score_is_optimal(w1: int, w2: int, len1: int, len2: int,
                           bump: int, place: int, go: int) -> None:
    """post: _"""
    if _realize is None:
        return
    len1 = _realize(len1) % 5 + 1
    len2 = _realize(len2) % 5 + 1
    seqj = _ints_to_dna(_realize(w1), len1)
    seqi = _ints_to_dna(_realize(w2), len2)
    mat = build_matrix(0, 5, 3)
    gi = build_gi(_realize(place) % 3, _realize(bump) % 5, len2 + 1)
    go = -(_realize(go) % 4)
    aj, ai, _score = CRISPResso2Align.global_align(seqj, seqi, mat, gi, go, -2)
    rescored = align_oracle.score_alignment(
        aj, ai, seqj, seqi, mat, gi, go, -2)
    optimum = align_oracle.dp_optimal_score(
        seqj, seqi, mat, gi, go, -2)
    assert rescored == optimum, (
        f"suboptimal: rescored {rescored} != optimum {optimum} for "
        f"seqj={seqj!r} seqi={seqi!r} gi={gi.tolist()} go={go}"
    )
    return True


# ---------------------------------------------------------------------------
# Concrete regression vectors (also run under plain pytest).
# Each vector is a historical failure of the diagonal-kernel port (P9) or an
# interesting boundary shape.
# ---------------------------------------------------------------------------

_P9_REGRESSIONS = [
    ("AAA", "AAAAAAA", [0, 0, 0, 0, 0, 0, 0, 6], -1, 0),   # gi bump at END
    ("AAAA", "AA", [0, 0, 2], -1, -2),                      # traceback tie
    ("AAA", "CCA", [0, 0, 0, 3], -1, -1),
]


def _simple_matrix(ident=5, mism=-4):
    m = np.full((256, 256), mism, dtype=np.int64)
    np.fill_diagonal(m, ident)
    return m


def test_p9_regression_vectors_match_pointers():
    for seqj, seqi, gi, go, ge in _P9_REGRESSIONS:
        mat = _simple_matrix()
        gi = np.array(gi, dtype=np.int64)
        got = CRISPResso2Align.global_align(seqj, seqi, mat, gi, go, ge)
        exp = CRISPResso2Align.global_align_pointers(seqj, seqi, mat, gi, go, ge)
        assert got == exp, f"{seqj!r} vs {seqi!r}: {got} != {exp}"


def test_p9_regression_vectors_are_optimal():
    for seqj, seqi, gi, go, ge in _P9_REGRESSIONS:
        mat = _simple_matrix()
        gi = np.array(gi, dtype=np.int64)
        aj, ai, _ = CRISPResso2Align.global_align(seqj, seqi, mat, gi, go, ge)
        rescored = align_oracle.score_alignment(
            aj, ai, seqj, seqi, mat, gi, go, ge)
        optimum = align_oracle.dp_optimal_score(seqj, seqi, mat, gi, go, ge)
        assert rescored == optimum


def test_symbolic_helpers_roundtrip():
    """The builders used inside symbolic properties behave as documented."""
    assert _ints_to_dna(0b11100100, 4) == "ACGT"  # LSB -> first char
    assert _ints_to_dna(0b00011011, 4) == "TGCA"
    assert _ints_to_dna(0, 3) == "AAA"
    assert _ints_to_dna(0b10101010, 4) == "GGGG"  # every digit 2
    gi = build_gi(2, 6, 8)
    assert gi.tolist() == [0, 0, 0, 0, 0, 0, 0, 6]
    assert build_gi(0, 3, 1).tolist() == [3]
