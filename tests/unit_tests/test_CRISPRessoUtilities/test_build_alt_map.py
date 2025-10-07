import os

import pandas as pd
import pytest

from CRISPResso2.CRISPRessoCOREResources import find_indels_substitutions
from CRISPResso2 import CRISPRessoUtilities as utilities

# ----------------------------- helpers -----------------------------

# This seq is altered for every test case; tests should not depend on its content.
REF_SEQ = "ATGCGTACGATCGTACGTAGCTAGCTAGCGTAGCTAGCTA"  # 40 bp
REF_POSITIONS = list(range(len(REF_SEQ)))

def _normalize_alt_map(alt_map):
    """Sort alt_edits so ordering does not affect equality checks."""
    norm = {}
    for key, val in alt_map.items():
        norm[key] = {
            "ref_seq": val["ref_seq"],
            "alt_edits": sorted(val["alt_edits"], key=lambda t: (t[0], t[1])),
        }
    return norm

def df_from_rows(*rows):
    return pd.DataFrame(list(rows))


# ----------------------------- alteration functions -----------------------------
# These functions are used throughout the tests to alter the base REF_SEQ before calling build_alt_map.

def make_unmodified(reads=1, ref_name="Reference"):
    """Unmodified read identical to REF_SEQ."""
    return {
        "#Reads": reads,
        "Aligned_Sequence": REF_SEQ,
        "Reference_Sequence": REF_SEQ,
        "n_inserted": 0,
        "n_deleted": 0,
        "n_mutated": 0,
        "Reference_Name": ref_name,
        "ref_positions": REF_POSITIONS,
        "insertion_coordinates": [],
        "deletion_coordinates": [],
        "substitution_positions": [],
        "insertion_sizes": [],
    }


def make_sub(position, alt_base, reads=1, ref_name="Reference"):
    """Single‑nucleotide substitution at 0‑based ref index 'position'."""
    aligned = REF_SEQ[:position] + alt_base + REF_SEQ[position + 1 :]
    return {
        "#Reads": reads,
        "Aligned_Sequence": aligned,
        "Reference_Sequence": REF_SEQ,
        "n_inserted": 0,
        "n_deleted": 0,
        "n_mutated": 1,
        "Reference_Name": ref_name,
        "ref_positions": REF_POSITIONS,
        "insertion_coordinates": [],
        "deletion_coordinates": [],
        "substitution_positions": [position],
        "insertion_sizes": [],
    }


def make_del(start, end, reads=1, ref_name="Reference"):
    """Deletion [start, end) in reference coordinates (end exclusive)."""
    deletion_len = end - start
    aligned = REF_SEQ[:start] + "-"*deletion_len + REF_SEQ[end:]  # aligned sequence without the deleted block
    return {
        "#Reads": reads,
        "Aligned_Sequence": aligned,
        "Reference_Sequence": REF_SEQ,
        "n_inserted": 0,
        "n_deleted": 1,
        "n_mutated": 0,
        "Reference_Name": ref_name,
        "ref_positions": REF_POSITIONS,
        "insertion_coordinates": [],
        "deletion_coordinates": [(start, end)],
        "substitution_positions": [],
        "insertion_sizes": [],
    }


def make_ins(after_index, ins_seq, reads=1, ref_name="Reference"):
    """
    Insertion of 'ins_seq' BETWEEN after_index and after_index+1 (0‑based).
    For build_alt_map, insertion_coordinates use the right‑anchor ref index (after_index+1),
    and the aligned start is exactly after_index+1 when there are no prior edits in the row.
    """
    aligned_start = after_index + 1
    aligned = REF_SEQ[:aligned_start] + ins_seq + REF_SEQ[aligned_start:]
    # updated_ref_seq = REF_SEQ [:aligned_start] + "-"*len(ins_seq) + REF_SEQ[aligned_start:]
    return {
        "#Reads": reads,
        "Aligned_Sequence": aligned,
        "Reference_Sequence": REF_SEQ,
        "n_inserted": len(ins_seq),
        "n_deleted": 0,
        "n_mutated": 0,
        "Reference_Name": ref_name,
        "ref_positions": REF_POSITIONS,
        "insertion_coordinates": [(after_index + 1, aligned_start)],
        "deletion_coordinates": [],
        "substitution_positions": [],
        "insertion_sizes": [len(ins_seq)],
    }



# ----------------------------- core mixed case (old test modernized) -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # Mixed: SNP at index 9, 1‑bp del at index 19, 10‑bp del at index 19, 2‑bp insertion after 30, plus an unmodified that must be skipped.
        (
            [
                make_sub(9, "G", reads=1),
                make_del(19, 20, reads=1),
                make_ins(30, "GG", reads=1),
                make_del(19, 29, reads=1),
                make_unmodified(reads=2),
            ],
            {"Reference": (1, 1)},
            {
            # Key is chrom, pos
            # Value is dict with ref_seq and list of alt_edits (each alt_edit is [type, alt_edit, reads])
                # substitution at 1 + 9 = 10
                (1, 10): {
                    "ref_seq": "A",
                    "alt_edits": [["sub", "G", 1]],
                },
                # deletions at 1 + 19 = 20; longest ref_seq retained
                (1, 20): {
                    "ref_seq": REF_SEQ[18:29],  # left flank + 10‑bp deleted block
                    "alt_edits": [
                        ["delete", REF_SEQ[19:20], 1],  # 1‑bp del
                        ["delete", REF_SEQ[19:29], 1],  # 10‑bp del
                    ],
                },
                # insertion after index 30 -> key (1, 32)
                (1, 32): {
                    "ref_seq": REF_SEQ[31],
                    "alt_edits": [["insert", "GG", 1]],
                },
            },
        ),
    ],
    ids=["mixed_edits_happy_path"],
)
def test_build_alt_map_mixed(rows, amplicon_positions, expected):
    df = df_from_rows(*rows)
    alt_map = utilities.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(alt_map) == _normalize_alt_map(expected)


# ----------------------------- substitutions -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # single SNP
        (
            [make_sub(9, "G", reads=1)],
            {"Reference": (1, 1)},
            {(1, 10): {"ref_seq": "A", "alt_edits": [["sub", "G", 1]]}},
        ),
        # merge same SNP (reads sum)
        (
            [make_sub(9, "G", reads=2), make_sub(9, "G", reads=3)],
            {"Reference": (1, 1)},
            {(1, 10): {"ref_seq": "A", "alt_edits": [["sub", "G", 5]]}},
        ),
        # split different alt bases
        (
            [make_sub(9, "G", reads=1), make_sub(9, "T", reads=2)],
            {"Reference": (1, 1)},
            {(1, 10): {"ref_seq": "A", "alt_edits": [["sub", "G", 1], ["sub", "T", 2]]}},
        ),
    ],
    ids=["sub_single", "sub_merge_same_base", "sub_split_diff_base"],
)
def test_build_alt_map_substitutions(rows, amplicon_positions, expected):
    df = df_from_rows(*rows)
    out = utilities.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)


# ----------------------------- deletions -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # single 1‑bp deletion at start=19
        (
            [make_del(19, 20, reads=1)],
            {"Reference": (1, 1)},
            {
                (1, 20): {
                    "ref_seq": REF_SEQ[18:20],  # left flank + deleted
                    "alt_edits": [["delete", REF_SEQ[19:20], 1]],
                }
            },
        ),
        # merge same‑length deletion at same key (reads sum)
        (
            [make_del(19, 20, reads=2), make_del(19, 20, reads=3)],
            {"Reference": (1, 1)},
            {
                (1, 20): {
                    "ref_seq": REF_SEQ[18:20],
                    "alt_edits": [["delete", REF_SEQ[19:20], 5]],
                }
            },
        ),
        # two different lengths at same key -> two entries; keep longest ref_seq
        (
            [make_del(19, 20, reads=1), make_del(19, 29, reads=1)],
            {"Reference": (1, 1)},
            {
                (1, 20): {
                    "ref_seq": REF_SEQ[18:29],  # longest window observed at this key
                    "alt_edits": [
                        ["delete", REF_SEQ[19:20], 1],
                        ["delete", REF_SEQ[19:29], 1],
                    ],
                }
            },
        ),
        # deletion at second element
        (
            [make_del(1, 2, reads=1)],
            {'Reference': (1, 1)},
            {
                (1, 2): {
                    "ref_seq": REF_SEQ[0:2],
                    "alt_edits": [["delete", REF_SEQ[1:2], 1]],
                }
            },
        ),
    ],
    ids=["del_single", "del_merge_same_len", "del_two_lengths_same_key", "del_second_element"],
)
def test_build_alt_map_deletions(rows, amplicon_positions, expected):
    df = df_from_rows(*rows)
    out = utilities.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        (
            [make_del(39, 40, reads=1)],
            {"Reference": (1, 1)},
            {
                (1, 40): {
                    "ref_seq": REF_SEQ[38:40],
                    "alt_edits": [["delete", REF_SEQ[39:40], 1]],
                }
            },
        ),
    ],
)
def test_build_alt_map_deletion_end_at_len_raises(rows, amplicon_positions, expected):
    # delete last base: (39, 40) for a 40‑bp ref
    df = df_from_rows(*rows)
    out = utilities.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)

def test_build_alt_map_deletion_start_at_zero_should_anchor_correctly():
    df = df_from_rows(make_del(0, 3, reads=1))  # delete first 3 bases
    amplicon_positions = {"Reference": (1, 1)}
    out = utilities.build_alt_map(df, amplicon_positions)
    # Desired behavior: ref_seq is just the deleted span (no left flank available)
    expected = {
        (1, 1): {
            "ref_seq": REF_SEQ[0:3],
            "alt_edits": [["delete", REF_SEQ[0:3], 1]],
        }
    }
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)


# ----------------------------- insertions -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # single insertion of "GG" after index 30 -> key (1, 32)
        (
            [make_ins(30, "GG", reads=1)],
            {"Reference": (1, 1)},
            {(1, 32): {"ref_seq": REF_SEQ[31], "alt_edits": [["insert", "GG", 1]]}},
        ),
        # merge same inserted string and key (reads sum)
        (
            [make_ins(30, "GG", reads=2), make_ins(30, "GG", reads=3)],
            {"Reference": (1, 1)},
            {(1, 32): {"ref_seq": REF_SEQ[31], "alt_edits": [["insert", "GG", 5]]}},
        ),
        # split different inserted strings at same key
        (
            [make_ins(30, "GG", reads=1), make_ins(30, "T", reads=1)],
            {"Reference": (1, 1)},
            {(1, 32): {"ref_seq": REF_SEQ[31], "alt_edits": [["insert", "GG", 1], ["insert", "T", 1]]}},
        ),
    ],
    ids=["ins_single", "ins_merge_same_seq", "ins_split_diff_seq"],
)
def test_build_alt_map_insertions(rows, amplicon_positions, expected):
    df = df_from_rows(*rows)
    out = utilities.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)

# ----------------------------- multiple amplicons & coordinate offsets -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        (
            # One SNP on amplicon A, one insertion on amplicon B, ensure independence and absolute coordinates.
            [make_sub(2, "T", reads=1, ref_name="A"),
             make_ins(0, "G", reads=2, ref_name="B")],
            {"A": (1, 1), "B": (2, 1000)},
            {
                # SNP: chrom 1, pos = 1 + 2 = 3
                (1, 3): {"ref_seq": REF_SEQ[2], "alt_edits": [["sub", "T", 1]]},
                # insertion on chrom 2: left_index = 1000 + (0+1) = 1001
                (2, 1001): {"ref_seq": REF_SEQ[1], "alt_edits": [["insert", "G", 2]]},
            },
        ),
    ],
    ids=["two_amplicons_offsets"],
)
def test_build_alt_map_multi_amplicon_and_offsets(rows, amplicon_positions, expected):
    df = df_from_rows(*rows)
    out = utilities.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)


# ----------------------------- unmodified rows are skipped -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected_size",
    [
        ([make_unmodified(reads=5)], {"Reference": (1, 1)}, 0),
        ([make_unmodified(), make_sub(0, "C", 1)], {"Reference": (1, 1)}, 1),
    ],
    ids=["only_unmodified", "mixed_with_unmodified"],
)
def test_build_alt_map_skips_unmodified(rows, amplicon_positions, expected_size):
    df = df_from_rows(*rows)
    out = utilities.build_alt_map(df, amplicon_positions)
    assert len(out) == expected_size

# ----------------------------- ensuring that insertions are keyed by right anchor -----------------------------

def test_build_alt_map_fidelity_like_real_row():

    # Short 40bp reference used across tests
    REF_SEQ = "ATGCGTACGATCGTACGTAGCTAGCTAGCGTAGCTAGCTA"

    # --- build a CRISPResso‑like "real" row ---
    left_pad = 4
    right_pad = 5
    ins_right_anchor = 11        # insertion between 10 and 11; right anchor == 11
    ins_seq = "GGG"
    reads = 7

    # Aligned includes inserted bases; Reference shows '-' at the inserted site
    aligned_seq = (
        "-" * left_pad
        + REF_SEQ[:ins_right_anchor]
        + ins_seq
        + REF_SEQ[ins_right_anchor:]
        + "-" * right_pad
    )
    reference_seq = (
        "-" * left_pad
        + REF_SEQ[:ins_right_anchor]
        + "-" * len(ins_seq)
        + REF_SEQ[ins_right_anchor:]
        + "-" * right_pad
    )

    # ref_positions spans the padded string length; negatives where there's no ref base
    ref_positions = (
        [-1] * left_pad
        + list(range(ins_right_anchor))        # 0..10 map one-to-one
        + [-999] * len(ins_seq)                # inserted segment has no ref coordinate
        + list(range(ins_right_anchor, len(REF_SEQ)))  # 11..39
        + [-1] * right_pad
    )

    row = {
        "#Reads": reads,
        "Aligned_Sequence": aligned_seq,
        "Reference_Sequence": reference_seq,
        "n_inserted": len(ins_seq),
        "n_deleted": 0,
        "n_mutated": 0,
        "Reference_Name": "Reference",
        "ref_positions": ref_positions,
        "insertion_coordinates": [(ins_right_anchor, left_pad + ins_right_anchor)],
        "deletion_coordinates": [],
        "substitution_positions": [],
        "insertion_sizes": [len(ins_seq)],
    }
    df = pd.DataFrame([row])

    # Amplicon starts at absolute position 500 on chromosome 1
    amplicon_positions = {"Reference": (1, 500)}

    # Expected: key uses the right anchor; ref_seq is the anchor base; alt carries inserted string & reads
    expected = {
        (1, 500 + ins_right_anchor): {
            "ref_seq": REF_SEQ[ins_right_anchor],
            "alt_edits": [["insert", ins_seq, reads]],
        }
    }

    out = utilities.build_alt_map(df, amplicon_positions)

    def _normalize(m):
        return {
            k: {"ref_seq": v["ref_seq"], "alt_edits": sorted(v["alt_edits"], key=lambda t: (t[0], t[1]))}
            for k, v in m.items()
        }

    assert _normalize(out) == _normalize(expected)


def test_aln_to_alt_map_ins_del_same_pos():
    ref1 = 'AATGCGTAC'
    aln1 = 'AA--CGTAC'
    payload1 = find_indels_substitutions(aln1, ref1, list(range(len(ref1)))).__dict__
    payload1['Reference_Sequence'] = ref1
    payload1['Aligned_Sequence'] = aln1

    ref2 = 'AA--TGCGTAC'
    aln2 = 'AATTTGCGTAC'
    payload2 = find_indels_substitutions(aln2, ref2, list(range(len(ref2)))).__dict__
    payload2['Reference_Sequence'] = ref2
    payload2['Aligned_Sequence'] = aln2

    rows = [
        payload1,
        payload2,
    ]

    for row in rows:
        row['#Reads'] = 1
        row['Reference_Name'] = 'Reference'
        row['n_inserted'] = row['insertion_n']
        row['n_deleted'] = row['deletion_n']
        row['n_mutated'] = row['substitution_n']

    df = pd.DataFrame(rows)
    amplicon_positions = {"Reference": (1, 1)}
    alt_map = utilities.build_alt_map(df, amplicon_positions)

    assert list(alt_map.keys()) == [(1, 2)]
    assert alt_map[(1, 2)] == {'ref_seq': 'ATG', 'alt_edits': [['delete', 'TG', 1], ['insert', 'TT', 1]]}


def test_aln_to_alt_map_to_vcf():
    ref1 = 'AATGCGTAC'
    aln1 = 'AATGCG-AC'
    #            ^ Interested in this deletion across each of these examples
    payload1 = find_indels_substitutions(aln1, ref1, list(range(len(ref1)))).__dict__
    payload1['Reference_Sequence'] = ref1
    payload1['Aligned_Sequence'] = aln1

    ref2 = 'AA--TGCGTAC'
    aln2 = 'AATTTGCG-AC'
    payload2 = find_indels_substitutions(aln2, ref2, list(range(len(ref2)))).__dict__
    payload2['Reference_Sequence'] = ref2
    payload2['Aligned_Sequence'] = aln2

    ref3 = 'AATGCGTAC'
    aln3 = 'AA-GCG-AC'
    payload3 = find_indels_substitutions(aln3, ref3, list(range(len(ref3)))).__dict__
    payload3['Reference_Sequence'] = ref3
    payload3['Aligned_Sequence'] = aln3

    ref4 = 'AATG--CGTAC'
    aln4 = 'AA-GAACG-AC'
    payload4 = find_indels_substitutions(aln4, ref4, list(range(len(ref4)))).__dict__
    payload4['Reference_Sequence'] = ref4
    payload4['Aligned_Sequence'] = aln4

    ref5 = 'AA--TGCGTAC'
    aln5 = 'AATTTGCG-GG'
    payload5 = find_indels_substitutions(aln5, ref5, list(range(len(ref5)))).__dict__
    payload5['Reference_Sequence'] = ref5
    payload5['Aligned_Sequence'] = aln5

    rows = [
        payload1,
        payload2,
        payload3,
        payload4,
        payload5,
    ]

    for row in rows:
        row['#Reads'] = 1
        row['Reference_Name'] = 'Reference'
        row['n_inserted'] = row['insertion_n']
        row['n_deleted'] = row['deletion_n']
        row['n_mutated'] = row['substitution_n']

    df = pd.DataFrame(rows)
    amplicon_positions = {"Reference": (1, 1)}
    alt_map = utilities.build_alt_map(df, amplicon_positions)

    # deletion at position 5 in each example above, should occur 5 times
    assert alt_map[(1, 7)] == {'ref_seq': 'GT', 'alt_edits': [['delete', 'T', 5]]}
    # insertion of TT occurs 2 times in 2, 5 and deletion of T occurs 2 times in 3, 4
    # I'm not entirely certain what the ref_seq should be in this case... but I do know that there should be a deletion in the alt_edits
    assert alt_map[(1, 2)] == {'ref_seq': 'A', 'alt_edits': [['insert', 'TT', 2], ['delete', 'T', 2]]}
    # insertion of AA occurs 1 time in 4
    assert alt_map[(1, 4)] == {'ref_seq': 'G', 'alt_edits': [['insert', 'AA', 1]]}
    # substitution of A -> G occurs 1 time in 5
    assert alt_map[(1, 8)] == {'ref_seq': 'A', 'alt_edits': [['sub', 'G', 1]]}
    # substitution of C -> G occurs 1 time in 5
    assert alt_map[(1, 9)] == {'ref_seq': 'C', 'alt_edits': [['sub', 'G', 1]]}

    temp_vcf_path = 'aln_to_alt_map_to_vcf.vcf'
    num_vcf_rows = utilities.vcf_lines_from_alt_map(alt_map, 5, ['Reference'], temp_vcf_path)
    assert num_vcf_rows == len(alt_map)

    # os.remove(temp_vcf_path)
