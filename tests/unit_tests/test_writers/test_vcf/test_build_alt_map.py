import os

import pandas as pd
import pytest

from CRISPResso2.CRISPRessoCOREResources import find_indels_substitutions
from CRISPResso2.writers import vcf

# ----------------------------- helpers -----------------------------

# This seq is altered for every test case; tests should not depend on its content.
REF_SEQ = "ATGCGTACGATCGTACGTAGCTAGCTAGCGTAGCTAGCTA"  # 40 bp
REF_POSITIONS = list(range(len(REF_SEQ)))

def _normalize_alt_map(alt_map):
    """Sort alt_edits so ordering does not affect equality checks."""
    norm = {}
    for key, val in alt_map.items():
        norm[key] = {
            "alt_edits": sorted(val["alt_edits"], key=lambda t: (t[0], t[1])),
        }
    return norm


def create_df_alleles(*refs_alns):
    payloads = []
    for ref, *aln in refs_alns:
        if isinstance(aln, list):
            if len(aln) == 1:
                aln = aln[0]
                num_reads = 1
                ref_name = 'Reference'
            elif len(aln) == 2:
                aln, num_reads = aln
                ref_name = 'Reference'
            elif len(aln) == 3:
                aln, num_reads, ref_name = aln
        payload = find_indels_substitutions(aln, ref, list(range(len(ref)))).__dict__
        payload['Reference_Sequence'] = ref
        payload['Aligned_Sequence'] = aln
        payload['#Reads'] = num_reads
        payload['Reference_Name'] = ref_name
        payloads += [payload]

    for payload in payloads:
        payload['n_inserted'] = payload['insertion_n']
        payload['n_deleted'] = payload['deletion_n']
        payload['n_mutated'] = payload['substitution_n']

    return pd.DataFrame(payloads)


# ----------------------------- alteration functions -----------------------------
# These functions are used throughout the tests to alter the base REF_SEQ before calling build_alt_map.

def make_unmodified(reads=1, ref_name="Reference"):
    """Unmodified read identical to REF_SEQ."""
    return (REF_SEQ, REF_SEQ, reads, ref_name)


def make_sub(position, alt_base, reads=1, ref_name="Reference"):
    """Single‑nucleotide substitution at 0‑based ref index 'position'."""
    aligned = REF_SEQ[:position] + alt_base + REF_SEQ[position + 1 :]
    return (REF_SEQ, aligned, reads, ref_name)


def make_del(start, end, reads=1, ref_name="Reference"):
    """Deletion [start, end) in reference coordinates (end exclusive)."""
    deletion_len = end - start
    aligned = REF_SEQ[:start] + "-"*deletion_len + REF_SEQ[end:]  # aligned sequence without the deleted block
    return (REF_SEQ, aligned, reads, ref_name)


def make_ins(after_index, ins_seq, reads=1, ref_name="Reference"):
    """
    Insertion of 'ins_seq' BETWEEN after_index and after_index+1 (0‑based).
    For build_alt_map, insertion_coordinates use the right‑anchor ref index (after_index+1),
    and the aligned start is exactly after_index+1 when there are no prior edits in the row.
    """
    aligned_start = after_index + 1
    aligned = REF_SEQ[:aligned_start] + ins_seq + REF_SEQ[aligned_start:]
    updated_ref_seq = REF_SEQ[:aligned_start] + "-"*len(ins_seq) + REF_SEQ[aligned_start:]
    return (updated_ref_seq, aligned, reads, ref_name)


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
            # Value is dict with list of alt_edits (each alt_edit is [type, alt_edit, reads, ref_seq])
                # substitution at 1 + 9 = 10
                (1, 10): {
                    "alt_edits": [["sub", "G", 1, "A"]],
                },
                # deletions at 19 (including the padding base, the deletion is at 20); each edit has its own ref_seq
                (1, 19): {
                    "alt_edits": [
                        ["delete", REF_SEQ[19:20], 1, REF_SEQ[18:20]],  # 1‑bp del
                        ["delete", REF_SEQ[19:29], 1, REF_SEQ[18:29]],  # 10‑bp del
                    ],
                },
                # insertion after index 30 -> key (1, 32)
                (1, 31): {
                    "alt_edits": [["insert", "GG", 1, REF_SEQ[30]]],
                },
            },
        ),
    ],
    ids=["mixed_edits_happy_path"],
)
def test_build_alt_map_mixed(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    alt_map = vcf.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(alt_map) == _normalize_alt_map(expected)


# ----------------------------- substitutions -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # single SNP
        (
            [make_sub(9, "G", reads=1)],
            {"Reference": (1, 1)},
            {(1, 10): {"alt_edits": [["sub", "G", 1, "A"]]}},
        ),
        # merge same SNP (reads sum)
        (
            [make_sub(9, "G", reads=2), make_sub(9, "G", reads=3)],
            {"Reference": (1, 1)},
            {(1, 10): {"alt_edits": [["sub", "G", 5, "A"]]}},
        ),
        # split different alt bases
        (
            [make_sub(9, "G", reads=1), make_sub(9, "T", reads=2)],
            {"Reference": (1, 1)},
            {(1, 10): {"alt_edits": [["sub", "G", 1, "A"], ["sub", "T", 2, "A"]]}},
        ),
    ],
    ids=["sub_single", "sub_merge_same_base", "sub_split_diff_base"],
)
def test_build_alt_map_substitutions(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_alt_map(df, amplicon_positions)
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
                (1, 19): {
                    "alt_edits": [["delete", REF_SEQ[19:20], 1, REF_SEQ[18:20]]],
                }
            },
        ),
        # merge same‑length deletion at same key (reads sum)
        (
            [make_del(19, 20, reads=2), make_del(19, 20, reads=3)],
            {"Reference": (1, 1)},
            {
                (1, 19): {
                    "alt_edits": [["delete", REF_SEQ[19:20], 5, REF_SEQ[18:20]]],
                }
            },
        ),
        # two different lengths at same key -> two entries; each has its own ref_seq
        (
            [make_del(19, 20, reads=1), make_del(19, 29, reads=1)],
            {"Reference": (1, 1)},
            {
                (1, 19): {
                    "alt_edits": [
                        ["delete", REF_SEQ[19:20], 1, REF_SEQ[18:20]],
                        ["delete", REF_SEQ[19:29], 1, REF_SEQ[18:29]],
                    ],
                }
            },
        ),
        # deletion at second element
        (
            [make_del(1, 2, reads=1)],
            {'Reference': (1, 1)},
            {
                (1, 1): {
                    "alt_edits": [["delete", REF_SEQ[1:2], 1, REF_SEQ[0:2]]],
                }
            },
        ),
    ],
    ids=["del_single", "del_merge_same_len", "del_two_lengths_same_key", "del_second_element"],
)
def test_build_alt_map_deletions(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)


@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        (
            [make_del(39, 40, reads=1)],
            {'Reference': (1, 1)},
            {
                (1, 39):{
                    'alt_edits': [['delete', REF_SEQ[39:49], 1, REF_SEQ[38:40]]],
                }
            },
        ),
        (
            [make_del(38, 40, reads=1)],
            {"Reference": (1, 1)},
            {
                (1, 38): {
                    "alt_edits": [["delete", REF_SEQ[38:40], 1, REF_SEQ[37:40]]],
                }
            },
        ),
    ],
    ids=['last_element', 'second_to_last_element']
)
def test_build_alt_map_deletion_end_at_len_raises(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)


def test_build_alt_map_deletion_start_at_zero_should_anchor_correctly():
    """When a deletion occurs at the start of a sequence, then you record the base after the deletion.

    Source: https://bioinformatics.stackexchange.com/questions/2476/how-to-represent-a-deletion-at-position-1-in-a-vcf-file
    """
    df = create_df_alleles(make_del(0, 3, reads=1))  # delete first 3 bases
    amplicon_positions = {"Reference": (1, 1)}
    out = vcf.build_alt_map(df, amplicon_positions)
    expected = {
        (1, 1): {
            "alt_edits": [["delete_start", REF_SEQ[0:3], 1, REF_SEQ[0:4]]],
        }
    }
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)


def test_build_alt_map_deletion_start():
    ref1 = 'GATTACA'
    aln1 = '-ATTACA'

    df = create_df_alleles(
        (ref1, aln1),
    )

    amplicon_positions = {"Reference": (1, 1)}
    alt_map = vcf.build_alt_map(df, amplicon_positions)

    assert alt_map[(1, 1)] == {
        'alt_edits': [['delete_start', 'G', 1, 'GA']],
    }

    vcf_lines = vcf._write_vcf_lines(1, 1, alt_map[(1, 1)]['alt_edits'], len(df))
    assert len(vcf_lines) == 1
    assert 'GA\tA' in vcf_lines[0]


def test_build_alt_map_multi_deletion_start():
    ref1 = 'AACCTTGG'
    aln1 = '-ACCTTGG'

    ref2 = 'AACCTTGG'
    aln2 = '--CCTTGG'

    ref3 = 'AACCTTGG'
    aln3 = '---CTTGG'

    ref4 = 'AACCTTGG'
    aln4 = '----TTGG'

    ref5 = 'AACCTTGG'
    aln5 = '-----TGG'

    ref6 = 'AACCTTGG'
    aln6 = '------GG'

    ref7 = 'AACCTTGG'
    aln7 = '-------G'

    df = create_df_alleles(
        (ref1, aln1),
        (ref2, aln2),
        (ref3, aln3),
        (ref4, aln4),
        (ref5, aln5),
        (ref6, aln6),
        (ref7, aln7),
    )

    amplicon_positions = {"Reference": (1, 1)}
    alt_map = vcf.build_alt_map(df, amplicon_positions)

    assert alt_map[(1, 1)] == {
        'alt_edits': [
            ['delete_start', 'A', 1, 'AA'],
            ['delete_start', 'AA', 1, 'AAC'],
            ['delete_start', 'AAC', 1, 'AACC'],
            ['delete_start', 'AACC', 1, 'AACCT'],
            ['delete_start', 'AACCT', 1, 'AACCTT'],
            ['delete_start', 'AACCTT', 1, 'AACCTTG'],
            ['delete_start', 'AACCTTG', 1, 'AACCTTGG'],
        ],
    }

    vcf_lines = vcf._write_vcf_lines(1, 1, alt_map[(1, 1)]['alt_edits'], len(df))
    # Biallelic: one line per deletion, each with its own REF
    assert len(vcf_lines) == 7
    # Check each line has the expected REF/ALT pair
    all_lines = '\n'.join(vcf_lines)
    assert 'AA\tA' in all_lines
    assert 'AAC\tC' in all_lines
    assert 'AACC\tC' in all_lines
    assert 'AACCT\tT' in all_lines
    assert 'AACCTT\tT' in all_lines
    assert 'AACCTTG\tG' in all_lines
    assert 'AACCTTGG\tG' in all_lines


def test_build_alt_map_multi_deletion_start_and_middle():
    ref1 = 'AACCTTGG'
    aln1 = '-ACCTTGG'

    ref2 = 'AACCTTGG'
    aln2 = '--CCTTGG'

    ref3 = 'AACCTTGG'
    aln3 = 'A--CTTGG'

    ref4 = 'AACCTTGG'
    aln4 = '----TTGG'

    ref5 = 'AACCTTGG'
    aln5 = '-----TGG'

    ref6 = 'AACCTTGG'
    aln6 = '------GG'

    ref7 = 'AACCTTGG'
    aln7 = '-------G'

    df = create_df_alleles(
        (ref1, aln1),
        (ref2, aln2),
        (ref3, aln3),
        (ref3, aln3),
        (ref4, aln4),
        (ref5, aln5),
        (ref6, aln6),
        (ref7, aln7),
        (ref7, aln7),
    )

    amplicon_positions = {"Reference": (1, 1)}
    alt_map = vcf.build_alt_map(df, amplicon_positions)

    assert alt_map[(1, 1)] == {
        'alt_edits': [
            ['delete_start', 'A', 1, 'AA'],
            ['delete_start', 'AA', 1, 'AAC'],
            ['delete', 'AC', 2, 'AAC'],
            ['delete_start', 'AACC', 1, 'AACCT'],
            ['delete_start', 'AACCT', 1, 'AACCTT'],
            ['delete_start', 'AACCTT', 1, 'AACCTTG'],
            ['delete_start', 'AACCTTG', 2, 'AACCTTGG'],
        ],
    }

    vcf_lines = vcf._write_vcf_lines(1, 1, alt_map[(1, 1)]['alt_edits'], len(df))
    # Biallelic: one line per edit
    assert len(vcf_lines) == 7
    all_lines = '\n'.join(vcf_lines)
    assert 'AA\tA' in all_lines       # delete_start 'A'
    assert 'AAC\tC' in all_lines       # delete_start 'AA'
    assert 'AAC\tA' in all_lines       # delete 'AC' (middle deletion)
    assert 'AACCT\tT' in all_lines     # delete_start 'AACC'
    assert 'AACCTT\tT' in all_lines    # delete_start 'AACCT'
    assert 'AACCTTG\tG' in all_lines   # delete_start 'AACCTT'
    assert 'AACCTTGG\tG' in all_lines  # delete_start 'AACCTTG'


# ----------------------------- insertions -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # single insertion of "GG" after index 30 -> key (1, 32)
        (
            [make_ins(30, "GG", reads=1)],
            {"Reference": (1, 1)},
            {(1, 31): {"alt_edits": [["insert", "GG", 1, REF_SEQ[30]]]}},
        ),
        # merge same inserted string and key (reads sum)
        (
            [make_ins(30, "GG", reads=2), make_ins(30, "GG", reads=3)],
            {"Reference": (1, 1)},
            {(1, 31): {"alt_edits": [["insert", "GG", 5, REF_SEQ[30]]]}},
        ),
        # split different inserted strings at same key
        (
            [make_ins(30, "GG", reads=1), make_ins(30, "T", reads=1)],
            {"Reference": (1, 1)},
            {(1, 31): {"alt_edits": [["insert", "GG", 1, REF_SEQ[30]], ["insert", "T", 1, REF_SEQ[30]]]}},
        ),
    ],
    ids=["ins_single", "ins_merge_same_seq", "ins_split_diff_seq"],
)
def test_build_alt_map_insertions(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_alt_map(df, amplicon_positions)
    assert _normalize_alt_map(out) == _normalize_alt_map(expected)


def test_upsert_edit_del_and_ins():
    alt_map = {
        ('chrX', 10): {
            'alt_edits': [['delete', 'T', 2, 'AT']]
        },
    }
    vcf._upsert_edit(alt_map, ('chrX', 10), 'AT', 'insert', 'TT', 2)
    assert alt_map[('chrX', 10)] == {'alt_edits': [['delete', 'T', 2, 'AT'], ['insert', 'TT', 2, 'AT']]}


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
                (1, 3): {"alt_edits": [["sub", "T", 1, REF_SEQ[2]]]},
                # insertion on chrom 2: left_index = 1000 + 0 = 1000
                (2, 1000): {"alt_edits": [["insert", "G", 2, REF_SEQ[0]]]},
            },
        ),
    ],
    ids=["two_amplicons_offsets"],
)
def test_build_alt_map_multi_amplicon_and_offsets(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_alt_map(df, amplicon_positions)
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
    df = create_df_alleles(*rows)
    out = vcf.build_alt_map(df, amplicon_positions)
    assert len(out) == expected_size

# ----------------------------- ensuring that insertions are keyed by right anchor -----------------------------

def test_build_alt_map_fidelity_like_real_row():
    # Short 40bp reference used across tests
    ref_seq = "ATGCGTACGA---TCGTACGTAGCTAGCTAGCGTAGCTAGCTA"
    aln_seq = "ATGCGTACGAGGGTCGTACGTAGCTAGCTAGCGTAGCTAGCTA"
    num_reads = 7
    ins_right_anchor = 9

    df = create_df_alleles((ref_seq, aln_seq, num_reads))

    # Amplicon starts at absolute position 500 on chromosome 1
    amplicon_positions = {"Reference": (1, 500)}

    # Expected: key uses the right anchor; ref_seq is the anchor base; alt carries inserted string & reads
    expected = {
        (1, 500 + ins_right_anchor): {
            "alt_edits": [["insert", "GGG", num_reads, ref_seq[ins_right_anchor]]],
        }
    }

    out = vcf.build_alt_map(df, amplicon_positions)

    def _normalize(m):
        return {
            k: {"alt_edits": sorted(v["alt_edits"], key=lambda t: (t[0], t[1]))}
            for k, v in m.items()
        }

    assert _normalize(out) == _normalize(expected)


def test_aln_to_alt_map_ins_del_same_pos():
    ref1 = 'AATGCGTAC'
    aln1 = 'AA--CGTAC'

    ref2 = 'AA--TGCGTAC'
    aln2 = 'AATTTGCGTAC'

    df = create_df_alleles((ref1, aln1), (ref2, aln2))
    amplicon_positions = {"Reference": (1, 1)}
    alt_map = vcf.build_alt_map(df, amplicon_positions)

    assert list(alt_map.keys()) == [(1, 2)]
    assert alt_map[(1, 2)] == {'alt_edits': [['delete', 'TG', 1, 'ATG'], ['insert', 'TT', 1, 'A']]}

    temp_vcf_path = 'aln_to_alt_map_ins_del_same_pos.vcf'
    num_reads = 5
    num_vcf_rows = vcf.vcf_lines_from_alt_map(alt_map, num_reads, {'Reference': 9}, temp_vcf_path)
    # Biallelic: 2 lines (one for deletion, one for insertion)
    assert num_vcf_rows == 2

    with open(temp_vcf_path) as fh:
        vcf_contents = fh.read()
    assert 'ATG\tA' in vcf_contents      # deletion
    assert 'A\tATT' in vcf_contents      # insertion

    os.remove(temp_vcf_path)


def test_aln_to_alt_map_ins_then_del():
    ref1 = 'AATTT---GCGTAC'
    aln1 = 'AATTTCCCGCG---'

    df = create_df_alleles((ref1, aln1))
    amplicon_positions = {'Reference': (1, 1)}
    alt_map = vcf.build_alt_map(df, amplicon_positions)

    assert list(alt_map.keys()) == [(1, 8), (1, 5)]
    assert alt_map[(1, 5)] == {'alt_edits': [['insert', 'CCC', 1, 'T']]}
    assert alt_map[(1, 8)] == {'alt_edits': [['delete', 'TAC', 1, 'GTAC']]}


def test_aln_to_alt_map_to_vcf():
    ref1 = 'AATGCGTAC'
    aln1 = 'AATGCG-AC'
    #             ^ Interested in this deletion across each of these examples

    ref2 = 'AA--TGCGTAC'
    aln2 = 'AAGGTGCG-AC'

    ref3 = 'AATGCGTAC'
    aln3 = 'AA-GCG-AC'

    ref4 = 'AATG--CGTAC'
    aln4 = 'AA-GAACG-AC'

    ref5 = 'AA--TGCGTAC'
    aln5 = 'AAGGTGCG-GG'

    ref6 = 'AATGCGTAC'
    aln6 = '-------AC'

    df = create_df_alleles(
        (ref1, aln1),
        (ref2, aln2),
        (ref3, aln3),
        (ref4, aln4),
        (ref5, aln5),
        (ref6, aln6),
    )

    amplicon_positions = {"Reference": (1, 1)}
    alt_map = vcf.build_alt_map(df, amplicon_positions)

    # deletion at position 7 (1-based) in each example above, should occur 5 times
    assert alt_map[(1, 6)] == {'alt_edits': [['delete', 'T', 5, 'GT']]}
    # insertion of GG occurs 2 times in 2, 5 and deletion of T occurs 2 times in 3, 4
    assert alt_map[(1, 2)] == {'alt_edits': [['insert', 'GG', 2, 'A'], ['delete', 'T', 2, 'AT']]}
    # insertion of AA occurs 1 time in 4
    assert alt_map[(1, 4)] == {'alt_edits': [['insert', 'AA', 1, 'G']]}
    # substitution of A -> G occurs 1 time in 5
    assert alt_map[(1, 8)] == {'alt_edits': [['sub', 'G', 1, 'A']]}
    # substitution of C -> G occurs 1 time in 5
    assert alt_map[(1, 9)] == {'alt_edits': [['sub', 'G', 1, 'C']]}
    # deletion at position 1 (1-based) occurs 1 time in 6
    assert alt_map[(1, 1)] == {'alt_edits': [['delete_start', 'AATGCGT', 1, 'AATGCGTA']]}

    temp_vcf_path = 'aln_to_alt_map_to_vcf.vcf'
    num_reads = 5
    num_vcf_rows = vcf.vcf_lines_from_alt_map(alt_map, num_reads, {'Reference': 9}, temp_vcf_path)
    # Biallelic: each edit is its own line. Position (1,2) has 2 edits, rest have 1 each = 7 total
    assert num_vcf_rows == 7

    with open(temp_vcf_path, 'r') as fh:
        vcf_contents = fh.read()

    # Biallelic: each edit on its own line with its own REF
    assert '\t'.join(('1', '6', '.', 'GT', 'G', '.', 'PASS', f'AF={5 / num_reads:.3f}')) in vcf_contents
    # Position 2: deletion and insertion are separate lines
    assert '\t'.join(('1', '2', '.', 'AT', 'A', '.', 'PASS', f'AF={2 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '2', '.', 'A', 'AGG', '.', 'PASS', f'AF={2 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '4', '.', 'G', 'GAA', '.', 'PASS', f'AF={1 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '8', '.', 'A', 'G', '.', 'PASS', f'AF={1 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '9', '.', 'C', 'G', '.', 'PASS', f'AF={1 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '1', '.', 'AATGCGTA', 'A', '.', 'PASS', f'AF={1 / num_reads:.3f}')) in vcf_contents

    os.remove(temp_vcf_path)
