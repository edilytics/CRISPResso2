import os
from types import SimpleNamespace

import pytest

# Adjust the import path if your module name is different
from CRISPResso2 import CRISPRessoShared
from CRISPResso2.writers import vcf

# ----------------------------- _alt_seq_from_edit -----------------------------

@pytest.mark.parametrize(
    "ref_seq,edit_type,alt_edit,expected",
    [
        # Insertions
        ("A",     "insert", "GG",  "AGG"),     # no trailing context
        ("AGC",   "insert", "TAGA",   "ATAGAGC"),    # trailing context preserved (ref_seq[1:])
        ("AGCGATTACA",   "insert", "T",   "ATGCGATTACA"),    # long ref_seq, trailing context preserved (ref_seq[1:])

        # Deletions (ref_seq = left_anchor + longest_deleted_span)
        ("AGCTAT",  "delete", "G",  "ACTAT"),      # shorter-than-longest deletion: keep remaining suffix
        ("AGCTAT",  "delete", "GCT", "AAT"),       # middle deletion at this locus: alt_edit removed center of ref_seq
        ("AGCTAT",  "delete", "GCTAT", "A"),       # longest deletion at this locus: ref is left anchor only

        # Substitutions
        ("A",     "sub",    "T",   "T"),
        ("AGC",   "sub",    "T",   "TGC"),
    ],
    ids=[
        "ins_one_base",
        "long_ins_trailing_context",
        "ins_long_ref_with_context",
        "del_one_base",
        "del_partial_removal",
        "del_long_ref_removal",
        "sub_single_base_no_context",
        "sub_single_base_with_context",
    ],
)
def test_alt_seq_from_edit_happy(ref_seq, edit_type, alt_edit, expected):
    """Tests whether _alt_seq_from_edit produces the expected ALT sequence given alt_edit from alt_map."""
    assert vcf._alt_seq_from_edit(ref_seq, edit_type, alt_edit) == expected


@pytest.mark.parametrize(
    "ref_seq,edit_type,alt_edit",
    [
        ("AGCT", "delete", "AGCT"),  # del length >= ref_len  -> error
        ("A",    "sub",    ""),      # sub must be length 1
        ("A",    "sub",    "GG"),    # sub can't be greater than 1
        ("A",    "foo",    "X"),     # unknown edit type
    ],
    ids=["del_too_long", "empty_sub", "too_long_sub", "unknown_type"],
)
def test_alt_seq_from_edit_errors(ref_seq, edit_type, alt_edit):
    """Tests whether _alt_seq_from_edit properly raises exceptions on invalid input."""
    with pytest.raises(CRISPRessoShared.BadParameterException):
        vcf._alt_seq_from_edit(ref_seq, edit_type, alt_edit)

# ----------------------------- _write_vcf_line -----------------------------

@pytest.mark.parametrize(
    "desc,chrom,pos,alt_edits,num_reads,expected_lines",
    [
        (
            "ins_and_sub_no_samples",
            "1", 100,
            [("insert", "GG", 2, "A"), ("sub", "T", 3, "A")],
            10,
            [
                "1\t100\t.\tA\tT\t.\tPASS\tAF=0.300",
                "1\t100\t.\tA\tAGG\t.\tPASS\tAF=0.200",
            ],
        ),
        (
            "mixed_with_deletion_context",
            "1", 200,
            [
                ("delete", "GCT", 5, "AGCT"),  # ALT=A
                ("delete", "GC",  2, "AGC"),    # ALT=A  (own ref_seq)
                ("insert", "T",   1, "A"),      # ALT=AT
                ("sub",    "T",   4, "A"),      # ALT=T
            ],
            10,
            [
                # Sorted by (len(alt), alt, len(ref), ref):
                "1\t200\t.\tAGC\tA\t.\tPASS\tAF=0.200",    # del GC: ref=AGC, alt=A
                "1\t200\t.\tAGCT\tA\t.\tPASS\tAF=0.500",   # del GCT: ref=AGCT, alt=A
                "1\t200\t.\tA\tT\t.\tPASS\tAF=0.400",      # sub T: ref=A, alt=T
                "1\t200\t.\tA\tAT\t.\tPASS\tAF=0.100",     # ins T: ref=A, alt=AT
            ],
        ),
        (
            "duplicate_alts_are_aggregated",
            "1", 150,
            [("insert", "GG", 2, "A"), ("insert", "GG", 3, "A")],  # same ALT twice -> 5 reads
            10,
            [
                "1\t150\t.\tA\tAGG\t.\tPASS\tAF=0.500",
            ],
        ),
    ],
    ids=lambda x: x if isinstance(x, str) else None,
)
def test_write_vcf_lines(desc, chrom, pos, alt_edits, num_reads, expected_lines):
    """Testing that _write_vcf_lines produces the expected biallelic VCF lines."""
    lines = vcf._write_vcf_lines(chrom, pos, alt_edits, num_reads)
    assert lines == expected_lines


@pytest.mark.parametrize("bad_num_reads", [0, None, -5])
def test_write_vcf_lines_bad_denominator_raises(bad_num_reads):
    with pytest.raises(CRISPRessoShared.BadParameterException):
        vcf._write_vcf_lines("1", 10, [("sub", "T", 1, "A")], bad_num_reads)


# ----------------------------- vcf_lines_from_alt_map -----------------------------

def test_vcf_lines_from_alt_map_header_and_ordering_no_samples():
    """Tests that vcf_lines_from_alt_map produces the expected header and ordering of records."""
    # Two edits on chrom "1", positions 10 then 20
    alt_map = {
        ("1", 20): {"alt_edits": [("delete", "G", 2, "AG")]},  # ALT=A
        ("1", 10): {"alt_edits": [("sub", "T", 3, "A")]},      # ALT=T
    }
    temp_vcf_path = os.path.join(os.path.dirname(__file__), "temp_test_no_samples.vcf")
    count = vcf.vcf_lines_from_alt_map(alt_map, num_reads=10, amplicon_lens={}, vcf_path=temp_vcf_path)
    # Header (4 lines prepended to all CRISPResso vcf files)
    lines = []
    with open(temp_vcf_path, "r") as f:
        for line in f:
            if line.strip() != "":
                lines.append(line.strip())

    assert lines[0].startswith("##fileformat=VCFv4.")
    assert lines[1] == "##source=CRISPResso2"
    assert lines[2].startswith("##INFO=<ID=AF")
    assert lines[3].startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    # Edits are listed in ascending chrom/pos order
    assert lines[4] == "1\t10\t.\tA\tT\t.\tPASS\tAF=0.300"
    assert lines[5] == "1\t20\t.\tAG\tA\t.\tPASS\tAF=0.200"

    # Total lines = 4 header + 2 records
    assert len(lines) == 6
    assert count == 2
    os.remove(temp_vcf_path)  # clean up


def test_vcf_lines_from_alt_map_header_and_ordering_with_contigs():
    """Test that vcf_lines_from_alt_map produces contig headers and the expected record."""
    alt_map = {
        ("X", 2): {"alt_edits": [("insert", "G", 2, "A")]},   # ALT=AG
    }
    temp_vcf_path = os.path.join(os.path.dirname(__file__), "temp_test_with_contigs.vcf")
    lines = []
    count = vcf.vcf_lines_from_alt_map(alt_map, num_reads=10, amplicon_lens={"Ref1": 100, "Ref2": 200}, vcf_path=temp_vcf_path)
    with open(temp_vcf_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    assert lines[3] == "##contig=<ID=Ref1,length=100>"
    assert lines[4] == "##contig=<ID=Ref2,length=200>"
    assert lines[5] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    assert lines[6] == "X\t2\t.\tA\tAG\t.\tPASS\tAF=0.200"
    assert len(lines) == 7
    assert count == 1
    os.remove(temp_vcf_path)


# ----------------------------- write_vcf_file (orchestrator) -----------------------------

def _stub_alt_map():
    # Two records to be emitted (biallelic: one line per edit)
    return {
        ("1", 10): {"alt_edits": [("sub", "T", 3, "A")]},
        ("1", 20): {"alt_edits": [("delete", "G", 2, "AG")]},
    }

@pytest.mark.parametrize(
    "ref_names,ref_lens,coords_str,expected_record_count,expected_contig_lines",
    [
        (["RefA"], {"RefA": 100}, "1:1", 2, 1),
        (["RefA", "RefB"], {"RefA": 100, "RefB": 200}, "1:1,2:1", 2, 2),  # 2 contig lines (different chroms)
    ],
)
def test_write_vcf_file_smoke(tmp_path, monkeypatch, ref_names, ref_lens, coords_str, expected_record_count, expected_contig_lines):
    # Monkeypatch build_alt_map so this test doesn't depend on upstream logic.
    monkeypatch.setattr(vcf, "build_alt_map", lambda df, amplicon_positions: _stub_alt_map())

    # Minimal df_alleles: only '#Reads' is used for denominator
    import pandas as pd
    df = pd.DataFrame([{"#Reads": 7}, {"#Reads": 3}])  # num_reads = 10

    args = SimpleNamespace(amplicon_coordinates=coords_str)
    out_path = tmp_path / "test.vcf"

    num_vcf_rows = vcf.write_vcf_file(df, ref_names, ref_lens, args, str(out_path))
    lines = []
    with open(out_path, "r") as f:
        for line in f:
            lines.append(line.strip())
    # Check: function returns number of lines written
    assert num_vcf_rows == 2

    # Expect 4 base header lines + contig lines + expected_record_count
    expected_header_lines = 4 + expected_contig_lines
    assert len(lines) == expected_header_lines + expected_record_count

    # First record is POS=10 (after header)
    assert lines[expected_header_lines].startswith("1\t10\t.\tA\tT\t.\tPASS\tAF=0.300")
    # Second record is POS=20
    assert lines[expected_header_lines + 1].startswith("1\t20\t.\tAG\tA\t.\tPASS\tAF=0.200")
