import os
from types import SimpleNamespace

import pytest

from CRISPResso2 import CRISPRessoShared
from CRISPResso2.writers import vcf


# ----------------------------- write_vcf_from_edits -----------------------------

def test_write_vcf_from_edits_header_and_ordering_no_samples():
    """Tests that write_vcf_from_edits produces the expected header and ordering of records."""
    edit_counts = {
        ("1", 20, "AG", "A"): 2,   # deletion
        ("1", 10, "A", "T"): 3,    # substitution
    }
    temp_vcf_path = os.path.join(os.path.dirname(__file__), "temp_test_no_samples.vcf")
    count = vcf.write_vcf_from_edits(edit_counts, num_reads=10, amplicon_lens={}, vcf_path=temp_vcf_path)

    with open(temp_vcf_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    assert lines[0].startswith("##fileformat=VCFv4.")
    assert lines[1] == "##source=CRISPResso2"
    assert lines[2].startswith("##INFO=<ID=AF")
    assert lines[3].startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    # Records in ascending chrom/pos order
    assert lines[4] == "1\t10\t.\tA\tT\t.\tPASS\tAF=0.300"
    assert lines[5] == "1\t20\t.\tAG\tA\t.\tPASS\tAF=0.200"

    assert len(lines) == 6
    assert count == 2
    os.remove(temp_vcf_path)


def test_write_vcf_from_edits_header_and_ordering_with_contigs():
    """Test that write_vcf_from_edits produces contig headers and the expected record."""
    edit_counts = {
        ("X", 2, "A", "AG"): 2,   # insertion
    }
    temp_vcf_path = os.path.join(os.path.dirname(__file__), "temp_test_with_contigs.vcf")
    count = vcf.write_vcf_from_edits(
        edit_counts, num_reads=10,
        amplicon_lens={"Ref1": 100, "Ref2": 200},
        vcf_path=temp_vcf_path,
    )

    with open(temp_vcf_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    assert lines[3] == "##contig=<ID=Ref1,length=100>"
    assert lines[4] == "##contig=<ID=Ref2,length=200>"
    assert lines[5] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    assert lines[6] == "X\t2\t.\tA\tAG\t.\tPASS\tAF=0.200"
    assert len(lines) == 7
    assert count == 1
    os.remove(temp_vcf_path)


def test_write_vcf_from_edits_bad_denominator_raises():
    for bad_num_reads in [0, None, -5]:
        with pytest.raises(CRISPRessoShared.BadParameterException):
            vcf.write_vcf_from_edits(
                {("1", 10, "A", "T"): 1},
                bad_num_reads,
                {},
                "/dev/null",
            )


# ----------------------------- write_vcf_file (orchestrator) -----------------------------

def _stub_edit_counts():
    return {
        ("1", 10, "A", "T"): 3,
        ("1", 20, "AG", "A"): 2,
    }

@pytest.mark.parametrize(
    "ref_names,ref_lens,coords_str,expected_record_count,expected_contig_lines",
    [
        (["RefA"], {"RefA": 100}, "1:1", 2, 1),
        (["RefA", "RefB"], {"RefA": 100, "RefB": 200}, "1:1,2:1", 2, 2),
    ],
)
def test_write_vcf_file_smoke(tmp_path, monkeypatch, ref_names, ref_lens, coords_str, expected_record_count, expected_contig_lines):
    monkeypatch.setattr(vcf, "build_edit_counts", lambda df, amplicon_positions: _stub_edit_counts())

    import pandas as pd
    df = pd.DataFrame([{"#Reads": 7}, {"#Reads": 3}])  # num_reads = 10

    args = SimpleNamespace(amplicon_coordinates=coords_str)
    out_path = tmp_path / "test.vcf"

    num_vcf_rows = vcf.write_vcf_file(df, ref_names, ref_lens, args, str(out_path))

    with open(out_path, "r") as f:
        lines = [line.strip() for line in f]

    assert num_vcf_rows == 2

    expected_header_lines = 4 + expected_contig_lines
    assert len(lines) == expected_header_lines + expected_record_count

    assert lines[expected_header_lines].startswith("1\t10\t.\tA\tT\t.\tPASS\tAF=0.300")
    assert lines[expected_header_lines + 1].startswith("1\t20\t.\tAG\tA\t.\tPASS\tAF=0.200")
