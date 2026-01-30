from CRISPResso2 import CRISPResso2Align, CRISPRessoShared
import tempfile
import os
import gzip

ALN_MATRIX = CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")


def test_get_mismatches():
    """Test get_mismatches function."""
    mismatch_cords = CRISPRessoShared.get_mismatches(
        "ATTA",
        "ATTA",
        ALN_MATRIX,
        -5,
        -3,
    )
    assert len(mismatch_cords) == 0

    mismatch_cords = CRISPRessoShared.get_mismatches(
        "GCAGTGGGCGCGCTA",
        "CCCACTGAAGGCCC",
        ALN_MATRIX,
        -5,
        -3,
    )
    expected_mismatch_count = 6
    assert len(mismatch_cords) == expected_mismatch_count


def test_get_relative_coordinates():
    """Test get_relative_coordinates function."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("ATCGT", "TTCGT")
    assert s1inds_gap_left == [0, 1, 2, 3, 4]
    assert s1inds_gap_right == [0, 1, 2, 3, 4]


def test_get_n_reads_fastq():
    """Test get_n_reads_fastq function."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA\n")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 1

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_n_reads_fastq_gzip():
    """Test get_n_reads_fastq function with gzipped file."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA\n")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # gzip file
    with open(f.name, "rb") as f_in, gzip.open(f.name + ".gz", "wb") as f_out:
        f_out.writelines(f_in)

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name + ".gz") == 1

    # Clean up: delete the file after the test
    os.remove(f.name)
    os.remove(f.name + ".gz")


def test_get_n_reads_fastq_three_extra_newlines():
    """Test get_n_reads_fastq with three extra newlines."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA\n")  # Ensure the file content is correct and ends with a newline
        f.write("\n\n")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 1

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_n_reads_fastq_four_extra_newlines():
    """Test get_n_reads_fastq with four extra newlines."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA\n")  # Ensure the file content is correct and ends with a newline
        f.write("\n\n\n\n\n\n\n\n")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 1

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_n_reads_fastq_100_reads():
    """Test get_n_reads_fastq with 100 reads."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        for _i in range(100):
            f.write("@SEQ_ID\n")
            f.write("GATTACA\n")
            f.write("+\n")
            f.write("AAAAAAA\n")
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    expected_read_count = 100
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == expected_read_count

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_n_reads_fastq_no_newline():
    """Test get_n_reads_fastq with no trailing newline."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        f.write("@SEQ_ID\n")
        f.write("GATTACA\n")
        f.write("+\n")
        f.write("AAAAAAA")  # Ensure the file content is correct and ends with a newline
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 1

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_n_reads_fastq_empty_file():
    """Test get_n_reads_fastq with empty file."""
    with tempfile.NamedTemporaryFile(encoding="utf-8", mode="w+", delete=False, suffix=".fastq") as f:
        f.flush()  # Flush the content to disk
        os.fsync(f.fileno())  # Ensure all internal buffers associated with the file are written to disk

    # Since the file needs to be read by another process, we ensure it's closed and written before passing it
    assert CRISPRessoShared.get_n_reads_fastq(f.name) == 0

    # Clean up: delete the file after the test
    os.remove(f.name)


def test_get_relative_coordinates_to_gap():
    """Test get_relative_coordinates with gap in to_sequence."""
    # unaligned sequences
    seq_1 = "TTCGT"
    seq_2 = "TTCT"

    # aligned_sequences
    to_sequence = "TTC-T"
    from_sequence = "TTCGT"

    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates(to_sequence, from_sequence)
    assert s1inds_gap_left == [0, 1, 2, 2, 3]
    assert s1inds_gap_right == [0, 1, 2, 3, 3]

    assert seq_1[0] == seq_2[s1inds_gap_left[0]]
    assert seq_1[1] == seq_2[s1inds_gap_left[1]]
    assert seq_1[2] == seq_2[s1inds_gap_left[2]]
    assert seq_1[4] == seq_2[s1inds_gap_left[4]]


def test_get_relative_coordinates_start_gap():
    """Test get_relative_coordinates with gap at start."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("--CGT", "TTCGT")
    assert s1inds_gap_left == [-1, -1, 0, 1, 2]
    assert s1inds_gap_right == [0, 0, 0, 1, 2]


def test_get_relative_coordinates_from_gap():
    """Test get_relative_coordinates with gap in from_sequence."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("ATCGT", "ATC-T")
    assert s1inds_gap_left == [0, 1, 2, 4]
    assert s1inds_gap_right == [0, 1, 2, 4]


def test_get_relative_coordinates_end_gap():
    """Test get_relative_coordinates with gap at end."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("ATC--", "ATCGT")
    assert s1inds_gap_left == [0, 1, 2, 2, 2]
    assert s1inds_gap_right == [0, 1, 2, 3, 3]


def test_get_relative_coordinates_multiple_gaps():
    """Test get_relative_coordinates with multiple gaps."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("AT--TC--G--CC", "ATCGTCGCGTTCC")
    assert s1inds_gap_left == [0, 1, 1, 1, 2, 3, 3, 3, 4, 4, 4, 5, 6]
    assert s1inds_gap_right == [0, 1, 2, 2, 2, 3, 4, 4, 4, 5, 5, 5, 6]


def test_get_relative_coordinates_ind_and_dels():
    """Test get_relative_coordinates with insertions and deletions."""
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates("ATG--C", "A-GCTC")
    assert s1inds_gap_left == [0, 2, 2, 2, 3]
    assert s1inds_gap_right == [0, 2, 3, 3, 3]


def test_get_quant_window_ranges_from_include_idxs():
    """Test get_quant_window_ranges_from_include_idxs function."""
    include_idxs = [0, 1, 2, 10, 11, 12]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(0, 2), (10, 12)]


def test_get_quant_window_ranges_from_include_idxs_empty():
    """Test get_quant_window_ranges_from_include_idxs with empty list."""
    include_idxs = []
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == []


def test_get_quant_window_ranges_from_include_idxs_single():
    """Test get_quant_window_ranges_from_include_idxs with single range."""
    include_idxs = [50, 51, 52, 53]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(50, 53)]


def test_get_quant_window_ranges_from_include_idxs_single_gap():
    """Test get_quant_window_ranges_from_include_idxs with single gap."""
    include_idxs = [50, 51, 52, 53, 55]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(50, 53), (55, 55)]


def test_get_quant_window_ranges_from_include_idxs_multiple_gaps():
    """Test get_quant_window_ranges_from_include_idxs with multiple gaps."""
    include_idxs = [50, 51, 52, 53, 55, 56, 57, 58, 60]
    assert CRISPRessoShared.get_quant_window_ranges_from_include_idxs(include_idxs) == [(50, 53), (55, 58), (60, 60)]


def test_get_silent_edits():
    """Test get_silent_edits function."""
    ref_seq = "AGS"
    seq = "AGS"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGC")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AGs"


def test_get_silent_edits_multiple_silent():
    """Test get_silent_edits with multiple silent edits."""
    ref_seq = "PHATKIDS"
    seq = "PHATKIDS"
    ref_codons = [("P", "CCT"), ("H", "CAT"), ("A", "GCT"), ("T", "ACT"), ("K", "AAA"), ("I", "ATT"), ("D", "GAT"), ("S", "AGC")]
    seq_codons = [("P", "CCC"), ("H", "CAC"), ("A", "GCT"), ("T", "ACT"), ("K", "AAA"), ("I", "ATT"), ("D", "GAC"), ("S", "AGT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "phATKIds"


def test_get_silent_edits_no_silent():
    """Test get_silent_edits with no silent edits."""
    ref_seq = "AGS"
    seq = "AGT"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT"), ("T", "ACT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AGT"


def test_get_silent_edits_deletion():
    """Test get_silent_edits with deletion."""
    ref_seq = "AGS"
    seq = "AG-"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AG-"


def test_get_silent_edits_multiple_deletions():
    """Test get_silent_edits with multiple deletions."""
    ref_seq = "AGS"
    seq = "A--"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "A--"


def test_get_silent_edits_insertion():
    """Test get_silent_edits with insertion."""
    ref_seq = "AGS-"
    seq = "AGST"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT"), ("T", "ACT")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AGST"


def test_get_silent_edits_middle_insertion():
    """Test get_silent_edits with insertion in middle."""
    ref_seq = "AG-S"
    seq = "AGTS"
    ref_codons = [("A", "GCT"), ("G", "GGT"), ("S", "AGT")]
    seq_codons = [("A", "GCT"), ("G", "GGT"), ("T", "ACT"), ("S", "AGC")]

    silent_edits = CRISPRessoShared.get_silent_edits(ref_seq, ref_codons, seq, seq_codons)
    assert silent_edits == "AGTs"
