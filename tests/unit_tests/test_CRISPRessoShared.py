from CRISPResso2 import CRISPResso2Align, CRISPRessoShared

ALN_MATRIX = CRISPResso2Align.read_matrix('./CRISPResso2/EDNAFULL')


def test_get_mismatches():
    mismatch_cords = CRISPRessoShared.get_mismatches(
        'ATTA',
        'ATTA',
        ALN_MATRIX,
        -5,
        -3,
    )
    assert len(mismatch_cords) == 0

    mismatch_cords = CRISPRessoShared.get_mismatches(
        'GCAGTGGGCGCGCTA',
        'CCCACTGAAGGCCC',
        ALN_MATRIX,
        -5,
        -3,
    )
    assert len(mismatch_cords) == 6


def test_get_relative_coordinates():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('ATCGT', 'TTCGT')
    assert s1inds_gap_left == [0, 1, 2, 3, 4]
    assert s1inds_gap_right == [0, 1, 2, 3, 4]


def test_get_relative_coordinates_to_gap():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('TTC-T', 'ATCGT')
    assert s1inds_gap_left == [0, 1, 2, 2, 3]
    assert s1inds_gap_right == [0, 1, 2, 3, 3]


def test_get_relative_coordinates_from_gap():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('ATCGT', 'TTC-T')
    assert s1inds_gap_left == [0, 1, 2, 4]
    assert s1inds_gap_right == [0, 1, 2, 4]
