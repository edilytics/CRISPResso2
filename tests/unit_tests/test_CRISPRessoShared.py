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
    to_sequence = 'TTC-T'
    from_sequence = 'TTCGT'
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates(to_sequence, from_sequence)
    assert s1inds_gap_left == [0, 1, 2, 2, 3] # should this be [0, 1, 2, 2, 4]?
    assert s1inds_gap_right == [0, 1, 2, 3, 3] # should this be [0, 1, 2, 4, 4]?
    assert from_sequence[0] == to_sequence[s1inds_gap_left[0]]
    assert from_sequence[1] == to_sequence[s1inds_gap_left[1]]
    assert from_sequence[2] == to_sequence[s1inds_gap_left[2]]
    assert from_sequence[2] == to_sequence[s1inds_gap_left[2]]
    # assert from_sequence[4] == to_sequence[s1inds_gap_left[4]] # this should accept, right?


def test_get_relative_coordinates_start_gap():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('--CGT', 'TTCGT')
    assert s1inds_gap_left == [-1, -1, 0, 1, 2] # should this be [-1, -1, 2, 3, 4]?
    assert s1inds_gap_right == [0, 0, 0, 1, 2] # should this be [2, 2, 2, 3, 4]?


def test_get_relative_coordinates_from_gap():
    s1inds_gap_left, s1inds_gap_right = CRISPRessoShared.get_relative_coordinates('ATCGT', 'ATC-T')
    assert s1inds_gap_left == [0, 1, 2, 4]
    assert s1inds_gap_right == [0, 1, 2, 4]
