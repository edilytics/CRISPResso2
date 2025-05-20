"""Unit tests for CRISPResso2CORE."""
import pytest
import pandas as pd


from CRISPResso2 import CRISPRessoCORE, CRISPRessoShared

def test_get_consensus_alignment_from_pairs():
    """Tests for generating consensus alignments from paired reads."""
    try:
        CRISPRessoCORE.get_consensus_alignment_from_pairs
    except AttributeError:
        pytest.xfail('get_consensus_alignment_from_pairs is not implemented yet!')

    print("testing Easy")

    #basic test
    qual1                =   "AAAA"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----TCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "AAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNCGATCGAT"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 80

    #test quality difference
    qual1                =   "AAAB"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----GCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "AAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNCGATCGAT"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 80

    #test quality difference
    qual1                =   "AAAA"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----GCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "BAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    print('got aln_seq ' + str(aln_seq))
    assert aln_seq ==      "NNCGAGCGAT"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 70

    #gaps between r1 and r2
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "--CGA-----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-------GA-"
    aln2_ref             = "ATCGATCGAT"
    qual2                = "AAAAAAAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNCGANNGAN"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 50

    print('Finished easy tests... now for the hard stuff')

    #insertion in r1
    qual1                =   "AAAA"
    aln1_seq             = "--CCGA-----".replace(" ","") #added replace for vertical alignment
    aln1_ref             = "ATC-GATCGAT".replace(" ","")
    aln2_seq             = "--- ----GA-".replace(" ","")
    aln2_ref             = "ATC GATCGAT".replace(" ","")
    qual2                =         "AA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNCCGANNGAN"
    assert ref_seq ==      "ATC-GATCGAT"
    assert score == 45 #double check this score... should be 5/11

    #deletion in r1
    qual1                =   "AA"
    aln1_seq             = "--C-A-----".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "-------GA-".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                =        "AA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    assert aln_seq ==      "NNC-ANNGAN"
    assert ref_seq ==      "ATCGATCGAT"
    assert score == 50 #double check this score... should be 5/10


def test_split_quant_window_coordinates_single():
    assert [(5, 10)] == CRISPRessoCORE.split_quant_window_coordinates('5-10')


def test_split_quant_window_coordinates_multiple():
    assert CRISPRessoCORE.split_quant_window_coordinates('2-5_10-12') == [(2, 5), (10, 12)]


def test_split_quant_window_coordinates_error():
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('a-5')


def test_split_quant_window_coordinates_empty():
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('_')


def test_split_quant_window_coordinates_partially_empty():
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('1-3_')


def test_split_quant_window_coordinates_blank():
    with pytest.raises(CRISPRessoShared.BadParameterException):
        CRISPRessoCORE.split_quant_window_coordinates('')


def test_get_include_idxs_from_quant_window_coordinates():
    quant_window_coordinates = '1-10_12-20'
    assert CRISPRessoCORE.get_include_idxs_from_quant_window_coordinates(quant_window_coordinates) == [*list(range(1, 11)), *list(range(12, 21))]


def test_get_cloned_include_idxs_from_quant_window_coordinates():
    quant_window_coordinates = '1-10_12-20'
    s1inds = list(range(22))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 11)), *list(range(12, 21))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_beginning():
    quant_window_coordinates = '1-10_12-20'
    # represents a 5bp insertion at the beginning (left)
    s1inds = list(range(5, 27))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(6, 16)), *list(range(17, 26))]

def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_beginning():
    quant_window_coordinates = '1-10_12-20'
    # represents a 5bp deletion at the beginning (left)
    s1inds = [-1, -1, -1, -1, -1 ] + list(range(26))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(0, 6)), *list(range(7, 16))]

def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion():
    quant_window_coordinates = '10-20_35-40'
    # represents a 7bp deletion in the middle
    s1inds = list(range(23)) + [22, 22, 22, 22, 22, 22, 22] + list(range(23, 34))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(10, 21)), *list(range(35-7, 41-7))]

def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_modified():
    quant_window_coordinates = '10-25_35-40'
    # represents a 7bp deletion in the middle, where part of the QW is deleted
    # [0, 1, 3, 4, ... , 21, 22, 22, 22, 22, 22, 22, 22, 22, 23, 24, ... , 33]
    s1inds = list(range(23)) + [22, 22, 22, 22, 22, 22, 22] + list(range(23, 34))
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(10, 23)), *list(range(35-7, 41-7))]


def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_end_modified(): 
    # 5 bp deletion at end of 20 bp sequence
    quant_window_coordinates = '1-5_10-20'
    s1inds = [*list(range(16)), *[15, 15, 15, 15, 15]]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 6)), *list(range(10, 16))]

def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_and_deletion():
    # 5 bp deletion and 5 bp insertion
    quant_window_coordinates = '1-5_10-20'
    s1inds = [0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 6, 7, 8, 9, 15, 16, 17, 18, 19, 20]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*list(range(1, 6)), *[6, 7, 8, 9, 15, 16, 17, 18, 19, 20]]

def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_and_deletion_modified():
    quant_window_coordinates = '1-5_10-20'
    s1inds = [0, 1, 2, 2, 4, 5, 6, 7, 7, 7, 7, 7, 7, 8, 9, 10, 15, 16, 17, 18, 19]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [*[1,2,4,5], *[8, 9, 10, 15, 16, 17, 18, 19]]

def test_get_cloned_include_idxs_from_quant_window_coordinates_insertion_across_qw():
    # 6 bp insertion in middle of 4 bp sequence
    quant_window_coordinates = '1-4'
    s1inds = [0,1,2,9,10]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1,2,9,10]

def test_get_cloned_include_idxs_from_quant_window_coordinates_deletion_entire_qw():
    # 5 bp deletion of entire qw
    quant_window_coordinates = '1-4_7-10'
    s1inds = [0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [1, 2, 3, 4]

def test_get_cloned_include_idxs_from_quant_window_coordinates_include_zero():
    quant_window_coordinates = '0-5'
    s1inds = [0, 1, 2, 3, 4, 5]
    assert CRISPRessoCORE.get_cloned_include_idxs_from_quant_window_coordinates(quant_window_coordinates, s1inds) == [0, 1, 2, 3, 4, 5]


# Testing parallelization functions
def test_regular_input():
    # Test with typical input
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(100, 4) == [0, 25, 50, 75, 100]

def test_remainder_input():
#     # Test with typical input
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(101, 4) == [0, 25, 50, 75, 101]

def test_similar_num_reads_input():
#     # Test with typical input
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(11, 10) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11]

def test_large_similar_num_reads_input():
#     # Test with typical input
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(101, 100) == list(range(0, 100)) + [101]

def test_more_processes_than_reads():
#     # Test with typical input
    # assert CRISPRessoCORE.get_variant_cache_equal_boundaries(3, 5) 
    # assert that an exception is raised
    with pytest.raises(Exception):
        CRISPRessoCORE.get_variant_cache_equal_boundaries(3, 5)

def test_single_process():
    # Test with a single process
    assert CRISPRessoCORE.get_variant_cache_equal_boundaries(50, 1) == [0, 50]

def test_zero_sequences():
    # Test with zero unique sequences
    with pytest.raises(Exception):
        CRISPRessoCORE.get_variant_cache_equal_boundaries(0, 3)

def test_large_numbers():
    # Test with large number of processes and sequences
    boundaries = CRISPRessoCORE.get_variant_cache_equal_boundaries(10000, 10)
    assert len(boundaries) == 11  # Check that there are 11 boundaries

def test_sublist_generation():
    n_processes = 4
    unique_reads = 100
    mock_variant_cache = [i for i in range(unique_reads)]
    assert len(mock_variant_cache) == 100
    boundaries = CRISPRessoCORE.get_variant_cache_equal_boundaries(unique_reads, n_processes) 
    assert boundaries == [0, 25, 50, 75, 100]
    sublists = []
    for i in range(n_processes):
        left_sublist_index = boundaries[i]
        right_sublist_index = boundaries[i+1]
        sublist = mock_variant_cache[left_sublist_index:right_sublist_index]
        sublists.append(sublist)
    assert [len(sublist) for sublist in sublists] == [25, 25, 25, 25]
    assert [s for sublist in sublists for s in sublist] == mock_variant_cache

def test_irregular_sublist_generation():
    n_processes = 4
    unique_reads = 113
    mock_variant_cache = [i for i in range(unique_reads)]
    assert len(mock_variant_cache) == 113
    boundaries = CRISPRessoCORE.get_variant_cache_equal_boundaries(unique_reads, n_processes) 
    # assert boundaries == [0, 25, 50, 75, 100]
    sublists = []
    for i in range(n_processes):
        left_sublist_index = boundaries[i]
        right_sublist_index = boundaries[i+1]
        sublist = mock_variant_cache[left_sublist_index:right_sublist_index]
        sublists.append(sublist)
    assert [len(sublist) for sublist in sublists] == [28,28,28,29]
    assert [s for sublist in sublists for s in sublist] == mock_variant_cache


def _make_df_alleles():
    """Construct a mock df_alleles DataFrame covering sub, del, ins edits."""
    ref_seq = "ATGCGTACGATCGTACGTAGCTAGCTAGCGTAGCTAGCTA"  # 40 bp
    ref_positions = list(range(len(ref_seq)))

    rows = []

    # ───────────── Unmodified allele (skipped by build_alt_map) ──────────────
    rows.append({
        "#Reads": 2,
        "Aligned_Sequence": ref_seq,
        "Reference_Sequence": ref_seq,
        "n_inserted": 0,
        "n_deleted": 0,
        "n_mutated": 0,
        "Reference_Name": "Reference",
        "ref_positions": ref_positions,
        "insertion_coordinates": [],
        "deletion_coordinates": [],
        "substitution_positions": [],
    })

    # ───────────── Single‑nucleotide substitution at index 9 (A→G) ──────────
    sub_seq = ref_seq[:9] + "G" + ref_seq[10:]
    rows.append({
        "#Reads": 1,
        "Aligned_Sequence": sub_seq,
        "Reference_Sequence": ref_seq,
        "n_inserted": 0,
        "n_deleted": 0,
        "n_mutated": 1,
        "Reference_Name": "Reference",
        "ref_positions": ref_positions,
        "insertion_coordinates": [],
        "deletion_coordinates": [],
        "substitution_positions": [9],
    })

    # ───────────── 1‑bp deletion removing G at index 19 ─────────────────────
    del_seq = ref_seq[:19] + ref_seq[20:]
    rows.append({
        "#Reads": 1,
        "Aligned_Sequence": del_seq,
        "Reference_Sequence": ref_seq,
        "n_inserted": 0,
        "n_deleted": 1,
        "n_mutated": 0,
        "Reference_Name": "Reference",
        "ref_positions": ref_positions,
        "insertion_coordinates": [],
        "deletion_coordinates": [(19, 20)],  # [start, end] (0‑based, end exclusive)
        "substitution_positions": [],
    })

    # ───────────── 2‑bp insertion "GG" after index 30 (between 30 and 31) ──
    ins_seq = ref_seq[:31] + "GG" + ref_seq[31:]
    rows.append({
        "#Reads": 1,
        "Aligned_Sequence": ins_seq,
        "Reference_Sequence": ref_seq,
        "n_inserted": 2,
        "n_deleted": 0,
        "n_mutated": 0,
        "Reference_Name": "Reference",
        "ref_positions": ref_positions,
        "insertion_coordinates": [(31, 33)],  # slice storing inserted sequence
        "deletion_coordinates": [],
        "substitution_positions": [],
    })
   # ───────────── 10‑bp deletion removing G... at index 19 ─────────────────────
    del_seq = ref_seq[:19] + ref_seq[29:]
    rows.append({
        "#Reads": 1,
        "Aligned_Sequence": del_seq,
        "Reference_Sequence": ref_seq,
        "n_inserted": 0,
        "n_deleted": 1,
        "n_mutated": 0,
        "Reference_Name": "Reference",
        "ref_positions": ref_positions,
        "insertion_coordinates": [],
        "deletion_coordinates": [(19, 29)],  # [start, end] (0‑based, end exclusive)
        "substitution_positions": [],
    })

    return pd.DataFrame(rows)


def _normalize_alt_map(alt_map):
    """Sort alt_seqs so order does not affect equality checks."""
    norm = {}
    for key, val in alt_map.items():
        norm[key] = {
            "ref_seq": val["ref_seq"],
            "alt_seqs": sorted(val["alt_seqs"], key=lambda t: (t[0], t[1])),
        }
    return norm


def test_build_alt_map():
    df_alleles = _make_df_alleles()

    # Chromosome 1, reference starts at genomic coordinate 1
    amplicon_positions = {"Reference": (1, 1)}

    alt_map = CRISPRessoCORE.build_alt_map(df_alleles, amplicon_positions)

    expected_alt_map = {
        # SNP at 1 + 9 = 10
        (1,10): {
            "ref_seq": "A",  # reference base
            "alt_seqs": [["sub", "G", 1]],
        },
        # 1‑bp deletion at 1 + 19 = 20
        (1,20): {
            "ref_seq": "AGCTAGCTAGC",  # flanking + deleted base (A|G)
            "alt_seqs": [["delete", "G", 1], ['delete', 'GCTAGCTAGC', 1]],  # content is not used beyond length
        },
        # 2‑bp insertion after index 30 -> coordinate 32
        (1,32): {
            "ref_seq": "G",  # base before insertion
            "alt_seqs": [["insert", "GG", 1]],
        },
    }
    if _normalize_alt_map(alt_map) != _normalize_alt_map(expected_alt_map):
        print(_normalize_alt_map(alt_map))
        print(_normalize_alt_map(expected_alt_map))
    assert _normalize_alt_map(alt_map) == _normalize_alt_map(expected_alt_map)


def test_write_vcf_from_alt_map():
    df_alleles = _make_df_alleles()
    amplicon_positions = {"Reference": (1, 1)}

    alt_map = CRISPRessoCORE.build_alt_map(df_alleles, amplicon_positions)
    num_reads = df_alleles["#Reads"].sum()

    vcf_text = CRISPRessoCORE.vcf_text_from_alt_map(alt_map, num_reads, ["Reference"])
    print("vcf_text:")
    print(vcf_text)

    expected_vcf = (
        "##fileformat=VCFv4.5\n"
        "##source=CRISPResso2\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tReference\n"
        "1\t10\t.\tA\tG\t.\tPASS\tAF=0.200\n"
        "1\t20\t.\tAGCTAGCTAGC\tA\t.\tPASS\tAF=0.200\n"
        "1\t32\t.\tG\tGGG\t.\tPASS\tAF=0.200\n"
    )

    assert vcf_text == expected_vcf
    
if __name__ == "__main__":
# execute only if run as a script
    test_get_consensus_alignment_from_pairs()
