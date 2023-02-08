"""Unit tests for CRISPResso2CORE."""

from pytest_check import check

from CRISPResso2 import CRISPRessoCORE

def test_get_consensus_alignment_from_pairs():
    """Tests for generating consensus alignments from paired reads."""

    print("testing Easy")

    #basic test
    qual1                =   "AAAA"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----TCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "AAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "NNCGATCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 80)

    #test quality difference
    qual1                =   "AAAB"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----GCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "AAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "NNCGATCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 80)

    #test quality difference
    qual1                =   "AAAA"
    aln1_seq             = "--CGAT----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-----GCGAT"
    aln2_ref             = "ATCGATCGAT"
    qual2                =      "BAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "NNCGAGCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 70)

    #gaps between r1 and r2
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "--CGA-----"
    aln1_ref             = "ATCGATCGAT"
    aln2_seq             = "-------GA-"
    aln2_ref             = "ATCGATCGAT"
    qual2                = "AAAAAAAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "NNCGANNGAN")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 50)

    print('Finished easy tests... now for the hard stuff')

    #insertion in r1
    qual1                =   "AAAA"
    aln1_seq             = "--CCGA-----".replace(" ","") #added replace for vertical alignment
    aln1_ref             = "ATC-GATCGAT".replace(" ","")
    aln2_seq             = "--- ----GA-".replace(" ","")
    aln2_ref             = "ATC GATCGAT".replace(" ","")
    qual2                =         "AA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "NNCCGANNGAN")
    check.equal(ref_seq, "ATC-GATCGAT")
    check.equal(score, 45) #double check this score... should be 5/11

    #deletion in r1
    qual1                =   "AA"
    aln1_seq             = "--C-A-----".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "-------GA-".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                =        "AA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "NNC-ANNGAN")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 40) #double check this score... should be 4/10

    # deletion in r2
    qual1                =   "AAA"
    aln1_seq             = "--CGA-----".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "-----T-GA-".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                =      "AAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "NNCGAT-GAN")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 60) #double check this score... should be 6/10

    # insertion at beginning of r1
    qual1                = "AAAA"
    aln1_seq             = "TA-CGA----- ".replace(" ","")
    aln1_ref             = "-ATCGATCGAT ".replace(" ","")
    aln2_seq             = " --------AT ".replace(" ","")
    aln2_ref             = " ATCGATCGAT".replace(" ","")
    qual2                =      "AAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "TA-CGANNNAT")
    check.equal(ref_seq, "-ATCGATCGAT")
    check.equal(score, 54) #double check this score... should be 6/11

    # alternating qualities
    qual1                = "BABABABABA"
    aln1_seq             = "ACCAACCAAT".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "TTGGTTGGTT".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                = "ABABABABAB"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "ATCGATCGAT")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 100)

    # large insertion in r1
    qual1                = "AAAAAA"
    aln1_seq             = "ACGTGA---------".replace(" ","")
    aln1_ref             = "A-----TCGATCGAT".replace(" ","")
    aln2_seq             = "------CGAT".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                =       "AAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "ACGTGANNNNNCGAT")
    check.equal(ref_seq, "A-----TCGATCGAT")
    check.equal(score, 33)

    # large insertion in r2
    qual1                = "AAAAA"
    aln1_seq             = "ATCGA-----".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "-----TTAGCT---".replace(" ","")
    aln2_ref             = "ATCGAT---C-GAT".replace(" ","")
    qual2                =      "AAAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "ATCGATTAGCTNNN")
    check.equal(ref_seq, "ATCGAT---C-GAT") #TODO: Is this right? ATCGATCCCCGGAT
    check.equal(score, 50)

    # Conflicts with reference
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "TAGCTAGCTA".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "TAGCTAGCTA".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                = "AAAAAAAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "TAGCTAGCTA")
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 0)

    # Conflicts between reads
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "TAGCTAGCTA".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "ATCGATCGAT".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                = "AAAAAAAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "TAGCTAGCTA") #Should it take r1 or the one that's like the reference sequence? - For Kendell
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 0)

    # Alternating reads
    qual1                = "AAAAAAAAAA"
    aln1_seq             = "AT--AT--AT".replace(" ","")
    aln1_ref             = "ATCGATCGAT".replace(" ","")
    aln2_seq             = "--CG--CG--".replace(" ","")
    aln2_ref             = "ATCGATCGAT".replace(" ","")
    qual2                = "AAAAAAAAAA"

    aln_seq, ref_seq, score = CRISPRessoCORE.get_consensus_alignment_from_pairs(aln1_seq, aln1_ref, qual1, aln2_seq, aln2_ref, qual2)
    check.equal(aln_seq, "ATCGATCGAT") #TODO: Failure returns AT-AT-AT
    check.equal(ref_seq, "ATCGATCGAT")
    check.equal(score, 100) #TODO: See above, score 60



if __name__ == "__main__":
# execute only if run as a script
    test_get_consensus_alignment_from_pairs()
