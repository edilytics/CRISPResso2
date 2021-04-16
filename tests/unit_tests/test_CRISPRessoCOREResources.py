from CRISPResso2 import CRISPRessoCOREResources

def test_find_indels_substitutions():
    payload=CRISPRessoCOREResources.find_indels_substitutions('CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG', 'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',[18,19])
    assert payload['insertion_n'] == 0
    assert payload['deletion_n'] == 0
    assert payload['substitution_n'] == 0

    #deletion outside of quantification window
    payload=CRISPRessoCOREResources.find_indels_substitutions('CATGGAATCCCTTCTGCA---CCTGGATCGCTTTTCCGAG', 'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',[21,22])
    assert payload['insertion_n'] == 0
    assert payload['substitution_n'] == 0
    assert payload['deletion_n'] == 0

    #deletion overlap quantification window
    payload=CRISPRessoCOREResources.find_indels_substitutions('CATGGAATCCCTTCTGCA---CCTGGATCGCTTTTCCGAG', 'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',[18,19])
    assert payload['insertion_n'] == 0
    assert payload['substitution_n'] == 0
    assert payload['deletion_n'] == 3


if __name__ == "__main__":
# execute only if run as a script
    test_find_indels_substitutions()
