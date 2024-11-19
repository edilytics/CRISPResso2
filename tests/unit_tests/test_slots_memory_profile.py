from CRISPResso2 import CRISPRessoCORE, CRISPResso2Align, CRISPRessoShared, CRISPRessoCOREResources

ALN_MATRIX = CRISPResso2Align.read_matrix('./CRISPResso2/EDNAFULL')

# ResultsSlotsDict is returned by CRISPRessoCOREResources.find_indels_substitutions()
# CRISPRessoCOREResources.find_indels_substitutions() is called in CRISPRessoCORE.get_new_variant_object()
# CRISPRessoCORE.get_new_variant_object() is called in CRISPRessoCORE.process_fastq() and CRISPRessoCORE.process_bam()
# I'll want to start with single-process runs

def test_process_fastq():
    inputs = {
        'fastq_filename': '',
        'variantCache': '', 
        'ref_names': '', 
        'refs': '', # building this will be a pain
        'args': '', 
        'files_to_remove': '', 
        'output_directory': '', 
        'fastq_write_out': False,
    }
    CRISPRessoCORE.process_fastq(**inputs)