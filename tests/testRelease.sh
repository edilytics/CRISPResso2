set -e
echo Running CRISPResso
CRISPResso -r1 FANC.Cas9.fastq -a CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG -g GGAATCCCTTCTGCAGCACC --debug &> CRISPResso_on_FANC.Cas9.log
echo Running CRISPResso with parameters
CRISPResso -r1 FANC.Cas9.fastq -a CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG -g GGAATCCCTTCTGCAGCACC -e CGGCCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCTGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG -c GGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTT --dump -qwc 20-30_45-50 -q 30 --default_min_aln_score 80 -an FANC -n params --debug &> CRISPResso_on_params.log
echo Running CRISPRessoBatch
CRISPRessoBatch -bs FANC.batch -a CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG -g GGAATCCCTTCTGCAGCACC -p 2 --debug --base_editor --debug &> CRISPRessoBatch_on_FANC.log
echo Running CRISPRessoPooled
CRISPRessoPooled -r1 Both.Cas9.fastq -f Cas9.amplicons.txt -p 2 --keep_intermediate --min_reads_to_use_region 100 --debug &> CRISPRessoPooled_on_Both.Cas9.log

echo TESTING CRISPRESSO2
diff CRISPResso_on_FANC.Cas9/Reference.nucleotide_frequency_table.txt expectedResults/CRISPResso_on_FANC.Cas9/Reference.nucleotide_frequency_table.txt
diff CRISPResso_on_FANC.Cas9/CRISPResso_quantification_of_editing_frequency.txt tests/expectedResults/CRISPResso_on_FANC.Cas9/CRISPResso_quantification_of_editing_frequency.txt

echo TESTING CRISPRESSO2 PARAMS
diff CRISPResso_on_params/FANC.nucleotide_frequency_table.txt expectedResults/CRISPResso_on_params/FANC.nucleotide_frequency_table.txt
diff CRISPResso_on_params/CRISPResso_quantification_of_editing_frequency.txt tests/expectedResults/CRISPResso_on_params/CRISPResso_quantification_of_editing_frequency.txt

echo TESTING BATCH
diff CRISPRessoBatch_on_FANC/Reference.MODIFICATION_FREQUENCY_SUMMARY.txt tests/expectedResults/CRISPRessoBatch_on_FANC/Reference.MODIFICATION_FREQUENCY_SUMMARY.txt

echo TESTING POOLED
diff CRISPRessoPooled_on_Both.Cas9/SAMPLES_QUANTIFICATION_SUMMARY.txt tests/expectedResults/CRISPRessoPooled_on_Both.Cas9/SAMPLES_QUANTIFICATION_SUMMARY.txt

echo Finished