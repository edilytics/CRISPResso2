set -e
echo Running CRISPRessoWGS
CRISPRessoWGS -b Both.Cas9.fastq.smallGenome.bam -r smallGenome/smallGenome.fa -f Cas9.regions.txt --debug &> CRISPRessoWGS_on_Both.Cas9.fastq.smallGenome.log

echo TESTING WGS
diff -w CRISPRessoWGS_on_Both.Cas9.fastq.smallGenome/SAMPLES_QUANTIFICATION_SUMMARY.txt expectedResults/CRISPRessoWGS_on_Both.Cas9.fastq.smallGenome/SAMPLES_QUANTIFICATION_SUMMARY.txt

echo Running CRISPRessoCompare
CRISPRessoCompare CRISPRessoBatch_on_FANC/CRISPResso_on_Cas9/ CRISPRessoBatch_on_FANC/CRISPResso_on_Untreated/ --debug &> CRISPRessoCompare.log

echo TESTING COMPARE
diff -w CRISPRessoCompare_on_Cas9_VS_Untreated/Deletions_quantification.txt expectedResults/CRISPRessoCompare_on_Cas9_VS_Untreated/Deletions_quantification.txt

echo Finished
