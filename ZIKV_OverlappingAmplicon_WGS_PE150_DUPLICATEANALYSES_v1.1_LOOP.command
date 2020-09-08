#!/bin/sh

outputdir=$1
sample=$2

reference="/PATH/TO/REFERENCE.fasta"

cd $outputdir

mkdir ./${sample}_combined

## Bfzip and index combined vcf
bgzip --threads 8 ./Duplicate_refalign_vcf/rep1and2_refalign${sample}.vcf
tabix -p vcf ./Duplicate_refalign_vcf/rep1and2_refalign${sample}.vcf.gz
echo "### COMPLETED BGZIP AND TABIX ###"

gunzip -k -d ./Duplicate_refalign_vcf/rep1and2_refalign${sample}.vcf.gz

## Extract consensus sequence for population
cat ${reference} | bcftools consensus -e "INFO/AF<0.50" ./Duplicate_refalign_vcf/rep1and2_refalign${sample}.vcf.gz > ./${sample}_combined/consensus_combined_${sample}.fasta
echo "### COMPLETED CONSENSUS EXTRACTION ###"

## Use SNPdat to annotate variants relative to reference
perl PATH/TO/SNPdat_v1.0.5.pl -i ./Duplicate_refalign_vcf/rep1and2_refalign${sample}.vcf -g /Volumes/HD1/GTF_files/ZIKVBC2_WGS_noBC.txt -f ${reference} -o ./${sample}_combined/snpdat_combined_refaligned_${sample}.tsv
echo "### COMPLETED SNPDAT ###"

## Use SNPdat to annotate variants relative to consensus
perl PATH/TO/SNPdat_v1.0.5.pl -i ./Duplicate_conalign_vcf/rep1and2_conalign${sample}.vcf -g /Volumes/HD1/GTF_files/ZIKVBC2_WGS_noBC.txt -f ./${sample}_combined/consensus_combined_${sample}.fasta -o ./${sample}_combined/snpdat_combined_conaligned_${sample}.tsv
echo "### COMPLETED SNPDAT ###"

## Use SNPgenie to characterize variants and calculate pi
PATH/TO/snpgenie.pl --minfreq=0.003 --snpreport=./Duplicate_conalign_vcf/rep1and2_conalign${sample}.vcf --vcfformat=2 --fastafile=./${sample}_combined/consensus_combined_${sample}.fasta --gtffile=PATH/TO/GTF_FILE.txt
echo "### COMPLETED SNPGENIE ###"

mv ./SNPGenie_Results ./${sample}_combined/SNPGenie_Results_${sample}_combined

## Calculate diversity metrics with R
Rscript PATH/TO/ZIKV_WGS_duplicatevcf_Diversity_dinuc_v1_1.R ${outputdir} ${sample}
echo "### COMPLETED DIVERSITY METRICS ###"

exit 1