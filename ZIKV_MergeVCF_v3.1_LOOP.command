#!/bin/sh

outputdir=$1
sample=$2

cd $outputdir/${sample}_rep1_v3.1/reference_aligned
gunzip -k -d ./lofreq_reference_${sample}_1.vcf.gz

cd $outputdir/${sample}_rep2_v3.1/reference_aligned
gunzip -k -d ./lofreq_reference_${sample}_2.vcf.gz

## Merge VCFs by taking average of each SNP
Rscript PATH/TO/VCFmerge_v1.1.R ${outputdir} ${sample}
echo "### COMPLETED VCF MERGE ###"

cd $outputdir

echo '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">' | cat - ./rep1and2_refalign${sample}.vcf > temp && mv temp ./rep1and2_refalign${sample}.vcf
echo '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">' | cat - ./rep1and2_conalign${sample}.vcf > temp && mv temp ./rep1and2_conalign${sample}.vcf

echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">' | cat - ./rep1and2_refalign${sample}.vcf > temp && mv temp ./rep1and2_refalign${sample}.vcf
echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">' | cat - ./rep1and2_conalign${sample}.vcf > temp && mv temp ./rep1and2_conalign${sample}.vcf

echo '##fileDate=20200224' | cat - ./rep1and2_refalign${sample}.vcf > temp && mv temp ./rep1and2_refalign${sample}.vcf
echo '##fileDate=20200224' | cat - ./rep1and2_conalign${sample}.vcf > temp && mv temp ./rep1and2_conalign${sample}.vcf

echo '##fileformat=VCFv4.0' | cat - ./rep1and2_refalign${sample}.vcf > temp && mv temp ./rep1and2_refalign${sample}.vcf
echo '##fileformat=VCFv4.0' | cat - ./rep1and2_conalign${sample}.vcf > temp && mv temp ./rep1and2_conalign${sample}.vcf

exit 1