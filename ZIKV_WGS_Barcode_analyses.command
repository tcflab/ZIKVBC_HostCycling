#!/bin/sh

## Set the sample name variable
echo "What is the name of the sample?"
read sample

## Change directory to the output directory for rep 1
cd PATH/TO/${sample}_rep1_v3.1

## Align processed paired-end reads to reference with default settings
bwa mem -t 4 PATH/TO/REFERENCE.fasta ./normalized_fastq_bam/norm_${sample}_1_r1.fastq ./normalized_fastq_bam/norm_${sample}_1_r2.fastq > relaxed_norm_${sample}_rep1.bam

## sort and index alignment
samtools sort ./relaxed_norm_${sample}_rep1.bam > sorted_relaxed_norm_${sample}_rep1.bam
samtools index ./sorted_relaxed_norm_${sample}_rep1.bam

## extract aligned reads overlapping the barcode region
samtools view -h ./sorted_relaxed_norm_${sample}_rep1.bam "KU501215_1:4007-4030" > BCregion_sorted_relaxed_norm_${sample}_rep1.bam

## convert aligned reads to fastq
bamToFastq -i ./BCregion_sorted_relaxed_norm_${sample}_rep1.bam -fq BCregion_sorted_relaxed_${sample}_rep1.fastq

## convert fastq to fasta
seqtk seq -a ./BCregion_sorted_relaxed_${sample}_rep1.fastq > BCregion_sorted_relaxed_${sample}_rep1.fasta

## extract just the barcode sequence from all the reads
sed -n '/^>/p; /CT.GC.GC.CT.AC.CC.CT.GC./p' ./BCregion_sorted_relaxed_${sample}_rep1.fasta > sed_BCregion_sorted_relaxed_${sample}_rep1.fasta
sed -e '$!N;/^>.*\n>/D' -e 'P;D' ./sed_BCregion_sorted_relaxed_${sample}_rep1.fasta > sedfilter_BCregion_sorted_relaxed_${sample}_rep1.fasta
sed 's/^.*\(CT.GC.GC.CT.AC.CC.CT.GC.\).*$/\1/' ./sedfilter_BCregion_sorted_relaxed_${sample}_rep1.fasta > BConly_${sample}_rep1.fasta

## Count kmers
PATH/TO/kmercountexact.sh in=./BConly_${sample}_rep1.fasta out=kmer_counts_${sample}_rep1.fasta k=24 khist=kmer_freq_hist_${sample}_rep1.tsv fastadump=t rcomp=f ow=t

## Generate tsv for barcode counts
seqkit fx2tab ./kmer_counts_${sample}_rep1.fasta -H > counts_BConly_${sample}_rep1.tsv

## Clean up
mkdir ./barcode_analyses
mv *relaxed*.fastq ./barcode_analyses
mv *.fasta ./barcode_analyses
mv *relaxed*.bam ./barcode_analyses
mv *relaxed*.bam.bai ./barcode_analyses
mv *.tsv ./barcode_analyses


## Change directory to the output directory for rep 2
cd PATH/TO/${sample}_rep2_v3.1

## Align processed paired-end reads to reference with default settings
bwa mem -t 4 PATH/TO/REFERENCE.fasta ./normalized_fastq_bam/norm_${sample}_2_r1.fastq ./normalized_fastq_bam/norm_${sample}_2_r2.fastq > relaxed_norm_${sample}_rep2.bam

## sort and index alignment
samtools sort ./relaxed_norm_${sample}_rep2.bam > sorted_relaxed_norm_${sample}_rep2.bam
samtools index ./sorted_relaxed_norm_${sample}_rep2.bam

## extract aligned reads overlapping the barcode region
samtools view -h ./sorted_relaxed_norm_${sample}_rep2.bam "KU501215_1:4007-4030" > BCregion_sorted_relaxed_norm_${sample}_rep2.bam

## convert aligned reads to fastq
bamToFastq -i ./BCregion_sorted_relaxed_norm_${sample}_rep2.bam -fq BCregion_sorted_relaxed_${sample}_rep2.fastq

## convert fastq to fasta
seqtk seq -a ./BCregion_sorted_relaxed_${sample}_rep2.fastq > BCregion_sorted_relaxed_${sample}_rep2.fasta

## extract just the barcode sequence from all the reads
sed -n '/^>/p; /CT.GC.GC.CT.AC.CC.CT.GC./p' ./BCregion_sorted_relaxed_${sample}_rep2.fasta > sed_BCregion_sorted_relaxed_${sample}_rep2.fasta
sed -e '$!N;/^>.*\n>/D' -e 'P;D' ./sed_BCregion_sorted_relaxed_${sample}_rep2.fasta > sedfilter_BCregion_sorted_relaxed_${sample}_rep2.fasta
sed 's/^.*\(CT.GC.GC.CT.AC.CC.CT.GC.\).*$/\1/' ./sedfilter_BCregion_sorted_relaxed_${sample}_rep2.fasta > BConly_${sample}_rep2.fasta

## Count kmers
PATH/TO/kmercountexact.sh in=./BConly_${sample}_rep2.fasta out=kmer_counts_${sample}_rep2.fasta k=24 khist=kmer_freq_hist_${sample}_rep2.tsv fastadump=t rcomp=f ow=t

## Generate tsv for barcode counts
seqkit fx2tab ./kmer_counts_${sample}_rep2.fasta -H > counts_BConly_${sample}_rep2.tsv

## Clean up
mkdir ./barcode_analyses
mv *relaxed*.fastq ./barcode_analyses
mv *.fasta ./barcode_analyses
mv *relaxed*.bam ./barcode_analyses
mv *relaxed*.bam.bai ./barcode_analyses
mv *.tsv ./barcode_analyses

exit 1