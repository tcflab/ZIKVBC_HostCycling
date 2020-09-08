#!/bin/sh

rawfastqdir=$1
outputdir=$2
sample=$3
rep=$4

cd $outputdir

read1lane1=`find ${rawfastqdir} -name "${sample}_${rep}*L001_R1_001.fastq.gz" ! -name '*._*'`
read2lane1=`find ${rawfastqdir} -name "${sample}_${rep}*L001_R2_001.fastq.gz" ! -name '*._*'`
read1lane2=`find ${rawfastqdir} -name "${sample}_${rep}*L002_R1_001.fastq.gz" ! -name '*._*'`
read2lane2=`find ${rawfastqdir} -name "${sample}_${rep}*L002_R2_001.fastq.gz" ! -name '*._*'`

reference="PATH/TO/REFERENCE.fasta"

mkdir ./${sample}_rep${rep}_v3.1
cd ./${sample}_rep${rep}_v3.1

## Quality trim and adapter trim 
cutadapt -j 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o at_${sample}_${rep}_r1_L1.fastq -p at_${sample}_${rep}_r2_L1.fastq ${read1lane1} ${read2lane1}

cutadapt -j 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o at_${sample}_${rep}_r1_L2.fastq -p at_${sample}_${rep}_r2_L2.fastq ${read1lane2} ${read2lane2}

## Trim non-overlapping regions, filter reads with mismatches, and merge reads
PATH/TO/bbmerge.sh \
	in1=at_${sample}_${rep}_r1_L1.fastq \
		in2=at_${sample}_${rep}_r2_L1.fastq \
			out=merged_${sample}_${rep}_r1_L1.fastq \
				out2=merged_${sample}_${rep}_r2_L1.fastq \
					outu=unmerged_${sample}_${rep}_r1_L1.fastq \
						outu2=unmerged_${sample}_${rep}_r2_L1.fastq \
							minoverlap=50 ecco=t mix=f

PATH/TO/bbmerge.sh \
	in1=at_${sample}_${rep}_r1_L2.fastq \
		in2=at_${sample}_${rep}_r2_L2.fastq \
			out=merged_${sample}_${rep}_r1_L2.fastq \
				out2=merged_${sample}_${rep}_r2_L2.fastq \
					outu=unmerged_${sample}_${rep}_r1_L2.fastq \
						outu2=unmerged_${sample}_${rep}_r2_L2.fastq \
							minoverlap=50 ecco=t mix=f

## Concatenate reads from lane 1 and lane 2
cat merged_${sample}_${rep}_r1_L1.fastq merged_${sample}_${rep}_r1_L2.fastq > cat_${sample}_${rep}_r1.fastq

cat merged_${sample}_${rep}_r2_L1.fastq merged_${sample}_${rep}_r2_L2.fastq > cat_${sample}_${rep}_r2.fastq

## Quality trim and adapter trim 
cutadapt -q 35,35 -m 50 --max-n 0 -j 4 -o qtatlf_${sample}_${rep}_r1.fastq -p qtatlf_${sample}_${rep}_r2.fastq ./cat_${sample}_${rep}_r1.fastq ./cat_${sample}_${rep}_r2.fastq

## Primer trim reads
cutadapt -g XGACAGTTCGAGTTTGAAGCGAAAG -g XAAGAAAGATCTGGCTGCCATGC -g XAGATGACGTCGATTGTTGGTGC \
	-g XTCAGGTGCATAGGAGTCAGCAA -g XAGAACGTTAGTGGACAGAGGCT -g XTTGATTGTGAACCGAGGACAGG -g XTGAAGGGCGTGTCATACTCCTT \
	-g XGGGAGAAGAAGATCACCCACCA -g XGCCTTAGGGGGAGTGTTGATCT -g XACGGTCGTTGTGGGATCTGTAA -g XCAGCCGTTATTGGAACAGCTGT \
	-g XCACTAAGGTCCACGTGGAGGAA -g XTGGCAGTGCTGGTAGCTATGAT -g XCAATGGTTTTGCTTTGGCCTGG -g XCCCTAGCGAAGTACTCACAGCT \
	-g XGTGGCATGAACCCAATAGCCAT -g XGTGGTCCATGGAAGCTAGATGC -g XCTGTTGAGTGCTTCGAGCCTTC -g XTATGGATGAGGCCCACTTCACA \
	-g XGGCTGGAAAACGGGTCATACAG -g XAGAGACTGACGAAGACCATGCA -g XTGGACCAGACACGGAGAGAAAA -g XCGTCTTGATGAGGAACAAGGGC \
	-g XTAATGGGAAGGAGAGAGGAGGG -g XCCCTGACCCTAATAGTGGCCAT -g XACTGGAACTCCTCTACAGCCAC -g XAGTGCAAAGCTGAGATGGTTGG \
	-g XGGTGGGGGATTGGCTTGAAAAA -g XAGGATGTGAATCTCGGCTCTGG -g XAAAAGTGGACACTAGGGTGCCA -g XACAAGGGGAATTTGGAAAGGCC \
	-g XAAATGGAAAAAGGGCACAGGGC -g XCAAACGAATGGCAGTCAGTGGA -g XATTTCCACAGAAGGGACCTCCG -g XACCACCTGGGCTGAGAACATTA \
	-G XAGTATGCACTCCCACGTCTAGT -G XTGATTCCAACCAGGTTTGCGAC -G XTACGGTGACACAACCTCCATGT -G XGGAGCCATGAACTGACAGCATT \
	-G XTGTGCGTCCTTGAACTCTACCA -G XCCATCTGTCCCTGCGTACTGTA -G XCGCCTCCAACTGATCCAAAGTC -G XTTGACTGCTGCTGCCAATCTAC \
	-G XGAGTGGGCATTCCTTCAGTGTG -G XGTGGGACTTTGGCCATTCACAT -G XCCTGGGCCTTATCTCCATTCCA -G XTATCAGCGCCAGATGAGCTACA \
	-G XAGAGAGAGGAGCATAAACCCCC -G XTTTCCCATGTGATGTCACCTGC -G XTACACTCCATCTGTGGTCTCCC -G XGCTCCAATGTCCCCATCCTTTG \
	-G XCCTCTAAGGGCCTCCTCCATTT -G XTGGTGAGTTGGAGTCCGGAAAT -G XGCCATCAAGTATGACCGGCTTT -G XCCTTTGCTCCGTCCTAAGCTTG \
	-G XCTCCAAAAGCCGCTCCTCTTTT -G XATTCTGGCTGGCTCAATTTCCG -G XAAGTGGTCACTGCATGTTGGAC -G XTCTCCACTTGGGGGTCAATTGT \
	-G XCCTTCCATTTCTCTCCCAGGGT -G XACCAGGGCCTCCTTTTGTGTAT -G XATGTGTAGAGTTGCGGGAGAGT -G XGGGCCTCATAGCTTCCATGGTA \
	-G XATGCTGCATTGCTACGAACCTT -G XTAATCCCAGCCCTTCAACACCA -G XCGTAAGTGACAACTTGTCCGCT -G XTGTCCCATCCAGTTGAGGGTTT \
	-G XATCCACACTCTGTTCCACACCA -G XTGACTAGCAGGCCTGACAACAT -G XACCACTAGTCCCTCTTCTGGAG \
	-a CTTTCGCTTCAAACTCGAACTGTCX -a GCATGGCAGCCAGATCTTTCTTX -a GCACCAACAATCGACGTCATCTX -a TTGCTGACTCCTATGCACCTGAX \
	-a AGCCTCTGTCCACTAACGTTCTX -a CCTGTCCTCGGTTCACAATCAAX -a AAGGAGTATGACACGCCCTTCAX -a TGGTGGGTGATCTTCTTCTCCCX \
	-a AGATCAACACTCCCCCTAAGGCX -a TTACAGATCCCACAACGACCGTX -a ACAGCTGTTCCAATAACGGCTGX -a TTCCTCCACGTGGACCTTAGTGX \
	-a ATCATAGCTACCAGCACTGCCAX -a CCAGGCCAAAGCAAAACCATTGX -a AGCTGTGAGTACTTCGCTAGGGX -a ATGGCTATTGGGTTCATGCCACX \
	-a GCATCTAGCTTCCATGGACCACX -a GAAGGCTCGAAGCACTCAACAGX -a TGTGAAGTGGGCCTCATCCATAX -a CTGTATGACCCGTTTTCCAGCCX \
	-a TGCATGGTCTTCGTCAGTCTCTX -a TTTTCTCTCCGTGTCTGGTCCAX -a GCCCTTGTTCCTCATCAAGACGX -a CCCTCCTCTCTCCTTCCCATTAX \
	-a ATGGCCACTATTAGGGTCAGGGX -a GTGGCTGTAGAGGAGTTCCAGTX -a CCAACCATCTCAGCTTTGCACTX -a TTTTTCAAGCCAATCCCCCACCX \
	-a CCAGAGCCGAGATTCACATCCTX -a TGGCACCCTAGTGTCCACTTTTX -a GGCCTTTCCAAATTCCCCTTGTX -a GCCCTGTGCCCTTTTTCCATTTX \
	-a TCCACTGACTGCCATTCGTTTGX -a CGGAGGTCCCTTCTGTGGAAATX -a TAATGTTCTCAGCCCAGGTGGTX -A ACTAGACGTGGGAGTGCATACTX \
	-A GTCGCAAACCTGGTTGGAATCAX -A ACATGGAGGTTGTGTCACCGTAX -A AATGCTGTCAGTTCATGGCTCCX -A TGGTAGAGTTCAAGGACGCACAX \
	-A TACAGTACGCAGGGACAGATGGX -A GACTTTGGATCAGTTGGAGGCGX -A GTAGATTGGCAGCAGCAGTCAAX -A CACACTGAAGGAATGCCCACTCX \
	-A ATGTGAATGGCCAAAGTCCCACX -A TGGAATGGAGATAAGGCCCAGGX -A TGTAGCTCATCTGGCGCTGATAX -A GGGGGTTTATGCTCCTCTCTCTX \
	-A GCAGGTGACATCACATGGGAAAX -A GGGAGACCACAGATGGAGTGTAX -A CAAAGGATGGGGACATTGGAGCX -A AAATGGAGGAGGCCCTTAGAGGX \
	-A ATTTCCGGACTCCAACTCACCAX -A AAAGCCGGTCATACTTGATGGCX -A CAAGCTTAGGACGGAGCAAAGGX -A AAAAGAGGAGCGGCTTTTGGAGX \
	-A CGGAAATTGAGCCAGCCAGAATX -A GTCCAACATGCAGTGACCACTTX -A ACAATTGACCCCCAAGTGGAGAX -A ACCCTGGGAGAGAAATGGAAGGX \
	-A ATACACAAAAGGAGGCCCTGGTX -A ACTCTCCCGCAACTCTACACATX -A TACCATGGAAGCTATGAGGCCCX -A AAGGTTCGTAGCAATGCAGCATX \
	-A TGGTGTTGAAGGGCTGGGATTAX -A AGCGGACAAGTTGTCACTTACGX -A AAACCCTCAACTGGATGGGACAX -A TGGTGTGGAACAGAGTGTGGATX \
	-A ATGTTGTCAGGCCTGCTAGTCAX -A CTCCAGAAGAGGGACTAGTGGTX \
		-e 0.1 -O 6 -j 4 -o ptqtatlf_${sample}_${rep}_r1.fastq -p ptqtatlf_${sample}_${rep}_r2.fastq \
			./qtatlf_${sample}_${rep}_r1.fastq ./qtatlf_${sample}_${rep}_r2.fastq

## Normalize depth with BBnorm
PATH/TO/bbnorm.sh in1=./ptqtatlf_${sample}_${rep}_r1.fastq in2=./ptqtatlf_${sample}_${rep}_r2.fastq out1=norm_${sample}_${rep}_r1.fastq out2=norm_${sample}_${rep}_r2.fastq target=2000

## Align reads with BWA MEM
bwa mem -t 8 -B 10 ${reference} ./norm_${sample}_${rep}_r1.fastq ./norm_${sample}_${rep}_r2.fastq > merged_norm_${sample}_${rep}.bam

## Sort sam by coordinate
samtools sort -o ./sorted_norm_${sample}_${rep}.bam ./merged_norm_${sample}_${rep}.bam

## Call variants with lofreq
lofreq call -l PATH/TO/GTF_FILE.txt -q 30 -Q 30 -C 300 -f ${reference} -o lofreq_reference_${sample}_${rep}.vcf ./sorted_norm_${sample}_${rep}.bam
echo "### COMPLETED LOFREQ VS REFERENCE ###"

## Use SNPdat to annotate variants
perl PATH/TO/SNPdat_v1.0.5.pl -i ./lofreq_reference_${sample}_${rep}.vcf -g PATH/TO/GTF_FILE.txt -f ${reference} -o snpdat_variants_norm_${sample}_${rep}.tsv
echo "### COMPLETED SNPDAT ###"

## Bfzip and index lofreq vcf
bgzip ./lofreq_reference_${sample}_${rep}.vcf
tabix -p vcf ./lofreq_reference_${sample}_${rep}.vcf.gz
echo "### COMPLETED BGZIP AND TABIX ###"

## Extract consensus sequence for population
cat ${reference} | bcftools consensus -e "INFO/AF<0.50" ./lofreq_reference_${sample}_${rep}.vcf.gz > consensus_${sample}_${rep}.fasta
echo "### COMPLETED CONSENSUS EXTRACTION ###"

## Call variants relative to consensus with lofreq
lofreq call -l PATH/TO/GTF_FILE.txt -q 30 -Q 30 -C 300 -f ./consensus_${sample}_${rep}.fasta -o lofreq_consensus_${sample}_${rep}.vcf ./sorted_norm_${sample}_${rep}.bam
echo "### COMPLETED LOFREQ VS CONSENSUS ###"

## Use SNPgenie to characterize variants and calculate pi
PATH/TO/snpgenie.pl --minfreq=0.003 --snpreport=./lofreq_consensus_${sample}_${rep}.vcf --vcfformat=2 --fastafile=./consensus_${sample}_${rep}.fasta --gtffile=PATH/TO/GTF_FILE.txt
echo "### COMPLETED SNPGENIE ###"

## Generate nucleotide counts per position
samtools index -b ./sorted_norm_${sample}_${rep}.bam

pysamstats -f ./consensus_${sample}_${rep}.fasta -t variation -D 1000000 --format \
	csv --output ./ntcounts_consensus_${sample}_${rep}.csv ./sorted_norm_${sample}_${rep}.bam

pysamstats -f ${reference} -t variation -D 1000000 --format \
	csv --output ./ntcounts_reference_${sample}_${rep}.csv ./sorted_norm_${sample}_${rep}.bam
echo "### COMPLETED PYSAMSTATS ###"

## Calculate diversity metrics with R
Rscript PATH/TO/ZIKV_WGS_pipeline_Diversity_dinuc_v1_1.R ${outputdir} ${sample} ${rep}
echo "### COMPLETED DIVERSITY METRICS ###"

### Organize and clean up files generated by pipeline ###
mkdir ./normalized_fastq_bam
mv *norm*.fastq ./normalized_fastq_bam
mv sorted_norm*.bam* ./normalized_fastq_bam
rm *lf_*.fastq
rm cat*.fastq
rm *norm*.bam
mkdir ./consensus_aligned
mkdir ./reference_aligned
mv *consensus_${sample}* ./consensus_aligned
mv *reference_${sample}* ./reference_aligned
mv *.tsv ./reference_aligned
mv *.vcf* ./reference_aligned
rm unmerged*.fastq
rm merged*.fastq
rm at*.fastq

exit 1