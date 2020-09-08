# ZIKVBC_HostCycling
Repository for scripts to process, analyze, and visualize NGS data from barcoded ZIKV. Scripts are listed in the order in which they should be run, since latter scripts will require outputs from preceding scripts.

### ZIKV Barcode Analyses
1. UMIconsensus_and_barcode_counts_2.0_12ntUMI.command

   Executable bash pipeline script to extract and count barcodes from ZIKV barcode amplicons from UMI sequencing.
   In addition to dependencies listed below, requires fastQ_to_UMI_list_12ntUMI.py, fastQ_to_UMI_consensus_2.0.py, and UMI_filter_list.R.
     
2. ZIKV_WGS_Barcode_analyses.command

   Executable bash pipeline script to extract and count barcodes from ZIKV barcodes sequenced by whole-genome sequencing.
   Required dependencies are listed below.
     
3. Mean_bc_frequencies.R

   R script to combine barcode frequencies from technical duplicate sequencing libraries.
     
4. Euclidean_distance_meanbcfreqs.R

   R script to calculate Euclidean distance between barcode populations of two samples.
     
### ZIKV Whole-Genome Sequencing Analyses
1. ZIKV_OverlappingAmplicon_WGS_PE150_v3.1_LOOP2.command

   Executable bash pipeline script to process paired-end Illumina reads, align to reference, call and annotate SNVs, and calculate diversity metrics from single technical replicates.
   In addition to dependencies listed below, requires ZIKV_WGS_pipeline_Diversity_dinuc_v1_1.R.
   Should be run with the loop command described below using the associated ArgList file.
     
2. ZIKV_MergeVCF_v3.1_LOOP.command

   Executable bash script to combine VCF files from technical duplicate sequencing libraries.
   Required dependencies are listed below.
   Should be run with the loop command described below using the associated ArgList file.
     
3. ZIKV_OverlappingAmplicon_WGS_PE150_DUPLICATEANALYSES_v1.1_LOOP.command

   Executable bash pipeline script to calculate diversity metrics from combined VCF files.
   In addition to dependencies listed below, requires ZIKV_WGS_duplicatevcf_Diversity_dinuc_v1_1.R.
   Should be run with the loop command described below using the associated ArgList file.

4. SNVs_overtime.R

   R script to collate and generate figures for SNVs frequencies over sequential passages.
     
5. Shared_SNVs.R
  
   R script to define and generate figures for convergent SNV trajectories over sequential passages.

6. Species_adaptive_SNVs.R

   R script to define and generate figures for species-specific, dual species, and alternate-specific SNV trajectories over sequential passages.
     
7. Coverage.R

     R script to collate and generate figures for read coverage depth.
     
### Loop Command
The bash scripts for WGS are written to be run as a loop on a list of samples. To run as a loop, the associated ArgList files will need to be updated with the arguments for each sample in the list.

The number of arguments (-n 4) should match the number of arguments on each line of the ArgList file.
```bash
Xargs -n 4 ~/Path/To/Pipeline.command <ArgList.txt
```

### Dependencies
Fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ \
Trimmomatic http://www.usadellab.org/cms/?page=trimmomatic \
BBTools Suite https://jgi.doe.gov/data-and-tools/bbtools/ \
Fastp https://github.com/OpenGene/fastp \
Python3 with the following packages https://www.python.org/downloads/
  * biopython https://biopython.org/
  * NumPy https://numpy.org/
  * pandas https://pandas.pydata.org/
  * seaborn https://seaborn.pydata.org/#
  * matplotlib https://matplotlib.org/

Umi-tools https://umi-tools.readthedocs.io/en/latest/ \
Burrows-Wheeler Aligner http://bio-bwa.sourceforge.net/ \
Samtools http://samtools.sourceforge.net/ \
R https://www.r-project.org/ \
seqtk https://github.com/lh3/seqtk \
seqkit https://bioinf.shenwei.me/seqkit/ \
cutadapt https://cutadapt.readthedocs.io/en/stable/ \
Bedtools https://bedtools.readthedocs.io/en/latest/index.html \
Lofreq https://csb5.github.io/lofreq/ \
Snpdat https://github.com/agdoran/snpdat \
SNPGenie https://github.com/chasewnelson/SNPGenie \
bcftools https://samtools.github.io/bcftools/ \
pysamstats https://github.com/alimanfoo/pysamstats
