#!/usr/bin/env python
# coding: utf-8

import sys

### Create a TSV file. Column 1 = UMI. Column 2 = fastqID. Quality_raw_reads = reads that were successfully paired, merged, and trimmed. ### 

# User must input location of UMI
# User must input the length of the expected amplicon. This should be the same number of base pairs that were input in cell #2
print('Generating UMI list for quality reads...')

from Bio import SeqIO
import csv

input_file = sys.argv[1] + '/lf_bbmerge_' + sys.argv[2] + '.fastq' 
fastq_sequences = SeqIO.parse(open(input_file),'fastq')

UMI_list = []

for fastq in fastq_sequences: 
    sequence = str(fastq.seq)
    fastqID = str(fastq.id)
    
    UMI = str(sequence[127:143])
    amplicon = str(sequence[:143])
    fastqID = str(fastqID[0:46])
    
    UMI_list.append([UMI, amplicon, fastqID])
    
with open('UMI_list_ZIKV_' + sys.argv[3] + '.tsv','w') as f:
    for i in UMI_list:
        f.write(i[0] + '\t' + i[1] + '\t' + i[2] + '\n')
        
with open('UMI_list_ZIKV_' + sys.argv[3] + '.tsv') as f: 
    quality_raw_reads = 0 
    for line in f: 
        quality_raw_reads += 1

print('UMI list of quality raw read generated...')        
print('Number of quality raw reads is:', quality_raw_reads)

### 3. Creating a dictionary. Key = UMI. Values = fastqIDs. ###
print('Creating UMI dictionary...')

from Bio import SeqIO

input_file = sys.argv[1] + '/lf_bbmerge_' + sys.argv[2] + '.fastq'
fastq_sequences = SeqIO.parse(open(input_file),'fastq')

UMI_dict = {}

for fastq in fastq_sequences: 
    sequence = str(fastq.seq)
    fastqID = str(fastq.id)
    
    fastqID = str(fastqID[0:46])
    UMI = str(sequence[127:143])
    
    if not UMI in UMI_dict:
        UMI_dict[UMI] = [fastqID]
    else:
        UMI_dict[UMI].append(fastqID)

print('UMI dictionary created...')
        
### Converts the UMI dictionary to a pandas dataframe. Index = UMI. Values = fastqID, which are listed in columns. ###
### The last column is equal to the number of reads per unique UMI. ###

print('Converting UMI dictionary to pandas dataframe...')

import pandas as pd
#from numpy import nan

UMI_DF = pd.DataFrame.from_dict(UMI_dict, orient='index')

UMI_DF['read_count'] = UMI_DF.count(axis=1)

UMI_DF

print('UMI dictionary converted...')

### Determine (m) maximum number of reads associated with any unique UMI. ###
### Apply m to read number cutoff model to determine (c) minimum number of reads per UMI in order for an UMI to be considered 'real' ### 

print('Calculating number of reads associated with UMI...')

import pandas as pd
from math import exp, expm1

UMI_DF = pd.DataFrame.from_dict(UMI_dict, orient='index')
UMI_DF['read_count'] = UMI_DF.count(axis=1)
m = UMI_DF['read_count'].max()

c = ((-1.24e-21)*(m**6)) + ((3.53e-17)*(m**5)) - ((3.9e-13)*(m**4)) + ((2.12e-9)*(m**3)) - ((6.06e-6)*(m**2)) + (.018*m) + 3.15

if c < 3:
    c2 = 102
else:
    c2 = c

print('Number of reads per UMI calculated...')

print('Maximum number of reads associated with any unique UMI, m =', m)
print('Minimum number of reads required per UMI to pass quality cutoff, c =', c2)

import csv

UMIreadsheader = ['Max_reads', 'Min_reads']
UMIreads = [m, c2]
csvfile = open ('Max_Min_reads_per_UMI.csv', 'w')
with csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(UMIreadsheader)
    writer.writerow(UMIreads)

### Generate distribution plot of reads per UMI ###
print('Generating distribution plot of reads per UMI...')

import numpy as np
import matplotlib.pyplot as pp
import pandas as pd
import seaborn

UMIs = pd.read_csv(('UMI_list_ZIKV_' + sys.argv[3] + '.tsv'), sep='\t', names = ['UMI', 'sequence', 'fastQ_ID'])

UMIs = UMIs.drop('fastQ_ID', axis=1) 
UMIs = UMIs.drop('sequence', axis=1)

UMIs['reads per unique UMI'] = UMIs.groupby('UMI')['UMI'].transform('count')

UMIs_sorted = UMIs.sort_values('reads per unique UMI')

UMIs_counts = UMIs_sorted.drop_duplicates('UMI')

UMIs_counts['UMIs per read count'] = UMIs_counts.groupby('reads per unique UMI')['reads per unique UMI'].transform('count')
UMIs_final = UMIs_counts.drop('UMI', axis=1)

dist_plot = UMIs_final.plot(x='reads per unique UMI', y='UMIs per read count', title="Read number distribution plot", kind='scatter', figsize=(15,5), legend=False, logy=True)

dist_plot.set_ylabel("Distinct UMIs per unique \n number of raw reads, logscale")
dist_plot.set_xlabel("Raw sequence reads per UMI")

fig = dist_plot.get_figure()
fig.savefig('UMI_list_ZIKV_' + sys.argv[3] + '_distribution_plot.pdf')

print('Distribution plot generated...')