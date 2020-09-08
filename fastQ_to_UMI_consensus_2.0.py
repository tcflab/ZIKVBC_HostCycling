#!/usr/bin/env python
# coding: utf-8

import sys

### Group fastQ sequences by header and create one consensus sequence for each quality UMI. ###
print('Grouping sequences by UMI in separate fasta files...')

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
import os 
import glob

directory = sys.argv[1] + '/grouped_reads_passing_cutoff_' + sys.argv[2] + '_' + sys.argv[3] + '.fasta.split/'

files_in_direct=glob.glob(directory+"*.fasta")

print('Sequences grouped...')
print("Script has detected "+str(len(files_in_direct))+" files in the directory "+directory)

print('Generating consensus sequences for each UMI group...')

for file in files_in_direct:
    align=AlignIO.read(file, "fasta")
    summary_align=AlignInfo.SummaryInfo(align)
    consensus=summary_align.dumb_consensus(threshold=0.5, ambiguous='N')
    str_con=str(consensus)
    
    filesplit = file.split("/")
    filesplitfurther=filesplit[-1].split("_")
    filesplitevenfurther=filesplitfurther[-1].split(".")
    ID=filesplitevenfurther[0]
    
    output_con= open(sys.argv[1] + '/consensus_' + sys.argv[2] + '_ZIKV_UMI.fasta',"a")
    
    output_con.write(">"+ID+"\n")
    output_con.write(str_con+"\n")
    
    output_con.close()
    
print('Consensus sequences generated in single fasta file...')