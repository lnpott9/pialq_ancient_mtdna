#!/usr/bin/env python3

# =============================================================================
# Title: 5a_merge_consensus.py
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: Python script to generate hapogroup calls from vcf files
# Usage: python 5a_merge_consensus.py <library1> <library2>
# run for 5 individuals with following commands:
# LQ02: python 5a_merge_consensus.py 83mt.fasta 84mt.fasta 
# LQ08: python 5a_merge_consensus.py 55mt.fasta 74mt.fasta
# LQ12: python 5a_merge_consensus.py 54mt.fasta 70mt.fasta
# LQ24: python 5a_merge_consensus.py 56mt_pmd3.fasta 78mt_pmd3.fasta
# LQ30: python 5a_merge_consensus.py 77mt.fasta LQ30a.1.fasta
# =============================================================================

import sys

# function to read in file & assign header & sequence
def read_fasta(filename):
    header = ""
    sequence = ""
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line
            else:
                sequence += line
    
    return header, sequence

# function to write new fasta file with same header & sequence
def write_fasta(filename, header, sequence):
    with open(filename, 'w') as f:
        f.write(header + '\n')
        f.write(sequence + '\n')

# function to loop through sequence 1 & compare to sequence 2
# if seq 1 has an N, use letter from other sequence
# if seq2 has an N, use letter from other sequence

def merge_sequences(seq1, seq2):
    merged_seq1 = []
    merged_seq2 = []
    
    for i in range(len(seq1)):
        base1 = seq1[i]
        base2 = seq2[i]
        
        if base1 == 'N' and base2 != 'N':
            merged_seq1.append(base2)
        else:
            merged_seq1.append(base1)
        
        if base2 == 'N' and base1 != 'N':
            merged_seq2.append(base1)
        else:
            merged_seq2.append(base2)
    
    return ''.join(merged_seq1), ''.join(merged_seq2)

file1 = sys.argv[1]
file2 = sys.argv[2]

# read in file 1
header1, seq1 = read_fasta(file1)

# read in file 2
header2, seq2 = read_fasta(file2)

merged_seq1, merged_seq2 = merge_sequences(seq1, seq2)

write_fasta(f"{file1.replace('.fasta', '')}_merged.fasta", header1, merged_seq1)
write_fasta(f"{file2.replace('.fasta', '')}_merged.fasta", header2, merged_seq2)