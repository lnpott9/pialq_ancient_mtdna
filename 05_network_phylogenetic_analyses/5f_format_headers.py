#!/usr/bin/env python3

# =============================================================================
# Title: 5f_format_headers.py
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: Python script to ensure fasta file headers match sequence name used in metadata
# Usage: python 5f_format_headers.py
# run like this:
########## sorted="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft"
########## for folder in "$sorted"/*/; do
########## python 5f_format_headers.py "$folder"
########## done
# =============================================================================

import os
import sys

# 1. if # arguments isn't two, fail
if len(sys.argv) != 2:
    print("Usage: 07_rename.py <fasta_directory>")
    sys.exit(1)

# 2. get directory containing fasta files from second argument
fasta_directory = sys.argv[1]

# 3. process each fasta file in directory
for fasta_file in os.listdir(fasta_directory):
    if fasta_file.endswith(".fasta"):
        fasta_path = os.path.join(fasta_directory, fasta_file)
        new_fasta = []
        
        # extract filename without the extension
        fasta_name = os.path.splitext(fasta_file)[0]
    
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    new_fasta.append(f">{fasta_name}\n")
                else:
                    new_fasta.append(line)
        
        # overwrite original fasta
        with open(fasta_path, 'w') as f:
            f.writelines(new_fasta)