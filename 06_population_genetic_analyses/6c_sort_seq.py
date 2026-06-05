#!/usr/bin/env python3

# =============================================================================
# Title: 6c_sort_seq.py
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: Python script to sort reference and PIALQ sequences into folders based on haplogroup
# Usage: python 6c_sort_seq.py
# =============================================================================

import os
import shutil
import pandas as pd

# 1. set paths
metadata_file = "diversity.meta"

ancient_folder = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/ancient_fasta"
modern_folder = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta"
destination_base = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/01_mafft"

# 2. load metadata
df = pd.read_csv(metadata_file, sep="\t", dtype=str).fillna("")

# 3.read in each row of file, strip whitespace
for idx, row in df.iterrows():
    file = row["file"].strip()
    paper = row["paper"].strip()
    ancient_or_modern = row["ancient_or_modern"].strip().lower()
    haplo = row["hap_2char"].strip()

    # Determine source path
    if ancient_or_modern == "ancient":
        folder = f"{paper}_ancient"
        src = os.path.join(ancient_folder, folder, file)
    elif ancient_or_modern == "modern":
        folder = f"{paper}_modern"
        src = os.path.join(modern_folder, folder, file)

    # Destination folder based on Macrohaplogroup
    dest_folder = os.path.join(destination_base, haplo)
    os.makedirs(dest_folder, exist_ok=True)

    dest_file = os.path.join(dest_folder, file)

    # Copy or warn
    if os.path.isfile(src):
        shutil.copy2(src, dest_file)
    else:
        print(f"missing file: {src}")