#!/usr/bin/env python3

# =============================================================================
# Title: 5l_add_traits_hap.py
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: Python script to add traits block ot nexus file based on haplogroup info from metadata: extra argument specifies which haplogroup column to use (1, 2, 3, 4, etc characters long)
# Usage: python 5l_add_traits_hap.py <input nexus file> <metadata file> <output nexus file> <haplogroup column name>
# =============================================================================

import pandas as pd
import re
import sys

# 1. get arguments from command line
nexus_input_path = sys.argv[1]
metadata_input_path = sys.argv[2]
nexus_output_path = sys.argv[3]
trait_column = sys.argv[4]

# 2. load metadata
metadata = pd.read_csv(metadata_input_path, sep="\t")

required_cols = ["file", trait_column]
for col in required_cols:
    if col not in metadata.columns:
        raise ValueError(f"Metadata file must contain a '{col}' column.")

# 3. clean sample names
metadata["file_clean"] = metadata["file"].str.replace(".fasta", "", regex=False)

# 4. read nexus and extract sample names
with open(nexus_input_path, "r") as f:
    nexus_text = f.read()

match = re.search(r"MATRIX\s*(.*?);\s*", nexus_text, re.IGNORECASE | re.DOTALL)
if not match:
    raise ValueError("Could not find MATRIX block in the NEXUS file.")

matrix_block = match.group(1)

dna_samples = []
for line in matrix_block.strip().splitlines():
    parts = line.strip().split()
    if parts:
        dna_samples.append(parts[0])

# 5. filter metadata to samples in matrix block
metadata = metadata[metadata["file_clean"].isin(dna_samples)]

missing_from_meta = set(dna_samples) - set(metadata["file_clean"])
if missing_from_meta:
    print("WARNING: DNA samples with no metadata:", missing_from_meta)

# 6. extract and clean trait labels
trait_labels = sorted(
    metadata[trait_column]
    .dropna()
    .unique()
)

trait_labels_clean = [
    lbl.replace(" ", "_") for lbl in trait_labels
]


# 7. build binary trait matrix
binary_traits = pd.DataFrame(
    0, index=metadata["file_clean"], columns=trait_labels_clean
)

for _, row in metadata.iterrows():
    raw = row[trait_column]
    if pd.isna(raw):
        continue
    cleaned = raw.replace(" ", "_")
    binary_traits.loc[row["file_clean"], cleaned] = 1

# 8. extract original nexus file up to first END;
parts = re.split(r"\bEND;\s*", nexus_text, flags=re.IGNORECASE)
nexus_data_block = parts[0].rstrip() + "\nEND;\n"

# 9. build TRAITS block
trait_block_lines = []
trait_block_lines.append("BEGIN TRAITS;")
trait_block_lines.append(f"\tDimensions NTRAITS={len(trait_labels_clean)};")
trait_block_lines.append("\tFormat labels=yes missing=? separator=Comma;")
trait_block_lines.append("\tTraitLabels " + " ".join(trait_labels_clean) + ";")
trait_block_lines.append("\tMatrix")

for sample in binary_traits.index:
    values = ",".join(str(v) for v in binary_traits.loc[sample])
    trait_block_lines.append(f"\t{sample}\t{values}")

trait_block_lines.append(";")
trait_block_lines.append("END;")

traits_block = "\n".join(trait_block_lines)

# 10. write output NEXUS with added TRAITS block
with open(nexus_output_path, "w") as f:
    f.write(nexus_data_block + "\n")
    f.write(traits_block + "\n")

print(f"Output written to: {nexus_output_path} 🧜")