# Ancient DNA PIALQ & Afro-descendant networks

**location of fasta files:**

```
base=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/ancient_fasta

"$base"/barquera_2020_ancient
"$base"/fleskes_2023_ancient
"$base"/harney_2023_ancient
"$base"/sandovalvelasco_2023_ancient
"$base"/schroeder_2015_ancient
"$base"/pialq_ancient
```

I copied all files to 

## 1. examine metadata

make new metadata for all afro-descendant ancient samples

```
srun --mem=5gb --time=5:00:00 --pty bash

cd $SHARED/projects/PIALQ/2025_ancient_lp/11_mjn

# 1. activate environment
module load conda/python3
source activate $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/msa_env

# 2. read in merged text file & reduce down to ancient samples from the 5 papers + pialq
python
import pandas as pd
all = pd.read_csv("all_merged.meta", sep="\t")

# 3. subset to pialq + afro-ancient papers of interest
pialq = all[all["paper"] == "pialq"]
ancient = all[all["paper"].isin(["barquera_2020", "fleskes_2023", "harney_2023", "sandovalvelasco_2023", "schroeder_2015", "pialq"])]

pialq.shape
(58, 22)
ancient.shape
(131, 22)

# 4. check lengths
pialq["length"].describe()
count       58.000000
mean     15065.017241
std       1412.672638
min      10915.000000
25%      13970.000000
50%      15565.500000
75%      16356.000000
max      16558.000000
Name: length, dtype: float64

pialq.sort_values(["length"])

ancient["length"].describe()
count      131.000000
mean     14593.587786
std       3065.717273
min        324.000000
25%      14621.500000
50%      15308.000000
75%      16360.000000
max      16569.000000
Name: length, dtype: float64

ancient.sort_values(["length"]).head(15)

ancient = ancient[ancient["length"] > 5000]

ancient.shape
(125, 22)

ancient["paper"].value_counts()
paper
pialq                   58
harney_2023             22
fleskes_2023            22
sandovalvelasco_2023    19
schroeder_2015           3
barquera_2020            1
Name: count, dtype: int64

ancient["country"].value_counts()
country
Peru             58
United States    44
St. Helena       19
St. Martin        3
Mexico            1
Name: count, dtype: int64

# no filtering based on haplogroup quality
pialq["Quality"].describe()
count    58.000000
mean      0.897045
std       0.059395
min       0.729700
25%       0.871675
50%       0.908300
75%       0.938600
max       1.000000
Name: Quality, dtype: float64

ancient["Quality"].describe()
count    125.000000
mean       0.917630
std        0.057244
min        0.653000
25%        0.890500
50%        0.929000
75%        0.951900
max        1.000000
Name: Quality, dtype: float64

# write new meta (samples of interest & outgroups)
pialq.to_csv("pialq_mjn.meta", sep="\t", index=False, header=True)
ancient.to_csv("ancient_mjn.meta", sep="\t", index=False, header=True)
```



## 2. sort into alignment folders

```
python 02_sort.py
```



```
# 02_sort.py

import os
import shutil
import pandas as pd

# 1. set paths
metadata_file = "ancient_mjn.meta"

ancient_folder = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/ancient_fasta"
destination = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/11_mjn/01_mafft"

# 2. load metadata
df = pd.read_csv(metadata_file, sep="\t", dtype=str).fillna("")

# 3.read in each row of file, strip whitespace
for idx, row in df.iterrows():
    file = row["file"].strip()
    paper = row["paper"].strip()
    folder = f"{paper}_ancient"
    src = os.path.join(ancient_folder, folder, file)

    dest_file = os.path.join(destination, file)

    # copy or warn if no file
    if os.path.isfile(src):
        shutil.copy2(src, dest_file)
    else:
        print(f"missing file: {src}")
```



## 03. reformat fasta headers

```
# set path to fasta folder
folder="$SHARED/projects/PIALQ/2025_ancient_lp/11_mjn/01_mafft"

# run on fastas
python 03_rename.py "$folder"
```



```
# 03_rename.py

import os
import sys

# 1. if # arguments isn't two, fail
if len(sys.argv) != 2:
    print("Usage: python rename.py <fasta_directory>")
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
                if line.startswith(">"):  # header
                    new_fasta.append(f">{fasta_name}\n")  # rename header using filename
                else:
                    new_fasta.append(line)
        
        # overwrite original fasta
        with open(fasta_path, 'w') as f:
            f.writelines(new_fasta)
        print(f"renamed header in {fasta_file}! 🎉")
```



## 04. align with mafft

```
#!/bin/bash -l
#SBATCH --job-name=linsi
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100gb
#SBATCH --tmp=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

module load mafft/7.475

BASE_DIR="$SHARED/projects/PIALQ/2025_ancient_lp/11_mjn"
INPUT="$BASE_DIR/01_mafft"
REF="$SHARED/ref_seqs"

# Build input file
cat "$REF/rCRS.fasta" $INPUT/*fasta | awk 'NF' > "$BASE_DIR/ancient_mafft_in.fasta"

# mafft with L-INS-i
if mafft --localpair --maxiterate 1000 \
    --thread 20 \
    "$BASE_DIR/ancient_mafft_in.fasta" \
    > "$BASE_DIR/ancient_mafft_out.fasta" ; then
    echo "brethren rejoice 💀"
else
	echo "once again we're fucked"
fi
```



```
#!/bin/bash -l
#SBATCH --job-name=pialq_linsi
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=150gb
#SBATCH --tmp=150gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

module load mafft/7.475

BASE_DIR="$SHARED/projects/PIALQ/2025_ancient_lp/11_mjn"
INPUT="$BASE_DIR/01_mafft"
REF="$SHARED/ref_seqs"

# Build input file
cat "$REF/rCRS.fasta" $INPUT/LQ*fasta | awk 'NF' > "$BASE_DIR/pialq_mafft_in.fasta"

# mafft with L-INS-i
if mafft --localpair --maxiterate 1000 \
    --thread 20 \
    "$BASE_DIR/pialq_mafft_in.fasta" \
    > "$BASE_DIR/pialq_mafft_out.fasta" ; then
    echo "brethren rejoice 💀"
else
        echo "once again we're fucked"
f
```



## 05. clean fastas

```
# clean. r
# Load libraries
library(pegas)
library(Biostrings)
library(ape)

# Set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/11_mjn/")

# 1. Paths to MSAs
pialq_path   <- "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/11_mjn/pialq_mafft_out.fasta"
ancient_path <- "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/11_mjn/ancient_mafft_out.fasta"

# 2. List of paths
path_list <- list(pialq_path, ancient_path)
names(path_list) <- c("pialq", "ancient")

# 3. Read sequences and convert to uppercase character matrices
aligned_seqs <- lapply(path_list, read.dna, format = "fasta")
aligned_seqs_char <- lapply(aligned_seqs, function(x) toupper(as.character(x)))

# 4. Identify rCRS sequence in each alignment
rCRS_seq <- lapply(aligned_seqs_char, function(x) {
  rcs_name <- grep("NC_012920", rownames(x), value = TRUE)
  if (length(rcs_name) == 0) stop("ERROR: rCRS not found in alignment")
  x[rcs_name, ]
})

# 5. Clean alignments — remove rCRS gaps, mask problem positions
problem_positions <- c(303:315, 515:522, 568:573, 3107, 16182:16194, 16519)

aligned_seqs_no_gaps <- lapply(seq_along(aligned_seqs_char), function(i) {
  alignment_mat <- aligned_seqs_char[[i]]
  rcs <- rCRS_seq[[i]]
  
  # Map rCRS alignment columns to rCRS coordinates
  rcs_coord <- 0
  mask_cols <- c()
  
  for (j in seq_along(rcs)) {
    if (rcs[j] != "-") {
      rcs_coord <- rcs_coord + 1
      if (rcs_coord %in% problem_positions) {
        mask_cols <- c(mask_cols, j)
      }
    }
  }
  
  if (length(mask_cols) > 0) {
    alignment_mat[, mask_cols] <- "N"
  }
  
  alignment_mat
})

# ✔ restore original names (“pialq” and “ancient”)
names(aligned_seqs_no_gaps) <- names(aligned_seqs_char)

# 6. Calculate missingness
count_missing <- function(seq_matrix) {
  apply(seq_matrix, 1, function(row) sum(row == "N"))
}

missing_counts_list <- lapply(aligned_seqs_no_gaps, count_missing)

missing_df_list <- lapply(names(aligned_seqs_no_gaps), function(name) {
  mat <- aligned_seqs_no_gaps[[name]]
  counts <- missing_counts_list[[name]]
  data.frame(
    Sample = rownames(mat),
    Alignment = name,
    Missing_Count = counts,
    Missing_Fraction = counts / ncol(mat)
  )
})

missing_df <- do.call(rbind, missing_df_list)

# 7. Convert back to DNAStringSet
clean_alignments_list <- lapply(aligned_seqs_no_gaps, function(mat) {
  dna_seq <- apply(mat, 1, function(row) DNAString(paste(row, collapse = "")))
  names(dna_seq) <- rownames(mat)
  DNAStringSet(dna_seq)
})

# 8. Remove rCRS
clean_alignments_list <- lapply(clean_alignments_list, function(dna_set) {
  drop <- grep("NC_012920", names(dna_set), value = TRUE)
  dna_set[!names(dna_set) %in% drop]
})

# 9. Write cleaned FASTA files with correct names
output_names <- c("clean_pialq.fasta", "clean_ancient.fasta")

for (i in seq_along(clean_alignments_list)) {
  writeXStringSet(
    clean_alignments_list[[i]],
    filepath = output_names[i],
    width = 18000
  )
}

```



## 5.5 split ancient MJN into 2

in r

```
# split_mjn.r
# load libraries
if (!require(pegas)) {
  install.packages("pegas")
  library(pegas)
}

if (!require(ape)) {
  install.packages("ape")
  library(ape)
}

if (!require(Biostrings)) {
  install.packages("Biostrings")
  library(Biostrings)
}

# 1. set wd
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/11_mjn/")

# 2. read in aligned fasta
fasta <- read.dna("clean_ancient.fasta", format = "fasta")
dim(fasta)

# 3. load metadata
meta <- read.csv("ancient_mjn.meta", sep="\t")

# 4. create objects for each subset
ooa <- meta[meta$macrohaplogroup %in% c("A", "B", "C", "D", "H", "J"), ]
final_ooa <- fasta[ooa$Sample, ]
dim(final_ooa) # 16 individuals

a <- meta[meta$macrohaplogroup %in% c("L0", "L1", "L2", "L3", "L4"), ]
final_a <- fasta[a$Sample, ]
dim(final_a) # 109 individuals

# 5. write out new fasta files
write.dna(ooa, file = "clean_ancient_ooa.fasta", format = "fasta")
write.dna(a, file = "clean_ancient_a.fasta", format = "fasta")
```



## 06. prepare for popart

a. **make nexus file**

```
# 1. activate environment
module load conda/python3
source activate $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/msa_env

python $SHARED/programs/AMAS/amas/AMAS.py convert -i clean*fasta -f fasta -u nexus -d dna
```

b. **replace Ns with ?s**

```
for nexus in clean*.fasta-out.nex; do
name=$(basename "$nexus" .fasta-out.nex)
python 06_replace.py "$nexus" "rep_${name}.nex"
done
```

```
for nexus in *.fasta-out.nex; do
name=$(basename "$nexus" .fasta-out.nex)
python 06_replace.py "$nexus" "rep_${name}.nex"
done
```



*06_replace.py*

```
#!/usr/bin/env python3
import argparse

def replace_N_with_missing(infile, outfile):
    with open(infile, "r") as f:
        lines = f.readlines()

    new_lines = []
    in_matrix = False

    for line in lines:
        stripped = line.strip()

        # Detect start/end of the MATRIX section
        if stripped.upper().startswith("MATRIX"):
            in_matrix = True
            new_lines.append(line)
            continue
        if in_matrix and stripped.upper().startswith(";"):
            in_matrix = False
            new_lines.append(line)
            continue

        if in_matrix and stripped:  # inside MATRIX
            parts = stripped.split()
            if len(parts) >= 2:
                name, seq = parts[0], parts[1]
                # Replace Ns (upper/lower) with ?
                seq = seq.replace("N", "?").replace("n", "?")
                new_line = f"{name}\t{seq}\n"
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    # Write out modified nexus file
    with open(outfile, "w") as f:
        f.writelines(new_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Replace Ns in NEXUS DNA sequences with ?"
    )
    parser.add_argument("infile", help="Input NEXUS file")
    parser.add_argument("outfile", help="Output NEXUS file with Ns replaced by ?")

    args = parser.parse_args()
    replace_N_with_missing(args.infile, args.outfile)


if __name__ == "__main__":
    main()
```



c. **add traits block for paper**

```
python 06_paper.py rep_clean_ancient.nex all_merged.meta popart_ancient_paper.nex
python 06_paper.py rep_clean_ancient_a.nex all_merged.meta popart_ancient_paper_a.nex
python 06_paper.py rep_clean_ancient_ooa.nex all_merged.meta popart_ancient_paper_ooa.nex

python 06_paper.py rep_clean_pialq.nex all_merged.meta popart_pialq_paper.nex
```



*06_paper.py*

```
#!/usr/bin/env python3
import pandas as pd
import re
import sys

# -----------------------------
# 1. Parse arguments
# -----------------------------
nexus_input_path = sys.argv[1]
metadata_input_path = sys.argv[2]
nexus_output_path = sys.argv[3]
trait_column = "paper"   # column in metadata containing trait labels


# -----------------------------
# 2. Load metadata
# -----------------------------
metadata = pd.read_csv(metadata_input_path, sep="\t")

required_cols = ["file", trait_column]
for col in required_cols:
    if col not in metadata.columns:
        raise ValueError(f"Metadata file must contain a '{col}' column.")


# -----------------------------
# 2b. Clean sample names
# -----------------------------
metadata["file_clean"] = metadata["file"].str.replace(".fasta", "", regex=False)


# -----------------------------
# 3. Read NEXUS and extract MATRIX samples
# -----------------------------
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


# -----------------------------
# 4. Filter metadata to samples in DNA block
# -----------------------------
metadata = metadata[metadata["file_clean"].isin(dna_samples)]

missing_from_meta = set(dna_samples) - set(metadata["file_clean"])
if missing_from_meta:
    print("WARNING: DNA samples with no metadata:", missing_from_meta)


# -----------------------------
# 5. Extract and clean trait labels
# -----------------------------
trait_labels = sorted(
    metadata[trait_column]
    .dropna()
    .unique()
)

trait_labels_clean = [
    lbl.replace(" ", "_") for lbl in trait_labels
]


# -----------------------------
# 6. Build binary trait matrix
# -----------------------------
binary_traits = pd.DataFrame(
    0, index=metadata["file_clean"], columns=trait_labels_clean
)

for _, row in metadata.iterrows():
    raw = row[trait_column]
    if pd.isna(raw):
        continue
    cleaned = raw.replace(" ", "_")
    binary_traits.loc[row["file_clean"], cleaned] = 1


# -----------------------------
# 7. Extract original NEXUS up to first END;
# -----------------------------
parts = re.split(r"\bEND;\s*", nexus_text, flags=re.IGNORECASE)
nexus_data_block = parts[0].rstrip() + "\nEND;\n"


# -----------------------------
# 8. Build TRAITS block
# -----------------------------
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


# -----------------------------
# 9. Write output NEXUS with added TRAITS block
# -----------------------------
with open(nexus_output_path, "w") as f:
    f.write(nexus_data_block + "\n")
    f.write(traits_block + "\n")

print(f"Output written to: {nexus_output_path} 🧜")
```



d. **add traits block of haplogroup**

run

```
python 06_hap.py rep_clean_ancient.nex all_merged.meta popart_ancient_hap.nex
python 06_hap.py rep_clean_ancient_a.nex all_merged.meta popart_ancient_hap_a.nex
python 06_hap.py rep_clean_ancient_ooa.nex all_merged.meta popart_ancient_hap_ooa.nex

python 06_hap.py rep_clean_pialq.nex all_merged.meta popart_pialq_hap.nex
```



*06_hap.py*

```
#!/usr/bin/env python3
import pandas as pd
import re
import sys

# -----------------------------
# 1. Parse arguments
# -----------------------------
nexus_input_path = sys.argv[1]
metadata_input_path = sys.argv[2]
nexus_output_path = sys.argv[3]
trait_column = "hap_2char"   # column in metadata containing trait labels


# -----------------------------
# 2. Load metadata
# -----------------------------
metadata = pd.read_csv(metadata_input_path, sep="\t")

required_cols = ["file", trait_column]
for col in required_cols:
    if col not in metadata.columns:
        raise ValueError(f"Metadata file must contain a '{col}' column.")


# -----------------------------
# 2b. Clean sample names
# -----------------------------
metadata["file_clean"] = metadata["file"].str.replace(".fasta", "", regex=False)


# -----------------------------
# 3. Read NEXUS and extract MATRIX samples
# -----------------------------
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


# -----------------------------
# 4. Filter metadata to samples in DNA block
# -----------------------------
metadata = metadata[metadata["file_clean"].isin(dna_samples)]

missing_from_meta = set(dna_samples) - set(metadata["file_clean"])
if missing_from_meta:
    print("WARNING: DNA samples with no metadata:", missing_from_meta)


# -----------------------------
# 5. Extract and clean trait labels
# -----------------------------
trait_labels = sorted(
    metadata[trait_column]
    .dropna()
    .unique()
)

trait_labels_clean = [
    lbl.replace(" ", "_") for lbl in trait_labels
]


# -----------------------------
# 6. Build binary trait matrix
# -----------------------------
binary_traits = pd.DataFrame(
    0, index=metadata["file_clean"], columns=trait_labels_clean
)

for _, row in metadata.iterrows():
    raw = row[trait_column]
    if pd.isna(raw):
        continue
    cleaned = raw.replace(" ", "_")
    binary_traits.loc[row["file_clean"], cleaned] = 1


# -----------------------------
# 7. Extract original NEXUS up to first END;
# -----------------------------
parts = re.split(r"\bEND;\s*", nexus_text, flags=re.IGNORECASE)
nexus_data_block = parts[0].rstrip() + "\nEND;\n"


# -----------------------------
# 8. Build TRAITS block
# -----------------------------
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


# -----------------------------
# 9. Write output NEXUS with added TRAITS block
# -----------------------------
with open(nexus_output_path, "w") as f:
    f.write(nexus_data_block + "\n")
    f.write(traits_block + "\n")

print(f"Output written to: {nexus_output_path} 🧜")
```



# 07. plot in popart

download popart from this website: [PopART version 1.7 for Windows](https://popart.maths.otago.ac.nz/wp-content/uploads/2022/08/popart-1.7.zip) from this website: https://popart.maths.otago.ac.nz/download/

manual: https://popart.maths.otago.ac.nz/wp-content/uploads/2022/08/pdf.zip

settings for each network:

```
popart_ancient_hap.nex
popart_ancient_paper.nex
popart_pialq_hap.nex
```

### ancient

popart_ancient_paper.nex stats:

```
Nucleotide diversity:	pi = 0.171011
Number of segregating sites:	406
Number of parsimony-informative sites:	228
Tajima's D statistic:	D = 61.3255
	p (D >= 61.3255) = 0
Node label	Matching Sequences
I8092	I8092
	I8093
LQ01	LQ01
	LQ58
LQ22	LQ22
	LQ29
LQ56	LQ56
	STH_347
	STH_441
LQ12	LQ12
	LQ17
	LQ31
	LQ55
LQ06	LQ06
	LQ23
	LQ62
I15334	I15334
	LQ21
	STH_358
LQ45	LQ45
	LQ48
STH_284	STH_284
	STH_499
LQ32	LQ32
	LQ36
I15335	I15335
	I15338
	I8089
	I8090
LQ16	LQ16
	LQ19
I15330	I15330
	I15331
	I15332
	I15340
LQ13	LQ13
	LQ50
	LQ57
I15336	I15336
	I15337
CHS36	CHS36
	CHS37
LQ28	LQ28
	LQ61

```

maria wanted them split into americas/eurasian + african so I did:

**ancient africa**

```
Nucleotide diversity:	pi = 0.00318632
Number of segregating sites:	362
Number of parsimony-informative sites:	190
Tajima's D statistic:	D = -2.01856
	p (D >= -2.01856) = 0.990196
Node label	Matching Sequences
I8092	I8092
	I8093
LQ22	LQ22
	LQ29
LQ56	LQ56
	STH_347
	STH_441
LQ12	LQ12
	LQ17
	LQ31
	LQ55
I15334	I15334
	LQ21
	STH_358
LQ45	LQ45
	LQ48
STH_284	STH_284
	STH_499
LQ32	LQ32
	LQ36
I15335	I15335
	I15338
	I8089
	I8090
LQ16	LQ16
	LQ19
I15330	I15330
	I15331
	I15332
	I15340
LQ13	LQ13
	LQ50
	LQ57
I15336	I15336
	I15337
CHS36	CHS36
	CHS37

```

**ancient out of africa**

```
Node label	Matching Sequences
LQ01	LQ01
	LQ58
LQ06	LQ06
	LQ62
LQ28	LQ28
	LQ61
Nucleotide diversity:	pi = 0.00219029
Number of segregating sites:	107
Number of parsimony-informative sites:	64
Tajima's D statistic:	D = -0.644698
	p (D >= -0.644698) = 0.716963
```





### pialq

```
popart_pialq_hap.nex
+ no sequences removed for having >5% undefined sites
+ log all stats to popart_pialq_hap.log
+ default vertex size: 50
+ edge color dark gray
+ label font: arial 26
+ legend font: arial 50
+ redraw network with 50 iterations

Node label	Matching Sequences
LQ01	LQ01
	LQ58
LQ22	LQ22
	LQ29
LQ12	LQ12
	LQ17
	LQ31
	LQ55
LQ06	LQ06
	LQ62
LQ45	LQ45
	LQ48
LQ32	LQ32
	LQ36
LQ16	LQ16
	LQ19
LQ13	LQ13
	LQ50
	LQ57
LQ28	LQ28
	LQ61
Nucleotide diversity:	pi = 0.00354345
Number of segregating sites:	326
Number of parsimony-informative sites:	196
Tajima's D statistic:	D = -1.75587
	p (D >= -1.75587) = 0.971415
```



After looking at this figure, I decided to use manually designated groupings based on what was present in this network:

pialq_final.meta

```
individual	network_grouping
LQ28	A2+(64)
LQ61	A2+(64)
LQ05	B2
LQ14	B2
LQ49	B2
LQ06	C1b
LQ23	C1b
LQ60	C1b
LQ62	C1b
LQ01	D1
LQ46	D1
LQ58	D1
LQ18	H1bw
LQ25	L0a2a
LQ44	L0a2a
LQ27	L0d1a
LQ08	L1c1a
LQ45	L1c1a
LQ48	L1c1a
LQ15	L1c1b
LQ41	L1c1d
LQ54	L1c1d
LQ39	L1c2'4
LQ22	L2a1
LQ29	L2a1
LQ56	L2a1
LQ59	L2a1
LQ63	L2a1
LQ02	L2b1a
LQ16	L2b1a
LQ19	L2b1a
LQ12	L2b1b
LQ17	L2b1b
LQ31	L2b1b
LQ55	L2b1b
LQ03	L2c
LQ30	L2c
LQ52	L2c
LQ11	L3b1a
LQ38	L3b1a
LQ09	L3d1
LQ24	L3d1
LQ26	L3d1
LQ33	L3d1
LQ64	L3d3
LQ40	L3d4
LQ10	L3e1a1a
LQ13	L3e1a1a
LQ50	L3e1a1a
LQ57	L3e1a1a
LQ32	L3e2
LQ35	L3e2
LQ36	L3e2
LQ42	L3e2
LQ43	L3e3b
LQ20	L3f1b
LQ21	L3f1b
LQ34	L3f1b
```



```
python 07_pialq_final.py rep_clean_pialq.nex pialq_final.meta popart_pialq_final.nex
```



07_pialq_final.py

```
#!/usr/bin/env python3
import pandas as pd
import re
import sys

# -----------------------------
# 1. Parse arguments
# -----------------------------
nexus_input_path = sys.argv[1]
metadata_input_path = sys.argv[2]
nexus_output_path = sys.argv[3]
trait_column = "network_grouping"   # column in metadata containing trait labels


# -----------------------------
# 2. Load metadata
# -----------------------------
metadata = pd.read_csv(metadata_input_path, sep="\t")

required_cols = ["individual", trait_column]
for col in required_cols:
    if col not in metadata.columns:
        raise ValueError(f"Metadata file must contain a '{col}' column.")

# -----------------------------
# 3. Read NEXUS and extract MATRIX samples
# -----------------------------
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


# -----------------------------
# 4. Filter metadata to samples in DNA block
# -----------------------------
metadata = metadata[metadata["individual"].isin(dna_samples)]

missing_from_meta = set(dna_samples) - set(metadata["individual"])
if missing_from_meta:
    print("WARNING: DNA samples with no metadata:", missing_from_meta)


# -----------------------------
# 5. Extract and clean trait labels
# -----------------------------
trait_labels = sorted(
    metadata[trait_column]
    .dropna()
    .unique()
)

trait_labels_clean = [
    lbl.replace(" ", "_") for lbl in trait_labels
]


# -----------------------------
# 6. Build binary trait matrix
# -----------------------------
binary_traits = pd.DataFrame(
    0, index=metadata["individual"], columns=trait_labels_clean
)

for _, row in metadata.iterrows():
    raw = row[trait_column]
    if pd.isna(raw):
        continue
    cleaned = raw.replace(" ", "_")
    binary_traits.loc[row["individual"], cleaned] = 1


# -----------------------------
# 7. Extract original NEXUS up to first END;
# -----------------------------
parts = re.split(r"\bEND;\s*", nexus_text, flags=re.IGNORECASE)
nexus_data_block = parts[0].rstrip() + "\nEND;\n"


# -----------------------------
# 8. Build TRAITS block
# -----------------------------
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


# -----------------------------
# 9. Write output NEXUS with added TRAITS block
# -----------------------------
with open(nexus_output_path, "w") as f:
    f.write(nexus_data_block + "\n")
    f.write(traits_block + "\n")

print(f"Output written to: {nexus_output_path} 🧜")
```

