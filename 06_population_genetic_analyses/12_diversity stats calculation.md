# 12. diversity stats calculation

## choose papers to use

```
# choose.r
# 1. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity")

# 2. load libraries, download if necessary
library(seqinr)
library(pegas)
library(dplyr)
library(ape)

# 3. load metadata & identify populations of interest
meta <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/all_merged.meta", sep="\t")
dim(meta)
# 18540    21

#### a. remove low quality ones, make this match tree stringency
meta <- meta[(meta$Quality > 0.8 & meta$length > 14000) | meta$paper == "pialq", ]

#### b. unique populations & papers
# look at unique combinations
c_unique <- meta %>% distinct(country, .keep_all = TRUE)
dim(c_unique)

# paper-country combo
pc_unique <- meta %>% distinct(paper, country, .keep_all = TRUE)
dim(pc_unique)
# 350  21

# count distinct options
pair_summary <- meta %>%
  count(paper, country)

big <- pair_summary[pair_summary$n > 60, ]


#### c. pick out populations to keep
keep <- meta$paper %in% c("pialq", "harney_2023", "sandovalvelasco_2023", "fleskes_2023", "1kg", "hgdp") | meta$paper %in% big$paper

meta_filtered <- meta[keep, ]
dim(meta_filtered)

# write df of sample sizes
meta_filt <- meta_filtered %>%
  count(paper, country)

unique_papers <- unique(meta_filt$paper)

write.table(meta_filt, "n_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```



## make new metadata 

```
# 1. activate environment
module load conda/python3
source activate $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/msa_env

cd /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity

python

# 2. import pandas
import pandas as pd

# 3. read in all_merged.meta as dataframe
meta = pd.read_csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/all_merged.meta", sep="\t")

# 4. only retain specific papers
papers = ["1kg", "arias_2017", "avila_2019", "barbieri_2012a", "barbieri_2013a", "barbieri_2013", "barbieri_2014a", "barbieri_2014b", "barbieri_2017", "bodner_2022", "brandini_2018", "brucato_2018", "chan_2019", "colombo_2025", "derenko_2007", "fleskes_2023", "garcia_2021", "garciaolivares_2022", "harney_2023", "hgdp", "huber_2025", "just_2008", "just_2015", "kumar_2011", "mccrow_2016", "oliveira_2018", "olivieri_2017", "pialq", "pierron_2017", "sandovalvelasco_2023", "silva_2021", "taylor_2020", "tito_unpub", "tribaldos_2021"]

subset = meta[meta["paper"].isin(papers)]
subset.shape # (17390, 21)

# 5. subset for length
ancient = subset[(subset["ancient_or_modern"] == "ancient")]
ancient = ancient.sort_values("length")
# remove everything < 10000 (4 harney samples)

quality = subset[(subset["length"] > 10000)]

quality.shape # (17384, 21)

# 6. write out new meta
quality.to_csv("diversity.meta", sep="\t", index=False, header=True)
```

## move

```
python sort.py
```



```
# sort.py

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
```

expected missing files:

```
missing file: /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta/tito_unpub_modern/JX669265.fasta
missing file: /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta/tito_unpub_modern/JX669269.fasta
missing file: /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta/tito_unpub_modern/JX669285.fasta
```



## rename

```
# set path to sorted fasta folder
sorted="$SHARED/projects/PIALQ/2025_ancient_lp/12_diversity/01_mafft"

for folder in "$sorted"/*/; do
python rename.py "$folder"
done
```



```
# rename.py

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
```



## align all:

```
#!/bin/bash -l
#SBATCH --job-name=01_m
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=80gb
#SBATCH --tmp=80gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-126

module load mafft/7.475

BASE_DIR="$SHARED/projects/PIALQ/2025_ancient_lp/12_diversity"
SORTED="$BASE_DIR/01_mafft"
REF="$SHARED/ref_seqs"

# Manual folder listing
FOLDERS=(
"A+"
"A1"
"A2"
"A3"
"A5"
"A6"
"A8"
"B2"
"B4"
"B5"
"B6"
"C"
"C1"
"C4"
"C5"
"C7"
"D1"
"D2"
"D3"
"D4"
"D5"
"D6"
"E1"
"E2"
"F1"
"F2"
"F3"
"F4"
"G1"
"G2"
"G3"
"H"
"H+"
"H1"
"H2"
"H3"
"H4"
"H5"
"H6"
"H7"
"H8"
"H9"
"HV"
"I"
"I1"
"I2"
"I3"
"I4"
"I5"
"I6"
"J1"
"J2"
"K1"
"K2"
"K3"
"L0"
"L1"
"L2"
"L3"
"L4"
"L5"
"M"
"M1"
"M2"
"M3"
"M4"
"M5"
"M6"
"M7"
"M8"
"M9"
"N1"
"N2"
"N5"
"N8"
"N9"
"P1"
"P9"
"Q1"
"R"
"R+"
"R0"
"R1"
"R2"
"R3"
"R5"
"R6"
"R7"
"R8"
"R9"
"T1"
"T2"
"U"
"U1"
"U2"
"U3"
"U4"
"U5"
"U6"
"U7"
"U8"
"U9"
"V"
"V+"
"V1"
"V2"
"V3"
"V7"
"V9"
"W"
"W+"
"W1"
"W3"
"W4"
"W5"
"W6"
"W7"
"W8"
"W9"
"X1"
"X2"
"X3"
"Y1"
"Y2"
"Z2"
"Z3"
"Z4"
)

# Folder for this task
FOLDER="${FOLDERS[$SLURM_ARRAY_TASK_ID]}"
FULL_PATH="$SORTED/$FOLDER"

echo "processing folder $FOLDER at $FULL_PATH"

# Build input file
cat "$REF/rCRS.fasta" $FULL_PATH/*fasta | awk 'NF' > "$SORTED/${FOLDER}_mafft_in_auto.fasta"
echo "rcrs added"

# real mafft with auto
if mafft --auto \
    --thread 10 \
    "$SORTED/${FOLDER}_mafft_in_auto.fasta" \
    > "$SORTED/${FOLDER}_mafft_out_auto.fasta" ; then
    echo "brethren rejoice 💀"
else
	echo "once again we're fucked"
fi
```





## clean & join alignments

```
# clean.r
# Load libraries
for(pkg in c("pegas", "ape", "Biostrings")) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# 1. Set directories
seq_dir   <- "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/01_mafft"
clean_dir <- "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/02_clean"

# 2. List input FASTAs
path_list <- list.files(seq_dir, pattern = "_mafft_out_auto\\.fasta$", full.names = TRUE)
names(path_list) <- sub("_mafft_out_auto\\.fasta$", "", basename(path_list))

# 3. Read sequences as DNAbin and convert to uppercase character matrices
aligned_seqs_char <- lapply(path_list, function(f) toupper(as.character(read.dna(f, format = "fasta"))))

# 4. Extract rCRS sequence from each alignment
rCRS_seq <- lapply(aligned_seqs_char, function(mat) mat[grep("NC_012920.1", rownames(mat)), ])

# 5. Mask problem positions in the original alignment
problem_positions <- c(303:315, 515:522, 568:573, 3107, 16182:16194, 16519)

aligned_seqs_masked <- mapply(function(mat, rcs) {
  rcs_coord <- 0
  mask_cols <- integer(0)
  
  for(j in seq_along(rcs)) {
    if(rcs[j] != "-") rcs_coord <- rcs_coord + 1
    if(rcs_coord %in% problem_positions) mask_cols <- c(mask_cols, j)
  }
  
  if(length(mask_cols) > 0) mat[, mask_cols] <- "N"
  mat
}, aligned_seqs_char, rCRS_seq, SIMPLIFY = FALSE)

# 6. Remove columns where rCRS has gaps
aligned_seqs_clean <- mapply(function(mat, rcs) {
  keep_cols <- rcs != "-"
  mat[, keep_cols, drop = FALSE]
}, aligned_seqs_masked, rCRS_seq, SIMPLIFY = FALSE)

# 7. Convert each cleaned alignment to DNAStringSet and remove rCRS
clean_alignments_list <- lapply(aligned_seqs_clean, function(mat) {
  dna_seq <- DNAStringSet(apply(mat, 1, paste, collapse = ""))
  names(dna_seq) <- rownames(mat)
  dna_seq[!grepl("NC_012920.1", names(dna_seq))]
})

# 8. Concatenate all DNAStringSets into one
all_seqs <- unlist(lapply(clean_alignments_list, as.character))
combined_alignment <- DNAStringSet(all_seqs)

# 9. Write single combined FASTA
writeXStringSet(combined_alignment, filepath = file.path(clean_dir, "clean_div.fasta"), width = 18000)
```



## diversity script

intrapop: pi, h

```
# INTRAPOPULATION CALCULATIONS: PI, HD

# 1. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/03_intra")

# 2. load libraries, download if necessary
library(seqinr)
library(pegas)
library(dplyr)
library(ape)
library(ggplot2)

# 3. load metadata
meta <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/diversity.meta", sep="\t")
dim(meta)

# 4. load aligned sequences
align <- read.dna("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/02_clean/clean_div.fasta", format = "fasta")

# Check dimensions
dim(align)  # rows = sequences, columns = alignment length
nrow(align)  # number of sequences
ncol(align)  # alignment length

# 4. change sequence names to accession numbers
rownames(align) <- sub("^[^.]+\\.", "", rownames(align))
head(rownames(align))

# drop duplicates
align_dedup <- align[!duplicated(rownames(align)), ]

# 5. make population metadata subsets!
# 1kg
onekg <- meta[(meta$paper == "1kg"),]
acb <- meta[(meta$paper == "1kg") & (meta$country == "ACB"),] # 98
asw <- meta[(meta$paper == "1kg") & (meta$country == "ASW"),] # 66
beb <- meta[(meta$paper == "1kg") & (meta$country == "BEB"),] # 104
cdx <- meta[(meta$paper == "1kg") & (meta$country == "CDX"),] # 107
ceu <- meta[(meta$paper == "1kg") & (meta$country == "CEU"),] # 101
chb <- meta[(meta$paper == "1kg") & (meta$country == "CHB"),] # 105
chs <- meta[(meta$paper == "1kg") & (meta$country == "CHS"),] # 112
clm <- meta[(meta$paper == "1kg") & (meta$country == "CLM"),] # 103
esn <- meta[(meta$paper == "1kg") & (meta$country == "ESN"),] # 111
fin <- meta[(meta$paper == "1kg") & (meta$country == "FIN"),] # 105
gbr <- meta[(meta$paper == "1kg") & (meta$country == "GBR"),] # 106
gih <- meta[(meta$paper == "1kg") & (meta$country == "GIH"),] # 109
gwd <- meta[(meta$paper == "1kg") & (meta$country == "GWD"),] # 120
ibs <- meta[(meta$paper == "1kg") & (meta$country == "IBS"),] # 108
itu <- meta[(meta$paper == "1kg") & (meta$country == "ITU"),] # 112
jpt <- meta[(meta$paper == "1kg") & (meta$country == "JPT"),] # 105
khv <- meta[(meta$paper == "1kg") & (meta$country == "KHV"),] # 103
lwk <- meta[(meta$paper == "1kg") & (meta$country == "LWK"),] # 112
msl <- meta[(meta$paper == "1kg") & (meta$country == "MSL"),] # 98
mxl <- meta[(meta$paper == "1kg") & (meta$country == "MXL"),] # 70
pjl <- meta[(meta$paper == "1kg") & (meta$country == "PJL"),] # 108
pel <- meta[(meta$paper == "1kg") & (meta$country == "PEL"),] # 90
pur <- meta[(meta$paper == "1kg") & (meta$country == "PUR"),] # 107
stu <- meta[(meta$paper == "1kg") & (meta$country == "STU"),] # 107
tsi <- meta[(meta$paper == "1kg") & (meta$country == "TSI"),] # 111
yri <- meta[(meta$paper == "1kg") & (meta$country == "YRI"),] # 109

# hgdp
# makrani_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MakraniHGDP"),] # 23
# mbuti_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MbutiHGDP"),] # 10
# biaka_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "BiakaHGDP"),] # 24
# french_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "FrenchHGDP"),] # 24
# sardinian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "SardinianHGDP"),] # 24
# colombian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "ColombianHGDP"),] # 5
# cambodian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "CambodianHGDP"),] # 8
# japanese_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "JapaneseHGDP"),] # 27
# han_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "HanHGDP"),] # 30
# surui_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "SuruiHGDP"),] # 6
# mandenka_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MandenkaHGDP"),] # 20
# yoruba_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "YorubaHGDP"),] # 18
# san_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "SanHGDP"),] # 2
# bantu_southafrica_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "BantuSouthAfricaHGDP"),] # 4
# karitiana_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "KaritianaHGDP"),] # 9
# tujia_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "TujiaHGDP"),] # 8
# bergamo_italian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "BergamoItalianHGDP"),] # 10
# tuscan_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "TuscanHGDP"),] # 6
# yi_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "YiHGDP"),] # 8
# miao_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MiaoHGDP"),] # 8
# oroqen_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "OroqenHGDP"),] # 7
# daur_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "DaurHGDP"),] # 9
# mongolian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MongolianHGDP"),] # 8
# xibo_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "XiboHGDP"),] # 7
# northern_han_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "NorthernHanHGDP"),] # 10
# uygur_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "UygurHGDP"),] # 8
# daihgdp_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "DaiHGDP"),] # 5
# lahu_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "LahuHGDP"),] # 6
# she_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "SheHGDP"),] # 8
# naxi_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "NaxiHGDP"),] # 6
# tu_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "TuHGDP"),] # 8
# bantu_kenya_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "BantuKenyaHGDP"),] # 10
mozabite_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MozabiteHGDP"),] # 26

# arias 2017
arias_2017_colombia <- meta[(meta$paper == "arias_2017") & (meta$country == "Colombia"),] # 436

# avila_2019
avila_2019_brazil <- meta[(meta$paper == "avila_2019") & (meta$country == "Brazil"),] # 96

# barbieri 2012a
barbieri_2012a_burkinafaso <- meta[(meta$paper == "barbieri_2012a") & (meta$country == "Burkina Faso"),] # 291

# barbieri 2013a
barbieri_2013a_zambia <- meta[(meta$paper == "barbieri_2013a") & (meta$country == "Zambia"),] # 169

# barbieri 2013
barbieri_2013_botswana <- meta[(meta$paper == "barbieri_2013") & (meta$country == "Botswana"),] # 306
barbieri_2013_namibia <- meta[(meta$paper == "barbieri_2013") & (meta$country == "Namibia"),] # 141

# barbieri 2014a
barbieri_2014a_angola <- meta[(meta$paper == "barbieri_2014a") & (meta$country == "Angola"),] # 170
barbieri_2014a_zambia <- meta[(meta$paper == "barbieri_2014a") & (meta$country == "Zambia"),] # 446

# barbieri 2014b
barbieri_2014b_botswana <- meta[(meta$paper == "barbieri_2014b") & (meta$country == "Botswana"),] # 74
barbieri_2014b_namibia <- meta[(meta$paper == "barbieri_2014b") & (meta$country == "Namibia"),] # 117

# bodner 2015
bodner_2015_italy <- meta[(meta$paper == "bodner_2015") & (meta$country == "Italy"),] # 216

# brandini 2018
brandini_2018_ecuador <- meta[(meta$paper == "brandini_2018") & (meta$country == "Ecuador"),] # 208

# brucato 2018
brucato_2018_comoros <- meta[(meta$paper == "brucato_2018") & (meta$country == "Comoros"),] # 48
brucato_2018_kenya <- meta[(meta$paper == "brucato_2018") & (meta$country == "Kenya"),] # 228

# Chan 2019
chan_2019_namibia <- meta[(meta$paper == "chan_2019") & (meta$country == "Namibia"),] # 77
chan_2019_southafrica <- meta[(meta$paper == "chan_2019") & (meta$country == "South Africa"),] # 121

# Colombo 2025
colombo_2025_algeria <- meta[(meta$paper == "colombo_2025") & (meta$country == "Algeria"),] # 49
colombo_2025_cameroon <- meta[(meta$paper == "colombo_2025") & (meta$country == "Cameroon"),] # 29
colombo_2025_chad <- meta[(meta$paper == "colombo_2025") & (meta$country == "Chad"),] # 65
colombo_2025_libya <- meta[(meta$paper == "colombo_2025") & (meta$country == "Libya"),] # 60
colombo_2025_morocco <- meta[(meta$paper == "colombo_2025") & (meta$country == "Morocco"),] # 95
colombo_2025_tunisia <- meta[(meta$paper == "colombo_2025") & (meta$country == "Tunisia"),] # 34

# Fleskes 2023
fleskes_2023_anson <- meta[(meta$paper == "fleskes_2023"),] # 22

# Garcia 2021
garcia_2021_argentina <- meta[(meta$paper == "garcia_2021"),] # 114

# Garcia Olivares 2022
garciaolivares_2022_canaryislands <- meta[(meta$paper == "garciaolivares_2022"),] # 896

# Harney 2023
harney_2023_catoctin <- meta[(meta$paper == "harney_2023"),] # 22

# Huber 2025
huber_2025_peru <- meta[(meta$paper == "huber_2025"),] # 143

# Just 2008 african & hispanic american
just_2008_usa <- meta[(meta$paper == "just_2008"),] # 265
just_2008_aa <- meta[(meta$paper == "just_2008") & (meta$extra_info == "African American"),] # 140
just_2008_h <- meta[(meta$paper == "just_2008") & (meta$extra_info == "Hispanic"),] # 125

# Just 2015 african & hispanic american
just_2015_usa <- meta[(meta$paper == "just_2015"),] # 588
just_2015_aa <- meta[(meta$paper == "just_2015") & (meta$extra_info == "African American"),] # 170
just_2015_h <- meta[(meta$paper == "just_2015") & (meta$extra_info == "Hispanic"),] # 155
just_2015_c <- meta[(meta$paper == "just_2015") & (meta$extra_info == "caucasian"),] # 263

# kumar 2011: Mexican American from San Antonio Family Heart Study
kumar_2011_usa <- meta[(meta$paper == "kumar_2011"),] # 215

# mccrow 2016
mccrow_2016_southafrica <- meta[(meta$paper == "mccrow_2016"),] # 87

# olivera 2018
oliveira_2018_angola <- meta[(meta$paper == "oliveira_2018"),] # 295

# olivieri 2017
olivieri_2017_italy <- meta[(meta$paper == "olivieri_2017"),] # 2157

# pialq
pialq_peru <- meta[(meta$paper == "pialq"),] # 58

# pierron 2017
pierron_2017_madagascar <- meta[(meta$paper == "pierron_2017"),] # 2849

# sandoval velasco 2023
sandovalvelasco_2023_sthelena <- meta[(meta$paper == "sandovalvelasco_2023"),] # 19

# silva et al
silva_2021_portugal <- meta[(meta$paper == "silva_2021") & (meta$country == "Portugal"),] #103
silva_2021_spain <- meta[(meta$paper == "silva_2021") & (meta$country == "Spain"),] #1023

# taylor et al
taylor_2020_usa <- meta[(meta$paper == "taylor_2020"),] # 1327

# tito un_pub
tito_unpub_peru <- meta[(meta$paper == "tito_unpub"),] # 192

# tribaldos 2021, may need to separate
tribaldos_2021_panama <- meta[(meta$paper == "tribaldos_2021"),] # 84

# 6. function to create subsets of fasta files by population
# function
subset_align_by_pop <- function(alignment, pop_meta) {
  samples <- pop_meta$Sample # get sample names
  matching_indices <- which(rownames(alignment) %in% samples) # get indices of alignment seqs
  pop_alignment <- alignment[matching_indices, ] # subset alignment
  return(pop_alignment)
}

# list of pops
pop_list <- list(
  # 1kg populations
  acb = acb,
  asw = asw,
  beb = beb,
  cdx = cdx,
  ceu = ceu,
  chb = chb,
  chs = chs,
  clm = clm,
  esn = esn,
  fin = fin,
  gbr = gbr,
  gih = gih,
  gwd = gwd,
  ibs = ibs,
  itu = itu,
  jpt = jpt,
  khv = khv,
  lwk = lwk,
  msl = msl,
  mxl = mxl,
  pel = pel,
  pjl = pjl,
  pur = pur,
  stu = stu,
  tsi = tsi,
  yri = yri,
  
  # other populations
  mozabite_hgdp = mozabite_hgdp,
  arias_2017_colombia = arias_2017_colombia,
  avila_2019_brazil = avila_2019_brazil,
  barbieri_2012a_burkinafaso = barbieri_2012a_burkinafaso,
  barbieri_2013a_zambia = barbieri_2013a_zambia,
  barbieri_2013_botswana = barbieri_2013_botswana,
  barbieri_2013_namibia = barbieri_2013_namibia,
  barbieri_2014a_angola = barbieri_2014a_angola,
  barbieri_2014a_zambia = barbieri_2014a_zambia,
  barbieri_2014b_botswana = barbieri_2014b_botswana,
  barbieri_2014b_namibia = barbieri_2014b_namibia,
  bodner_2015_italy = bodner_2015_italy,
  brandini_2018_ecuador = brandini_2018_ecuador,
  brucato_2018_comoros = brucato_2018_comoros,
  brucato_2018_kenya = brucato_2018_kenya,
  chan_2019_namibia = chan_2019_namibia,
  chan_2019_southafrica = chan_2019_southafrica,
  colombo_2025_algeria = colombo_2025_algeria,
  colombo_2025_chad = colombo_2025_chad,
  colombo_2025_libya = colombo_2025_libya,
  colombo_2025_tunisia = colombo_2025_tunisia,
  fleskes_2023_anson = fleskes_2023_anson,
  garcia_2021_argentina = garcia_2021_argentina,
  garciaolivares_2022_canaryislands = garciaolivares_2022_canaryislands,
  harney_2023_catoctin = harney_2023_catoctin,
  huber_2025_peru = huber_2025_peru,
  just_2008_usa = just_2008_usa,
  just_2008_aa = just_2008_aa,
  just_2008_h = just_2008_h,
  just_2015_usa = just_2015_usa,
  just_2015_aa = just_2015_aa,
  just_2015_h = just_2015_h,
  just_2015_c = just_2015_c,
  kumar_2011_usa = kumar_2011_usa,
  mccrow_2016_southafrica = mccrow_2016_southafrica,
  oliveira_2018_angola = oliveira_2018_angola,
  olivieri_2017_italy = olivieri_2017_italy,
  pialq_peru = pialq_peru,
  pierron_2017_madagascar = pierron_2017_madagascar,
  sandovalvelasco_2023_sthelena = sandovalvelasco_2023_sthelena,
  silva_2021_portugal = silva_2021_portugal,
  silva_2021_spain = silva_2021_spain,
  taylor_2020_usa = taylor_2020_usa,
  tito_unpub_peru = tito_unpub_peru,
  tribaldos_2021_panama = tribaldos_2021_panama
)

# apply function to pop list
pop_alignments <- lapply(pop_list, function(pop_meta) {
  subset_align_by_pop(align_dedup, pop_meta)
})

# check sizes
pop_sizes <- lapply(pop_alignments, nrow)

# 7. get number of segregating sites per population with ape 
# seg.sites(pop, strict = FALSE, trailingGapsAsN = TRUE)
seg_sites <- lapply(pop_alignments, function(x) {
  seg.sites(x, strict = TRUE, trailingGapsAsN = TRUE)
})

seg_sites_count <- sapply(seg_sites, length)

# 8. get number of unique haplotypes per population with pegas
# haplotype(x, labels = NULL, strict = FALSE, trailingGapsAsN = TRUE, ...)
haps <- lapply(pop_alignments, function(x) {
  pegas::haplotype(x, strict = TRUE, trailingGapsAsN = TRUE)
})

# Apply to all populations
hap_stats <- lapply(names(pop_alignments), function(pop_name) {
  x <- pop_alignments[[pop_name]]
  h <- pegas::haplotype(x, strict = TRUE, trailingGapsAsN = TRUE)
  
  list(
    haplotypes = h,
    n_haplotypes = nrow(h)
  )
})
names(hap_stats) <- names(pop_alignments)

# extract n & freqs of haplotypes
n_haplotypes <- sapply(hap_stats, function(x) x$n_haplotypes)

# build summary table
hap_summary_table <- data.frame(
  population = names(hap_stats),
  n_sequences = sapply(pop_alignments, nrow),
  n_haplotypes = sapply(hap_stats, function(x) x$n_haplotypes)
)

write.table(hap_summary_table, "all_hap_count.csv", row.names = FALSE, sep = "\t", quote = FALSE)

# 9. get haplotype diversity for each population with pegas
# hap.div(x, variance = FALSE, method = "Nei", ...)
hap_diversity <- lapply(pop_alignments, function(x) {
  hap.div(x, variance = TRUE)
})

# 10. get nucleotide diversity with variance w pairwise deletion with ape
#nuc_diversity <- lapply(pop_alignments, function(x) {
#  nuc.div(x, variance = TRUE, pairwise.deletion = TRUE)
#})

nuc_diversity <- list()
for(pop_name in names(pop_alignments)) {
  cat("Processing:", pop_name, "- n =", nrow(pop_alignments[[pop_name]]), "\n")
  x <- pop_alignments[[pop_name]]
  
  nuc_diversity[[pop_name]] <- tryCatch({
    result <- pegas::nuc.div(x, variance = TRUE, pairwise.deletion = FALSE)
    cat("  Done! Pi =", result[1], "\n")
    result
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(c(NA, NA))
  })
  
  # Save progress after each population
  saveRDS(nuc_diversity, "pi_progress.rds")
}

# 11. tajima's d for each population
#taj_d <- lapply(pop_alignments, function(x) {
#  tajima.test(x)
#})

taj_d <- list()
for(pop_name in names(pop_alignments)) {
  cat("Processing:", pop_name, "- n =", nrow(pop_alignments[[pop_name]]), "\n")
  x <- pop_alignments[[pop_name]]
  
  taj_d[[pop_name]] <- tryCatch({
    result <- pegas::tajima.test(x)
    cat("  Done! D =", result$D, "\n")
    result
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(list(D = NA, Pval.normal = NA, Pval.beta = NA))
  })
  
  # Save progress after each population
  saveRDS(taj_d, "taj_progress.rds")
}

# Create a summary table of all results
summary_stats <- data.frame(
  population = names(pop_alignments),
  n_sequences = sapply(pop_alignments, nrow),
  segregating_sites = seg_sites_count,  # Now this is 80 values, not 95,444
  n_haplotypes = sapply(hap_stats, function(x) x$n_haplotypes),
  haplotype_diversity = sapply(hap_diversity, function(x) x[1]),
  hap_div_variance = sapply(hap_diversity, function(x) x[2]), 
  nucleotide_diversity = sapply(nuc_diversity, function(x) x[1]),
  nuc_div_variance = sapply(nuc_diversity, function(x) x[2]),
  tajimas_d = sapply(taj_d, function(x) x$D),
  tajimas_d_pval_normal = sapply(taj_d, function(x) x$Pval.normal),
  tajimas_d_pval_beta = sapply(taj_d, function(x) x$Pval.beta)
)

write.table(summary_stats, "intrapop_summ.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# plotting
# pi vs sample size
pi_n <- ggplot(summary_stats,
               aes(x = n_sequences, y = nucleotide_diversity)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    x = "Sample size (n)",
    y = expression("Nucleotide diversity (" * pi * ")")
  ) +
  theme_classic()

# plot of pi per population with a bar plot
pi_bar <- ggplot(summary_stats, aes(x = reorder(population,nucleotide_diversity), y = nucleotide_diversity)) +
  geom_col() +
  labs(
    x = "Population",
    y = expression("Nucleotide diversity (" * pi * ")")
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# bixplot
ggplot(summary_stats, aes(x = factor(1), y = nucleotide_diversity)) +
  geom_boxplot(width = 0.4, fill = "white") +
  geom_point(data = summary_stats[summary_stats$population == "pialq_peru", ],
             aes(x = factor(1), y = nucleotide_diversity), 
             color = "black", size = 2.5) +  # Highlight pialq
  theme_bw() + 
  labs(y = expression("Nucleotide diversity (" * pi * ")"))

# pi violin plot
pi_violin <- ggplot(summary_stats, aes(x = factor(1), y = nucleotide_diversity)) +
  geom_violin(fill = "lightsteelblue3", alpha = 0.8) +
  geom_point(data = summary_stats[summary_stats$population == "pialq_peru", ],
             aes(x = factor(1), y = nucleotide_diversity), 
             color = "black", size = 2.5) +  # Highlight pialq
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  labs(y = expression("Nucleotide diversity (" * pi * ")"))

# plot of h vs sample size
h_n <- ggplot(summary_stats,
              aes(x = n_sequences, y = haplotype_diversity)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    x = "Sample size (n)",
    y = expression("Haplotype diversity (Hd)")
  ) +
  theme_classic()

# h barplot
h_bar <- ggplot(summary_stats, aes(x = reorder(population, haplotype_diversity), y = haplotype_diversity)) +
  geom_col() +
  labs(
    x = "Population",
    y = expression("Haplotype diversity (Hd)")
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# plot of h per population with a boxplot
ggplot(summary_stats, aes(x = factor(1), y = haplotype_diversity)) +
  geom_boxplot(width = 0.4, fill = "white") +
  geom_point(data = summary_stats[summary_stats$population == "pialq_peru", ],
             aes(x = factor(1), y = haplotype_diversity), 
             color = "black", size = 2.5) +  # Highlight pialq
  theme_bw() + 
  labs(y = expression("Haplotype diversity (Hd)"))

# H violin plot
h_violin <- ggplot(summary_stats, aes(x = factor(1), y = haplotype_diversity)) +
  geom_violin(fill = "lightsteelblue3", alpha = 0.8) +
  geom_point(data = summary_stats[summary_stats$population == "pialq_peru", ],
             aes(x = factor(1), y = haplotype_diversity), 
             color = "black", size = 2.5) +  # Highlight pialq
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  labs(y = expression("Haplotype diversity (Hd)"))

# save to pdfs
pdf("intrapop_plots.pdf", width = 8, height = 11)
pi_n
pi_bar
pi_violin
h_n
h_bar
h_violin
dev.off()

# combine all sequences used in calculations
all_pop_sequences <- do.call(rbind, pop_alignments)
all_pop_seq <- all_pop_sequences[!duplicated(rownames(all_pop_sequences)), ]
nrow(all_pop_seq)

# write combined sequences to FASTA
write.dna(all_pop_seq, 
          file = "div_final.fasta", 
          format = "fasta",
          nbcol = -1)

# also save metadata for these sequences
all_pop_samples <- rownames(all_pop_sequences)
all_pop_meta <- meta[meta$Sample %in% all_pop_samples, ]
write.table(all_pop_meta, 
            "div_final.meta", 
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

# make text file that has population and file
# save population assigments
pop_assign_list <- lapply(names(pop_list), function(pop_name) {
  data.frame(
    Sample = pop_list[[pop_name]]$Sample,
    Population = pop_name,
    stringsAsFactors = FALSE
  )
})

pop_assign_df <- do.call(rbind, pop_assign_list)
rownames(pop_assign_df) <- NULL

write.table(
  pop_assign_df,
  "pop_assign.txt",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)
```



**interpop:**

create environment for r packages:

```
module load miniforge

conda create --copy -p $SHARED/projects/PIALQ/2025_ancient_lp/12_diversity/st_env r-base=4.3 r-ape r-pegas r-hierfstat r-codetools -c conda-forge -c bioconda

# activate env
source activate /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/st_env

# install haplotypes package afterwards
R --vanilla
install.packages("haplotypes")
```

**interpop**

```
#!/bin/bash
#SBATCH --job-name=interpop
#SBATCH --time=70:00:00
#SBATCH --mem=80gb
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

# Load conda
module load conda

# Activate your environment
source activate $SHARED/projects/PIALQ/2025_ancient_lp/12_diversity/st_env

# go to wd
cd /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/04_inter

# run r script
Rscript --vanilla pval.r

```

pval.r

```
# pval.r
# Pairwise PhiST + MDS
# with permutation testing (nperm = 1000)

# 1. load libraries
library(ape)
library(haplotypes)

# 2. logging helper
log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
  flush.console()
}

# 3. read deta: alignment & metadata
log_msg("loading data...")

tryCatch({
  align <- read.dna("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/03_intra/div_final.fasta", 
                    format = "fasta")
  meta  <- read.delim("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/03_intra/div_final.meta", 
                      stringsAsFactors = FALSE)
  log_msg(paste("Loaded", nrow(align), "sequences with", ncol(align), "sites"))
}, error = function(e) {
  stop(paste("ERROR loading data:", e$message))
})

# 4. load population assignments from file
log_msg("Loading population assignments from pop_assign.txt...")
if (!file.exists("pop_assign.txt")) {
  stop("ERROR: pop_assign.txt not found in working directory. Please run intrapop.r first to generate this file.")
}

pop_assignments <- tryCatch({
  read.delim("pop_assign.txt", stringsAsFactors = FALSE, header = FALSE)
}, error = function(e) {
  stop(paste("ERROR reading pop_assign.txt:", e$message))
})

names(pop_assignments) <- c("sample", "population")
log_msg(paste("Loaded", nrow(pop_assignments), "population assignments"))

# create lookup: sample name -> population
sample_to_pop <- setNames(pop_assignments$population, pop_assignments$sample)

# match populations to alignment samples
log_msg("matching samples to pops")
align_sample_names <- rownames(align)
align_to_pop <- sample_to_pop[align_sample_names]

# keep only samples with population assignments
valid <- !is.na(align_to_pop)
if (any(!valid)) {
  log_msg(paste("Dropping", sum(!valid), "samples with no population assignment:"))
  invalid_samples <- align_sample_names[!valid]
  for (s in head(invalid_samples, 10)) {
    cat("  ", s, "\n")
  }
  if (length(invalid_samples) > 10) {
    cat("  ... and", length(invalid_samples) - 10, "more\n")
  }
}

align <- align[valid, ]
align_to_pop <- align_to_pop[valid]
grouping <- factor(align_to_pop)

log_msg(paste("Using", nrow(align), "samples in", nlevels(grouping), "populations"))

# Check population sizes
pop_sizes <- table(grouping)
log_msg("\nPopulation sizes:")
for (pop_name in names(pop_sizes)) {
  log_msg(paste("  ", pop_name, ":", pop_sizes[pop_name], "samples"))
}

# 5. clean alignment to handle ambiguity codes
log_msg("\ncleaning alignment and handling ambiguity codes")

# convert DNAbin to character matrix
align_char <- as.character(align)

# replace ambiguity codes with N
# IUPAC ambiguity codes: R, Y, S, W, K, M, B, D, H, V
ambig_codes <- c("r", "y", "s", "w", "k", "m", "b", "d", "h", "v",
                 "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")

n_ambig_replaced <- 0
for (code in ambig_codes) {
  count <- sum(align_char == code, na.rm = TRUE)
  if (count > 0) {
    align_char[align_char == code] <- "n"
    n_ambig_replaced <- n_ambig_replaced + count
  }
}

if (n_ambig_replaced > 0) {
  log_msg(paste("Replaced", n_ambig_replaced, "ambiguity codes with 'n'"))
} else {
  log_msg("No ambiguity codes found")
}

# remove sites that are all gaps or all missing
log_msg("Filtering uninformative sites")
n_sites_orig <- ncol(align_char)

valid_sites <- apply(align_char, 2, function(col) {
  # Count actual bases (not gaps or missing)
  bases <- col[!is.na(col) & col != "-" & col != "n" & col != "?"]
  length(bases) > 0
})

align_char_clean <- align_char[, valid_sites]
n_sites_removed <- n_sites_orig - ncol(align_char_clean)
log_msg(paste("Retained", ncol(align_char_clean), "of", n_sites_orig, "sites"))
if (n_sites_removed > 0) {
  log_msg(paste("  (removed", n_sites_removed, "uninformative sites)"))
}

# 6. write cleaned alignment to temporary file
temp_fasta <- tempfile(fileext = ".fasta")
log_msg(paste("Writing cleaned alignment to temporary file"))

fasta_lines <- character()
for (i in 1:nrow(align_char_clean)) {
  fasta_lines <- c(fasta_lines, 
                   paste0(">", rownames(align_char_clean)[i]),
                   paste(align_char_clean[i, ], collapse = ""))
}
writeLines(fasta_lines, temp_fasta)

# 7. read cleaned fasta with haplotypes package
dna_hap <- tryCatch({
  read.fas(temp_fasta)
}, error = function(e) {
  unlink(temp_fasta)
  stop(paste("ERROR reading fasta with haplotypes package:", e$message))
})

log_msg(paste("Loaded", nrow(dna_hap), "sequences with", ncol(dna_hap), "sites"))

# Verify sample order matches
if (!all(rownames(dna_hap) == names(grouping))) {
  unlink(temp_fasta)
  stop("ERROR: Sample order mismatch between alignment and grouping")
}

# 7. Pairwise PhiST (with 1000 permutations)
log_msg("\n=== calculating pairwise PhiST ===")

phist_result <- tryCatch({
  pairPhiST(
    x = dna_hap,
    populations = grouping,
    nperm = 1000
  )
}, error = function(e) {
  unlink(temp_fasta)
  stop(paste("ERROR in pairPhiST calculation:", e$message))
})

log_msg("PhiST calculation complete!")

# clean up temporary file
unlink(temp_fasta)

# 8. extract and validate results
log_msg("\n=== Processing PhiST Results ===")

# extract distance matrix with validation
phist_dist <- NULL
if (is.list(phist_result) && "distance.matrix" %in% names(phist_result)) {
  log_msg("Extracting distance matrix from list element 'distance.matrix'")
  phist_dist <- phist_result$distance.matrix
} else if (is.matrix(phist_result)) {
  log_msg("Using phist_result directly as matrix")
  phist_dist <- phist_result
} else if (is.list(phist_result) && "PhiST" %in% names(phist_result)) {
  log_msg("Extracting distance matrix from list element 'PhiST'")
  phist_dist <- phist_result$PhiST
} else {
  log_msg("Attempting to convert result to matrix")
  phist_dist <- tryCatch({
    as.matrix(phist_result)
  }, error = function(e) {
    log_msg(paste("ERROR converting to matrix:", e$message))
    log_msg("Available result elements:")
    print(str(phist_result))
    stop("ERROR: Could not extract distance matrix from PhiST result")
  })
}

# validate distance matrix exists
if (is.null(phist_dist)) {
  stop("ERROR: Failed to extract PhiST distance matrix")
}

# Check for and report NA values
if (any(is.na(phist_dist))) {
  log_msg("WARNING: PhiST matrix contains NA values")
  na_count <- sum(is.na(phist_dist))
  log_msg(paste("  Number of NA values:", na_count))
}

log_msg("PhiST distance matrix validated successfully")

# 9. save results
write.table(phist_dist, "phist_matrix.txt", quote = FALSE, sep = "\t")
log_msg("Saved PhiST matrix to phist_matrix.txt")

# save p-values
if ("p" %in% names(phist_result)) {
  pval_matrix <- phist_result$p
  write.table(pval_matrix, "phist_pvalues.txt",
              quote = FALSE, sep = "\t")
  log_msg("Saved p-values to phist_pvalues.txt")

  # Summary of significance
  sig_count <- sum(pval_matrix[upper.tri(pval_matrix)] < 0.05, na.rm = TRUE)
  total_pairs <- sum(upper.tri(pval_matrix))
  log_msg(paste("  Significant comparisons (p < 0.05):", sig_count, "of", total_pairs))
} else {
  log_msg("WARNING: No p-values found in PhiST results")
}

# 10. MDS calculation
log_msg("\n=== performing metric MDS ===")

# Ensure matrix is symmetric and has no negative values
phist_dist_clean <- phist_dist
phist_dist_clean[is.na(phist_dist_clean)] <- 0
phist_dist_clean[phist_dist_clean < 0] <- 0
phist_dist_clean <- (phist_dist_clean + t(phist_dist_clean)) / 2

if (any(is.na(phist_dist))) {
  log_msg("NOTE: NA values in PhiST matrix were replaced with 0 for MDS")
}

mds_result <- tryCatch({
  cmdscale(phist_dist_clean, k = min(2, nlevels(grouping) - 1), eig = TRUE)
}, error = function(e) {
  stop(paste("ERROR in MDS calculation:", e$message))
})

# calculate variance explained
eig_vals <- mds_result$eig[mds_result$eig > 0]
if (length(eig_vals) == 0) {
  log_msg("WARNING: No positive eigenvalues in MDS - results may not be meaningful")
  var_exp <- c(0, 0)
} else {
  var_exp <- 100 * eig_vals / sum(eig_vals)
}

log_msg(paste("MDS Axis 1 explains", round(var_exp[1], 2), "% of variance"))
if (length(var_exp) > 1) {
  log_msg(paste("MDS Axis 2 explains", round(var_exp[2], 2), "% of variance"))
  log_msg(paste("Total variance explained:", round(sum(var_exp[1:2]), 2), "%"))
}

# 11. save MDS coordinates
n_dims <- ncol(mds_result$points)

mds_output <- data.frame(
    Population = rownames(mds_result$points),
    MDS1 = mds_result$points[, 1],
    MDS2 = mds_result$points[, 2],
    Variance1 = var_exp[1],
    Variance2 = ifelse(length(var_exp) > 1, var_exp[2], 0)
  )

write.table(mds_output, "mds_coordinates.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)
log_msg("Saved MDS coordinates to mds_coordinates.txt")

# 12. summary statistics
log_msg("\n=== Summary Statistics ===")
log_msg(paste("Total samples analyzed:", nrow(dna_hap)))
log_msg(paste("Number of populations:", nlevels(grouping)))
log_msg(paste("Informative sites used:", ncol(dna_hap)))
log_msg(paste("Number of pairwise comparisons:", 
              nlevels(grouping) * (nlevels(grouping) - 1) / 2))

# 13. final output summary
log_msg("\nAnalysis finished successfully!")


```





MDS script:

```
# plot MDS

# 1. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/04_inter")

# 2. load libraries, download if necessary
library(dplyr)
library(ggplot2)

# 3. plot mds using coordinates already calculated
# read in coordinates
mds_data <- read.table("all_mds_coordinates.txt", header = TRUE)

# MDS plot
ggplot(mds_data, aes(x = MDS1, y = MDS2, label = Population)) +
  geom_point(size = 1) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "MDS Plot",
       x = "MDS1",
       y = "MDS2")

# avoid overlaps
library(ggrepel)

final_pop <- ggplot(mds_data, aes(x = MDS1, y = MDS2, label = Population)) +
  geom_point(size = 1, color = "black") +
  geom_text_repel(size = 1, max.overlaps=Inf) +
  theme_bw() +
  labs(x = "MDS1",
       y = "MDS2")

# plot colored by region
# create new column for regions
unique(mds_data$Population)

# Create a mapping of populations to regions
# based on UN GeoScheme: https://www.emiw.org/fileadmin/emiw/UserActivityDocs/Geograph.Representation/Geographic-Representation-Appendix_1.pdf
region_mapping <- c(
  "acb" = "African Diaspora in Americas: Caribbean",
  "asw" = "African Diaspora in Americas: Northern",
  "beb" = "Asia: Southern",
  "cdx" = "Asia: Eastern",
  "ceu" = "Europe: Central",
  "chb" = "Asia: Eastern",
  "chs" = "Asia: Eastern",
  "clm" = "Americas: Southern",
  "esn" = "Africa: Western",
  "fin" = "Europe: Northern",
  "gbr" = "Europe: Northern",
  "gih" = "Asia: Southern",
  "gwd" = "Africa: Western",
  "ibs" = "Europe: Southern",
  "itu" = "Asia: Southern",
  "jpt" = "Asia: Eastern",
  "khv" = "Asia: Southeastern",
  "lwk" = "Africa: Eastern",
  "msl" = "Africa: Western",
  "mxl" = "Americas: Central",
  "pel" = "Americas: Southern",
  "pjl" = "Asia: Southern",
  "pur" = "Americas: Northern",
  "stu" = "Asia: Southern",
  "tsi" = "Europe: Southern",
  "yri" = "Africa: Western",
  "mozabite_hgdp" = "Africa: Northern",
  "arias_2017_colombia" = "Americas: Southern",
  "avila_2019_brazil" = "Americas: Southern",
  "barbieri_2012a_burkinafaso" = "Africa: Western",
  "barbieri_2012b_zambia" = "Africa: Eastern",
  "barbieri_2013_botswana" = "Africa: Southern",
  "barbieri_2013_namibia" = "Africa: Southern",
  "barbieri_2014a_zambia" = "Africa: Eastern",
  "barbieri_2014a_angola" = "Africa: Southern",
  "barbieri_2014b_botswana" = "Africa: Southern",
  "barbieri_2014b_namibia" = "Africa: Southern",
  "bodner_2022_italy" = "Europe: Southern",
  "brandini_2018_ecuador" = "Americas: Southern",
  "brucato_2018_kenya" = "Africa: Eastern",
  "brucato_2018_comoros" = "Africa: Eastern",
  "chan_2019_southafrica" = "Africa: Southern",
  "chan_2019_namibia" = "Africa: Southern",
  "colombo_2025_algeria" = "Africa: Northern",
  "colombo_2025_chad" = "Africa: Central",
  "colombo_2025_libya" = "Africa: Northern",
  "colombo_2025_tunisia" = "Africa: Northern",
  "fleskes_2023_anson" = "historical Afro-diasporic",
  "garcia_2021_argentina" = "Americas: Southern",
  "garciaolivares_2023_canaryislands" = "Europe: Southern",
  "harney_2023_catoctin" = "historical Afro-diasporic",
  "huber_2025_peru" = "Americas: Southern",
  "just_2008_usa" = "Americas: Northern",
  "just_2015_usa" = "Americas: Northern",
  "kumar_2011_usa" = "Americas: Northern",
  "mccrow_2016_southafrica" = "Africa: Southern",
  "oliveira_2018_angola" = "Africa: Southern",
  "olivieri_2017_italy" = "Europe: Southern",
  "pialq_peru" = "historical Afro-diasporic",
  "pierron_2017_madagascar" = "Africa: Eastern",
  "sandovalvelasco_2023_sthelena" = "historical Afro-diasporic",
  "silva_2021_spain" = "Europe: Southern",
  "silva_2021_portugal" = "Europe: Southern",
  "taylor_2020_usa" = "Americas: Northern",
  "tito_unpub_peru" = "Americas: Southern",
  "tribaldos_2021_panama" = "Americas: Central"
)

# Add the region column to data
mds_data <- mds_data %>%
  mutate(region = region_mapping[Population])

# plot colored by region
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw() +
  labs(title = "MDS Plot by Region",
       x = "MDS1",
       y = "MDS2")


# custom colors
custom_colors <- c(
  "Africa: Central" = "#006400",
  "Africa: Eastern" = "#9ACD32",
  "Africa: Northern" = "#71BC78",
  "Africa: Southern" = "#4CBB17",  
  "Africa: Western" = "#4B5320",
  "African Diaspora in Americas: Caribbean" = "#39FF14",
  "African Diaspora in Americas: Northern" = "#D0F0C0",
  "Americas: Central" = "#EE82EE",
  "Americas: Northern" = "#4B0082",
  "Americas: Southern" = "#800080",
  "Europe: Southern" = "#1E90FF",
  "Europe: Northern" = "#00008B",
  "Europe: Central" = "#AFEEEE",
  "Asia: Eastern" = "#F9629F",
  "Asia: Southern" = "#eec0c8",
  "historical Afro-diasporic" = "black"
)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  labs(x = "MDS1",
       y = "MDS2")

# size points by sample size
size_mapping <- c(
  "acb" = 98,
  "asw" = 66,
  "beb" = 104,
  "cdx" = 107,
  "ceu" = 101,
  "chb" = 105,
  "chs" = 112,
  "clm" = 103,
  "esn" = 111,
  "fin" = 105,
  "gbr" = 106,
  "gih" = 109,
  "gwd" = 120,
  "ibs" = 108,
  "itu" = 112,
  "jpt" = 105,
  "khv" = 103,
  "lwk" = 112,
  "msl" = 98,
  "mxl" = 70,
  "pel" = 90,
  "pjl" = 108,
  "pur" = 107,
  "stu" = 111,
  "tsi" = 111,
  "yri" = 109,
  "mozabite_hgdp" = 26,
  "arias_2017_colombia" = 436,
  "avila_2019_brazil" = 96,
  "barbieri_2012a_burkinafaso" = 291,
  "barbieri_2012b_zambia" = 169,
  "barbieri_2013_botswana" = 306,
  "barbieri_2013_namibia" = 141,
  "barbieri_2014a_zambia" = 446,
  "barbieri_2014a_angola" = 170,
  "barbieri_2014b_botswana" = 74,
  "barbieri_2014b_namibia" = 117,
  "bodner_2022_italy" = 216,
  "brandini_2018_ecuador" = 208,
  "brucato_2018_kenya" = 228,
  "brucato_2018_comoros" = 48,
  "chan_2019_southafrica" = 121,
  "chan_2019_namibia" = 77,
  "colombo_2025_algeria" = 49,
  "colombo_2025_chad" = 65,
  "colombo_2025_libya" = 60,
  "colombo_2025_tunisia" = 34,
  "fleskes_2023_anson" = 22,
  "garcia_2021_argentina" = 114,
  "garciaolivares_2023_canaryislands" = 896,
  "harney_2023_catoctin" = 22,
  "huber_2025_peru" = 143,
  "just_2008_usa" = 265,
  "just_2015_usa" = 588,
  "kumar_2011_usa" = 215,
  "mccrow_2016_southafrica" = 87,
  "oliveira_2018_angola" = 295,
  "olivieri_2017_italy" = 2157,
  "pialq_peru" = 58,
  "pierron_2017_madagascar" = 2849,
  "sandovalvelasco_2023_sthelena" = 19,
  "silva_2021_spain" = 1023,
  "silva_2021_portugal" = 103,
  "taylor_2020_usa" = 1327,
  "tito_unpub_peru" = 189,
  "tribaldos_2021_panama" = 84
)

# add sample size
mds_data <- mds_data %>%
  mutate(n = size_mapping[Population])

# replot
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region, size = n, label = Population)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  labs(x = "MDS1 (61.9% of variance)",
       y = "MDS2 (17.9% of variance)")

# w text labels
final_region <- ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region, size = n, label = Population)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(size = 1, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  labs(x = "MDS1 (61.9% of variance)",
       y = "MDS2 (17.9% of variance)")

# save final mds plot
pdf("all_mds.pdf", width = 11, height = 8)
final_pop
final_region
dev.off()

```

script for more chill mds:

```
# plot MDS with continents

# 1. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/04_inter")

# 2. load libraries, download if necessary
library(dplyr)
library(ggplot2)

# 3. plot mds using coordinates already calculated
# read in coordinates
mds_data <- read.table("all_mds_coordinates.txt", header = TRUE)

# MDS plot
ggplot(mds_data, aes(x = MDS1, y = MDS2, label = Population)) +
  geom_point(size = 1) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "MDS Plot",
       x = "MDS1",
       y = "MDS2")

# avoid overlaps
library(ggrepel)
ggplot(mds_data, aes(x = MDS1, y = MDS2, label = Population)) +
  geom_point(size = 1, color = "black") +
  geom_text_repel(size = 1, max.overlaps=Inf) +
  theme_bw() +
  labs(x = "MDS1",
       y = "MDS2")

# plot colored by region
# create new column for regions
unique(mds_data$Population)

# Create a mapping of populations to regions
# based on UN GeoScheme: https://www.emiw.org/fileadmin/emiw/UserActivityDocs/Geograph.Representation/Geographic-Representation-Appendix_1.pdf
region_mapping <- c(
  "acb" = "present-day African Diaspora",
  "asw" = "present-day African Diaspora",
  "beb" = "Asia",
  "cdx" = "Asia",
  "ceu" = "Europe",
  "chb" = "Asia",
  "chs" = "Asia",
  "clm" = "Americas",
  "esn" = "Africa",
  "fin" = "Europe",
  "gbr" = "Europe",
  "gih" = "Asia",
  "gwd" = "Africa",
  "ibs" = "Europe",
  "itu" = "Asia",
  "jpt" = "Asia",
  "khv" = "Asia",
  "lwk" = "Africa",
  "msl" = "Africa",
  "mxl" = "Americas",
  "pel" = "Americas",
  "pjl" = "Asia",
  "pur" = "Americas",
  "stu" = "Asia",
  "tsi" = "Europe",
  "yri" = "Africa",
  "mozabite_hgdp" = "Africa",
  "arias_2017_colombia" = "Americas",
  "avila_2019_brazil" = "Americas",
  "barbieri_2012a_burkinafaso" = "Africa",
  "barbieri_2012b_zambia" = "Africa",
  "barbieri_2013_botswana" = "Africa",
  "barbieri_2013_namibia" = "Africa",
  "barbieri_2014a_zambia" = "Africa",
  "barbieri_2014a_angola" = "Africa",
  "barbieri_2014b_botswana" = "Africa",
  "barbieri_2014b_namibia" = "Africa",
  "bodner_2022_italy" = "Europe",
  "brandini_2018_ecuador" = "Americas",
  "brucato_2018_kenya" = "Africa",
  "brucato_2018_comoros" = "Africa",
  "chan_2019_southafrica" = "Africa",
  "chan_2019_namibia" = "Africa",
  "colombo_2025_algeria" = "Africa",
  "colombo_2025_chad" = "Africa",
  "colombo_2025_libya" = "Africa",
  "colombo_2025_tunisia" = "Africa",
  "fleskes_2023_anson" = "historical African Diaspora",
  "garcia_2021_argentina" = "Americas",
  "garciaolivares_2023_canaryislands" = "Europe",
  "harney_2023_catoctin" = "historical African Diaspora",
  "huber_2025_peru" = "Americas",
  "just_2008_usa" = "Americas",
  "just_2015_usa" = "Americas",
  "kumar_2011_usa" = "Americas",
  "mccrow_2016_southafrica" = "Africa",
  "oliveira_2018_angola" = "Africa",
  "olivieri_2017_italy" = "Europe",
  "pialq_peru" = "THIS STUDY: Hacienda La Quebrada",
  "pierron_2017_madagascar" = "Africa",
  "sandovalvelasco_2023_sthelena" = "historical African Diaspora",
  "silva_2021_spain" = "Europe",
  "silva_2021_portugal" = "Europe",
  "taylor_2020_usa" = "Americas",
  "tito_unpub_peru" = "Americas",
  "tribaldos_2021_panama" = "Americas"
)

# Add the region column to data
mds_data <- mds_data %>%
  mutate(region = region_mapping[Population])

# plot colored by region
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region)) +
  geom_point(size = 2, alpha = 1) +
  theme_bw() +
  labs(title = "MDS Plot by Region",
       x = "MDS1",
       y = "MDS2")

# custom colors & shapes
custom_colors <- c(
  "Africa" = "#006400",
  "Americas" = "#800080",
  "Asia" = "#F9629F",
  "Europe" = "#1E90FF",
  "present-day African Diaspora" = "#39FF14",
  "historical African Diaspora" = "black",
  "THIS STUDY: Hacienda La Quebrada" = "red"
)

custom_shapes <- c(
  "Africa" = 20,
  "Americas" = 20,
  "Asia" = 20,
  "Europe" = 20,
  "present-day African Diaspora" = 20,
  "historical African Diaspora" = 20,
  "THIS STUDY: Hacienda La Quebrada" = 4
)

# size points by sample size
size_mapping <- c(
  "acb" = 98,
  "asw" = 66,
  "beb" = 104,
  "cdx" = 107,
  "ceu" = 101,
  "chb" = 105,
  "chs" = 112,
  "clm" = 103,
  "esn" = 111,
  "fin" = 105,
  "gbr" = 106,
  "gih" = 109,
  "gwd" = 120,
  "ibs" = 108,
  "itu" = 112,
  "jpt" = 105,
  "khv" = 103,
  "lwk" = 112,
  "msl" = 98,
  "mxl" = 70,
  "pel" = 90,
  "pjl" = 108,
  "pur" = 107,
  "stu" = 111,
  "tsi" = 111,
  "yri" = 109,
  "mozabite_hgdp" = 26,
  "arias_2017_colombia" = 436,
  "avila_2019_brazil" = 96,
  "barbieri_2012a_burkinafaso" = 291,
  "barbieri_2012b_zambia" = 169,
  "barbieri_2013_botswana" = 306,
  "barbieri_2013_namibia" = 141,
  "barbieri_2014a_zambia" = 446,
  "barbieri_2014a_angola" = 170,
  "barbieri_2014b_botswana" = 74,
  "barbieri_2014b_namibia" = 117,
  "bodner_2022_italy" = 216,
  "brandini_2018_ecuador" = 208,
  "brucato_2018_kenya" = 228,
  "brucato_2018_comoros" = 48,
  "chan_2019_southafrica" = 121,
  "chan_2019_namibia" = 77,
  "colombo_2025_algeria" = 49,
  "colombo_2025_chad" = 65,
  "colombo_2025_libya" = 60,
  "colombo_2025_tunisia" = 34,
  "fleskes_2023_anson" = 22,
  "garcia_2021_argentina" = 114,
  "garciaolivares_2023_canaryislands" = 896,
  "harney_2023_catoctin" = 22,
  "huber_2025_peru" = 143,
  "just_2008_usa" = 265,
  "just_2015_usa" = 588,
  "kumar_2011_usa" = 215,
  "mccrow_2016_southafrica" = 87,
  "oliveira_2018_angola" = 295,
  "olivieri_2017_italy" = 2157,
  "pialq_peru" = 58,
  "pierron_2017_madagascar" = 2849,
  "sandovalvelasco_2023_sthelena" = 19,
  "silva_2021_spain" = 1023,
  "silva_2021_portugal" = 103,
  "taylor_2020_usa" = 1327,
  "tito_unpub_peru" = 189,
  "tribaldos_2021_panama" = 84
)

# add sample size
mds_data <- mds_data %>%
  mutate(n = size_mapping[Population])

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region, shape = region)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  theme_bw() +
  labs(x = "MDS1",
       y = "MDS2")

# replot
simple <- ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region, shape = region, size = n, label = Population)) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  theme_bw() +
  labs(x = "MDS1 (61.9% of variance)",
       y = "MDS2 (17.9% of variance)",
       color = "Continent", 
       shape = "Continent")


# w text labels --> don't like this, very crazy
final <- simple +
  geom_text_repel(
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.2,
    segment.alpha = 0.7,
    box.padding = 0.5,
    point.padding = 0.5,
    force = 3,
    force_pull = 0.2,
    seed = 42
  )

# clean up
# fix legend order
region_order <- c(
  "Africa",
  "Americas", 
  "Asia",
  "Europe",
  "present-day African Diaspora",
  "historical African Diaspora",
  "THIS STUDY: Hacienda La Quebrada"
)

mds_data$region <- factor(mds_data$region, levels = region_order)

# label only a few reference pops
label_pops <- c("acb", "asw", "yri", "gwd", "msl", "esn", "lwk", "pel",
                "tito_unpub_peru", "huber_2025_peru",
                "barbieri_2014b_namibia", "colombo_2025_chad", "brucato_2018_comoros",
                "barbieri_2012a_burkinafaso", "pialq_peru",
                "fleskes_2023_anson", "harney_2023_catoctin",
                "sandovalvelasco_2023_sthelena")

label_map <- c(
  "acb" = "ACB",
  "asw" = "ASW",
  "yri" = "YRI",
  "gwd" = "GWD",
  "msl" = "MSL",
  "esn" = "ESN",
  "lwk" = "LWK",
  "pel" = "PEL",
  "tito_unpub_peru" = "Tito unpub Peru",
  "huber_2025_peru" = "Huber 2025\nPeru",
  "barbieri_2014b_namibia" = "Barbieri 2014\nNamibia",
  "colombo_2025_chad" = "Colombo 2025\nChad",
  "brucato_2018_comoros" = "Brucato 2018\nComoros",
  "barbieri_2012a_burkinafaso" = "Barbieri 2012\nBurkina Faso",
  "fleskes_2023_anson" = "Fleskes 2023\nAnson Street",
  "harney_2023_catoctin" = "Harney 2023\nCatoctin Furnace",
  "sandovalvelasco_2023_sthelena" = "Sandoval-Velasco\n2023 St. Helena",
  "pialq_peru" = "THIS STUDY"
)

mds_data <- mds_data %>%
  mutate(label_wrapped = label_map[Population])

final <- ggplot(mds_data, 
                aes(x = MDS1, y = MDS2, 
                    color = region, shape = region,
                    size = n, label = label_wrapped)) +
                geom_point(alpha = 0.9) +
                geom_text_repel(
                  data = subset(mds_data, Population %in% label_pops),
                  segment.color = "black",
                  size = 3,
                  direction = "both",
                  show.legend = FALSE,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.15,
                  segment.alpha = 1,
                  box.padding = 0.4,
                  point.padding = 0.01,
                  force = 10,
                  seed = 42
                ) +
                scale_color_manual(values = custom_colors, breaks = region_order) +
                scale_shape_manual(values = custom_shapes, breaks = region_order) +
                scale_size_continuous(range = c(1, 5), guide = "legend") +
                theme_bw() +
                theme(
                  legend.title = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 14),
                  legend.key.size = unit(0.4, "cm"),
                  axis.title = element_text(size = 16),
                  axis.text = element_text(size = 14),
                  panel.grid.minor = element_blank()
                ) +
                labs(
                  x = "MDS1 (61.9% of variance)",
                  y = "MDS2 (17.9% of variance)",
                  color = "Continent",
                  shape = "Continent",
                  size = "Paper Sample Size"
                )

# save final mds plot
pdf("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/mds_paper.pdf", width = 11, height = 8)
final
dev.off()


```



phist script:

```
# load libraries
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

if (!requireNamespace("viridis", quietly = TRUE)) {
  install.packages("viridis")
}

if (!requireNamespace("grid", quietly = TRUE)) {
  install.packages("grid")
}

library(pheatmap)
library(viridis)
library(grid)

# 1. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/04_inter")

# 2. load phist matrix & p value smatrix
phist_matrix <- as.matrix(read.table("all_phist_matrix.txt", header = TRUE, row.names = 1))
pval_matrix <- as.matrix(read.table("all_phist_pvalues.txt", header = TRUE, row.names = 1))

# 3. visualize heatmap of phist values
phist <- pheatmap(phist_matrix, 
                  cluster_rows = FALSE, 
                  cluster_cols = FALSE, 
                  border_color = "darkgray", 
                  na_col = "white",
                  color = viridis(100, option = "viridis", direction = -1),
                  main = "Pairwise PhiST values"
)

# 4. visualize heatmap of p values
pval <- pheatmap(pval_matrix, 
                 cluster_rows = FALSE, 
                 cluster_cols = FALSE,
                 border_color = "darkgray", 
                 na_col = "white",
                 color = viridis(100, option = "viridis", direction = -1),
                 main = "P-values of pairwise PhiST comparisons"
)

# 5. visualize heatmap of significant phist values
# 16525 individuals in 66 populations
# bonferroni corrected  66 × 65 / 2 = 2145 comparison

# bonferroni correction for 1,225 comparisons
bonf_threshold <- 0.05 / 2145
bonf_threshold

# create binary matrix for yes/no significance
sig_binary <- ifelse(pval_matrix < bonf_threshold, 1, 0)

pval_bin <- pheatmap(sig_binary, 
                     cluster_rows = FALSE, 
                     cluster_cols = FALSE,
                     border_color = "black", 
                     na_col = "black",
                     color = c("lightgray", "red"),  # Not sig = gray, Sig = red
                     legend_breaks = c(0, 1),
                     legend_labels = c("Not Significant", "Significant"),
                     main = "Pairwise PhiST Significance (a = 2.33e-05)"
)

# harney & fleskes not significant
# harney p
pval_matrix["pialq_peru", "harney_2023_catoctin"]
# 0.02197802

# fleskes p
pval_matrix["pialq_peru", "fleskes_2023_anson"]
# 0.6643357

# sandoval velasco p
pval_matrix["sandovalvelasco_2023_sthelena", "pialq_peru"]
# 0.3026973

# save final plots
pdf("all_phist_heatmap.pdf", width = 12, height = 12)
print(phist)
grid.newpage()
print(pval)
grid.newpage()
print(pval_bin) 
dev.off()


```





counting Ls in south america

```
import pandas as pd

df = pd.read_csv("all_merged.meta", sep = "\t")

# create county subsets
peru = df[(df["country"] == "Peru") & (df["hap_1char"] == "L") & (df["paper"] != "pialq")] # 10

col = df[(df["country"] == "Colombia") & (df["hap_1char"] == "L")] # 0
ven = df[(df["country"] == "Venezuela") & (df["hap_1char"] == "L")] # 0
ecua = df[(df["country"] == "Ecuador") & (df["hap_1char"] == "L")] # 0
arg = df[(df["country"] == "Argentina") & (df["hap_1char"] == "L")] # 0
chile = df[(df["country"] == "Chile") & (df["hap_1char"] == "L")] # 0
brazil = df[(df["country"] == "Brazil") & (df["hap_1char"] == "L")] # 0

# results
>>> peru
                 file          paper ancient_or_modern original_format  ... hap_6char hap_7char hap_8char length
4519   PQ827220.fasta     huber_2025            modern           fasta  ...    L1c2b1   L1c2b1b   L1c2b1b  16570
4525   PQ827226.fasta     huber_2025            modern           fasta  ...    L2b1a3    L2b1a3    L2b1a3  16568
4529   PQ827230.fasta     huber_2025            modern           fasta  ...    L3d3a1   L3d3a1a   L3d3a1a  16568
4550   PQ827251.fasta     huber_2025            modern           fasta  ...    L1c1d1    L1c1d1    L1c1d1  16570
4576   PQ827277.fasta     huber_2025            modern           fasta  ...    L3d1a1   L3d1a1a   L3d1a1a  16568
4615   PQ827316.fasta     huber_2025            modern           fasta  ...    L3b1a1   L3b1a1a   L3b1a1a  16568
4619   PQ827320.fasta     huber_2025            modern           fasta  ...    L3e1a2    L3e1a2    L3e1a2  16568
4620   PQ827321.fasta     huber_2025            modern           fasta  ...    L3e1a2    L3e1a2    L3e1a2  16568
4637   PQ827338.fasta     huber_2025            modern           fasta  ...     L2a1f     L2a1f     L2a1f  16568
10189  MG571168.fasta  barbieri_2017            modern           fasta  ...    L3e1a1   L3e1a1a   L3e1a1a  16569
```

goodness of fit of mds:

```
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/04_inter")


phist_dist_clean <- as.matrix(read.table("phist_matrix.txt", header = TRUE, sep = "\t"))

# make symmetric and remove negatives (same as original script)
phist_dist_clean[is.na(phist_dist_clean)] <- 0
phist_dist_clean[phist_dist_clean < 0] <- 0
phist_dist_clean <- (phist_dist_clean + t(phist_dist_clean)) / 2

# re-run MDS — takes seconds
mds_result <- cmdscale(phist_dist_clean, k = 2, eig = TRUE)

# get GOF
mds_result$GOF

# and variance explained
eig_vals <- mds_result$eig[mds_result$eig > 0]
var_exp <- 100 * eig_vals / sum(eig_vals)
var_exp[1:2]
```

