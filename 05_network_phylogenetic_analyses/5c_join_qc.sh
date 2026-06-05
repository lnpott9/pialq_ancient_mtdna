# =============================================================================
# Title: 5c_join_qc.sh
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: Python script to join reference metadata with haplogroup results and QC
# Usage: run interactively
# =============================================================================

srun --mem=10gb --time=5:00:00 --pty bash

# 1. activate environment
module load conda/python3
source activate $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/msa_env

cd $SHARED/projects/PIALQ/2025_ancient_lp/09_msa/

python

# 2. import pandas
import pandas as pd

# 3. read in short_all_haplogrep_summary.txt as dataframe
hg = pd.read_csv("short_all_hap_summ.txt", sep="\t")

# 4. look at head to make sure all is well
hg.head()
hg.shape #(18540, 5)

# 5. read in metadata file
md = pd.read_csv("all.meta", sep="\t")
md.shape #(18566, 7)

# 6. add new column to md df of sample name without fasta
md["Sample"] = md["file"].str.replace(".fasta", "", regex=False)

# 7. ensure haplogrep & metadata are same length
assert hg.shape[0] == md.shape[0]

# Samples in md that don't have haplogroup data
missing_in_hg = set(md["Sample"]) - set(hg["Sample"])
# remove these samples from md
md = md[~md["Sample"].isin(missing_in_hg)]

# 8. read in fasta files
from pathlib import Path
from Bio import SeqIO

# directories w fasta files
ancient_dir = Path("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/ancient_fasta")
modern_dir = Path("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta")

# find all fasta files in both directories
fasta_files = list(ancient_dir.rglob("*.fasta")) + list(modern_dir.rglob("*.fasta"))

# read all sequences into a dataframe
fasta_list = []
for fasta_file in fasta_files:
    df = pd.DataFrame([{
        "sequence_id": record.id,
        "description": record.description,
        "sequence": str(record.seq),
        "file": fasta_file.name
    } for record in SeqIO.parse(fasta_file, "fasta")])
    fasta_list.append(df)

# 9. combine all fasta dataframes
all_fasta = pd.concat(fasta_list, ignore_index=True)

# 10. add sample name
all_fasta["Sample"] = all_fasta["file"].str.replace(".fasta", "", regex=False)

# 11. count number of characters (ATCG)
all_fasta["length"] = all_fasta["sequence"].str.upper().str.count(r"[ACGT]")
# look at distribution
all_fasta["length"].describe()
count    18541.000000
mean     16534.996171
std        518.195750
min         14.000000
25%      16564.000000
50%      16568.000000
75%      16569.000000
max      16659.000000
Name: length, dtype: float64

fasta_lengths = all_fasta[["Sample", "length"]]

# 12. join cleaned haplogrep & metadata text files based on sample column
merged = pd.merge(md, hg, on="Sample", how="inner")

merged.shape
(18540, 12)

# 13. replace Haplogroup for PIALQ C samples that were manually confirmed to be C1b
samples_to_update = ["LQ06", "LQ23", "LQ60", "LQ62"]

# update Haplogroup to C1b
merged.loc[merged["Sample"].isin(samples_to_update), "Haplogroup"] = "C1b"

# 14. add row to describe multiple character length haplogroups
def get_nchar(haplo, n):
    if pd.isna(haplo):
        return None
    return haplo[:n]
    
for n in [1, 2, 3, 4, 5, 6, 7, 8]:
    merged[f"hap_{n}char"] = merged["Haplogroup"].apply(lambda x: get_nchar(x, n))

# 15. add fasta lengths to df
# join fasta length with merged
missing_in_hg_id = set(merged["Sample"]) - set(all_fasta["Sample"])
missing_in_sq_id = set(all_fasta["Sample"]) - set(merged["Sample"])

combo = pd.merge(merged, fasta_lengths, on="Sample", how="inner")

# 16. fix some of the incorrect metadata
combo['paper'] = combo['paper'].replace({
    'bodner_2023': 'bodner_2012',
    'valverde_2016': 'llamas_2016',
    'garciaolivares_2023': 'garciaolivares_2022',
    'barbieri_2012b': 'barbieri_2013a'
})

# 17. add country labels & join
for country in sorted(one['country'].unique()):
    print(country)
    
labels = pd.read_csv('country_labels.txt', sep='\t')
new = combo.merge(labels, on='country', how='left')

# 18. change all that start with LQ as "Hacienda La Quebrada"
new.loc[new['file'].str.startswith('LQ'), 'country_label'] = 'Hacienda La Quebrada'
final = new

# 19. save all merged data
final.to_csv("all_merged.meta", sep="\t", index=False, header=True)

# 20. remove anything that isnt A, B, C, D, H, L
subset = final[final["hap_2char"].isin(["A2", "B2", "C1", "D1", "H1", "L0", "L1", "L2", "L3"])].copy()

subset.shape # (10065, 21)

# 21. only keep samples longer than 14000 characters AND quality over 80%, or if it's a pialq sequence
quality = subset[((subset["length"] > 14000) & (subset["Quality"] > 0.8)) | (subset["paper"] == "pialq")].copy()
remove = subset[((subset["length"] < 14000) | (subset["Quality"] < 0.8) & (subset["paper"] != "pialq"))].copy()
remove.shape # (93, 21)

# 22. only keep samples that are in target haplogroups
# A (target: A2 and A2+(64))
# B (target: B2 and B2b)
# C (target: C1 and C1b)
# D (target: D1)
# H (target: H1 and H1b)
# L0 (targets: L0a2a, L0d1a)
# L1 (targets: L1c, L1c1, L1c2'4)
# L2 (targets: L2a1, L2a1a3, L2a1d, L2a1f, L2a1q, L2b1, L2c
# L3 (targets: L3b1a, L3d1, L3d3, L3d4, L3d3 & L3d4, L3e1a1a, L3e2, L3e3, L3f1b)

# define keep_list
keep_2 = ["A2", "B2", "C1", "D1", "H1", "L0", "L1", "L2", "L3"]
keep_3 = ["B2b", "C1b", "L1c", "L2c"]
keep_4 = ["H1bw", "L1c1", "L1c2", "L1c4", "L2a1", "L2b1", "L3d1", "L3d3", "L3d4", "L3e2", "L3e3"]
keep_5 = ["L0a2a", "L0d1a", "L2a1d", "L2a1f", "L2a1q", "L3b1a", "L3f1b"]
keep_6 = ["L2a1a3"]
keep_7 = ["A2+(64)", "L3e1a1a"]

# 23. subset all by samples that have a desired haplogroup
small = quality[(quality["hap_2char"].isin(keep_2)) | (quality["hap_3char"].isin(keep_3)) | (quality["hap_4char"].isin(keep_4)) | (quality["hap_5char"].isin(keep_5)) | (quality["hap_6char"].isin(keep_6)) | (quality["hap_7char"].isin(keep_7))]

small.shape #(9987, 21)

# 24. deduplicate
small_final = small.drop_duplicates("Sample", keep="first")

# 25. write new meta (samples of interest & outgroups)
small_final.to_csv("final_cut.meta", sep="\t", index=False, header=True)