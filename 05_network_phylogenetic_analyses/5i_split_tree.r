# =============================================================================
# Title: 5i_split_tree.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to create input files for iqtree within each broad haplogroup (phylogenetic analyses)
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. load libraries & install if necessary
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

if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

# 2. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/03_iqtree")

# 3. set paths
a2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_A2.fasta")
b2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_B2.fasta")
c1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_C1.fasta")
d1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_D1.fasta")
h1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_H1.fasta")
l0_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L0.fasta")
l1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L1.fasta")
l2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L2.fasta")
l3_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean/clean_L3.fasta")

# 4. read fastas as dnabin objects
read_in <- function(path) {
  fasta <- read.dna(path, format = "fasta")
  return(fasta)
}

seqs_a2 <- read_in(a2_path)
seqs_b2 <- read_in(b2_path)
seqs_c1 <- read_in(c1_path)
seqs_d1 <- read_in(d1_path)
seqs_h1 <- read_in(h1_path)
seqs_l0 <- read_in(l0_path)
seqs_l1 <- read_in(l1_path)
seqs_l2 <- read_in(l2_path)
seqs_l3 <- read_in(l3_path)

# 5. load metadata
meta <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/all_merged.meta", sep="\t")

# 6. add neanderthal root to the metadata so it's included in all subsets
ref_rows <- data.frame(
  file = "AM948965.fasta",
  paper = "REF",
  ancient_or_modern = NA,
  original_format = NA,
  country = NA,
  extra_info = "REF",
  mtdna_type = NA,
  Sample = "AM948965",
  HaplogrepSampleID = "REF",
  Haplogroup = "REF",
  Rank = NA,
  Quality = NA,
  hap_1char = "REF",
  hap_2char = "REF",
  hap_3char = "REF",
  hap_4char = "REF",
  hap_5char = "REF",
  hap_6char = "REF",
  hap_7char = "REF",
  hap_8char = "REF",
  length = NA,
  country_label = "REF"
)

meta <- rbind(meta, ref_rows)

# 7. create objects for each tree subset
subset_dna <- function(dna_obj, meta_df, filter_expr) {
  rows <- meta_df %>% filter(!!enquo(filter_expr))
  idx <- match(rows$Sample, rownames(dna_obj))
  if(any(is.na(idx))) {
    warning("sequences in metadata but missing in tree: ", paste(rows$Sample[is.na(idx)], collapse = ", "))
  }
  # Subset and return DNAbin (only existing sequences)
  dna_obj[idx[!is.na(idx)], ]
}

# A2: A2, A2+(64)
final_a2 <- subset_dna(
  seqs_a2,
  meta,
  (hap_2char == "A2" & !hap_3char %in% c("A20","A21","A22","A23","A24","A25","A26")) |
    hap_2char == "REF"
)

final_a264 <- subset_dna(seqs_a2, meta,(hap_7char == "A2+(64)" | hap_2char == "REF"))

# B2: B2, B2b
final_b2 <- subset_dna(seqs_b2, meta, hap_2char == "B2" | hap_2char == "REF")
final_b2b <- subset_dna(seqs_b2, meta, hap_3char == "B2b" | hap_2char == "REF" | Sample == "LQ14")

# C1: C1b
final_c1 <- subset_dna(seqs_c1, meta, hap_2char == "C1" | hap_2char == "REF")
final_c1b <- subset_dna(seqs_c1, meta, hap_3char == "C1b" | hap_2char == "REF")

# D1: D1
final_d1 <- subset_dna(seqs_d1, meta, hap_2char == "D1" | hap_2char == "REF")

# H1: H1, H1bw
final_h1 <- subset_dna(seqs_h1, meta, hap_3char %in% c("REF", "H1a", "H1b", "H1c", "H1e", "H1f", "H1g", "H1h", "H1i", "H1j", "H1k", "H1l", "H1m", "H1n", "H1o", "H1p", "H1q", "H1r", "H1s", "H1t", "H1u", "H1v", "H1w", "H1x", "H1y", "H1z"))
final_h1bw <- subset_dna(seqs_h1, meta, hap_4char == "H1bw" | hap_2char == "REF")

# L0: L0a2a, L0d1a
final_l0 <- subset_dna(seqs_l0, meta, hap_2char == "L0" | hap_2char == "REF")
final_l0a2a <- subset_dna(seqs_l0, meta, hap_5char == "L0a2a" | hap_2char == "REF")
final_l0d1a <- subset_dna(seqs_l0, meta, hap_5char == "L0d1a" | hap_2char == "REF")

# L1: L1c, L1c1, L1c2'4
final_l1 <- subset_dna(seqs_l1, meta, hap_2char == "L1" | hap_2char == "REF")
final_l1c1 <- subset_dna(seqs_l1, meta, hap_4char == "L1c1" | hap_2char == "REF")
final_l1c1a <- subset_dna(seqs_l1, meta, hap_5char == "L1c1a" | hap_2char == "REF")
final_l1c1b <- subset_dna(seqs_l1, meta, hap_5char == "L1c1b" | hap_2char == "REF")
final_l1c1d <- subset_dna(seqs_l1, meta, hap_5char == "L1c1d" | hap_2char == "REF")
final_l1c2 <- subset_dna(seqs_l1, meta, (hap_4char == "L1c2" & hap_5char != "L1c2'") | hap_2char == "REF")
final_l1c4 <- subset_dna(seqs_l1, meta, hap_4char == "L1c4" | hap_2char == "REF")
final_l1c24 <- subset_dna(seqs_l1, meta, (hap_4char == "L1c2" | hap_4char == "L1c4") | hap_2char == "REF")

# L2: L2a1, L2a1a, L2a1d, L2a1f, L2a1q, L2b1, L2c
final_l2 <- subset_dna(seqs_l2, meta, hap_2char == "L2" | hap_2char == "REF")
final_l2a1a <- subset_dna(seqs_l2, meta, hap_5char == "L2a1a" | hap_2char == "REF")
final_l2a1d <- subset_dna(seqs_l2, meta, hap_5char == "L2a1d" | hap_2char == "REF")
final_l2a1f <- subset_dna(seqs_l2, meta, hap_5char == "L2a1f" | hap_2char == "REF")
final_l2a1q <- subset_dna(seqs_l2, meta, hap_5char == "L2a1q" | hap_2char == "REF")
final_l2b1 <- subset_dna(seqs_l2, meta, hap_4char == "L2b1" | hap_2char == "REF")
final_l2c <- subset_dna(seqs_l2, meta, hap_3char == "L2c" | hap_2char == "REF")

# L3: L3b1a, L3d1, L3d3, L3d4, L3d3 & L3d4, L3e1a1a, L3e2, L3e3, L3f1b
final_l3 <- subset_dna(seqs_l3, meta, hap_2char == "L3" | hap_2char == "REF")
final_l3b1a <- subset_dna(seqs_l3, meta, hap_5char == "L3b1a" | hap_2char == "REF")
final_l3d1 <- subset_dna(seqs_l3, meta, hap_4char == "L3d1" | hap_2char == "REF")
final_l3d3 <- subset_dna(seqs_l3, meta, hap_4char == "L3d3" | hap_2char == "REF")
final_l3d4 <- subset_dna(seqs_l3, meta, hap_4char == "L3d4" | hap_2char == "REF")
final_l3e1a1a <- subset_dna(seqs_l3, meta, hap_7char == "L3e1a1a" | hap_2char == "REF")
final_l3e1 <- subset_dna(seqs_l3, meta, hap_4char == "L3e1" | hap_2char == "REF")
final_l3e2 <- subset_dna(seqs_l3, meta, hap_4char == "L3e2" | hap_2char == "REF")
final_l3e3b <- subset_dna(seqs_l3, meta, hap_5char == "L3e3b" | hap_2char == "REF")
final_l3f1b <- subset_dna(seqs_l3, meta, hap_5char == "L3f1b" | hap_2char == "REF")
final_l3f1 <- subset_dna(seqs_l3, meta, hap_4char == "L3f1" | hap_2char == "REF")

# 5. write fasta files for iqtree
# A2
write.dna(final_a2, file = "iq_a2.fasta", format = "fasta")
write.dna(final_a264, file = "iq_a264.fasta", format = "fasta")

# B2
write.dna(final_b2, file = "iq_b2.fasta", format = "fasta")
write.dna(final_b2b, file = "iq_b2b.fasta", format = "fasta")

# C1
write.dna(final_c1, file = "iq_c1.fasta", format = "fasta")
write.dna(final_c1b, file = "iq_c1b.fasta", format = "fasta")

# D1
write.dna(final_d1, file = "iq_d1.fasta", format = "fasta")

# H1
write.dna(final_h1, file = "iq_h1.fasta", format = "fasta")
write.dna(final_h1bw, file = "iq_h1bw.fasta", format = "fasta")

# L0
write.dna(final_l0, file = "iq_l0.fasta", format = "fasta")
write.dna(final_l0a2a, file = "iq_l0a2a.fasta", format = "fasta")
write.dna(final_l0d1a, file = "iq_l0d1a.fasta", format = "fasta")

# L1
write.dna(final_l1, file = "iq_l1.fasta", format = "fasta")
write.dna(final_l1c1, file = "iq_l1c1.fasta", format = "fasta")
write.dna(final_l1c1a, file = "iq_l1c1a.fasta", format = "fasta")
write.dna(final_l1c1b, file = "iq_l1c1b.fasta", format = "fasta")
write.dna(final_l1c1d, file = "iq_l1c1d.fasta", format = "fasta")
write.dna(final_l1c2, file = "iq_l1c2.fasta", format = "fasta")
write.dna(final_l1c4, file = "iq_l1c4.fasta", format = "fasta")
write.dna(final_l1c24, file = "iq_l1c24.fasta", format = "fasta")

# L2
write.dna(final_l2, file = "iq_l2.fasta", format = "fasta")
write.dna(final_l2a1a, file = "iq_l2a1a.fasta", format = "fasta")
write.dna(final_l2a1d, file = "iq_l2a1d.fasta", format = "fasta")
write.dna(final_l2a1f, file = "iq_l2a1f.fasta", format = "fasta")
write.dna(final_l2a1q, file = "iq_l2a1q.fasta", format = "fasta")
write.dna(final_l2b1, file = "iq_l2b1.fasta", format = "fasta")
write.dna(final_l2c, file = "iq_l2c.fasta", format = "fasta")

# L3
write.dna(final_l3, file = "iq_l3.fasta", format = "fasta")
write.dna(final_l3b1a, file = "iq_l3b1a.fasta", format = "fasta")
write.dna(final_l3d1, file = "iq_l3d1.fasta", format = "fasta")
write.dna(final_l3d3, file = "iq_l3d3.fasta", format = "fasta")
write.dna(final_l3d4, file = "iq_l3d4.fasta", format = "fasta")
write.dna(final_l3e1a1a, file = "iq_l3e1a1a.fasta", format = "fasta")
write.dna(final_l3e1, file = "iq_l3e1.fasta", format = "fasta")
write.dna(final_l3e2, file = "iq_l3e2.fasta", format = "fasta")
write.dna(final_l3e3b, file = "iq_l3e3.fasta", format = "fasta")
write.dna(final_l3f1b, file = "iq_l3f1b.fasta", format = "fasta")
write.dna(final_l3f1, file = "iq_l3f1.fasta", format = "fasta")

# 8. create metadata for iqtree samples
tree_samples <- unique(c(
  rownames(final_a2),
  rownames(final_b2),
  rownames(final_c1),
  rownames(final_d1),
  rownames(final_h1),
  rownames(final_l0),
  rownames(final_l1),
  rownames(final_l2),
  rownames(final_l3),
  rownames(final_a264),
  rownames(final_b2b),
  rownames(final_c1b),
  rownames(final_h1bw),
  rownames(final_l0a2a),
  rownames(final_l0d1a),
  rownames(final_l1c1a),
  rownames(final_l1c1b),
  rownames(final_l1c1d),
  rownames(final_l1c2),
  rownames(final_l1c24),
  rownames(final_l1c4),
  rownames(final_l2a1a),
  rownames(final_l2a1d),
  rownames(final_l2a1f),
  rownames(final_l2a1q),
  rownames(final_l2b1),
  rownames(final_l2c),
  rownames(final_l3b1a),
  rownames(final_l3d1),
  rownames(final_l3d3),
  rownames(final_l3d4),
  rownames(final_l3e1),
  rownames(final_l3e1a1a),
  rownames(final_l3e2),
  rownames(final_l3e3b),
  rownames(final_l3f1b),
  rownames(final_l3f1)
))

# 9. subset metadata to only those samples
meta_iqtree <- meta %>%
  filter(Sample %in% tree_samples)

# 10. write metadata
write.table(
  meta_iqtree,
  file = "iqtree.meta",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)