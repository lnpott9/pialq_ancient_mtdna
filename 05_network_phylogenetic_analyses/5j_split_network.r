# =============================================================================
# Title: 5j_split_network.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to create input files for popart within each broad haplogroup (network analyses)
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. load libraries
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

# 2. set wd
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn")

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

# 4. load metadata
meta <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/final_cut.meta", sep="\t")

# 5. create objects for each tree subset
subset_dna <- function(dna_obj, meta_df, filter_expr) {
  rows <- meta_df %>% filter(!!enquo(filter_expr))
  idx <- match(rows$Sample, rownames(dna_obj))
  if(any(is.na(idx))) {
    warning("sequences in metadata but missing in tree: ", paste(rows$Sample[is.na(idx)], collapse = ", "))
  }
  dna_obj[idx[!is.na(idx)], ]
}

# A2: A2+(64)
final_a264 <- subset_dna(seqs_a2, meta, hap_7char == "A2+(64)")

# B2: B2, B2b
final_b2b <- subset_dna(seqs_b2, meta, hap_3char == "B2b")

# C1: C1b
final_c1b <- subset_dna(seqs_c1, meta, hap_3char == "C1b")

# D1: D1
final_d1 <- subset_dna(seqs_d1, meta, hap_2char == "D1")

# H1: H1, H1bw
final_h1bw <- subset_dna(seqs_h1, meta, hap_4char == "H1bw")

# L0: L0a2a, L0d1a
final_l0a2a <- subset_dna(seqs_l0, meta, hap_5char == "L0a2a")
final_l0d1a <- subset_dna(seqs_l0, meta, hap_5char == "L0d1a")

# L1: L1c, L1c1, L1c2'4
final_l1c <- subset_dna(seqs_l1, meta, hap_3char == "L1c")
final_l1c1 <- subset_dna(seqs_l1, meta, hap_4char == "L1c1")
final_l1c24 <- subset_dna(seqs_l1, meta, hap_4char %in% c("L1c2", "L1c4"))

# L2: L2a1, L2a1a, L2a1d, L2a1f, L2a1q, L2b1, L2c
final_l2a1a   <- subset_dna(seqs_l2, meta, hap_5char == "L2a1a")
final_l2a1d   <- subset_dna(seqs_l2, meta, hap_5char == "L2a1d")
final_l2a1f   <- subset_dna(seqs_l2, meta, hap_5char == "L2a1f")
final_l2a1q   <- subset_dna(seqs_l2, meta, hap_5char == "L2a1q")
final_l2b1    <- subset_dna(seqs_l2, meta, hap_4char == "L2b1")
final_l2c     <- subset_dna(seqs_l2, meta, hap_3char == "L2c")

# L3: L3b1a, L3d1, L3d3, L3d4, L3d3 & L3d4, L3e1a1a, L3e2, L3e3, L3f1b
final_l3b1a <- subset_dna(seqs_l3, meta, hap_5char == "L3b1a")
final_l3d1 <- subset_dna(seqs_l3, meta, hap_4char == "L3d1")
final_l3d3 <- subset_dna(seqs_l3, meta, hap_4char == "L3d3")
final_l3d4 <- subset_dna(seqs_l3, meta, hap_4char == "L3d4")
final_l3e1a1a <- subset_dna(seqs_l3, meta, hap_7char == "L3e1a1a")
final_l3e2 <- subset_dna(seqs_l3, meta, hap_4char == "L3e2")
final_l3e3 <- subset_dna(seqs_l3, meta, hap_4char == "L3e3")
final_l3f1b <- subset_dna(seqs_l3, meta, hap_5char == "L3f1b")

# 5. write fasta files for iqtree
# A2
write.dna(final_a264, file = "mjn_a2_64.fasta", format = "fasta")

# B2
write.dna(final_b2b, file = "mjn_b2b.fasta", format = "fasta")

# C1
write.dna(final_c1b, file = "mjn_c1b.fasta", format = "fasta")

# D1
write.dna(final_d1, file = "mjn_d1.fasta", format = "fasta")

# H1
write.dna(final_h1bw, file = "mjn_h1bw.fasta", format = "fasta")

# L0
write.dna(final_l0a2a, file = "mjn_l0a2a.fasta", format = "fasta")
write.dna(final_l0d1a, file = "mjn_l0d1a.fasta", format = "fasta")

# L1
write.dna(final_l1c1, file = "mjn_l1c1.fasta", format = "fasta")
write.dna(final_l1c24, file = "mjn_l1c24.fasta", format = "fasta")

# L2
write.dna(final_l2a1a, file = "mjn_l2a1a.fasta", format = "fasta")
write.dna(final_l2a1d, file = "mjn_l2a1d.fasta", format = "fasta")
write.dna(final_l2a1f, file = "mjn_l2a1f.fasta", format = "fasta")
write.dna(final_l2a1q, file = "mjn_l2a1q.fasta", format = "fasta")
write.dna(final_l2b1, file = "mjn_l2b1.fasta", format = "fasta")
write.dna(final_l2c, file = "mjn_l2c.fasta", format = "fasta")

# L3
write.dna(final_l3b1a, file = "mjn_l3b1a.fasta", format = "fasta")
write.dna(final_l3d1, file = "mjn_l3d1.fasta", format = "fasta")
write.dna(final_l3d3, file = "mjn_l3d3.fasta", format = "fasta")
write.dna(final_l3d4, file = "mjn_l3d4.fasta", format = "fasta")
write.dna(final_l3e1a1a, file = "mjn_l3e1a1a.fasta", format = "fasta")
write.dna(final_l3e2, file = "mjn_l3e2.fasta", format = "fasta")
write.dna(final_l3e3, file = "mjn_l3e3.fasta", format = "fasta")
write.dna(final_l3f1b, file = "mjn_l3f1b.fasta", format = "fasta")

# 6. create metadata for iqtree samples
tree_samples <- unique(c(
  rownames(final_a264),
  rownames(final_b2b),
  rownames(final_c1b),
  rownames(final_d1),
  rownames(final_h1bw),
  rownames(final_l0a2a),
  rownames(final_l0d1a),
  rownames(final_l1c1),
  rownames(final_l1c24),
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
  rownames(final_l3e1a1a),
  rownames(final_l3e2),
  rownames(final_l3e3),
  rownames(final_l3f1b)
))

# 7. subset metadata to only those samples
meta_iqtree <- meta %>%
  filter(Sample %in% tree_samples)

# 8. write metadata
write.table(
  meta_iqtree,
  file = "mjn.meta",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)