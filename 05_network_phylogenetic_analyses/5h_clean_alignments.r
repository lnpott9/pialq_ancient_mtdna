# =============================================================================
# Title: 5h_clean_alignments.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to clean & standardize sequences after alignment
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

# 2. set wd
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/02_clean")

# 3. set paths to MSA output files
a2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/A2_mafft_out.fasta")
b2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/B2_mafft_out.fasta")
c1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/C1_mafft_out.fasta")
d1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/D1_mafft_out.fasta")
h1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/H1_mafft_out.fasta")
l0_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/L0_mafft_out.fasta")
l1_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/L1_mafft_out.fasta")
l2_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/L2_mafft_out.fasta")
l3_path <- ("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft/L3_mafft_out.fasta")

# 4. read paths into a list
path_list <- list(a2_path, b2_path, c1_path, d1_path, h1_path, l0_path, l1_path, l2_path, l3_path)

# 5. read sequences into DNAbin objects, convert all to upper case characters
aligned_seqs <- lapply(path_list, read.dna, format = "fasta")

aligned_seqs_char <- lapply(aligned_seqs, function(x) toupper(as.character(x)))

# length/width should be pretty close to 16569 (length of rCRS)
lapply (aligned_seqs_char, function (x) dim(x))

# 6. identify rCRS sequence in each alignment
rCRS_seq <- lapply(aligned_seqs_char, function(x) {
  rcs_name <- grep("NC_012920.1", rownames(x), value = TRUE)
  x[rcs_name, ]
})

# 7. clean alignments: remove rCRS gaps, mask problem positions, replace ambiguous codes
problem_positions <- c(303:315, 515:522, 568:573, 3107, 16182:16194, 16519)
haplogroups <- c("A2","B2","C1","D1","H1","L0","L1","L2","L3")

aligned_seqs_no_gaps <- lapply(seq_along(aligned_seqs_char), function(i) {
  alignment_mat <- aligned_seqs_char[[i]]
  rcs <- rCRS_seq[[i]]
  
  # mask problem positions
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

names(aligned_seqs_no_gaps) <- haplogroups

# 8. Calculate missingness after cleaning
count_missing <- function(seq_matrix) {
  apply(seq_matrix, 1, function(row) sum(row == "N"))
}

missing_counts_list <- lapply(aligned_seqs_no_gaps, count_missing)

missing_df_list <- lapply(seq_along(aligned_seqs_no_gaps), function(i) {
  data.frame(
    Haplogroup = names(aligned_seqs_no_gaps)[i],
    Sample = rownames(aligned_seqs_no_gaps[[i]]),
    Missing_Count = missing_counts_list[[i]],
    Missing_Fraction = missing_counts_list[[i]] / ncol(aligned_seqs_no_gaps[[i]])
  )
})
missing_df <- do.call(rbind, missing_df_list)

# 9. examine missing_df manually --> did not remove any based on missingness

# 10. convert back to DNAStringSet
clean_alignments_list <- lapply(aligned_seqs_no_gaps, function(mat) {
  # convert each row to a DNAString
  dna_seq <- apply(mat, 1, function(seq_row) {
    DNAString(paste(seq_row, collapse = ""))
  })
  
  # assign rownames as sequence names
  names(dna_seq) <- rownames(mat)
  
  # wrap as DNAStringSet
  DNAStringSet(dna_seq)
})

# 11. Remove the rCRS reference (NC_012920.1) from each
clean_alignments_list <- lapply(clean_alignments_list, function(dna_set) {
  names_to_drop <- grep("NC_012920.1", names(dna_set), value = TRUE)
  dna_set[!names(dna_set) %in% names_to_drop]
})

# 12. write cleaned alignments
haplogroups <- names(aligned_seqs_no_gaps)
for (i in seq_along(clean_alignments_list)) {
  writeXStringSet(
    clean_alignments_list[[i]],
    filepath = paste0("clean_", haplogroups[i], ".fasta"),
    width = 18000
  )
}
