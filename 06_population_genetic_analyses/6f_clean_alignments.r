# =============================================================================
# Title: 6f_clean_alignments.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to clean & standardize sequences after alignment
# Usage: run interactively
# =============================================================================

# 1. load libraries
for(pkg in c("pegas", "ape", "Biostrings")) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# 2. set directories
seq_dir   <- "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/01_mafft"
clean_dir <- "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/02_clean"

# 3. list input FASTAs
path_list <- list.files(seq_dir, pattern = "_mafft_out_auto\\.fasta$", full.names = TRUE)
names(path_list) <- sub("_mafft_out_auto\\.fasta$", "", basename(path_list))

# 4. read sequences as DNAbin and convert to uppercase character matrices
aligned_seqs_char <- lapply(path_list, function(f) toupper(as.character(read.dna(f, format = "fasta"))))

# 5. extract rCRS sequence from each alignment
rCRS_seq <- lapply(aligned_seqs_char, function(mat) mat[grep("NC_012920.1", rownames(mat)), ])

# 6. mask problem positions in the original alignment
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

# 7. remove columns where rCRS has gaps
aligned_seqs_clean <- mapply(function(mat, rcs) {
  keep_cols <- rcs != "-"
  mat[, keep_cols, drop = FALSE]
}, aligned_seqs_masked, rCRS_seq, SIMPLIFY = FALSE)

# 8. convert each cleaned alignment to DNAStringSet and remove rCRS
clean_alignments_list <- lapply(aligned_seqs_clean, function(mat) {
  dna_seq <- DNAStringSet(apply(mat, 1, paste, collapse = ""))
  names(dna_seq) <- rownames(mat)
  dna_seq[!grepl("NC_012920.1", names(dna_seq))]
})

# 9. concatenate all DNAStringSets into one
all_seqs <- unlist(lapply(clean_alignments_list, as.character))
combined_alignment <- DNAStringSet(all_seqs)

# 10. Write single combined FASTA
writeXStringSet(combined_alignment, filepath = file.path(clean_dir, "clean_div.fasta"), width = 18000)