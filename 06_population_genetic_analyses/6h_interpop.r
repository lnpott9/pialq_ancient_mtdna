# =============================================================================
# Title: 6h_interpop.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to run calculation of pairwise PhiST and calculate MDS with permutation testing
# Usage: n/a, see wrapper script 6h_interpop_wrapper.sh
# =============================================================================

# 1. load libraries
library(ape)
library(haplotypes)

# 2. logging helper
log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
  flush.console()
}

# 3. read data: alignment & metadata
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
  stop("ERROR: pop_assign.txt not found")
}

pop_assignments <- tryCatch({
  read.delim("pop_assign.txt", stringsAsFactors = FALSE, header = FALSE)
}, error = function(e) {
  stop(paste("ERROR reading pop_assign.txt:", e$message))
})

names(pop_assignments) <- c("sample", "population")
log_msg(paste("Loaded", nrow(pop_assignments), "population assignments"))

# 5. create lookup: sample name -> population
sample_to_pop <- setNames(pop_assignments$population, pop_assignments$sample)

# 6. match populations to alignment samples
log_msg("matching samples to pops")
align_sample_names <- rownames(align)
align_to_pop <- sample_to_pop[align_sample_names]

# 7. keep only samples with population assignments
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

# 8. check population sizes
pop_sizes <- table(grouping)
log_msg("\nPopulation sizes:")
for (pop_name in names(pop_sizes)) {
  log_msg(paste("  ", pop_name, ":", pop_sizes[pop_name], "samples"))
}

# 9. clean alignment to handle ambiguity codes
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

# 10. write cleaned alignment to temporary file
temp_fasta <- tempfile(fileext = ".fasta")
log_msg(paste("writing cleaned alignment to temporary file"))

fasta_lines <- character()
for (i in 1:nrow(align_char_clean)) {
  fasta_lines <- c(fasta_lines, 
                   paste0(">", rownames(align_char_clean)[i]),
                   paste(align_char_clean[i, ], collapse = ""))
}
writeLines(fasta_lines, temp_fasta)

# 11. read cleaned fasta with haplotypes package
dna_hap <- tryCatch({
  read.fas(temp_fasta)
}, error = function(e) {
  unlink(temp_fasta)
  stop(paste("error reading fasta with haplotypes package:", e$message))
})

log_msg(paste("Loaded", nrow(dna_hap), "sequences with", ncol(dna_hap), "sites"))

# 12. verify sample order matches
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

# 13. extract and validate results
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

# 14. validate distance matrix exists
if (is.null(phist_dist)) {
  stop("ERROR: Failed to extract PhiST distance matrix")
}

# 15. check for NA values
if (any(is.na(phist_dist))) {
  log_msg("WARNING: PhiST matrix contains NA values")
  na_count <- sum(is.na(phist_dist))
  log_msg(paste("  Number of NA values:", na_count))
}

log_msg("PhiST distance matrix validated successfully")

# 16. save results
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

# 17. MDS calculation
log_msg("\n=== performing metric MDS ===")

# ensure matrix is symmetric and has no negative values
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

# get GOF
mds_gof <- mds_result$GOF

# 18. calculate variance explained
eig_vals <- mds_result$eig[mds_result$eig > 0]
if (length(eig_vals) == 0) {
  log_msg("no positive eigenvalues")
  var_exp <- c(0, 0)
} else {
  var_exp <- 100 * eig_vals / sum(eig_vals)
}

log_msg(paste("MDS Axis 1 explains", round(var_exp[1], 2), "% of variance"))
if (length(var_exp) > 1) {
  log_msg(paste("MDS Axis 2 explains", round(var_exp[2], 2), "% of variance"))
  log_msg(paste("Total variance explained:", round(sum(var_exp[1:2]), 2), "%"))
  log_msg(paste("Goodness of fit:", mds_gof)))
}

# 19. save MDS coordinates
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

# 20. summary statistics
log_msg("\n=== Summary Statistics ===")
log_msg(paste("Total samples analyzed:", nrow(dna_hap)))
log_msg(paste("Number of populations:", nlevels(grouping)))
log_msg(paste("Informative sites used:", ncol(dna_hap)))
log_msg(paste("Number of pairwise comparisons:", 
              nlevels(grouping) * (nlevels(grouping) - 1) / 2))

# 21. final output summary
log_msg("\nAnalysis finished successfully!")

