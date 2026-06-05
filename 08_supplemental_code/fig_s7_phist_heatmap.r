# =============================================================================
# Title: fig_s7_phist_heatmap.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to create heatmap showing PhiST valuesand statistical significance
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. load libraries
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

# 2. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/04_inter")

# 3. load phist matrix & p value smatrix
phist_matrix <- as.matrix(read.table("all_phist_matrix.txt", header = TRUE, row.names = 1))
pval_matrix <- as.matrix(read.table("all_phist_pvalues.txt", header = TRUE, row.names = 1))

# 4. visualize heatmap of phist values
phist <- pheatmap(phist_matrix, 
                  cluster_rows = FALSE, 
                  cluster_cols = FALSE, 
                  border_color = "darkgray", 
                  na_col = "white",
                  color = viridis(100, option = "viridis", direction = -1),
                  main = "Pairwise PhiST values"
)

# 5. visualize heatmap of p values
pval <- pheatmap(pval_matrix, 
                 cluster_rows = FALSE, 
                 cluster_cols = FALSE,
                 border_color = "darkgray", 
                 na_col = "white",
                 color = viridis(100, option = "viridis", direction = -1),
                 main = "P-values of pairwise PhiST comparisons"
)

# 6. visualize heatmap of significant phist values
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

# 7. save final plots
pdf("all_phist_heatmap.pdf", width = 12, height = 12)
print(phist)
grid.newpage()
print(pval)
grid.newpage()
print(pval_bin) 
dev.off()

