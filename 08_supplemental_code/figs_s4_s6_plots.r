# =============================================================================
# Title: figs_s4_s6_plots.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to create a depth histogram & plot of haplogroup classification vs depth
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. load libraries
library(ggplot2)

# 2. setwd
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/02_circmap/seq_stats")

# 3. read in text file
stats <- read.table("seq_stats.txt", 
                    header = TRUE,
                    sep = "\t", # tab
                    quote = "", # Ignore quotes
                    fill = TRUE, # Fill missing columns
                    stringsAsFactors = FALSE)

# 4. histogram
hist <- ggplot(stats, aes(x=cov_mean)) + 
  geom_histogram(binwidth = 20, color = "black", fill = "green") +
  scale_x_continuous(breaks = seq(0,350, by=20)) +
  xlab("Mean Depth (X)") +
  ylab("Library Count")
  
# 5. graph of haplogrep quality vs coverage
scatter <- ggplot(stats, aes(x=cov_mean, y=Quality_5x)) +
  geom_point() +
  xlab("Mean Depth (X)") +
  ylab("Haplogrep Quality (based on 5X VCF)")

# 6. print to pdf
pdf("fig_s4_hist.pdf")
hist
dev.off()

pdf("fig_s6_scatter.pdf")
scatter
dev.off()

