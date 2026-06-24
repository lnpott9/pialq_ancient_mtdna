# =============================================================================
# Title: figs_s4_s6_plots.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to create a depth histogram & plot of haplogroup classification vs depth
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. load libraries
library(ggplot2)
library(dplyr)

# 2. setwd
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/02_circmap/seq_stats")

# 3. read in text file
stats <- read.table("table_s8.txt", 
                    header = TRUE,
                    sep = "\t", # tab
                    quote = "", # Ignore quotes
                    fill = TRUE, # Fill missing columns
                    stringsAsFactors = FALSE)

# 4. histogram
stats$Mean.cov..X. <- as.numeric(gsub(",", "", stats$Mean.cov..X.))

hist <- ggplot(stats, aes(x=Mean.cov..X.)) + 
  geom_histogram(binwidth = 100, color = "black", fill = "green") +
  scale_x_continuous(breaks = seq(0,3800, by=100)) +
  xlab("Mean Depth (X)") +
  ylab("Library Count")

# 5. graph of haplogrep quality vs coverage
graph <- read.table("table_s10_5x.txt", 
                    header = TRUE,
                    sep = "\t", # tab
                    quote = "", # Ignore quotes
                    fill = TRUE, # Fill missing columns
                    stringsAsFactors = FALSE)

stats$Sample.Name <- gsub("\\.mt$", "", stats$Sample.Name)

join <- left_join(stats, graph, by = c("Sample.Name" = "Library"))

scatter <- ggplot(join, aes(x=Mean.cov..X., y=Quality)) +
  geom_point() +
  xlab("Mean Depth (X)") +
  ylab("Haplogrep Quality (based on 5X VCF)")

# 6. print to svg
svg("fig_s4_hist.svg")
hist
dev.off()

svg("fig_s6_scatter.svg")
scatter
dev.off()

