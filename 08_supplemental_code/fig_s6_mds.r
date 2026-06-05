# =============================================================================
# Title: fig_s6_mds.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to create MDS with every point labeled
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/04_inter")

# 2. load libraries, download if necessary
library(dplyr)
library(ggplot2)

# 3. plot mds using coordinates already calculated
# read in coordinates
mds_data <- read.table("all_mds_coordinates.txt", header = TRUE)

# 4. avoid overlaps
library(ggrepel)

final_pop <- ggplot(mds_data, aes(x = MDS1, y = MDS2, label = Population)) +
  geom_point(size = 1, color = "black") +
  geom_text_repel(size = 1, max.overlaps=Inf) +
  theme_bw() +
  labs(x = "MDS1",
       y = "MDS2")

# save mds plot
pdf("all_mds.pdf", width = 11, height = 8)
final_pop
dev.off()
