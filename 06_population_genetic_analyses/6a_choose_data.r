# =============================================================================
# Title: 6a_choose_data.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to choose papers/datasets to include in population genetics analyses
# Usage: run interactively
# =============================================================================

# 1. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity")

# 2. load libraries, download if necessary
library(seqinr)
library(pegas)
library(dplyr)
library(ape)

# 3. load metadata & identify populations of interest
meta <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/all_merged.meta", sep="\t")
dim(meta)
# 18540    21

#### a. remove low quality ones, make this match tree stringency
meta <- meta[(meta$Quality > 0.8 & meta$length > 14000) | meta$paper == "pialq", ]

#### b. unique populations & papers
# look at unique combinations
c_unique <- meta %>% distinct(country, .keep_all = TRUE)
dim(c_unique)

# paper-country combo
pc_unique <- meta %>% distinct(paper, country, .keep_all = TRUE)
dim(pc_unique)
# 350  21

# count distinct options
pair_summary <- meta %>%
  count(paper, country)

big <- pair_summary[pair_summary$n > 60, ]

#### c. pick out populations to keep
keep <- meta$paper %in% c("pialq", "harney_2023", "sandovalvelasco_2023", "fleskes_2023", "1kg", "hgdp") | meta$paper %in% big$paper

meta_filtered <- meta[keep, ]
dim(meta_filtered)

# 4. write df of sample sizes
meta_filt <- meta_filtered %>%
  count(paper, country)

unique_papers <- unique(meta_filt$paper)

write.table(meta_filt, "n_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)