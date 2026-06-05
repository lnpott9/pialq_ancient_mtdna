#!/bin/bash

# =============================================================================
# Title: 6h_interpop_wrapper.sh
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: Wrapper script to run 6h_interpop.r
# Usage: sbatch 6h_interpop_wrapper.sh
# =============================================================================

#SBATCH --job-name=interpop
#SBATCH --time=70:00:00
#SBATCH --mem=80gb
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

# 1. load conda
module load conda

# 2. activate environment
source activate $SHARED/projects/PIALQ/2025_ancient_lp/12_diversity/st_env

# 3. go to wd
cd /projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/04_inter

# 4. run r script
Rscript --vanilla pval.r
