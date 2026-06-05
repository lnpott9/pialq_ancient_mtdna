#!/bin/bash -l

# =============================================================================
# Title: 5g_mafft.sh
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: SLURM job array to align sequences in broad haplogroups before network/phylogenetic analyses: A2, B2, C1, D1, L0, L1, L2, L3
# Usage: sbatch 5g_mafft.sh
# =============================================================================

#SBATCH --job-name=auto
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80gb
#SBATCH --tmp=80gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-8

module load mafft/7.475

BASE_DIR="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa"
SORTED="$BASE_DIR/01_mafft"
REF="$SHARED/ref_seqs"

# 1. folder listing
FOLDERS=(
    "A2"
    "B2"
    "C1"
    "D1"
    "H1"
    "L0"
    "L1"
    "L2"
    "L3"
)

# 2. choose folder for this task
FOLDER="${FOLDERS[$SLURM_ARRAY_TASK_ID]}"
FULL_PATH="$SORTED/$FOLDER"

echo "processing folder $FOLDER at $FULL_PATH"

# 3. build input file from fasta files in input folder, adding rCRS as well
cat "$REF/rCRS.fasta" $FULL_PATH/*fasta | awk 'NF' > "$SORTED/${FOLDER}_mafft_in_auto.fasta"
echo "rcrs added"

# 4. run mafft with auto-selection of algorithm
if mafft --auto \
    --thread 10 \
    "$SORTED/${FOLDER}_mafft_in_auto.fasta" \
    > "$SORTED/${FOLDER}_mafft_out_auto.fasta" ; then
    echo "alignment succeeded"
else
	echo "alignment failed"
fi