#!/bin/bash -l

# =============================================================================
# Title: 6e_mafft.sh
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: SLURM job array to align sequences in broad haplogroups before popgen analyses
# Usage: sbatch 6e_mafft.sh
# =============================================================================

#SBATCH --job-name=01_m
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=80gb
#SBATCH --tmp=80gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-126

module load mafft/7.475

BASE_DIR="$SHARED/projects/PIALQ/2025_ancient_lp/12_diversity"
SORTED="$BASE_DIR/01_mafft"
REF="$SHARED/ref_seqs"

# 1. folder listing
FOLDERS=(
"A+"
"A1"
"A2"
"A3"
"A5"
"A6"
"A8"
"B2"
"B4"
"B5"
"B6"
"C"
"C1"
"C4"
"C5"
"C7"
"D1"
"D2"
"D3"
"D4"
"D5"
"D6"
"E1"
"E2"
"F1"
"F2"
"F3"
"F4"
"G1"
"G2"
"G3"
"H"
"H+"
"H1"
"H2"
"H3"
"H4"
"H5"
"H6"
"H7"
"H8"
"H9"
"HV"
"I"
"I1"
"I2"
"I3"
"I4"
"I5"
"I6"
"J1"
"J2"
"K1"
"K2"
"K3"
"L0"
"L1"
"L2"
"L3"
"L4"
"L5"
"M"
"M1"
"M2"
"M3"
"M4"
"M5"
"M6"
"M7"
"M8"
"M9"
"N1"
"N2"
"N5"
"N8"
"N9"
"P1"
"P9"
"Q1"
"R"
"R+"
"R0"
"R1"
"R2"
"R3"
"R5"
"R6"
"R7"
"R8"
"R9"
"T1"
"T2"
"U"
"U1"
"U2"
"U3"
"U4"
"U5"
"U6"
"U7"
"U8"
"U9"
"V"
"V+"
"V1"
"V2"
"V3"
"V7"
"V9"
"W"
"W+"
"W1"
"W3"
"W4"
"W5"
"W6"
"W7"
"W8"
"W9"
"X1"
"X2"
"X3"
"Y1"
"Y2"
"Z2"
"Z3"
"Z4"
)

# 2. pick folder for this task
FOLDER="${FOLDERS[$SLURM_ARRAY_TASK_ID]}"
FULL_PATH="$SORTED/$FOLDER"

echo "processing folder $FOLDER at $FULL_PATH"

# 3. build input file
cat "$REF/rCRS.fasta" $FULL_PATH/*fasta | awk 'NF' > "$SORTED/${FOLDER}_mafft_in_auto.fasta"
echo "rcrs added"

# 4. real mafft with auto
if mafft --auto \
    --thread 10 \
    "$SORTED/${FOLDER}_mafft_in_auto.fasta" \
    > "$SORTED/${FOLDER}_mafft_out_auto.fasta" ; then
    echo "alignment succeeded"
else
	echo "alignment failed"
fi