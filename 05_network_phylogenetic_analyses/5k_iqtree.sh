#!/bin/bash -l
# =============================================================================
# Title: 5k_iqtree.sh
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: SLURM job array to run IQTREE
# Usage: sbatch 5k_iqtree.sh
# =============================================================================

#SBATCH --job-name=10k
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --tmp=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-27

module load conda

cd $SHARED/projects/PIALQ/2025_ancient_lp/10_iqtree/

source activate $SHARED/projects/PIALQ/2025_ancient_lp/10_iqtree/tree_env

# 1. set names to files
TREES=(
"iq_a2.fasta"
"iq_a264.fasta"
"iq_b2.fasta"
"iq_b2b.fasta"
"iq_c1.fasta"
"iq_c1b.fasta"
"iq_d1.fasta"
"iq_h1.fasta"
"iq_h1bw.fasta"
"iq_l0a2a.fasta"
"iq_l0d1a.fasta"
"iq_l1c1.fasta"
"iq_l1c24.fasta"
"iq_l2a1a.fasta"
"iq_l2a1d.fasta"
"iq_l2a1f.fasta"
"iq_l2a1q.fasta"
"iq_l2b1.fasta"
"iq_l2c.fasta"
"iq_l3b1a.fasta"
"iq_l3d1.fasta"
"iq_l3d3.fasta"
"iq_l3d4.fasta"
"iq_l3e1.fasta"
"iq_l3e1a1a.fasta"
"iq_l3e2.fasta"
"iq_l3e3.fasta"
"iq_l3e3b.fasta"
"iq_l3f1b.fasta"
)

base="/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/03_iqtree"

# 2. select file for this task
FILE="${TREES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing $FILE"

# 3. infer rooted trees with outgroup of neanderthal mtDNA
PREFIX="${FILE%.fasta}"
iqtree -s "$base/$FILE" -m MFP -bb 10000 -alrt 1000 -T 1 -safe --prefix "$PREFIX" -o AM948965