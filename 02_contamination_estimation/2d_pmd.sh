# =============================================================================
# Title: 2d_pmd.sh
# Author: Laura N. Pott
# Date: 2026-06-03
# Description: script to separate ancient DNA molecules from others based on PMD-scores
# Usage: run interactively
# =============================================================================

# 1. start interactive session
srun --mem=8gb --time=60:00 --pty bash

# 2. load modules
module load python2/2.7.12_anaconda4.2
module load samtools/1.21

# 3. set paths
working=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/04_contam/02_pmd
pmdtools=/home/mnievesc/shared/programs/PMDtools/pmdtools.0.60.py

# 4. make array of intput files
files=("78mt_sorted.md.bam" "56mt_sorted.md.bam")

# decontaminate at different thresholds
cd "$working"

for bam in "${files[@]}"; do
	filename=$(basename "$bam" _sorted.md.bam)
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 2 --header | samtools view -Sb - > "${filename}_pmd2_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 2.5 --header | samtools view -Sb - > "${filename}_pmd2.5_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 3 --header | samtools view -Sb - > "${filename}_pmd3_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 3.5 --header | samtools view -Sb - > "${filename}_pmd3.5_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 4 --header | samtools view -Sb - > "${filename}_pmd4_sorted.md.bam"
	samtools view -h "$bam" | python2 "$pmdtools" --threshold 4.5 --header | samtools view -Sb - > "${filename}_pmd4.5_sorted.md.bam"
done