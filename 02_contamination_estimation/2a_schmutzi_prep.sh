# =============================================================================
# Title: 2a_schmutzi_prep.sh
# Author: Laura N. Pott
# Date: 2026-06-03
# Description: script to process fastq reads: trim adapters, merge reads, map reads to GRCh37 using EAGER v2.5.0
# Usage: run interactively
# =============================================================================

# 1. set paths to rCRS and deduplicated BAM files (separated by library treatment)
ref=$SHARED/ref_seqs/rCRS.fasta
working=$SHARED/projects/PIALQ/2025_ancient_lp/03_schmutzi/filt/dsdna
working=$SHARED/projects/PIALQ/2025_ancient_lp/03_schmutzi/filt/ssdna

# 2. load samtools
module load samtools/1.21

# 3. sort all BAM files, add MD tags, then index
for bam_file in "$working"/*rmdup.bam; do
	filename=$(basename "$bam_file" .mt_rmdup.bam)
	samtools sort "$bam_file" > "${filename}_sorted.bam"
	samtools calmd -b "${filename}_sorted.bam" $ref > "${filename}_sorted.md.bam"
	samtools index "${filename}_sorted.md.bam"
done