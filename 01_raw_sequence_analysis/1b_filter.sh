#!/bin/bash -l

# =============================================================================
# Title: 1b_filter.sh
# Author: Laura N. Pott
# Date: 2026-06-02
# Description: SLURM script to filter primary alignments to mtDNA after mapping to GRCh37
# Usage: sbatch 1b_filter.sh
# =============================================================================

#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

cd $SHARED/projects/PIALQ/2025_remap_lp

module load samtools 
module load bedtools

bams=$SHARED/projects/PIALQ/2025_ancient_lp/01_wg_bwa/results/mapping/bwa
output_dir=$SHARED/projects/PIALQ/2025_ancient_lp/02_mtdna_circularmapper

# Summary variables
summary_file="$output_dir"/mtdna_filtering_summary.txt
echo -e "Sample\tRemoved_XA_X1\tTotal_Mapped_MT\tPercent_Removed" > "$summary_file"

for bam in "$bams"/*bam; do
    basename=$(basename "$bam" _PE.mapped.bam)

    # filter primary alignments to mt only, MQ ≥ 20
    samtools view -bh -F 4 -F 256 -F 2048 -q 20 "$bam" MT > "${basename}.mt.primary.bam"

    # convert to sam
    samtools view -h "${basename}.mt.primary.bam" > "${basename}.mt.primary.sam"

    # count reads with XA or X1 tags
    num_removed=$(grep -v '^@' "${basename}.mt.primary.sam" | grep -E 'XA:Z:|X1:Z:' | wc -l)
    total_reads=$(grep -v '^@' "${basename}.mt.primary.sam" | wc -l)
    percent_removed=$(awk "BEGIN {printf \"%.2f\", ($num_removed/$total_reads)*100}")

    # Add to summary
    echo -e "$basename\t$num_removed\t$total_reads\t$percent_removed" >> "$summary_file"

    # remove reads with XA or X1 tags
    grep -v -e 'XA:Z:' -e 'X1:Z:' "${basename}.mt.primary.sam" > "${basename}.mt.clean.sam"

    # convert back to BAM
    samtools view -bh "${basename}.mt.clean.sam" > "${basename}.mt.clean.bam"

    # convert to FASTQ
    bedtools bamtofastq -i "${basename}.mt.clean.bam" -fq "${basename}.mt.fastq"

    # clean up
    rm "${basename}.mt.primary.bam" "${basename}.mt.primary.sam" "${basename}.mt.clean.sam" "${basename}.mt.clean.bam"
    gzip "${basename}.mt.fastq"
    mv "${basename}.mt.fastq.gz" "$output_dir"

    echo "Done with $basename"
    echo
done