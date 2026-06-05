#!/bin/bash -l

# =============================================================================
# Title: 1a_wg.sh
# Author: Laura N. Pott
# Date: 2026-06-02
# Description: SLURM script to process fastq reads: trim adapters, merge reads, map reads to GRCh37 using EAGER v2.5.0
# Usage: sbatch 1a_wg.sh
# =============================================================================

#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=160g
#SBATCH --tmp=160g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

# dependencies
nextflow=/home/mnievesc/shared/programs/nextflow/nextflow
module load singularity/current
module load java/openjdk-17.0.2

# set working directory
working=/home/mnievesc/shared/projects/PIALQ/2025_remap_lp/01_wg_bwa
cd "$working"

# run eager script through mapping
$nextflow pull nf-core/eager

NXF_VER='22.10.6' $nextflow run nf-core/eager \
-r 2.5.0 \
-c /home/mnievesc/shared/projects/PIALQ/2025_sequencing5_analyses_lp/nextflow.config \
-profile singularity \
-resume \
--input '*{R1,R2}_001.fastq.gz' \
--fasta /home/mnievesc/shared/ref_seqs/human_g1k_v37.fasta \
--fasta_index /home/mnievesc/shared/ref_seqs/human_g1k_v37.fasta.fai \
--outdir './results' \
--preserve5p \
--mergedonly \
--run_bam_filtering \
--bam_mapping_quality_threshold 30 \
--bam_filter_minreadlength 30 \
--bam_unmapped_type 'discard'

# clean intermediate files
$nextflow clean -f -k