#!/bin/bash -l

# =============================================================================
# Title: 2g_vcf_decontam.sh
# Author: Laura N. Pott
# Date: 2026-06-03
# Description: SLURM script to regenerate vcfs for 56mt and 78mt based on PMD-filtered BAM files at PMD=3
# Usage: sbatch 2g_vcf_decontam.sh
# =============================================================================

#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

# dependencies
nextflow=/home/mnievesc/shared/programs/nextflow/nextflow
module load singularity/current
module load java/openjdk-17.0.2

# set working directory
working=$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/02_pmd
cd "$working"

# pull most recent version of eager
NXF_VER='22.10.6' $nextflow pull nf-core/eager

# run EAGER script
NXF_VER='22.10.6' $nextflow run nf-core/eager \
-r 2.5.0 \
-profile singularity \
--input 'f*pmd3_sorted.md.bam' \
--single_end \
--bam \
--fasta $SHARED/ref_seqs/rCRS.fasta \
--fasta_index $SHARED/ref_seqs/rCRS.fasta.fai \
--outdir './05_vcf' \
--skip_fastqc \
--skip_adapterremoval \
--skip_preseq \
--skip_deduplication \
--skip_qualimap \
--damage_calculation_tool 'mapdamage' \
--run_mapdamage_rescaling \
--run_genotyping \
--genotyping_tool 'ug' \
--genotyping_source 'rescaled' \
--gatk_ploidy 1 \
--gatk_ug_out_mode 'EMIT_ALL_CONFIDENT_SITES' \
--gatk_ug_genotype_model 'BOTH' \
--gatk_ug_keep_realign_bam

# clean 
$nextflow clean -f -k