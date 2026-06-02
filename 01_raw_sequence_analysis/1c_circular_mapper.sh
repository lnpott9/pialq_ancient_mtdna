# =============================================================================
# Title: 1c_circular_mapper.sh
# Author: Laura N. Pott
# Date: 2026-06-02
# Description: SLURM script to map reads to rCRS using EAGER v2.5.0 (circular mapper), deduplicate, and call variants
# Usage: sbatch 1c_circular_mapper.sh
# =============================================================================

#!/bin/bash -l
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128g
#SBATCH --tmp=128g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

# dependencies
nextflow=/home/mnievesc/shared/programs/nextflow/nextflow
module load singularity/current
module load java/openjdk-17.0.2

# set working directory
working=$SHARED/projects/PIALQ/2025_ancient_lp/02_circmap
cd "$working"

# pull most recent version of eager
NXF_VER='22.10.6' $nextflow pull nf-core/eager

# run EAGER script
NXF_VER='22.10.6' $nextflow run nf-core/eager \
-r 2.5.0 \
-profile singularity \
-c $working/nextflow.config \
--resume \
--input '*fastq.gz' \
--single_end \
--fasta $SHARED/ref_seqs/rCRS.fasta \
--fasta_index $SHARED/ref_seqs/rCRS.fasta.fai \
--outdir './results' \
--skip_fastqc \
--skip_adapterremoval \
--preserve5p \
--mergedonly \
--mapper 'circularmapper' \
--circulartarget 'gi\|251831106\|ref\|NC_012920.1\|' \
--run_bam_filtering \
--bam_mapping_quality_threshold 30 \
--bam_filter_minreadlength 30 \
--bam_unmapped_type 'discard' \
--dedupper 'dedup' \
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