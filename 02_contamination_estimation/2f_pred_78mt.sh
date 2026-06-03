# =============================================================================
# Title: 2f_pred_78mt.sh
# Author: Laura N. Pott
# Date: 2026-06-03
# Description: SLURM job array to estimate contamination levels in decontaminated dsDNA library 78mt with contaminant prediction, worldwide database
# Usage: sbatch 2f_pred_78mt.sh
# =============================================================================

#!/bin/bash -l
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-5

# 1. set path to reference sequences, contam database, & go to working directory
working=$SHARED/projects/PIALQ/2025_ancient_lp/04_contam/04_wpred
ref=$SHARED/ref_seqs/rCRS.fasta
database=$SHARED/programs/schmutzi/share/schmutzi/alleleFreqMT/197/freqs
cd "$working"

# 2. load modules
module load samtools/1.21
module load R/4.4.0-openblas-rocky8

# 3. set R library path
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"

# 4. create a file list of sorted, MD-tagged bams and select one
bam_list=(
    "78mt_pmd2.5_sorted.md.bam"
    "78mt_pmd2_sorted.md.bam"
    "78mt_pmd3.5_sorted.md.bam"
    "78mt_pmd3_sorted.md.bam"
    "78mt_pmd4.5_sorted.md.bam"
    "78mt_pmd4_sorted.md.bam"
)
bam_file=${bam_list[$SLURM_ARRAY_TASK_ID]}

# 5. pick current bam file
filename=$(basename "$bam_file" _sorted.md.bam)

# 6. Run schmutzi pipeline
# a. run contDeam for initial contamination estimate based on deamination patterns
$SHARED/programs/schmutzi/src/contDeam.pl --lengthDeam 20 --library double --out "${filename}_contdeam" "$ref" "$bam_file"

# b. schmutzi with prediction
$SHARED/programs/schmutzi/src/schmutzi.pl --uselength --ref "$ref" --out "${filename}_schmutzi_wpred" "${filename}_contdeam" "$database" "$bam_file"

# 7. organize outputs
mkdir -p "${filename}_schmutzi"
mv "${filename}"_contdeam* "${filename}_schmutzi/"
mv "${filename}"_schmutzi_wpred* "${filename}_schmutzi/"
mv "${filename}"*bam "${filename}_schmutzi/"
mv "${filename}"*bam.bai "${filename}_schmutzi/"