# =============================================================================
# Title: 2b_nopred_dsdna.sh
# Author: Laura N. Pott
# Date: 2026-06-03
# Description: SLURM job array to estimate contamination levels in dsDNA libraries with no contaminant prediction
# Usage: sbatch 2b_nopred_dsdna.sh
# =============================================================================


#!/bin/bash -l
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-57

# 1. set path to reference sequences, contam database, & go to working directory
working=$SHARED/projects/PIALQ/2025_ancient_lp/03_schmutzi/filt/dsdna
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
    "106mt_sorted.md.bam"
    "107mt_sorted.md.bam"
    "108mt_sorted.md.bam"
    "109mt_sorted.md.bam"
    "110mt_sorted.md.bam"
    "111mt_sorted.md.bam"
    "112mt_sorted.md.bam"
    "113mt_sorted.md.bam"
    "114mt_sorted.md.bam"
    "115mt_sorted.md.bam"
    "116mt_sorted.md.bam"
    "117mt_sorted.md.bam"
    "118mt_sorted.md.bam"
    "119mt_sorted.md.bam"
    "120mt_sorted.md.bam"
    "121mt_sorted.md.bam"
    "122mt_sorted.md.bam"
    "123mt_sorted.md.bam"
    "124mt_sorted.md.bam"
    "125mt_sorted.md.bam"
    "126mt_sorted.md.bam"
    "127mt_sorted.md.bam"
    "128mt_sorted.md.bam"
    "129mt_sorted.md.bam"
    "130mt_sorted.md.bam"
    "131mt_sorted.md.bam"
    "132mt_sorted.md.bam"
    "133mt_sorted.md.bam"
    "134mt_sorted.md.bam"
    "135mt_sorted.md.bam"
    "136mt_sorted.md.bam"
    "19dsV2_sorted.md.bam"
    "33dsV2_sorted.md.bam"
    "59mt_sorted.md.bam"
    "60mt_sorted.md.bam"
    "61mt_sorted.md.bam"
    "62mt_sorted.md.bam"
    "63mt_sorted.md.bam"
    "64mt_sorted.md.bam"
    "66mt_sorted.md.bam"
    "67mt_sorted.md.bam"
    "68mt_sorted.md.bam"
    "69mt_sorted.md.bam"
    "70mt_sorted.md.bam"
    "71mt_sorted.md.bam"
    "72mt_sorted.md.bam"
    "73mt_sorted.md.bam"
    "74mt_sorted.md.bam"
    "75mt_sorted.md.bam"
    "76mt_sorted.md.bam"
    "77mt_sorted.md.bam"
    "78mt_sorted.md.bam"
    "79mt_sorted.md.bam"
    "80mt_sorted.md.bam"
    "81mt_sorted.md.bam"
    "82mt_sorted.md.bam"
    "84mt_sorted.md.bam"
    "85mt_sorted.md.bam"
)
bam_file=${bam_list[$SLURM_ARRAY_TASK_ID]}

# 5. pick current bam file
filename=$(basename "$bam_file" _sorted.md.bam)

# 6. Run schmutzi pipeline
# a. run contDeam for initial contamination estimate based on deamination patterns
$SHARED/programs/schmutzi/src/contDeam.pl --lengthDeam 20 --library double --out "${filename}_contdeam" "$ref" "$bam_file"

# b. schmutzi with no prediction of contaminant
$SHARED/programs/schmutzi/src/schmutzi.pl --notusepredC --uselength --ref "$ref" --out "${filename}_schmutzi_npred" "${filename}_contdeam" "$database" "$bam_file"

# 7. organize outputs
mkdir -p "${filename}_schmutzi"
mv "${filename}"_contdeam* "${filename}_schmutzi/"
mv "${filename}"_schmutzi_npred* "${filename}_schmutzi/"
mv "${filename}"*bam "${filename}_schmutzi/"
mv "${filename}"*bam.bai "${filename}_schmutzi/"