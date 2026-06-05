#!/bin/bash -l

# =============================================================================
# Title: 5b_ref_haplogroup.sh
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: SLURM script to generate hapogroup calls from fasta files
# Usage: sbatch 5b_ref_haplogroup.sh
# =============================================================================

#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --tmp=2g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu

# 1. set paths
# parent folder is analysis folder
parent_folder="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/"
temp="$parent_folder/tmp"
mkdir -p "$temp"

haplogrep="$SHARED/programs/haplogrep3"

# set paths to input fasta files
a_folder="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/ancient_fasta"
m_folder="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/00_ref/03_organized_ref_fastas/modern_fasta"

# set paths to output files & clear
summary_file="$parent_folder/all_hap_summ.txt"
short_summary_file="$parent_folder/short_all_hap_summ.txt"
> "$summary_file"
> "$short_summary_file"

# 2. add headings to full summary file
echo -e "Sample\tHaplogrepSampleID\tHaplogroup\tRank\tQuality\tRange\tNot_Found_Polys\tFound_Polys\tRemaining_Polys\tAAC_In_Remainings\tInput_Sample" > "$summary_file"

# 3. run haplogrep to generate the top hit
for folder in "$a_folder"/* "$m_folder"/*; do
    for fasta in "$folder"/*.fasta; do
        sample=$(basename "$fasta" .fasta)      
        output_tmp="$temp/${sample}.haplogrep.txt"
        "$haplogrep" classify \
            --tree phylotree-rcrs@17.2 \
            --in "$fasta" \
            --extend-report \
            --hits 1 \
            --out "$output_tmp"
        # extract second line from output and append to summary file
        if [ -s "$output_tmp" ]; then
            second_line=$(awk 'NR==2 {gsub(/"/, ""); print}' "$output_tmp")
            echo -e "$sample\t${second_line}" >> "$summary_file"
        else
            echo "⚠️⚠️ no output for $sample ⚠️⚠️"
        fi
    done
done

# 4. create short version of summary file
cut -f1-5 "$summary_file" > "$short_summary_file"