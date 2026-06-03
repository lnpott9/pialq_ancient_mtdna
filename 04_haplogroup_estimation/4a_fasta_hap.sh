# =============================================================================
# Title: 4a_fasta_hap.sh
# Author: Laura N. Pott
# Date: 2026-06-03
# Description: script to generate hapogroup calls from consensus fasta sequences
# Usage: run interactively
# =============================================================================

# 1. define folder variables
output_folder=$SHARED/projects/PIALQ/2025_ancient_lp/06_haplogrep/01_fasta
fasta_folder=$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus/02_pialq
haplo=$SHARED/programs/haplogrep3

# 2. run haplogrep loop to get individual haplogrep results (top 3 hits per individual)
for fasta in "$fasta_folder"/*fasta; do
	basename=$(basename "$fasta" .fasta)
	"$haplo" classify --tree phylotree-rcrs@17.2 --in "$fasta" --extend-report --hits 3 --out "$output_folder/$basename.haplogrep.txt"
done

# 3. generate one summary file with top hit for each individual
echo -e "Sample\tHaplogrepSampleID\tHaplogroup\tRank\tQuality\tRange\tNot_Found_Polys\tFound_Polys\tRemaining_Polys\tAAC_In_Remainings\tInput_Sample" > "$output_folder"/pialq_fasta_haplogrep.txt

for file in "$output_folder"/*.haplogrep.txt; do 
    basename=$(basename "$file" .haplogrep.txt)
    second_line=$(sed -n '2p' "$file" | sed 's/"//g')
    echo -e "$basename\t${second_line}" >> "$output_folder"/pialq_fasta_haplogrep.txt
done