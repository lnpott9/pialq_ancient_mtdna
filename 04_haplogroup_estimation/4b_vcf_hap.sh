# =============================================================================
# Title: 4b_vcf_hap.sh
# Author: Laura N. Pott
# Date: 2026-06-03
# Description: script to generate hapogroup calls from vcf files
# Usage: run interactively
# =============================================================================

# 1. define folder variables
output_folder=$SHARED/projects/PIALQ/2025_ancient_lp/08_haplogrep/01_all
vcf_folder=$SHARED/projects/PIALQ/2025_ancient_lp/02_circmap/results/genotyping
haplo=$SHARED/programs/haplogrep3

# 2. run haplogrep loop to get individual haplogrep results (top 3 hits)
for vcf in "$vcf_folder"/*3x.vcf.gz; do
	basename=$(basename "$vcf" _3x.vcf.gz)
	"$haplo" classify --tree phylotree-rcrs@17.2 --in "$vcf" --extend-report --hits 3 --out "$output_folder/${basename}_3x.haplogrep.txt"
done

for vcf in "$vcf_folder"/*5x.vcf.gz; do
	basename=$(basename "$vcf" _5x.vcf.gz)
	"$haplo" classify --tree phylotree-rcrs@17.2 --in "$vcf" --extend-report --hits 3 --out "$output_folder/${basename}_5x.haplogrep.txt"
done

for vcf in "$vcf_folder"/*10x.vcf.gz; do
	basename=$(basename "$vcf" _10x.vcf.gz)
	"$haplo" classify --tree phylotree-rcrs@17.2 --in "$vcf" --extend-report --hits 3 --out "$output_folder/${basename}_10x.haplogrep.txt"
done

# 3. create summary file
output_haplogrep_file="$SHARED/projects/PIALQ/2025_ancient_lp/08_haplogrep/pialq_vcf_haplogrep.txt"
echo -e "Sample\tHaplogrepSampleID\tHaplogroup\tRank\tQuality\tRange\tNot_Found_Polys\tFound_Polys\tRemaining_Polys\tAAC_In_Remainings\tInput_Sample" > "$output_haplogrep_file"

# 3. generate one summary file
for file in "$output_folder"/*.haplogrep.txt; do 
    basename=$(basename "$file" .haplogrep.txt)
    second_line=$(sed -n '2p' "$file" | sed 's/"//g')
    echo -e "$basename\t${second_line}" >> "$output_haplogrep_file"
done