# 08_haplogrep

*$SHARED/projects/PIALQ/2025_ancient_lp/08_haplogrep/*

*finished 7/14/25*

```
# 1. define folder variables
output_folder=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/08_haplogrep/01_all
vcf_folder=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/07_mutect2/03_diroma

# 2. run haplogrep loop to get individual haplogrep results (top 3 hits)
for vcf in "$vcf_folder"/*norm.vcf.gz; do
	basename=$(basename "$vcf" _norm.vcf.gz)
	haplogrep3 classify --tree phylotree-rcrs@17.2 --in "$vcf" --extend-report --hits 3 --out "$output_folder/$basename.haplogrep.txt"
done
```



```
# 1. set output file name for haplogrep results
parent_folder=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/08_haplogrep/
output_folder=/home/mnievesc/shared/projects/PIALQ/2025_ancient_lp/08_haplogrep/01_all
output_haplogrep_file="$parent_folder/pialq_haplogrep_summary_vcf.txt"

# 2. clear or create output file
echo -e "Sample\tHaplogrepSampleID\tHaplogroup\tRank\tQuality\tRange\tNot_Found_Polys\tFound_Polys\tRemaining_Polys\tAAC_In_Remainings\tInput_Sample" > "$output_haplogrep_file"

# 3. generate one summary file
for file in "$output_folder"/*.haplogrep.txt; do 
    basename=$(basename "$file" .haplogrep.txt)
    second_line=$(sed -n '2p' "$file" | sed 's/"//g')
    echo -e "$basename\t${second_line}" >> "$output_haplogrep_file"
done
```



run on consensus fasta files

```
# 3. generate one summary file
output_folder=/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/06_haplogrep/01_fasta

echo -e "Sample\tHaplogrepSampleID\tHaplogroup\tRank\tQuality\tRange\tNot_Found_Polys\tFound_Polys\tRemaining_Polys\tAAC_In_Remainings\tInput_Sample" > pialq_fasta_haplogrep.txt

for file in "$output_folder"/f_*.haplogrep.txt; do 
    basename=$(basename "$file" .haplogrep.txt)
    second_line=$(sed -n '2p' "$file" | sed 's/"//g')
    echo -e "$basename\t${second_line}" >> pialq_fasta_haplogrep.txt
done

cut -f1,3,11 pialq_fasta_haplogrep.txt > short_pialq_fasta_haplogrep.txt
```

