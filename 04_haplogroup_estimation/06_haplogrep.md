# 06_haplogrep

*$SHARED/projects/PIALQ/2025_ancient_lp/06_haplogrep/*

**for fastas (these are all decontaminated)**

*finished 7/8/25*

```
# 1. define folder variables
output_folder=$SHARED/projects/PIALQ/2025_ancient_lp/06_haplogrep/01_all
fasta_folder=$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus/02_pialq
haplogrep=$SHARED/programs/haplogrep3

# 2. run haplogrep loop to get individual haplogrep results (top 3 hits)
for fasta in "$fasta_folder"/*fasta; do
	basename=$(basename "$fasta" .fasta)
	"$haplogrep" classify --tree phylotree-rcrs@17.2 --in "$fasta" --extend-report --hits 3 --out "$output_folder/$basename.haplogrep.txt"
done
```



```
# 1. set output file name for haplogrep results
parent_folder=$SHARED/projects/PIALQ/2025_ancient_lp/06_haplogrep/
output_folder=$SHARED/projects/PIALQ/2025_ancient_lp/06_haplogrep/01_all
output_haplogrep_file="$parent_folder/pialq_haplogrep_summary_fasta.txt"

# 2. clear or create output file
echo -e "Sample\tHaplogrepSampleID\tHaplogroup\tRank\tQuality\tRange\tNot_Found_Polys\tFound_Polys\tRemaining_Polys\tAAC_In_Remainings\tInput_Sample" > "$output_haplogrep_file"

# 3. generate one summary file
for file in "$output_folder"/f*.haplogrep.txt; do 
	basename=$(basename "$file" .haplogrep.txt)
    second_line=$(sed -n '2p' "$file" | sed 's/"//g')
    third_line=$(sed -n '3p' "$file" | sed 's/"//g')
    fourth_line=$(sed -n '4p' "$file" | sed 's/"//g')
    echo -e "$basename\t${second_line}" >> "$output_haplogrep_file"
done
```



**for vcfs (mixed contam/not contam)**

in $SHARED/projects/PIALQ/2025_ancient_lp/02_circmap/results/genotyping

```
# index
for vcf in *.mt.unifiedgenotyper.vcf.gz; do
bcftools index $vcf
done	

# make summary files at 3 filtering levels

# 3
for vcf in *.mt.unifiedgenotyper.vcf.gz; do
base=$(basename "$vcf" .mt.unifiedgenotyper.vcf.gz)
bcftools view -i 'QUAL>=30 && INFO/DP>=3' "$vcf" -Oz -o "${base}_3x.vcf.gz"
tabix -p vcf "${base}_3x.vcf.gz"
done

# 5
for vcf in *.mt.unifiedgenotyper.vcf.gz; do
base=$(basename "$vcf" .mt.unifiedgenotyper.vcf.gz)
bcftools view -i 'QUAL>=30 && INFO/DP>=5' "$vcf" -Oz -o "${base}_5x.vcf.gz"
tabix -p vcf "${base}_5x.vcf.gz"
done

# 10
for vcf in *.mt.unifiedgenotyper.vcf.gz; do
base=$(basename "$vcf" .mt.unifiedgenotyper.vcf.gz)
bcftools view -i 'QUAL>=30 && INFO/DP>=10' "$vcf" -Oz -o "${base}_10x.vcf.gz"
tabix -p vcf "${base}_10x.vcf.gz"
done

# make combined files
bcftools merge *_3x.vcf.gz -Oz -o pialq_3x.vcf.gz
tabix -p vcf pialq_3x.vcf.gz

bcftools merge *_10x.vcf.gz -Oz -o pialq_10x.vcf.gz
tabix -p vcf pialq_10x.vcf.gz

bcftools merge *_5x.vcf.gz -Oz -o pialq_5x.vcf.gz
tabix -p vcf pialq_5x.vcf.gz
```



**for 56mt_pmd3 and 78mt_pmd3 vcfs **

in $SHARED/projects/PIALQ/2025_ancient_lp/04_contam/02_pmd/05_vcf/genotyping

```
# index
for vcf in *pmd3.unifiedgenotyper.vcf.gz; do
bcftools index $vcf
done	

# make summary files at 3 filtering levels

# 3
for vcf in *pmd3.unifiedgenotyper.vcf.gz; do
base=$(basename "$vcf" .mt.unifiedgenotyper.vcf.gz)
bcftools view -i 'QUAL>=30 && INFO/DP>=3' "$vcf" -Oz -o "${base}_3x.vcf.gz"
tabix -p vcf "${base}_3x.vcf.gz"
done

# 5
for vcf in *pmd3.unifiedgenotyper.vcf.gz; do
base=$(basename "$vcf" .mt.unifiedgenotyper.vcf.gz)
bcftools view -i 'QUAL>=30 && INFO/DP>=5' "$vcf" -Oz -o "${base}_5x.vcf.gz"
tabix -p vcf "${base}_5x.vcf.gz"
done

# 10
for vcf in *pmd3.unifiedgenotyper.vcf.gz; do
base=$(basename "$vcf" .mt.unifiedgenotyper.vcf.gz)
bcftools view -i 'QUAL>=30 && INFO/DP>=10' "$vcf" -Oz -o "${base}_10x.vcf.gz"
tabix -p vcf "${base}_10x.vcf.gz"
done

# make combined files
bcftools merge *_3x.vcf.gz -Oz -o decontam_3x.vcf.gz
tabix -p vcf decontam_3x.vcf.gz

bcftools merge *_10x.vcf.gz -Oz -o decontam_10x.vcf.gz
tabix -p vcf decontam_10x.vcf.gz

bcftools merge *_5x.vcf.gz -Oz -o decontam_5x.vcf.gz
tabix -p vcf decontam_5x.vcf.gz
```







variant depth check on realigned bam files:

```
# libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("Rsamtools")
install.packages("readr")

library("Rsamtools")
library("GenomicRanges")
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(readr)

# 1. read in text file
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/06_haplogrep")
data <- read_tsv("data.txt", na = "-", quote = "")

# drop empty paths
data <- data %>%
  filter(!is.na(path) & path != "")

# 2. function to get depth at positions using samtools
get_depth_samtools <- function(bam_file, positions, 
                               chrom = "gi|251831106|ref|NC_012920.1|",
                               samtools_path = "samtools") {
  # --- 1. Input checks ---
  if (!file.exists(bam_file)) {
    stop("BAM file not found: ", bam_file)
  }
  if (!length(positions)) {
    stop("No positions provided.")
  }
  
  # --- 2. Create temporary BED file ---
  bed_file <- tempfile(fileext = ".bed")
  bed_data <- data.frame(
    chrom = chrom,
    start = as.integer(positions) - 1,  # BED is 0-based
    end   = as.integer(positions)
  )
  write.table(
    bed_data, bed_file, sep = "\t", quote = FALSE, 
    row.names = FALSE, col.names = FALSE
  )
  
  # --- 3. Build and run samtools command ---
  cmd <- sprintf('%s depth -a -b "%s" "%s"', samtools_path, bed_file, bam_file)
  depth_output <- tryCatch(system(cmd, intern = TRUE), 
                           warning = function(w) character(0),
                           error   = function(e) character(0))
  
  # --- 4. Parse output ---
  if (length(depth_output) > 0) {
    depth_df <- read.table(text = depth_output, 
                           col.names = c("chrom", "pos", "depth"))
  } else {
    depth_df <- data.frame(chrom = chrom, pos = positions, depth = 0)
  }
  
  # --- 5. Clean up temporary file ---
  unlink(bed_file)
  
  # --- 6. Return result ---
  return(depth_df)
}

# 3. change mutation list to numbers not strings
data <- data %>%
  mutate(expected_mutations = str_split(expected_mutations, ",\\s*")) %>%
  mutate(expected_mutations = lapply(expected_mutations, function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    x
  }))

data <- data %>%
  mutate(bam_path = file.path(path, file))

# 4. Calculate depth for each sample at their haplogroup's diagnostic positions
data <- data %>%
  mutate(depths = map2(bam_path, expected_mutations, 
                       ~ get_depth_samtools(bam_file = .x, positions = .y)))

depths_long <- data %>%
  select(individual, file, reported_haplogroup, depths) %>%
  unnest(depths)

# 5. summary statistics
summary_library <- depths_long %>%
  group_by(file) %>%
  summarise(
    mean_depth = mean(depth, na.rm = TRUE),
    median_depth = median(depth, na.rm = TRUE),
    min_depth = min(depth, na.rm = TRUE),
    max_depth = max(depth, na.rm = TRUE),
    n_positions = n()
  )

summary_ind <- depths_long %>%
  group_by(individual) %>%
  summarise(
    mean_depth = mean(depth, na.rm = TRUE),
    median_depth = median(depth, na.rm = TRUE),
    min_depth = min(depth, na.rm = TRUE),
    max_depth = max(depth, na.rm = TRUE),
    n_positions = n()
  )

# Save results
write.table(summary_ind, "expected_mutations_summary_individual.txt", sep = "\t", row.names = FALSE)
write.table(summary_library, "expected_mutations_summary_library.txt", sep = "\t", row.names = FALSE)
write.table(depths_long, "expected_mutations_all.txt", sep = "\t", row.names = FALSE)

```

