# =============================================================================
# Title: 1d_depth_vcf.sh
# Author: Laura N. Pott
# Date: 2026-06-02
# Description: Bash script to filter variant in VCFs by depth (3X, 5X, 10X)
# Usage: commands entered manually
# =============================================================================

module load bcftools/1.16-gcc-8.2.0-5d4xg4y

# 1. index vcf files
for vcf in *.mt.unifiedgenotyper.vcf.gz; do
bcftools index $vcf
done	

# 2. make summary files at 3 filtering levels for depth: 3X, 5X, 10X

### a. 3X
for vcf in *.mt.unifiedgenotyper.vcf.gz; do
base=$(basename "$vcf" .mt.unifiedgenotyper.vcf.gz)
bcftools view -i 'QUAL>=30 && INFO/DP>=3' "$vcf" -Oz -o "${base}_3x.vcf.gz"
tabix -p vcf "${base}_3x.vcf.gz"
done

### b. 5X
for vcf in *.mt.unifiedgenotyper.vcf.gz; do
base=$(basename "$vcf" .mt.unifiedgenotyper.vcf.gz)
bcftools view -i 'QUAL>=30 && INFO/DP>=5' "$vcf" -Oz -o "${base}_5x.vcf.gz"
tabix -p vcf "${base}_5x.vcf.gz"
done

### c. 10X
for vcf in *.mt.unifiedgenotyper.vcf.gz; do
base=$(basename "$vcf" .mt.unifiedgenotyper.vcf.gz)
bcftools view -i 'QUAL>=30 && INFO/DP>=10' "$vcf" -Oz -o "${base}_10x.vcf.gz"
tabix -p vcf "${base}_10x.vcf.gz"
done

# 3. make combined files
bcftools merge *_3x.vcf.gz -Oz -o pialq_3x.vcf.gz
tabix -p vcf pialq_3x.vcf.gz

bcftools merge *_10x.vcf.gz -Oz -o pialq_10x.vcf.gz
tabix -p vcf pialq_10x.vcf.gz

bcftools merge *_5x.vcf.gz -Oz -o pialq_5x.vcf.gz
tabix -p vcf pialq_5x.vcf.gz