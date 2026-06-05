# =============================================================================
# Title: 5l_network.sh
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: script to create NEXUS files for popART
# Usage: ru ninteractively
# =============================================================================

# 1. create nexus files using AMAS
module load python3/3.8.3_anaconda2020.07_mamba

base="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn"
python $SHARED/programs/AMAS/amas/AMAS.py convert -i "$base"/mjn*fasta -f fasta -u nexus -d dna

# 2. replace all N nucleotides with "?"
for nexus in mjn*.fasta-out.nex; do
name=$(basename "$nexus" .fasta-out.nex)
python 5l_replace_n.py "$nexus" "rep_${name}.nex"
done

# 3. create nexus files with trait information from different metadata variables
meta="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn/mjn.meta"

# a. paper sequence came from (author_year)
for nexus in rep_mjn*.nex; do
name=$(basename "$nexus" .nex)
python 5l_add_traits_paper.py "$nexus" "$meta" popart/"${name}_paper.nex"
done

# b. country of individual who gave dna sample
for nexus in rep_mjn*.nex; do
name=$(basename "$nexus" .nex)
python 5l_add_traits_country.py "$nexus" "$meta" popart/"${name}_country.nex"
done

# c. haplogroup (different depths L3 vs L3e etc
three=("rep_mjn_d1.nex")
four=("rep_mjn_l2c.nex")
five=("rep_mjn_b2b.nex" "rep_mjn_c1b.nex" "rep_mjn_h1bw.nex" "rep_mjn_l1c1.nex" "rep_mjn_l1c24.nex" "rep_mjn_l2b1.nex" "rep_mjn_l3d1.nex" "rep_mjn_l3d3.nex" "rep_mjn_l3d4.nex" "rep_mjn_l3e1.nex" "rep_mjn_l3e2.nex" "rep_mjn_l3e3.nex")
six=("rep_mjn_l0a2a.nex" "rep_mjn_l0d1a.nex" "rep_mjn_l2a1a.nex" "rep_mjn_l2a1d.nex" "rep_mjn_l2a1f.nex" "rep_mjn_l2a1q.nex" "rep_mjn_l3b1a.nex" "rep_mjn_l3e3b.nex" "rep_mjn_l3f1b.nex")
eight=("rep_mjn_a2_64.nex" "rep_mjn_l3e1a1a.nex")

for nexus in "${three[@]}"; do
name=$(basename "$nexus" .nex)
python 5l_add_traits_hap.py "$nexus" "$meta" popart/"${name}_hap.nex" hap_3char
done

for nexus in "${four[@]}"; do
name=$(basename "$nexus" .nex)
python 5l_add_traits_hap.py "$nexus" "$meta" popart/"${name}_hap.nex" hap_4char
done

for nexus in "${five[@]}"; do
name=$(basename "$nexus" .nex)
python 5l_add_traits_hap.py "$nexus" "$meta" popart/"${name}_hap.nex" hap_5char
done

for nexus in "${six[@]}"; do
name=$(basename "$nexus" .nex)
python 5l_add_traits_hap.py "$nexus" "$meta" popart/"${name}_hap.nex" hap_6char
done

for nexus in "${eight[@]}"; do
name=$(basename "$nexus" .nex)
python 5l_add_traits_hap.py "$nexus" "$meta" popart/"${name}_hap.nex" Haplogroup
done