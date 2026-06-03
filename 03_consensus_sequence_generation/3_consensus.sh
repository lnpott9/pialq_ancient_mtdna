# =============================================================================
# Title: 3_consensus.sh
# Author: Laura N. Pott
# Date: 2026-06-03
# Description: script to generate consensus fasta sequences
# Usage: run interactively
# =============================================================================

# 1. start interactive session
srun --mem=4gb --time=01:00:00 --pty bash

# s2. et paths to executable and output directory
log2fasta=$SHARED/programs/schmutzi/src/log2fasta
output=$SHARED/projects/PIALQ/2025_ancient_lp/05_consensus/02_pialq

# 3. loop through files and call consensus fasta with minimum of base quality 30
# npred: 56mt_pmd3, 78mt_pmd3, all other libraries
for log in *_schmutzi_npred_final_endo.log; do
	basename=$(basename "$log" _schmutzi_npred_final_endo.log)
	"$log2fasta" -q 30 -name "${basename}" -indel 20 "$log" > "$output/${basename}.fasta"
done

# wpred: 19dsV2
for log in *_schmutzi_wpred_final_endo.log; do
	basename=$(basename "$log" _schmutzi_wpred_final_endo.log)
	"$log2fasta" -q 30 -name "${basename}" -indel 20 "$log" > "$output/${basename}.fasta"
done