# =============================================================================
# Title: 5e_outgroup.sh
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: script to download neanderthal reference mtDNA genome using NCBI's "E-Utilities"
# Usage: ru ninteractively
# =============================================================================

# neanderthal ref
efetch -db nucleotide -id AM948965 -format fasta > AM948965.fasta

# copy ref to each alignment folder
BASE="$SHARED/projects/PIALQ/2025_ancient_lp/09_msa/01_mafft"
for dir in $BASE/*/; do
cp $SHARED/ref_seqs/AM948965.fasta "$dir"
done