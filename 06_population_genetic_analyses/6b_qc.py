# =============================================================================
# Title: 6b_qc.py
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: Python script to create new metadata file for seqs included in papers selected for pop gen analyses
# Usage: run interactively
# =============================================================================

module load python3/3.8.3_anaconda2020.07_mamba

python

# 2. import pandas
import pandas as pd

# 3. read in all_merged.meta as dataframe
meta = pd.read_csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/all_merged.meta", sep="\t")

# 4. only retain specific papers
papers = ["1kg", "arias_2017", "avila_2019", "barbieri_2012a", "barbieri_2013a", "barbieri_2013", "barbieri_2014a", "barbieri_2014b", "barbieri_2017", "bodner_2022", "brandini_2018", "brucato_2018", "chan_2019", "colombo_2025", "derenko_2007", "fleskes_2023", "garcia_2021", "garciaolivares_2022", "harney_2023", "hgdp", "huber_2025", "just_2008", "just_2015", "kumar_2011", "mccrow_2016", "oliveira_2018", "olivieri_2017", "pialq", "pierron_2017", "sandovalvelasco_2023", "silva_2021", "taylor_2020", "tito_unpub", "tribaldos_2021"]

subset = meta[meta["paper"].isin(papers)]
subset.shape # (17390, 21)

# 5. subset for length
ancient = subset[(subset["ancient_or_modern"] == "ancient")]
ancient = ancient.sort_values("length")
# remove everything < 10000 (4 harney samples)
quality = subset[(subset["length"] > 10000)]
quality.shape # (17384, 21)

# 6. write out new meta
quality.to_csv("diversity.meta", sep="\t", index=False, header=True)