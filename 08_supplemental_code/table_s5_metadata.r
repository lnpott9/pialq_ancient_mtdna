# =============================================================================
# Title: table_s5_metadata.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to create Supplemental Table S5
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. set wd
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa")
library(dplyr)

# 2. getting all papers
tree <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/03_iqtree/iqtree.meta", sep = "\t")
network <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn/mjn.meta", sep = "\t")
diversity <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/diversity.meta", sep = "\t")

# 3. bind together
all_meta <- bind_rows(tree, network, diversity)
all_meta <- unique(all_meta)
meta <- all_meta[!(all_meta$file %in% c("AM948965.fasta")), ]

# 4. summarize using dplyr
all_summary <- meta %>%
  group_by(ancient_or_modern, paper, country, hap_2char) %>%
  summarize(.groups = "keep") %>%
  count()

tree_summary <- tree %>%
  group_by(ancient_or_modern, paper, country, hap_2char) %>%
  summarize(.groups = "keep") %>%
  count()

network_summary <- network %>%
  group_by(ancient_or_modern, paper, country, hap_2char) %>%
  summarize(.groups = "keep") %>%
  count()

diversity_summary <- diversity %>%
  group_by(ancient_or_modern, paper, country, hap_2char) %>%
  summarize(.groups = "keep") %>%
  count()

# 5. save separate summaries
write.table(all_summary, file = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/summ_all.txt", quote = FALSE, sep = "\t", col.names = TRUE)
write.table(tree_summary, file = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/summ_tree.txt", quote = FALSE, sep = "\t", col.names = TRUE)
write.table(network_summary, file = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/summ_network.txt", quote = FALSE, sep = "\t", col.names = TRUE)
write.table(diversity_summary, file = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/summ_diversity.txt", quote = FALSE, sep = "\t", col.names = TRUE)

# 6. iff accesion is in tree, add columns that says yes
# if accession is in network, add column that says yes
# if accession is in diversity, add column that says yes
meta <- meta %>%
  mutate(
    used_in_trees = ifelse(file %in% tree$file, "yes", "no"),
    used_in_network = ifelse(file %in% network$file, "yes", "no"),
    used_in_comparative_analyses = ifelse(file %in% diversity$file, "yes", "no")
  )

# get DOIs
unique(meta$paper)

# 7. create mapping for new column of doi
doi <- c(
    "arencibia_2023" = "10.1002/ajpa.24822",
    "bodner_2012" = "10.1101/gr.131722.111",
    "bodner_2015" = "10.1016/j.fsigen.2014.09.012 ",
    "delgado_2021" = "10.1016/j.quaint.2020.08.031",
    "derenko_2007" = "10.1086/522933",
    "fehrenschmitz_2015" = "10.1371/journal.pone.0127141",
    "garcia_2021" = "10.1093/hmg/ddab105",      
    "harney_2023" = "10.1126/science.ade4995",
    "heinz_2015" = "10.1371/journal.pone.0134129",
    "hgdp" = "10.1126/science.aay501",            
    "just_2008" = "10.1016/j.fsigen.2007.12.001",
    "llamas_2016" = "10.1126/sciadv.1501385",
    "nakatsuka_2020" = "10.1016/j.cell.2020.04.015 ",   
    "posth_2018" = "10.1016/j.cell.2018.10.027 ",
    "rocarada_2021" = "10.1016/j.isci.2021.102553",
    "russo_2023" = "10.1002/ajpa.24795",          
    "sandovalvelasco_2023" = "10.1016/j.ajhg.2023.08.001",
    "schroeder_2015" = "10.1073/pnas.1421784112",
    "silva_2015" = "10.1038/srep12526",
    "sgdp" = "10.1038/nature18964",
    "fleskes_2023" = "10.1073/pnas.2201620120",        
    "1kg" = "10.1038/nature15393",
    "brandini_2018" = "10.1093/molbev/msx267",
    "cardoso_2012" = "10.1038/hdy.2011.131",        
    "gonder_2007" = "10.1093/molbev/msl209",
    "tamm_2007" = "10.1371/journal.pone.0000829",
    "huber_2025" = "10.1038/s41598-025-08241-6",          
    "garciaolivares_2022" = "10.1016/j.isci.2022.105907",
    "pialq" = "this study",
    "arias_2017" = " 10.1002/ajpa.23345",          
    "avila_2019" = "10.1016/j.fsigen.2019.07.004",
    "barbieri_2012a" = "10.1093/molbev/msr291 ",
    "barbieri_2013a" = "10.1038/ejhg.2012.192",  
    "barbieri_2013" = "10.1016/j.ajhg.2012.12.010 ",
    "barbieri_2014a" = " 10.1371/journal.pone.0099117",
    "barbieri_2014b" = "10.1002/ajpa.22441",
    "cabrera_2018" = "10.1186/s12862-018-1211-4",        
    "brucato_2018" = "10.1016/j.ajhg.2017.11.011 ",
    "chan_2019" = " 10.1038/s41586-019-1714-1",
    "colombo_2025" = "10.1038/s41598-025-12209-x",        
    "hernandez_2015" = "10.1371/journal.pone.0139784",
    "just_2015" = "10.1016/j.fsigen.2014.09.021",
    "mccrow_2016" = " 10.1002/pros.23126",         
    "oliveira_2018" = "10.1002/ajpa.23378",
    "tito_unpub" = "unpublished, GenBank accessions JX669136-327",
    "kumar_2011" = "10.1186/1471-2148-11-293",          
    "barbieri_2017" = "10.1038/s41598-017-17728-w",
    "olivieri_2017" = "10.1093/molbev/msx082",       
    "pierron_2017" = "10.1073/pnas.1704906114",
    "rito_2013" = "10.1371/journal.pone.0080031",
    "silva_2021" = "10.1038/s41598-021-95996-3",          
    "simao_2019" = "10.1016/j.fsigen.2018.12.007",
    "soares_2012" = "10.1093/molbev/msr245",
    "taylor_2020" = " 10.3390/genes11111290",         
    "torroni_2006" = "10.1016/j.tig.2006.04.001",
    "tribaldos_2021" = "10.1016/j.cell.2021.02.040",
    "pardinas_2014" = "10.1002/ajhb.22601",       
    "gomezcarballa_2018" = "10.1101/gr.234674.118",
    "lee_2015" = "10.13110/humanbiology.87.1.0029",
    "ottoni_2010" = "10.1371/journal.pone.0013378",         
    "costa_2009" = "10.1016/j.mad.2008.12.001"
)

# 8. join based on paper name
meta$doi <- doi[meta$paper]
sum(is.na(meta$doi))

# find missing w/ unique(meta$paper[is.na(meta$doi)])

# 9. output large text file with accession number, paper, Haplogroup, country, extra_info, and doi
final <- meta %>%
  select(Sample, paper, doi, ancient_or_modern, original_format, country, extra_info, 
         Haplogroup, Quality, hap_2char, used_in_trees, used_in_network,
         used_in_comparative_analyses)

# 10. write out table to final metadata file
write.table(final, file = "/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/supplemental_metadata.txt", quote = FALSE, sep = "\t", col.names = TRUE)
