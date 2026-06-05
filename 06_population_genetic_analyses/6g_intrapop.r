# =============================================================================
# Title: 6g_intrapop.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to calculate measures of intrapopulation diversity
# Usage: run interactively
# =============================================================================

# 1. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/03_intra")

# 2. load libraries, download if necessary
library(seqinr)
library(pegas)
library(dplyr)
library(ape)
library(ggplot2)

# 3. load metadata
meta <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/diversity.meta", sep="\t")
dim(meta)

# 4. load aligned sequences
align <- read.dna("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/02_clean/clean_div.fasta", format = "fasta")

# check dimensions
dim(align) # rows = sequences, columns = alignment length
nrow(align) # number of sequences
ncol(align) # alignment length

# 4. change sequence names to accession numbers
rownames(align) <- sub("^[^.]+\\.", "", rownames(align))
head(rownames(align))

# drop duplicates
align_dedup <- align[!duplicated(rownames(align)), ]

# 5. make population metadata subsets
# 1kg
onekg <- meta[(meta$paper == "1kg"),]
acb <- meta[(meta$paper == "1kg") & (meta$country == "ACB"),] # 98
asw <- meta[(meta$paper == "1kg") & (meta$country == "ASW"),] # 66
beb <- meta[(meta$paper == "1kg") & (meta$country == "BEB"),] # 104
cdx <- meta[(meta$paper == "1kg") & (meta$country == "CDX"),] # 107
ceu <- meta[(meta$paper == "1kg") & (meta$country == "CEU"),] # 101
chb <- meta[(meta$paper == "1kg") & (meta$country == "CHB"),] # 105
chs <- meta[(meta$paper == "1kg") & (meta$country == "CHS"),] # 112
clm <- meta[(meta$paper == "1kg") & (meta$country == "CLM"),] # 103
esn <- meta[(meta$paper == "1kg") & (meta$country == "ESN"),] # 111
fin <- meta[(meta$paper == "1kg") & (meta$country == "FIN"),] # 105
gbr <- meta[(meta$paper == "1kg") & (meta$country == "GBR"),] # 106
gih <- meta[(meta$paper == "1kg") & (meta$country == "GIH"),] # 109
gwd <- meta[(meta$paper == "1kg") & (meta$country == "GWD"),] # 120
ibs <- meta[(meta$paper == "1kg") & (meta$country == "IBS"),] # 108
itu <- meta[(meta$paper == "1kg") & (meta$country == "ITU"),] # 112
jpt <- meta[(meta$paper == "1kg") & (meta$country == "JPT"),] # 105
khv <- meta[(meta$paper == "1kg") & (meta$country == "KHV"),] # 103
lwk <- meta[(meta$paper == "1kg") & (meta$country == "LWK"),] # 112
msl <- meta[(meta$paper == "1kg") & (meta$country == "MSL"),] # 98
mxl <- meta[(meta$paper == "1kg") & (meta$country == "MXL"),] # 70
pjl <- meta[(meta$paper == "1kg") & (meta$country == "PJL"),] # 108
pel <- meta[(meta$paper == "1kg") & (meta$country == "PEL"),] # 90
pur <- meta[(meta$paper == "1kg") & (meta$country == "PUR"),] # 107
stu <- meta[(meta$paper == "1kg") & (meta$country == "STU"),] # 107
tsi <- meta[(meta$paper == "1kg") & (meta$country == "TSI"),] # 111
yri <- meta[(meta$paper == "1kg") & (meta$country == "YRI"),] # 109

# hgdp
# makrani_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MakraniHGDP"),] # 23
# mbuti_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MbutiHGDP"),] # 10
# biaka_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "BiakaHGDP"),] # 24
# french_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "FrenchHGDP"),] # 24
# sardinian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "SardinianHGDP"),] # 24
# colombian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "ColombianHGDP"),] # 5
# cambodian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "CambodianHGDP"),] # 8
# japanese_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "JapaneseHGDP"),] # 27
# han_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "HanHGDP"),] # 30
# surui_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "SuruiHGDP"),] # 6
# mandenka_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MandenkaHGDP"),] # 20
# yoruba_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "YorubaHGDP"),] # 18
# san_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "SanHGDP"),] # 2
# bantu_southafrica_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "BantuSouthAfricaHGDP"),] # 4
# karitiana_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "KaritianaHGDP"),] # 9
# tujia_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "TujiaHGDP"),] # 8
# bergamo_italian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "BergamoItalianHGDP"),] # 10
# tuscan_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "TuscanHGDP"),] # 6
# yi_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "YiHGDP"),] # 8
# miao_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MiaoHGDP"),] # 8
# oroqen_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "OroqenHGDP"),] # 7
# daur_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "DaurHGDP"),] # 9
# mongolian_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MongolianHGDP"),] # 8
# xibo_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "XiboHGDP"),] # 7
# northern_han_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "NorthernHanHGDP"),] # 10
# uygur_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "UygurHGDP"),] # 8
# daihgdp_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "DaiHGDP"),] # 5
# lahu_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "LahuHGDP"),] # 6
# she_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "SheHGDP"),] # 8
# naxi_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "NaxiHGDP"),] # 6
# tu_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "TuHGDP"),] # 8
# bantu_kenya_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "BantuKenyaHGDP"),] # 10
mozabite_hgdp <- meta[(meta$paper == "hgdp") & (meta$country == "MozabiteHGDP"),] # 26

# arias 2017
arias_2017_colombia <- meta[(meta$paper == "arias_2017") & (meta$country == "Colombia"),] # 436

# avila_2019
avila_2019_brazil <- meta[(meta$paper == "avila_2019") & (meta$country == "Brazil"),] # 96

# barbieri 2012a
barbieri_2012a_burkinafaso <- meta[(meta$paper == "barbieri_2012a") & (meta$country == "Burkina Faso"),] # 291

# barbieri 2013a
barbieri_2013a_zambia <- meta[(meta$paper == "barbieri_2013a") & (meta$country == "Zambia"),] # 169

# barbieri 2013
barbieri_2013_botswana <- meta[(meta$paper == "barbieri_2013") & (meta$country == "Botswana"),] # 306
barbieri_2013_namibia <- meta[(meta$paper == "barbieri_2013") & (meta$country == "Namibia"),] # 141

# barbieri 2014a
barbieri_2014a_angola <- meta[(meta$paper == "barbieri_2014a") & (meta$country == "Angola"),] # 170
barbieri_2014a_zambia <- meta[(meta$paper == "barbieri_2014a") & (meta$country == "Zambia"),] # 446

# barbieri 2014b
barbieri_2014b_botswana <- meta[(meta$paper == "barbieri_2014b") & (meta$country == "Botswana"),] # 74
barbieri_2014b_namibia <- meta[(meta$paper == "barbieri_2014b") & (meta$country == "Namibia"),] # 117

# bodner 2015
bodner_2015_italy <- meta[(meta$paper == "bodner_2015") & (meta$country == "Italy"),] # 216

# brandini 2018
brandini_2018_ecuador <- meta[(meta$paper == "brandini_2018") & (meta$country == "Ecuador"),] # 208

# brucato 2018
brucato_2018_comoros <- meta[(meta$paper == "brucato_2018") & (meta$country == "Comoros"),] # 48
brucato_2018_kenya <- meta[(meta$paper == "brucato_2018") & (meta$country == "Kenya"),] # 228

# Chan 2019
chan_2019_namibia <- meta[(meta$paper == "chan_2019") & (meta$country == "Namibia"),] # 77
chan_2019_southafrica <- meta[(meta$paper == "chan_2019") & (meta$country == "South Africa"),] # 121

# Colombo 2025
colombo_2025_algeria <- meta[(meta$paper == "colombo_2025") & (meta$country == "Algeria"),] # 49
colombo_2025_cameroon <- meta[(meta$paper == "colombo_2025") & (meta$country == "Cameroon"),] # 29
colombo_2025_chad <- meta[(meta$paper == "colombo_2025") & (meta$country == "Chad"),] # 65
colombo_2025_libya <- meta[(meta$paper == "colombo_2025") & (meta$country == "Libya"),] # 60
colombo_2025_morocco <- meta[(meta$paper == "colombo_2025") & (meta$country == "Morocco"),] # 95
colombo_2025_tunisia <- meta[(meta$paper == "colombo_2025") & (meta$country == "Tunisia"),] # 34

# Fleskes 2023
fleskes_2023_anson <- meta[(meta$paper == "fleskes_2023"),] # 22

# Garcia 2021
garcia_2021_argentina <- meta[(meta$paper == "garcia_2021"),] # 114

# Garcia Olivares 2022
garciaolivares_2022_canaryislands <- meta[(meta$paper == "garciaolivares_2022"),] # 896

# Harney 2023
harney_2023_catoctin <- meta[(meta$paper == "harney_2023"),] # 22

# Huber 2025
huber_2025_peru <- meta[(meta$paper == "huber_2025"),] # 143

# Just 2008 african & hispanic american
just_2008_usa <- meta[(meta$paper == "just_2008"),] # 265
just_2008_aa <- meta[(meta$paper == "just_2008") & (meta$extra_info == "African American"),] # 140
just_2008_h <- meta[(meta$paper == "just_2008") & (meta$extra_info == "Hispanic"),] # 125

# Just 2015 african & hispanic american
just_2015_usa <- meta[(meta$paper == "just_2015"),] # 588
just_2015_aa <- meta[(meta$paper == "just_2015") & (meta$extra_info == "African American"),] # 170
just_2015_h <- meta[(meta$paper == "just_2015") & (meta$extra_info == "Hispanic"),] # 155
just_2015_c <- meta[(meta$paper == "just_2015") & (meta$extra_info == "caucasian"),] # 263

# kumar 2011: Mexican American from San Antonio Family Heart Study
kumar_2011_usa <- meta[(meta$paper == "kumar_2011"),] # 215

# mccrow 2016
mccrow_2016_southafrica <- meta[(meta$paper == "mccrow_2016"),] # 87

# olivera 2018
oliveira_2018_angola <- meta[(meta$paper == "oliveira_2018"),] # 295

# olivieri 2017
olivieri_2017_italy <- meta[(meta$paper == "olivieri_2017"),] # 2157

# pialq
pialq_peru <- meta[(meta$paper == "pialq"),] # 58

# pierron 2017
pierron_2017_madagascar <- meta[(meta$paper == "pierron_2017"),] # 2849

# sandoval velasco 2023
sandovalvelasco_2023_sthelena <- meta[(meta$paper == "sandovalvelasco_2023"),] # 19

# silva et al
silva_2021_portugal <- meta[(meta$paper == "silva_2021") & (meta$country == "Portugal"),] #103
silva_2021_spain <- meta[(meta$paper == "silva_2021") & (meta$country == "Spain"),] #1023

# taylor et al
taylor_2020_usa <- meta[(meta$paper == "taylor_2020"),] # 1327

# tito un_pub
tito_unpub_peru <- meta[(meta$paper == "tito_unpub"),] # 192

# tribaldos 2021, may need to separate
tribaldos_2021_panama <- meta[(meta$paper == "tribaldos_2021"),] # 84

# 6. function to create subsets of fasta files by population
# function
subset_align_by_pop <- function(alignment, pop_meta) {
  samples <- pop_meta$Sample # get sample names
  matching_indices <- which(rownames(alignment) %in% samples) # get indices of alignment seqs
  pop_alignment <- alignment[matching_indices, ] # subset alignment
  return(pop_alignment)
}

# list of pops
pop_list <- list(
  # 1kg populations
  acb = acb,
  asw = asw,
  beb = beb,
  cdx = cdx,
  ceu = ceu,
  chb = chb,
  chs = chs,
  clm = clm,
  esn = esn,
  fin = fin,
  gbr = gbr,
  gih = gih,
  gwd = gwd,
  ibs = ibs,
  itu = itu,
  jpt = jpt,
  khv = khv,
  lwk = lwk,
  msl = msl,
  mxl = mxl,
  pel = pel,
  pjl = pjl,
  pur = pur,
  stu = stu,
  tsi = tsi,
  yri = yri,
  
  # other populations
  mozabite_hgdp = mozabite_hgdp,
  arias_2017_colombia = arias_2017_colombia,
  avila_2019_brazil = avila_2019_brazil,
  barbieri_2012a_burkinafaso = barbieri_2012a_burkinafaso,
  barbieri_2013a_zambia = barbieri_2013a_zambia,
  barbieri_2013_botswana = barbieri_2013_botswana,
  barbieri_2013_namibia = barbieri_2013_namibia,
  barbieri_2014a_angola = barbieri_2014a_angola,
  barbieri_2014a_zambia = barbieri_2014a_zambia,
  barbieri_2014b_botswana = barbieri_2014b_botswana,
  barbieri_2014b_namibia = barbieri_2014b_namibia,
  bodner_2015_italy = bodner_2015_italy,
  brandini_2018_ecuador = brandini_2018_ecuador,
  brucato_2018_comoros = brucato_2018_comoros,
  brucato_2018_kenya = brucato_2018_kenya,
  chan_2019_namibia = chan_2019_namibia,
  chan_2019_southafrica = chan_2019_southafrica,
  colombo_2025_algeria = colombo_2025_algeria,
  colombo_2025_chad = colombo_2025_chad,
  colombo_2025_libya = colombo_2025_libya,
  colombo_2025_tunisia = colombo_2025_tunisia,
  fleskes_2023_anson = fleskes_2023_anson,
  garcia_2021_argentina = garcia_2021_argentina,
  garciaolivares_2022_canaryislands = garciaolivares_2022_canaryislands,
  harney_2023_catoctin = harney_2023_catoctin,
  huber_2025_peru = huber_2025_peru,
  just_2008_usa = just_2008_usa,
  just_2008_aa = just_2008_aa,
  just_2008_h = just_2008_h,
  just_2015_usa = just_2015_usa,
  just_2015_aa = just_2015_aa,
  just_2015_h = just_2015_h,
  just_2015_c = just_2015_c,
  kumar_2011_usa = kumar_2011_usa,
  mccrow_2016_southafrica = mccrow_2016_southafrica,
  oliveira_2018_angola = oliveira_2018_angola,
  olivieri_2017_italy = olivieri_2017_italy,
  pialq_peru = pialq_peru,
  pierron_2017_madagascar = pierron_2017_madagascar,
  sandovalvelasco_2023_sthelena = sandovalvelasco_2023_sthelena,
  silva_2021_portugal = silva_2021_portugal,
  silva_2021_spain = silva_2021_spain,
  taylor_2020_usa = taylor_2020_usa,
  tito_unpub_peru = tito_unpub_peru,
  tribaldos_2021_panama = tribaldos_2021_panama
)

# apply function to pop list
pop_alignments <- lapply(pop_list, function(pop_meta) {
  subset_align_by_pop(align_dedup, pop_meta)
})

# check sizes
pop_sizes <- lapply(pop_alignments, nrow)

# 7. get number of segregating sites per population with ape 
# seg.sites(pop, strict = FALSE, trailingGapsAsN = TRUE)
seg_sites <- lapply(pop_alignments, function(x) {
  seg.sites(x, strict = TRUE, trailingGapsAsN = TRUE)
})

seg_sites_count <- sapply(seg_sites, length)

# 8. get number of unique haplotypes per population with pegas
# haplotype(x, labels = NULL, strict = FALSE, trailingGapsAsN = TRUE, ...)
haps <- lapply(pop_alignments, function(x) {
  pegas::haplotype(x, strict = TRUE, trailingGapsAsN = TRUE)
})

# Apply to all populations
hap_stats <- lapply(names(pop_alignments), function(pop_name) {
  x <- pop_alignments[[pop_name]]
  h <- pegas::haplotype(x, strict = TRUE, trailingGapsAsN = TRUE)
  
  list(
    haplotypes = h,
    n_haplotypes = nrow(h)
  )
})
names(hap_stats) <- names(pop_alignments)

# extract n & freqs of haplotypes
n_haplotypes <- sapply(hap_stats, function(x) x$n_haplotypes)

# build summary table
hap_summary_table <- data.frame(
  population = names(hap_stats),
  n_sequences = sapply(pop_alignments, nrow),
  n_haplotypes = sapply(hap_stats, function(x) x$n_haplotypes)
)

write.table(hap_summary_table, "all_hap_count.csv", row.names = FALSE, sep = "\t", quote = FALSE)

# 9. get haplotype diversity for each population with pegas
# hap.div(x, variance = FALSE, method = "Nei", ...)
hap_diversity <- lapply(pop_alignments, function(x) {
  hap.div(x, variance = TRUE)
})

# 10. get nucleotide diversity with variance w pairwise deletion with ape
#nuc_diversity <- lapply(pop_alignments, function(x) {
#  nuc.div(x, variance = TRUE, pairwise.deletion = TRUE)
#})

nuc_diversity <- list()
for(pop_name in names(pop_alignments)) {
  cat("Processing:", pop_name, "- n =", nrow(pop_alignments[[pop_name]]), "\n")
  x <- pop_alignments[[pop_name]]
  
  nuc_diversity[[pop_name]] <- tryCatch({
    result <- pegas::nuc.div(x, variance = TRUE, pairwise.deletion = FALSE)
    cat("  Done! Pi =", result[1], "\n")
    result
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(c(NA, NA))
  })
  
  # Save progress after each population
  saveRDS(nuc_diversity, "pi_progress.rds")
}

# 11. tajima's d for each population
#taj_d <- lapply(pop_alignments, function(x) {
#  tajima.test(x)
#})

taj_d <- list()
for(pop_name in names(pop_alignments)) {
  cat("Processing:", pop_name, "- n =", nrow(pop_alignments[[pop_name]]), "\n")
  x <- pop_alignments[[pop_name]]
  
  taj_d[[pop_name]] <- tryCatch({
    result <- pegas::tajima.test(x)
    cat("  Done! D =", result$D, "\n")
    result
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(list(D = NA, Pval.normal = NA, Pval.beta = NA))
  })
  
  # Save progress after each population
  saveRDS(taj_d, "taj_progress.rds")
}

# 12. create a summary table of all results
summary_stats <- data.frame(
  population = names(pop_alignments),
  n_sequences = sapply(pop_alignments, nrow),
  segregating_sites = seg_sites_count,  # Now this is 80 values, not 95,444
  n_haplotypes = sapply(hap_stats, function(x) x$n_haplotypes),
  haplotype_diversity = sapply(hap_diversity, function(x) x[1]),
  hap_div_variance = sapply(hap_diversity, function(x) x[2]), 
  nucleotide_diversity = sapply(nuc_diversity, function(x) x[1]),
  nuc_div_variance = sapply(nuc_diversity, function(x) x[2]),
  tajimas_d = sapply(taj_d, function(x) x$D),
  tajimas_d_pval_normal = sapply(taj_d, function(x) x$Pval.normal),
  tajimas_d_pval_beta = sapply(taj_d, function(x) x$Pval.beta)
)

write.table(summary_stats, "intrapop_summ.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# 13. plotting
# pi vs sample size
pi_n <- ggplot(summary_stats,
               aes(x = n_sequences, y = nucleotide_diversity)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    x = "Sample size (n)",
    y = expression("Nucleotide diversity (" * pi * ")")
  ) +
  theme_classic()

# plot of pi per population with a bar plot
pi_bar <- ggplot(summary_stats, aes(x = reorder(population,nucleotide_diversity), y = nucleotide_diversity)) +
  geom_col() +
  labs(
    x = "Population",
    y = expression("Nucleotide diversity (" * pi * ")")
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# bixplot
ggplot(summary_stats, aes(x = factor(1), y = nucleotide_diversity)) +
  geom_boxplot(width = 0.4, fill = "white") +
  geom_point(data = summary_stats[summary_stats$population == "pialq_peru", ],
             aes(x = factor(1), y = nucleotide_diversity), 
             color = "black", size = 2.5) +  # Highlight pialq
  theme_bw() + 
  labs(y = expression("Nucleotide diversity (" * pi * ")"))

# pi violin plot
pi_violin <- ggplot(summary_stats, aes(x = factor(1), y = nucleotide_diversity)) +
  geom_violin(fill = "lightsteelblue3", alpha = 0.8) +
  geom_point(data = summary_stats[summary_stats$population == "pialq_peru", ],
             aes(x = factor(1), y = nucleotide_diversity), 
             color = "black", size = 2.5) +  # Highlight pialq
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  labs(y = expression("Nucleotide diversity (" * pi * ")"))

# plot of h vs sample size
h_n <- ggplot(summary_stats,
              aes(x = n_sequences, y = haplotype_diversity)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    x = "Sample size (n)",
    y = expression("Haplotype diversity (Hd)")
  ) +
  theme_classic()

# h barplot
h_bar <- ggplot(summary_stats, aes(x = reorder(population, haplotype_diversity), y = haplotype_diversity)) +
  geom_col() +
  labs(
    x = "Population",
    y = expression("Haplotype diversity (Hd)")
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# plot of h per population with a boxplot
ggplot(summary_stats, aes(x = factor(1), y = haplotype_diversity)) +
  geom_boxplot(width = 0.4, fill = "white") +
  geom_point(data = summary_stats[summary_stats$population == "pialq_peru", ],
             aes(x = factor(1), y = haplotype_diversity), 
             color = "black", size = 2.5) +  # Highlight pialq
  theme_bw() + 
  labs(y = expression("Haplotype diversity (Hd)"))

# H violin plot
h_violin <- ggplot(summary_stats, aes(x = factor(1), y = haplotype_diversity)) +
  geom_violin(fill = "lightsteelblue3", alpha = 0.8) +
  geom_point(data = summary_stats[summary_stats$population == "pialq_peru", ],
             aes(x = factor(1), y = haplotype_diversity), 
             color = "black", size = 2.5) +  # Highlight pialq
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  labs(y = expression("Haplotype diversity (Hd)"))

# 14. save to pdfs
pdf("intrapop_plots.pdf", width = 8, height = 11)
pi_n
pi_bar
pi_violin
h_n
h_bar
h_violin
dev.off()

# 15. combine all sequences used in calculations
all_pop_sequences <- do.call(rbind, pop_alignments)
all_pop_seq <- all_pop_sequences[!duplicated(rownames(all_pop_sequences)), ]
nrow(all_pop_seq)

# 16. write combined sequences to FASTA
write.dna(all_pop_seq, 
          file = "div_final.fasta", 
          format = "fasta",
          nbcol = -1)

# 17. also save metadata for these sequences
all_pop_samples <- rownames(all_pop_sequences)
all_pop_meta <- meta[meta$Sample %in% all_pop_samples, ]
write.table(all_pop_meta, 
            "div_final.meta", 
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

# make text file that has population and file for use in 6h_interpop.r
# save population assigments
pop_assign_list <- lapply(names(pop_list), function(pop_name) {
  data.frame(
    Sample = pop_list[[pop_name]]$Sample,
    Population = pop_name,
    stringsAsFactors = FALSE
  )
})

pop_assign_df <- do.call(rbind, pop_assign_list)
rownames(pop_assign_df) <- NULL

write.table(
  pop_assign_df,
  "pop_assign.txt",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)