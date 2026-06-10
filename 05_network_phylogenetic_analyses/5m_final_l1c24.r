# =============================================================================
# Title: 5m_final_l1c24.r
# Author: Laura N. Pott
# Date: 2026-06-10
# Description: R script to create phylogenetic tree using ggtree for all L1c2'4 sequences
# Usage: run interactively in RStudio IDE
# =============================================================================

# https://arftrhmn.net/creating-a-publication-quality-phylogeny-using-ggtree/
# https://luisdva.github.io/rstats/reduce-tree/
# https://yulab-smu.top/treedata-book/chapter6.html

# 1. WORKING DIRECTORY
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/10_iqtree/")

# 2. LIBRARIES
library(ggtree)
library(treeio)
library(ape)
library(ggplot2)
library(dplyr)
library(tidytree)
library(viridis)

# 3. READ IN TREEFILE WITH TREEIO
tree <- read.iqtree("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/10_iqtree/iq_l1c24.treefile")
tree # check it was read in correctly

# 4. READ IN METADATA
meta <- read.table("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/03_iqtree/iqtree.meta",
                   sep = "\t",
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   fill = TRUE,
                   quote = "",
                   comment.char = "")

# 5. READ IN THE TIP/SAMPLE ORDER FROM LOG FILE
iqtree_lines <- readLines("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/10_iqtree/iq_l1c24.log")
iqtree_lines <- iqtree_lines[-(1:20)] # skip first 20 lines
total_line <- grep("^\\*\\*\\*\\*\\s+TOTAL", iqtree_lines) # everything after numbers
iqtree_lines <- iqtree_lines[1:(total_line - 1)]

# split and convert to data frame
tip_split <- strsplit(iqtree_lines, "\\s+")
tip_map <- do.call(rbind, lapply(tip_split, function(x) x[2:3]))
tip_map <- as.data.frame(tip_map, stringsAsFactors = FALSE)
colnames(tip_map) <- c("tip_id", "sample_name")

# 6. ROOT TREE WITH APE ON NEANDERTHAL MTDNA & CHECK IT'S ROOTED
rooted_tree <- root(tree, outgroup = "AM948965", resolve.root = TRUE)
is.rooted(rooted_tree)

# 7. REPLACE TIP LABELS WITH SAMPLE NAMES
# extract the current tip labels (numbers)
current_tips <- rooted_tree@phylo$tip.label

# convert tip ids to characters
tip_map$tip_id <- as.character(tip_map$tip_id)

# match each tree tip to its row in tip_map
idx <- match(current_tips, tip_map$tip_id)

# warn if any tips were not matched
if(any(is.na(idx))) {
  warning("Some tree tip labels were not found in tip_map")
}

# replace the labels
rooted_tree@phylo$tip.label <- tip_map$sample_name[idx]

# inspect the first 20 updated labels
rooted_tree@phylo$tip.label[1:20]

# 8. CONNECT TREE DATA WITH METADATA
# subset metadata to only rows where 'Sample' is in tree tips
fixed_tips <- rooted_tree@phylo$tip.label
meta_filt <- meta[meta$Sample %in% fixed_tips, ]

# base tree
full <- ggtree(rooted_tree, layout = "circular", linewidth = 0.1)

# join metadata Sample to tree label
full$data <- full$data %>%
  left_join(meta_filt, by = c("label" = "Sample"))

# 9. ADD CUSTOM PALETTE 
# maximally different color list from https://godsnotwheregodsnot.blogspot.com/2013/11/kmeans-color-quantization-seeding.html
my_colors <- c(
  "USA" = "#FFFF00",
  "Madagascar" = "#1CE6FF",
  "Peru" = "#FF34FF",
  "Colombia" = "#FF4A46",
  "Zambia" = "#008941",
  "Angola" = "#006FA6",
  "Burkina Faso" = "#A30059",
  "Ecuador" = "#FFDBE5",
  "Argentina" = "#7A4900",
  "Spain" = "#0000A6",
  "Italy" = "#63FFAC",
  "Kenya" = "#B79762",
  "Canary Islands" = "#004D43",
  "Namibia" = "#8FB0FF",
  "Brazil" = "#997D87",
  "Colombian in Medellín, Colombia (CLM)" = "#5A0007",
  "South Africa" = "#809693",
  "Panama" = "#1B4400",
  "Puerto Rican in Puerto Rico (PUR)" = "#4FC601",
  "Peruvian in Lima, Peru (PEL)" = "#3B5DFF",
  "Esan in Nigeria (ESN)" = "#4A3B53",
  "Yoruba in Ibadan, Nigeria (YRI)" = "#FF2F80",
  "Gambian in Western Division – Mandinka (GWD)" = "#61615A",
  "African Caribbean in Barbados (ACB)" = "#BA0900",
  "Botswana" = "#6B7900",
  "Mende in Sierra Leone (MSL)" = "#00C2A0",
  "Mexican Ancestry in Los Angeles CA USA (MXL)" = "#FFAA92",
  "Bolivia" = "#FF90C9",
  "Luhya in Webuye, Kenya (LWK)" = "#B903AA",
  "Chad" = "#D16100",
  "African Ancestry in SW USA (ASW)" = "#DDEFFF",
  "Nigeria" = "#000035",
  "Comoros" = "#7B4F4B",
  "Senegal" = "#A1C299",
  "Central African Republic" = "#300018",
  "Morocco" = "#0AA6D8",
  "Cameroon" = "#013349",
  "Paraguay" = "#00846F",
  "Chile" = "#372101",
  "Somalia" = "#FFB500",
  "Algeria" = "#C2FFED",
  "Sudan" = "#A079BF",
  "Finnish in Finland (FIN)" = "#CC0744",
  "Mozambique" = "#C0B9B2",
  "Iberian populations in Spain (IBS)" = "#C2FF99",
  "Toscani in Italia (TSI)" = "#001E09",
  "Tunisia" = "#00489C",
  "Venezuela" = "#6F0062",
  "Mexico" = "#0CBD66",
  "unknown" = "#EEC3FF",
  "Russia" = "#456D75",
  "Libya" = "#B77B68",
  "Portugal" = "#7A87A1",
  "St. Helena" = "#788D66",
  "Utah residents with Northern and Western European ancestry (CEU)" = "#885578",
  "British from England and Scotland (GBR)" = "#FAD09F",
  "Niger" = "#FF8A9A",
  "Tanzania" = "#D157A0",
  "Ethiopia" = "#BEC459",
  "Pakistan" = "#456648",
  "France" = "#0086ED",
  "São Tomé and Príncipe" = "#886F4C",
  "Dominican Republic" = "#34362D",
  "Ghana" = "#B4A8BD",
  "Saudi Arabia" = "#00A6AA",
  "St. Martin" = "#452C2C",
  "Botswana or Namibia" = "#636375",
  "Sierra Leone" = "#A3C8C9",
  "Canada" = "#FF913F",
  "China" = "#938A81",
  "Equatorial Guinea" = "#575329",
  "Estonia" = "#00FECF",
  "Gujarati Indians in Houston, Texas, USA (GIH)" = "#B05B6F",
  "Jordan" = "#8CD0FF",
  "Mali" = "#3B9700",
  "Mauritania" = "#04F757",
  "Uruguay" = "#C8A1A1",
  "Western Sahara" = "#1E6E00"
)

# 14. DRAFT FIGURE
# full tree, colored by UFBoot with clade highlighted
final_hap <- full +
  geom_tippoint(aes(color = hap_5char), size = 0.5) +
  geom_tiplab(size = 1.5) +
  theme(
    legend.position = "right",
    plot.margin = unit(c(2, 1, 2, 1), "cm"))

# clade tree, UFBoot labels, tips are countries (need to add margin space around the plot to make it look normal)
final_clade <- full +
  geom_tippoint(aes(color = country_label), size = 0.2) +
  scale_color_manual(values = my_colors, name = "Country") +
  geom_tiplab(size = 1.5) +
  theme(
    legend.position = "right",
    plot.margin = unit(c(2, 1, 2, 1), "cm")) +
  geom_text2(aes(label = UFboot, subset = !isTip), 
             size = 1, color = "red", hjust = 0.0, vjust = 0.0)

# save plots
library(gridExtra)

pdf("final_l1c24.pdf", width = 20, height = 10)
grid.arrange(final_hap, final_clade, ncol = 2)
dev.off()
