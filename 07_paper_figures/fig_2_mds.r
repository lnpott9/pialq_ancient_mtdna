# =============================================================================
# Title: fig_2_mds.R
# Author: Laura N. Pott
# Date: 2026-06-02
# Description: R script to create MDS plot from pre-calculated points
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. set working directory
setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/04_inter")

# 2. load libraries, download if necessary
library(dplyr)
library(ggplot2)

# 3. plot mds using coordinates already calculated
# read in coordinates
mds_data <- read.table("all_mds_coordinates.txt", header = TRUE)

# MDS plot
ggplot(mds_data, aes(x = MDS1, y = MDS2, label = Population)) +
  geom_point(size = 1) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  theme_bw() +
  labs(title = "MDS Plot",
       x = "MDS1",
       y = "MDS2")

# avoid overlaps
library(ggrepel)
ggplot(mds_data, aes(x = MDS1, y = MDS2, label = Population)) +
  geom_point(size = 1, color = "black") +
  geom_text_repel(size = 1, max.overlaps=Inf) +
  theme_bw() +
  labs(x = "MDS1",
       y = "MDS2")

# plot colored by region
# create new column for regions
unique(mds_data$Population)

# Create a mapping of populations to regions
# based on UN GeoScheme: https://www.emiw.org/fileadmin/emiw/UserActivityDocs/Geograph.Representation/Geographic-Representation-Appendix_1.pdf
region_mapping <- c(
  "acb" = "present-day African Diaspora",
  "asw" = "present-day African Diaspora",
  "beb" = "Asia",
  "cdx" = "Asia",
  "ceu" = "Europe",
  "chb" = "Asia",
  "chs" = "Asia",
  "clm" = "Americas",
  "esn" = "Africa",
  "fin" = "Europe",
  "gbr" = "Europe",
  "gih" = "Asia",
  "gwd" = "Africa",
  "ibs" = "Europe",
  "itu" = "Asia",
  "jpt" = "Asia",
  "khv" = "Asia",
  "lwk" = "Africa",
  "msl" = "Africa",
  "mxl" = "Americas",
  "pel" = "Americas",
  "pjl" = "Asia",
  "pur" = "Americas",
  "stu" = "Asia",
  "tsi" = "Europe",
  "yri" = "Africa",
  "mozabite_hgdp" = "Africa",
  "arias_2017_colombia" = "Americas",
  "avila_2019_brazil" = "Americas",
  "barbieri_2012a_burkinafaso" = "Africa",
  "barbieri_2012b_zambia" = "Africa",
  "barbieri_2013_botswana" = "Africa",
  "barbieri_2013_namibia" = "Africa",
  "barbieri_2014a_zambia" = "Africa",
  "barbieri_2014a_angola" = "Africa",
  "barbieri_2014b_botswana" = "Africa",
  "barbieri_2014b_namibia" = "Africa",
  "bodner_2022_italy" = "Europe",
  "brandini_2018_ecuador" = "Americas",
  "brucato_2018_kenya" = "Africa",
  "brucato_2018_comoros" = "Africa",
  "chan_2019_southafrica" = "Africa",
  "chan_2019_namibia" = "Africa",
  "colombo_2025_algeria" = "Africa",
  "colombo_2025_chad" = "Africa",
  "colombo_2025_libya" = "Africa",
  "colombo_2025_tunisia" = "Africa",
  "fleskes_2023_anson" = "historical African Diaspora",
  "garcia_2021_argentina" = "Americas",
  "garciaolivares_2023_canaryislands" = "Europe",
  "harney_2023_catoctin" = "historical African Diaspora",
  "huber_2025_peru" = "Americas",
  "just_2008_usa" = "Americas",
  "just_2015_usa" = "Americas",
  "kumar_2011_usa" = "Americas",
  "mccrow_2016_southafrica" = "Africa",
  "oliveira_2018_angola" = "Africa",
  "olivieri_2017_italy" = "Europe",
  "pialq_peru" = "THIS STUDY: Hacienda La Quebrada",
  "pierron_2017_madagascar" = "Africa",
  "sandovalvelasco_2023_sthelena" = "historical African Diaspora",
  "silva_2021_spain" = "Europe",
  "silva_2021_portugal" = "Europe",
  "taylor_2020_usa" = "Americas",
  "tito_unpub_peru" = "Americas",
  "tribaldos_2021_panama" = "Americas"
)

# Add the region column to data
mds_data <- mds_data %>%
  mutate(region = region_mapping[Population])

# plot colored by region
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region)) +
  geom_point(size = 2, alpha = 1) +
  theme_bw() +
  labs(title = "MDS Plot by Region",
       x = "MDS1",
       y = "MDS2")

# custom colors & shapes
custom_colors <- c(
  "Africa" = "#006400",
  "Americas" = "#800080",
  "Asia" = "#F9629F",
  "Europe" = "#1E90FF",
  "present-day African Diaspora" = "#39FF14",
  "historical African Diaspora" = "black",
  "THIS STUDY: Hacienda La Quebrada" = "red"
)

custom_shapes <- c(
  "Africa" = 20,
  "Americas" = 20,
  "Asia" = 20,
  "Europe" = 20,
  "present-day African Diaspora" = 20,
  "historical African Diaspora" = 20,
  "THIS STUDY: Hacienda La Quebrada" = 4
)

# size points by sample size
size_mapping <- c(
  "acb" = 98,
  "asw" = 66,
  "beb" = 104,
  "cdx" = 107,
  "ceu" = 101,
  "chb" = 105,
  "chs" = 112,
  "clm" = 103,
  "esn" = 111,
  "fin" = 105,
  "gbr" = 106,
  "gih" = 109,
  "gwd" = 120,
  "ibs" = 108,
  "itu" = 112,
  "jpt" = 105,
  "khv" = 103,
  "lwk" = 112,
  "msl" = 98,
  "mxl" = 70,
  "pel" = 90,
  "pjl" = 108,
  "pur" = 107,
  "stu" = 111,
  "tsi" = 111,
  "yri" = 109,
  "mozabite_hgdp" = 26,
  "arias_2017_colombia" = 436,
  "avila_2019_brazil" = 96,
  "barbieri_2012a_burkinafaso" = 291,
  "barbieri_2012b_zambia" = 169,
  "barbieri_2013_botswana" = 306,
  "barbieri_2013_namibia" = 141,
  "barbieri_2014a_zambia" = 446,
  "barbieri_2014a_angola" = 170,
  "barbieri_2014b_botswana" = 74,
  "barbieri_2014b_namibia" = 117,
  "bodner_2022_italy" = 216,
  "brandini_2018_ecuador" = 208,
  "brucato_2018_kenya" = 228,
  "brucato_2018_comoros" = 48,
  "chan_2019_southafrica" = 121,
  "chan_2019_namibia" = 77,
  "colombo_2025_algeria" = 49,
  "colombo_2025_chad" = 65,
  "colombo_2025_libya" = 60,
  "colombo_2025_tunisia" = 34,
  "fleskes_2023_anson" = 22,
  "garcia_2021_argentina" = 114,
  "garciaolivares_2023_canaryislands" = 896,
  "harney_2023_catoctin" = 22,
  "huber_2025_peru" = 143,
  "just_2008_usa" = 265,
  "just_2015_usa" = 588,
  "kumar_2011_usa" = 215,
  "mccrow_2016_southafrica" = 87,
  "oliveira_2018_angola" = 295,
  "olivieri_2017_italy" = 2157,
  "pialq_peru" = 58,
  "pierron_2017_madagascar" = 2849,
  "sandovalvelasco_2023_sthelena" = 19,
  "silva_2021_spain" = 1023,
  "silva_2021_portugal" = 103,
  "taylor_2020_usa" = 1327,
  "tito_unpub_peru" = 189,
  "tribaldos_2021_panama" = 84
)

# add sample size
mds_data <- mds_data %>%
  mutate(n = size_mapping[Population])

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region, shape = region)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  theme_bw() +
  labs(x = "MDS1",
       y = "MDS2")

# replot
simple <- ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region, shape = region, size = n, label = Population)) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  theme_bw() +
  labs(x = "MDS1 (61.9% of variance)",
       y = "MDS2 (17.9% of variance)",
       color = "Continent", 
       shape = "Continent")


# w text labels --> don't like this, very crazy
final <- simple +
  geom_text_repel(
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.size = 0.2,
    segment.alpha = 0.7,
    box.padding = 0.5,
    point.padding = 0.5,
    force = 3,
    force_pull = 0.2,
    seed = 42
  )

# clean up
# fix legend order
region_order <- c(
  "Africa",
  "Americas", 
  "Asia",
  "Europe",
  "present-day African Diaspora",
  "historical African Diaspora",
  "THIS STUDY: Hacienda La Quebrada"
)

mds_data$region <- factor(mds_data$region, levels = region_order)

# label only a few reference pops
label_pops <- c("acb", "asw", "yri", "gwd", "msl", "esn", "lwk", "pel",
                "tito_unpub_peru", "huber_2025_peru",
                "barbieri_2014b_namibia", "colombo_2025_chad", "brucato_2018_comoros",
                "barbieri_2012a_burkinafaso", "pialq_peru",
                "fleskes_2023_anson", "harney_2023_catoctin",
                "sandovalvelasco_2023_sthelena")

label_map <- c(
  "acb" = "ACB",
  "asw" = "ASW",
  "yri" = "YRI",
  "gwd" = "GWD",
  "msl" = "MSL",
  "esn" = "ESN",
  "lwk" = "LWK",
  "pel" = "PEL",
  "tito_unpub_peru" = "Tito unpub Peru",
  "huber_2025_peru" = "Huber 2025\nPeru",
  "barbieri_2014b_namibia" = "Barbieri 2014\nNamibia",
  "colombo_2025_chad" = "Colombo 2025\nChad",
  "brucato_2018_comoros" = "Brucato 2018\nComoros",
  "barbieri_2012a_burkinafaso" = "Barbieri 2012\nBurkina Faso",
  "fleskes_2023_anson" = "Fleskes 2023\nAnson Street",
  "harney_2023_catoctin" = "Harney 2023\nCatoctin Furnace",
  "sandovalvelasco_2023_sthelena" = "Sandoval-Velasco\n2023 St. Helena",
  "pialq_peru" = "THIS STUDY"
)

mds_data <- mds_data %>%
  mutate(label_wrapped = label_map[Population])

final <- ggplot(mds_data, 
                aes(x = MDS1, y = MDS2, 
                    color = region, shape = region,
                    size = n, label = label_wrapped)) +
                geom_point(alpha = 0.9) +
                geom_text_repel(
                  data = subset(mds_data, Population %in% label_pops),
                  segment.color = "black",
                  size = 3,
                  direction = "both",
                  show.legend = FALSE,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.15,
                  segment.alpha = 1,
                  box.padding = 0.4,
                  point.padding = 0.01,
                  force = 10,
                  seed = 42
                ) +
                scale_color_manual(values = custom_colors, breaks = region_order) +
                scale_shape_manual(values = custom_shapes, breaks = region_order) +
                scale_size_continuous(range = c(1, 5), guide = "legend") +
                theme_bw() +
                theme(
                  legend.title = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 14),
                  legend.key.size = unit(0.4, "cm"),
                  axis.title = element_text(size = 16),
                  axis.text = element_text(size = 14),
                  panel.grid.minor = element_blank()
                ) +
                labs(
                  x = "MDS1 (61.9% of variance)",
                  y = "MDS2 (17.9% of variance)",
                  color = "Continent",
                  shape = "Continent",
                  size = "Paper Sample Size"
                )

# save final mds plot
pdf("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/12_diversity/mds_paper.pdf", width = 11, height = 8)
final
dev.off()