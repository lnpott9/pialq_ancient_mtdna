# =============================================================================
# Title: fig_3.R
# Author: Laura N. Pott
# Date: 2026-06-02
# Description: R script to create figure of geographical affinity based on network & phylogenetic analyses
# Usage: run interactively in RStudio IDE
# =============================================================================

# Load libraries
library(ComplexHeatmap)
library(tidyr)
library(dplyr)
library(tibble)

setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp")

# 1. Load the data
raw_data <- read.table("tidy_geography.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

# 2. Standardize the countries
unique(raw_data$country)

data <- raw_data %>%
  mutate(country = case_when(
    country == "PEL" ~ "Peru",
    country == "ACB" ~ "Barbados",
    country == "ASW" ~ "United States",
    country == "LWK" ~ "Kenya",
    country == "ESN" ~ "Nigeria",
    country == "YRI" ~ "Nigeria",
    country == "GWD" ~ "The Gambia",
    country == "MandenkaHGDP" ~ "Senegal",
    country == "PUR" ~ "Puerto Rico",
    country == "MSL" ~ "Sierra Leone",
    country == "IBS" ~ "Spain",
    country == "YorubaHGDP" ~ "Nigeria",
    country == "MozabiteHGDP" ~ "Algeria",
    country == "St. Helena" ~ "St. Helena: 1800s",
    country == "St. Martin" ~ "St. Martin: 1600s",
    country == "Catoctin Furnace" ~ "Catoctin Furnace: 1700-1800s",
    TRUE ~ country
  ))

unique(data$country)

# 3. Define countries and their specific regions
region_lookup <- c(
  "United States" = "North America",
  "Barbados" = "Caribbean",
  "Puerto Rico" = "Caribbean",
  "Dominican Republic" = "Caribbean",
  "Peru" = "South America", 
  "Colombia" = "South America", 
  "Ecuador" = "South America",
  "Brazil" = "South America",
  "Paraguay" = "South America",
  "Spain" = "Southern Europe",
  "Portugal" = "Southern Europe",
  "Canary Islands" = "Northern Africa",
  "Algeria" = "Northern Africa",
  "Libya" = "Northern Africa",
  "Morocco" = "Northern Africa",
  "The Gambia" = "Western Africa",
  "Senegal" = "Western Africa",
  "Sierra Leone" = "Western Africa",
  "Burkina Faso" = "Western Africa",
  "Nigeria" = "Western Africa",
  "Cameroon" = "Central Africa",
  "Sao Tome & Principe" = "Central Africa",
  "Angola" = "Central Africa",
  "Zambia" = "Southern Africa",
  "South Africa" = "Southern Africa",
  "Namibia" = "Southern Africa",
  "Botswana" = "Southern Africa",
  "Mozambique" = "Southern Africa",
  "Somalia" = "Eastern Africa",
  "Kenya" = "Eastern Africa",
  "Comoros" = "Eastern Africa",
  "Tanzania" = "Eastern Africa",
  "Madagascar" = "Eastern Africa",
  "St. Helena: 1800s" = "historical diaspora",
  "St. Martin: 1600s" = "historical diaspora",
  "Catoctin Furnace: 1700-1800s" = "historical diaspora"
)

country_order <- names(region_lookup)

# 4. Collapse multiple types per individual-country pair
collapsed_data <- data %>%
  group_by(individual, country) %>%
  summarise(
    type = case_when(
      all(type == "network") ~ "network",
      all(type == "tree")    ~ "tree",
      TRUE ~ "both"   # mixed or already "both"
    ),
    .groups = "drop"
  )

# 5. Process into wide matrix (values = type string)
processed_data <- collapsed_data %>%
  mutate(country = factor(country, levels = country_order)) %>%
  arrange(country) %>%
  pivot_wider(names_from = country, values_from = type, values_fill = NA) %>%
  column_to_rownames("individual")

# 6. Reorder columns to match country_order (only keep countries present in data)
col_order <- intersect(country_order, colnames(processed_data))
processed_data <- processed_data[, col_order]

heatmap_matrix <- as.matrix(processed_data)

# 7. Keep only individuals present in the data, in a specified order based on haplogroup
individual_order <- c(
  "LQ28",
  "LQ61",
  "LQ14",
  "LQ05",
  "LQ49",
  "LQ06",
  "LQ23",
  "LQ60",
  "LQ62",
  "LQ01",
  "LQ46",
  "LQ58",
  "LQ18",
  "LQ44",
  "LQ25",
  "LQ27",
  "LQ08",
  "LQ15",
  "LQ41",
  "LQ45",
  "LQ48",
  "LQ54",
  "LQ39",
  "LQ22",
  "LQ29",
  "LQ56",
  "LQ59",
  "LQ63",
  "LQ02",
  "LQ12",
  "LQ16",
  "LQ17",
  "LQ19",
  "LQ31",
  "LQ55",
  "LQ03",
  "LQ30",
  "LQ52",
  "LQ11",
  "LQ38",
  "LQ09",
  "LQ24",
  "LQ26",
  "LQ33",
  "LQ64",
  "LQ40",
  "LQ10",
  "LQ13",
  "LQ50",
  "LQ57",
  "LQ32",
  "LQ35",
  "LQ36",
  "LQ42",
  "LQ43",
  "LQ20",
  "LQ21",
  "LQ34"
)


row_order <- intersect(individual_order, rownames(heatmap_matrix))
heatmap_matrix <- heatmap_matrix[row_order, ]

# 8. Create region vector aligned with matrix columns
column_groups <- region_lookup[colnames(heatmap_matrix)]

# 9. Define colors for regions
region_colors <- c(
  "North America" = "#4e148c",
  "Caribbean" = "#990066",
  "South America" = "#ff2e8c", 
  "Southern Europe" = "#006699",
  "Northern Africa" = "#3399cc", 
  "Western Africa" = "darkgreen",
  "Central Africa" = "#669900",
  "Southern Africa" = "#d9d22e", 
  "Eastern Africa" = "#EC9F05",
  "historical diaspora" = "#a50104"
)

# 10. Create column annotation
col_ann <- HeatmapAnnotation(
  Region = column_groups,
  col = list(Region = region_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# 11. Cell function: fill by region color, label by type
col_fun <- function(j, i, x, y, w, h, fill) {
  current_region <- column_groups[j]
  val <- heatmap_matrix[i, j]
  
  if (!is.na(val)) {
    grid.rect(x, y, w, h, gp = gpar(fill = region_colors[current_region], col = "grey80"))
    
    label <- switch(val,
                    "network" = "N",
                    "tree" = "T",
                    "both" = "B",
                    ""
    )
    grid.text(label, x, y, gp = gpar(fontsize = 7, col = "white", fontface = "bold"))
    
  } else {
    grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey80"))
  }
}

# 12. Extract haplogroup labels aligned to matrix rows
row_labels <- data %>%
  select(individual, haplogroup) %>%
  distinct() %>%
  mutate(haplogroup = substr(haplogroup, 1, 2)) %>%
  filter(individual %in% rownames(heatmap_matrix)) %>%
  deframe() %>%                    
  .[rownames(heatmap_matrix)]

# 13. Build type legend
type_legend <- Legend(
  labels = c("Network", "Tree", "Both"),
  title = "Type",
  legend_gp = gpar(fill = "grey40"),
  labels_gp = gpar(fontsize = 9),
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  graphics = list(
    function(x, y, w, h) {
      grid.rect(x, y, w, h, gp = gpar(fill = "grey40", col = NA))
      grid.text("N", x, y, gp = gpar(fontsize = 7, col = "white", fontface = "bold"))
    },
    function(x, y, w, h) {
      grid.rect(x, y, w, h, gp = gpar(fill = "grey40", col = NA))
      grid.text("T", x, y, gp = gpar(fontsize = 7, col = "white", fontface = "bold"))
    },
    function(x, y, w, h) {
      grid.rect(x, y, w, h, gp = gpar(fill = "grey40", col = NA))
      grid.text("B", x, y, gp = gpar(fontsize = 7, col = "white", fontface = "bold"))
    }
  )
)

# 14. Generate heatmap
row_ann_right <- rowAnnotation(
  Individual = anno_text(rownames(heatmap_matrix), 
                         gp = gpar(fontsize = 8))
)

ht <- Heatmap(heatmap_matrix, 
              name = "Match",
              rect_gp = gpar(type = "none"), 
              cell_fun = col_fun,            
              top_annotation = col_ann, 
              right_annotation = row_ann_right,
              
              # COLUMN SETTINGS
              column_split = factor(column_groups, levels = unique(region_lookup)),
              column_title_side = "top",
              column_title_rot = 45,
              column_title_gp = gpar(fontsize = 10, fontface = "bold", hjust = 0),
              column_gap = unit(1, "mm"),
              
              # ROW SETTINGS
              row_split = factor(row_labels, levels = unique(row_labels)), 
              row_gap = unit(1, "mm"),    
              row_title_rot = 0, 
              row_title_gp = gpar(fontsize = 10, fontface = "bold"),
              row_names_gp = gpar(fontsize = 8),
              
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              row_order = seq_len(nrow(heatmap_matrix)),
              show_heatmap_legend = FALSE,
              row_names_side = "left",
              column_names_gp = gpar(fontsize = 9),
              column_names_rot = 45
)

# 15. Draw heatmap with type legend
pdf("test.pdf", width = 14, height = 10)
draw(ht, annotation_legend_list = list(type_legend))
dev.off()