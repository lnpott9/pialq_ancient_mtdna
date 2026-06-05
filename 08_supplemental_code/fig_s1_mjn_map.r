# =============================================================================
# Title: fig_s1_mjn_map.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to create map showing geographic metadata for reference database of mtDNA sequences used for networks
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. load necessary libraries
library(ggplot2)
library(maps)
library(dplyr)
library(stringr)

setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn")

# 2. read into data
meta <- read.csv("mjn.meta", sep = "\t")

# 3. clean country names & standardize
sort(unique(meta$country))

# 4. remove empty ones
meta <- meta[!(meta$country %in% c("unknown", "")), ]

# 5. remove parentheses
meta <- meta %>%
  mutate(country_clean = str_trim(gsub("\\s*\\([^\\)]*\\)", "", country)))

# 6. fix 3 code 1kg & weird ones
meta <- meta %>%
  mutate(country_new = recode(country_clean,
                              "ACB" = "Barbados",
                              "ASW" = "USA",
                              "BEB" = "Bangladesh",
                              "CDX" = "China",
                              "CEU" = "USA",
                              "CHB" = "China",
                              "CHS" = "China",
                              "CLM" = "Colombia",
                              "ESN" = "Nigeria",
                              "FIN" = "Finland",
                              "GBR" = "UK",
                              "GIH" = "USA",
                              "GWD" = "Gambia",
                              "IBS" = "Spain",
                              "ITU" = "UK",
                              "JPT" = "Japan",
                              "KHV" = "Vietnam",
                              "LWK" = "Kenya",
                              "MSL" = "Sierra Leone",
                              "MXL" = "USA",
                              "PEL" = "Peru",
                              "PJL" = "Pakistan",
                              "PUR" = "Puerto Rico",
                              "STU" = "UK",
                              "TSI" = "Italy",
                              "YRI" = "Nigeria",
                              "Czechia" = "Czech Republic",
                              "United States" = "USA",
                              "Colombian" = "Colombia",
                              "Cambodian" = "Cambodia",
                              "PapuaNewGuinea" = "Papua New Guinea",     
                              "BotswanaOrNamibia" = "Namibia",
                              "SierraLeone" = "Sierra Leone",
                              "Abkhazia" = "Georgia",
                              "Congo" = "Democratic Republic of the Congo",              
                              "SouthAfrica" = "South Africa",         
                              "Kyrgyzystan" = "Kyrgyzstan",
                              "Korea" = "South Korea",
                              "Czechoslovia" = "Czech Republic",    
                              "Guinea Equatorial" = "Equatorial Guinea",
                              "Tenerife" = "Canary Islands",
                              "La Gomera" = "Canary Islands",
                              "Lanzarote" = "Canary Islands",
                              "Oromo" = "Ethiopia",
                              "São Tomé and Príncipe" = "Sao Tome and Principe",
                              "Sao Tome e Principe" = "Sao Tome and Principe",
                              "St. Helena" = "Saint Helena",
                              "St. Martin" = "Saint Martin",
                              "Scotland" = "UK",
                              "Great Britain" = "UK",
                              "United Kingdom" = "UK",
                              "England" = "UK",
                              "OrkneyIslands" = "UK",
                              "MakraniHGDP" = "Pakistan",        
                              "MbutiHGDP" = "Democratic Republic of the Congo",
                              "BiakaHGDP" = "Central African Republic",
                              "FrenchHGDP" = "France",       
                              "ColombianHGDP" = "Colombia",     
                              "SuruiHGDP" = "Brazil",
                              "MandenkaHGDP" = "Senegal",
                              "YorubaHGDP" = "Nigeria",    
                              "SanHGDP" = "Namibia",
                              "BantuSouthAfricaHGDP" = "South Africa",
                              "KaritianaHGDP" = "Brazil",
                              "BergamoItalianHGDP" = "Italy",
                              "DaurHGDP" = "China",
                              "BantuKenyaHGDP" = "Kenya",
                              "MozabiteHGDP" = "Algeria",
                              "BantuTswanaSGDP" = "Botswana",
                              "BantuHereroSGDP" = "Namibia",
                              "BantuKenyaSGDP" = "Kenya",
                              "BasqueSGDP" = "France",
                              "BiakaSGDP" = "Central African Republic",        
                              "BrahuiSGDP" = "Pakistan",  
                              "ChaneSGDP" = "Argentina",     
                              "ChukchiSGDP" = "Russia",    
                              "DinkaSGDP" = "Sudan",        
                              "EsanSGDP" = "Nigeria",        
                              "EskimoChaplinSGDP" = "Russia",   
                              "EskimoNaukanSGDP" = "Russia",  
                              "EskimoSirenikiSGDP" = "Russia", 
                              "EstonianSGDP" = "Estonia",      
                              "FrenchSGDP" = "France",         
                              "GambianSGDP" = "Gambia",       
                              "JordanianSGDP" = "Jordan",   
                              "JuhoanNorthSGDP" = "Namibia",     
                              "KaritianaSGDP" = "Brazil",      
                              "KhomaniSanSGDP" = "South Africa",    
                              "LuhyaSGDP" = "Kenya",         
                              "LuoSGDP" = "Kenya",             
                              "MakraniSGDP" = "Pakistan",      
                              "MansiSGDP" = "Russia",     
                              "MandenkaSGDP" = "Senegal",   
                              "MasaiSGDP" = "Kenya",        
                              "MayanSGDP" = "Mexico",       
                              "MbutiSGDP" = "Democratic Republic of the Congo",       
                              "MendeSGDP" = "Sierra Leone",        
                              "MixeSGDP" = "Mexico",       
                              "MixtecSGDP" = "Mexico",          
                              "MozabiteSGDP" = "Algeria",       
                              "PiapocoSGDP" = "Colombia",      
                              "PimaSGDP" = "Mexico",         
                              "QuechuaSGDP" = "Peru",       
                              "RussianSGDP" = "Russia",      
                              "SaharawiSGDP" = "Western Sahara",    
                              "SardinianSGDP" = "Italy",  
                              "SpanishSGDP" = "Spain",         
                              "SuruiSGDP" = "Brazil",           
                              "YorubaSGDP" = "Nigeria",      
                              "ZapotecSGDP" = "Mexico", 
                              .default = country_clean))  # keep any names not in codes


# 7. compare country names between map data and mine & standardize if needed
# base world map
map_data <- map_data("world")
setdiff(meta$country_new, map_data$region)

# 8. make new dataframe for sample size by country
country_counts <- meta %>%
  group_by(country_new) %>%
  summarise(sample_size = n())

# 9. merge my dataset w map data
merged_map_data <- map_data %>%
  left_join(country_counts, by = c("region" = "country_new"))

full <- ggplot(merged_map_data, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "black", linewidth = 0.1) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1)  +
  ggtitle("sample sizes of all used in networks")

# 10. function to plot
plot_haplotype_map <- function(meta,
                               map_data,
                               hap_col,
                               hap_value,
                               title_prefix = "mjn sample sizes by country") {
  
  hap_data <- meta %>%
    filter(grepl(hap_value, .data[[hap_col]]))
  
  hap_counts <- hap_data %>%
    group_by(country_new) %>%
    summarise(sample_size = n(), .groups = "drop")
  
  merged_data <- map_data %>%
    left_join(hap_counts, by = c("region" = "country_new"))
  
  plot <- ggplot(merged_data,
                 aes(x = long, y = lat, group = group, fill = sample_size)) +
    geom_polygon(color = "black", linewidth = 0.1) +
    coord_fixed(1.1) +
    theme_void() +
    theme(legend.position = "right") +
    scale_fill_viridis_c(option = "B", direction = -1) +
    ggtitle(paste(paste(hap_value, collapse = ", "), title_prefix))
  
  list(
    plot = plot,
    map_data = merged_data,
    meta_data = hap_data
  )
}


# 11. plot A2+(64)
a264 <- plot_haplotype_map(
  meta = meta,
  map_data = map_data,
  hap_col = "hap_7char",
  hap_value = "A2\\+\\(64\\)"
)

# 12. plot B2b
b2b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value  = "B2b"
)

# 13. plot C1b
c1b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value  = "C1b"
)

# 14. plot D1
d1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_2char",
  hap_value  = "D1"
)

# 15. plot H1bw
h1bw <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "H1bw"
)


# 16. plot L0
#l0a2a
l0a2a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L0a2a"
)

l0d1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L0d1a"
)

# 17. plot L1
l1c1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L1c1"
)

l1c24 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = c("L1c2", "L1c4")
)

# 18. plot L2
l2a1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1a"
)

l2a1d <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1d"
)

l2a1f <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1f"
)

l2a1q <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L2a1q"
)

l2b1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L2b1"
)

l2c <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_3char",
  hap_value  = "L2c"
)

# 19. plot L3
l3b1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L3b1a"
)

l3d1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3d1"
)

l3d3 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3d3"
)

l3d4 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3d4"
)

l3e1 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3e1"
)

l3e1a1a <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_7char",
  hap_value  = "L3e1a1a"
)

l3e2 <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_4char",
  hap_value  = "L3e2"
)

l3e3b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L3e3b"
)

l3f1b <- plot_haplotype_map(
  meta       = meta,
  map_data   = map_data,
  hap_col    = "hap_5char",
  hap_value  = "L3f1b"
)

# 20. save
library(gridExtra)

pdf("maps_mjn.pdf", width = 10, height = 8)
grid.arrange(full, ncol = 1)
#grid.arrange(a264$plot, b2b$plot, ncol = 1)
#grid.arrange(c1b$plot, d1$plot, ncol = 1)
#grid.arrange(h1bw$plot, l0a2a$plot, ncol = 1)
#grid.arrange(l0d1a$plot, l1c1$plot, ncol = 1)
#grid.arrange(l1c24$plot, l2a1a$plot, ncol = 1)
#grid.arrange(l2a1d$plot, l2a1f$plot, ncol = 1)
#grid.arrange(l2a1q$plot, l2b1$plot, ncol = 1)
#grid.arrange(l2c$plot, l3b1a$plot, ncol = 1)
#grid.arrange(l3d1$plot, l3d3$plot, ncol = 1)
#grid.arrange(l3d4$plot, l3e1a1a$plot, ncol = 1)
#grid.arrange(l3e2$plot, l3e3b$plot, ncol = 1)
#grid.arrange(l3f1b$plot, ncol = 1)
dev.off()

# 21. also save metadata summary as text file
meta_list <- list(
  a264 = a264$meta_data,
  b2b = b2b$meta_data,
  c1b = c1b$meta_data,
  d1 = d1$meta_data,
  h1bw = h1bw$meta_data,
  l0a2a = l0a2a$meta_data,
  l0d1a = l0d1a$meta_data,
  l1c1 = l1c1$meta_data,
  l1c24 = l1c24$meta_data,
  l2a1a = l2a1a$meta_data,
  l2a1d = l2a1d$meta_data,
  l2a1f = l2a1f$meta_data,
  l2a1q = l2a1q$meta_data,
  l2b1  = l2b1$meta_data,
  l2c = l2c$meta_data,
  l3b1a = l3b1a$meta_data,
  l3d1  = l3d1$meta_data,
  l3d3 = l3d3$meta_data,
  l3d4 = l3d4$meta_data,
  l3e1 = l3e1$meta_data,
  l3e1a1a = l3e1a1a$meta_data,
  l3e2 = l3e2$meta_data,
  l3e3b = l3e3b$meta_data,
  l3f1b = l3f1b$meta_data
)

for (name in names(meta_list)) {
  paper <- meta_list[[name]] %>%
    count(paper, sort = TRUE)
  country <- meta_list[[name]] %>%
    count(country, sort = TRUE)
  extra_info <- meta_list[[name]] %>%
    count(extra_info, sort = TRUE)
  hap_summ <- meta_list[[name]] %>%
    count(hap_3char, sort = TRUE)
  hap <- meta_list[[name]] %>%
    count(Haplogroup, sort = TRUE)
  filename <- paste0(name, "_summ.txt")
  
  capture.output({
    cat("paper\n")
    print(paper)
    cat("country\n")
    print(country)
    cat("extra info from metadata\n")
    print(extra_info)
    cat("hap_3char\n")
    print(hap_summ)
    cat("haplogroup\n")
    print(hap)
  }, file = filename)
}
