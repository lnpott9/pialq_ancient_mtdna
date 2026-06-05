# =============================================================================
# Title: fig_s3_map.r
# Author: Laura N. Pott
# Date: 2026-06-05
# Description: R script to create map showing under/overrepresentation by country in reference database of mtDNA sequences used for phylogenetic trees
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. load necessary libraries
library(ggplot2)
library(maps)
library(dplyr)
library(stringr)

setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/")

# 2. read into data
# getting all papers
tree <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/03_iqtree/iqtree.meta", sep = "\t", header = TRUE)
network <- read.csv("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/04_mjn/mjn.meta", sep = "\t", header = TRUE)

# 3. bind together & drop duplicates
all_meta <- bind_rows(tree, network)
uniq_meta <- unique(all_meta)

# 4. look at country labels & sort by frequency
sort(unique(uniq_meta$country_label))

uniq_meta %>%
  count(country_label) %>%
  arrange(desc(n))

# 5. check for na values
sum(is.na(uniq_meta$country)) # 1
uniq_meta %>% filter(is.na(country))

# 6. remove empty ones & ref
meta <- uniq_meta[!(uniq_meta$country %in% c("unknown", "")), ]
meta <- meta[!(meta$file %in% c("AM948965.fasta")), ]

# 7. remove parentheses
meta <- meta %>%
  mutate(country_clean = str_trim(gsub("\\s*\\([^\\)]*\\)", "", country)))

# 8. fix 3 code 1kg & weird ones
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


# 9. compare country names between map data and mine
# base world map
map_data <- map_data("world")
setdiff(meta$country_new, map_data$region)

# 10. make new dataframe for sample size by country
country_counts <- meta %>%
  group_by(country_new) %>%
  summarise(sample_size = n())

# 11. calculate each country's percent of total sample
country_counts <- country_counts %>%
  mutate(sample_perc = sample_size / sum(sample_size))

# 12. add populations of each country
# install.packages("wbstats")
library(wbstats)

pop <- wb_data("SP.POP.TOTL", start_date = 2022, end_date = 2022) %>%
  select(country, SP.POP.TOTL) %>%
  rename(country_new = country, population = SP.POP.TOTL) %>%
  mutate(country_new = case_when(
    country_new == "United Kingdom" ~ "UK",
    country_new == "United States" ~ "USA",
    country_new == "Gambia, The" ~ "Gambia",
    country_new == "Congo, Dem. Rep." ~ "Democratic Republic of the Congo",
    country_new == "Russian Federation" ~ "Russia",
    country_new == "Venezuela, RB" ~ "Venezuela",
    country_new == "Yemen, Rep." ~ "Yemen",
    country_new == "Puerto Rico (US)" ~ "Puerto Rico",
    country_new == "Somalia, Fed. Rep." ~ "Somalia",
    TRUE ~ country_new
  ))

country_counts <- country_counts %>%
  left_join(pop, by = "country_new")

sum(is.na(country_counts$population))
country_counts %>% filter(is.na(population))

# add manually
country_counts <- country_counts %>%
  mutate(population = case_when(
    country_new == "Canary Islands" ~ 2000000,
    country_new == "Saint Helena"   ~ 4439,
    country_new == "Saint Martin"   ~ 73777,
    country_new == "Western Sahara" ~ 600000,
    TRUE ~ population
  ))

# 13. calculate each country's percent of total sample
country_counts <- country_counts %>%
  mutate(pop_perc = population / sum(population))

# 14. calculate % diff
country_counts <- country_counts %>%
  mutate(diff = (sample_perc - pop_perc) / pop_perc)

# 15. binary coding
country_counts <- country_counts %>%
  mutate(diff_sign = ifelse(diff > 0, "overrepresented", "underrepresented"))

# 16. get world map
ggplot(map_data, aes(x = long, y = lat, group = group, fill = region)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "none")

# 17. merge my dataset w map data
merged_map_data <- map_data %>%
  left_join(country_counts, by = c("region" = "country_new"))

# 18. full map
map <- ggplot(merged_map_data, aes(x = long, y = lat, group = group, fill = diff_sign)) +
  geom_polygon(color = "black", linewidth = 0.1) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "bottom") +
  ggtitle("Percent change between database frequency vs country population")

# 19. save
pdf("representation.pdf", width = 11, height = 8)
map
dev.off()
