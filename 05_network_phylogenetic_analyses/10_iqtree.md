# 10: ML phylogenetic tree with iqtree



in $SHARED/projects/PIALQ/2025_ancient_lp/10_iqtree:

````
module load miniforge

conda create --copy -p $SHARED/projects/PIALQ/2025_ancient_lp/10_iqtree/tree_env python numpy pandas iqtree -c conda-forge -c bioconda

source activate $SHARED/projects/PIALQ/2025_ancient_lp/10_iqtree/tree_env
````



run iqtree! based on this tutorial: https://iqtree.github.io/doc/Tutorial#input-data



```
iqtree -s A.fasta -m MFP -B 1000 -T AUTO
```

outputs:

+ .iqtree: the main report file that is self-readable. You should look at this file to see the computational results. It also contains a textual representation of the final tree
+ .treefile: the ML tree in NEWICK format, which can be visualized by any supported tree viewer programs like FigTree or iTOL
+ .log: log file of the entire run (also printed on the screen). To report bugs, please send this log file and the original alignment file to the authors
+ example.phy.model: log-likelihoods for all models tested. It serves as a checkpoint file to recover an interrupted model selection



http://www.iqtree.org/doc/iqtree-doc.pdf



(ran single-threaded to avoid numerical underflow issues)



**10k:**

```
#!/bin/bash -l
#SBATCH --job-name=10k
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --tmp=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pott0195@umn.edu
#SBATCH --array=0-27

module load conda

cd $SHARED/projects/PIALQ/2025_ancient_lp/10_iqtree/

source activate $SHARED/projects/PIALQ/2025_ancient_lp/10_iqtree/tree_env

# files
TREES=(
"iq_a2.fasta"
"iq_a264.fasta"
"iq_b2.fasta"
"iq_b2b.fasta"
"iq_c1.fasta"
"iq_c1b.fasta"
"iq_d1.fasta"
"iq_h1.fasta"
"iq_h1bw.fasta"
"iq_l0a2a.fasta"
"iq_l0d1a.fasta"
"iq_l1c1.fasta"
"iq_l1c24.fasta"
"iq_l2a1.fasta"
"iq_l2a1a.fasta"
"iq_l2a1f.fasta"
"iq_l2b1.fasta"
"iq_l2c.fasta"
"iq_l3b1a.fasta"
"iq_l3d1.fasta"
"iq_l3d3.fasta"
"iq_l3d4.fasta"
"iq_l3e1.fasta"
"iq_l3e1a1a.fasta"
"iq_l3e2.fasta"
"iq_l3e3.fasta"
"iq_l3e3b.fasta"
"iq_l3f1b.fasta"
)

base="/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/03_iqtree"

# select file for this task
FILE="${TREES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing $FILE"

# Infer rooted trees with outgroup of neanderthal mtDNA
PREFIX="${FILE%.fasta}"
iqtree -s "$base/$FILE" -m MFP -bb 10000 -alrt 1000 -T 1 -safe --prefix "$PREFIX" -o AM948965
```



ultrafsst bootstrap info: https://iqtree.github.io/doc/Frequently-Asked-Questions#how-do-i-interpret-ultrafast-bootstrap-ufboot-support-values



maps of all sample sizes are here:

```
# iqtree sample sizes

# load necessary libraries
library(ggplot2)
library(maps)
library(dplyr)
library(stringr)

setwd("/projects/standard/mnievesc/shared/projects/PIALQ/2025_ancient_lp/09_msa/03_iqtree")

# read into data
meta <- read.csv("iqtree.meta", sep = "\t")

# clean country names & standardize
sort(unique(meta$country))

# remove empty ones
meta <- meta[!(meta$country %in% c("unknown", "")), ]

# remove parentheses
meta <- meta %>%
  mutate(country_clean = str_trim(gsub("\\s*\\([^\\)]*\\)", "", country)))

# fix 3 code 1kg & weird ones
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

# compare between map data and mine
# base world map
map_data <- map_data("world")
setdiff(meta$country_new, map_data$region)

# grep("Congo", map_data$region, value = TRUE)

# make new dataframe for sample size by country
country_counts <- meta %>%
  group_by(country_new) %>%
  summarise(sample_size = n())

# get world map
ggplot(map_data, aes(x = long, y = lat, group = group, fill = region)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "none")

# merge my dataset w map data
merged_map_data <- map_data %>%
  left_join(country_counts, by = c("region" = "country_new"))

full <- ggplot(merged_map_data, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1)  +
  ggtitle("sample sizes of all used in trees")

# plot A2 countries
a2_data <- meta[meta$hap_2char == "A2", ] # subset a2 data

a2_counts <- a2_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())

merged_a2 <- map_data %>%
  left_join(a2_counts, by = c("region" = "country_new"))

a2_plot <- ggplot(merged_a2, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("A2 iqtree sample sizes by country")

# plot B2 countries
b2_data <- meta[meta$hap_2char == "B2", ]
b2_counts <- b2_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_b2 <- map_data %>%
  left_join(b2_counts, by = c("region" = "country_new"))
b2_plot <- ggplot(merged_b2, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("B2 iqtree sample sizes by country")

# plot C1 countries
c1_data <- meta[meta$hap_2char == "C1", ]
c1_counts <- c1_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_c1 <- map_data %>%
  left_join(c1_counts, by = c("region" = "country_new"))
c1_plot <- ggplot(merged_c1, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("C1 iqtree sample sizes by country")

# plot D1 countries
d1_data <- meta[meta$hap_2char == "D1", ]
d1_counts <- d1_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_d1 <- map_data %>%
  left_join(d1_counts, by = c("region" = "country_new"))
d1_plot <- ggplot(merged_d1, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", size = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("D1 iqtree sample sizes by country")

# plot H1 countries
h1_data <- meta[meta$hap_2char == "H1", ]
h1_counts <- h1_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_h1 <- map_data %>%
  left_join(h1_counts, by = c("region" = "country_new"))
h1_plot <- ggplot(merged_h1, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("H1 iqtree sample sizes by country")


# plot L0 countries
l0_data <- meta[meta$hap_2char == "L0", ]
l0_counts <- l0_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l0 <- map_data %>%
  left_join(l0_counts, by = c("region" = "country_new"))
l0_plot <- ggplot(merged_l0, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L0 iqtree sample sizes by country")

# plot L1 countries
l1_data <- meta[meta$hap_2char == "L1", ]
l1_counts <- l1_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l1 <- map_data %>%
  left_join(l1_counts, by = c("region" = "country_new"))
l1_plot <- ggplot(merged_l1, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L1 iqtree sample sizes by country")

# plot L2 countries
l2_data <- meta[meta$hap_2char == "L2", ]
l2_counts <- l2_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l2 <- map_data %>%
  left_join(l2_counts, by = c("region" = "country_new"))
l2_plot <- ggplot(merged_l2, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L2 iqtree sample sizes by country")

# plot L3 countries
l3_data <- meta[meta$hap_2char == "L3", ]
l3_counts <- l3_data %>% # get counts per country
  group_by(country_new) %>%
  summarise(sample_size = n())
merged_l3 <- map_data %>%
  left_join(l3_counts, by = c("region" = "country_new"))
l3_plot <- ggplot(merged_l3, aes(x = long, y = lat, group = group, fill = sample_size)) +
  geom_polygon(color = "white", linewidth = 0.2) +
  coord_fixed(1.1) +
  theme_void() +
  theme(legend.position = "right") +
  scale_fill_viridis_c(option = "C", direction = -1) +
  ggtitle("L3 iqtree sample sizes by country")

# save
library(gridExtra)

pdf("sample_maps_fine.pdf", width = 10, height = 8)
grid.arrange(full, ncol = 1)
grid.arrange(a2_plot, b2_plot, ncol = 1)
grid.arrange(c1_plot, d1_plot, ncol = 1)
grid.arrange(h1_plot, l0_plot, ncol = 1)
grid.arrange(l1_plot, l2_plot, ncol = 1)
grid.arrange(l3_plot, ncol = 1)
dev.off()

# also save metadata summary as text file
meta_list = list(a2_data, b2_data, c1_data, d1_data, h1_data, l0_data, l1_data, l2_data, l3_data)
names(meta_list) <- c("a2", "b2", "c1", "d1", "h1", "l0", "l1", "l2", "l3")

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
```

