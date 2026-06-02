# =============================================================================
# Title: fig_1_map.R
# Author: Laura N. Pott
# Date: 2026-06-02
# Description: R script to create map of Peru
# Usage: run interactively in RStudio IDE
# =============================================================================

# 1. install & load libraries
install.packages("sf")
install.packages("geoperu")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("ggplot2")
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")
install.packages("ggspatial")
install.packages("cowplot")

library(sf)
library(geoperu)
library(rnaturalearth)
library(dplyr)
library(ggplot2)
library(rnaturalearthhires)
library(ggspatial)
library(cowplot)


# 2. get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# 3. center on Peru
# Source - https://stackoverflow.com/a
# Posted by Jindra Lacko, modified by community. See post 'Timeline' for change history
# Retrieved 2025-12-10, License - CC BY-SA 4.0

### projection string for globe centered on Peru
crs_string <- "+proj=ortho +lon_0=-75 +lat_0=-10"

### background for the globe - center buffered by earth radius
ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = crs_string)

# 4. country polygons, cut to visible hemisphere and reproject
world_globe <- world %>% 
  st_intersection(ocean %>% st_transform(4326)) %>%
  st_transform(crs = crs_string)

# peru highlighted on globe
fig_globe <- ggplot() +
  geom_sf(data = ocean, fill = "aliceblue", color = NA) +
  geom_sf(data = world_globe, fill = "gray90", color = "gray70", size = 0.15) +
  geom_sf(data = world_globe[world_globe$name == "Peru", ], 
          fill = "black", color = "black", size = 0.1) +
  theme_void() + 
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )


# 5. canete province data
canete <- geoperu::get_geo_peru(geography = "CAÑETE", 
                                level = "prov",
                                simplified = FALSE)

# 6. create point for La Quebrada (defined before use in district_zoom)
la_quebrada <- data.frame(
  name = "La Quebrada",
  lon = -76.39449422059920,
  lat = -13.037168770199518
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(st_crs(canete))

# 7. zoom to San Luis district
san_luis <- canete[canete$distrito == "SAN LUIS", ]

district_zoom <- ggplot() +
  geom_sf(data = san_luis, fill = "gray90", color = "#9D0208", size = 0.3) +
  geom_sf(data = la_quebrada, color = "black", size = 3, shape = 19) +
  coord_sf(xlim = st_bbox(san_luis)[c(1,3)], 
           ylim = st_bbox(san_luis)[c(2,4)]) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# 8. canete province map with La Quebrada point
fig_lq <- ggplot() +
  geom_sf(data = canete, fill = "#EDABB2", color = "#9D0208", size = 0.3) +
  geom_sf(data = canete[canete$distrito == "SAN LUIS", ], 
          fill = "#9D0208", color = "#9D0208", size = 0.3) +
  geom_sf(data = la_quebrada, color = "black", size = 3, shape = 20) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  annotation_scale(location = "bl") +
  geom_sf_text(data = la_quebrada, aes(label = name), 
               nudge_x = -0.23, nudge_y = -0.05, size = 2, fontface = "bold")

# 9. Peru country map with Canete highlighted
peru_data <- geoperu::get_geo_peru(level = "all", simplified = FALSE)
canete_highlight <- peru_data[peru_data$provincia == "CAÑETE", ]

peru_deps <- ne_states("Peru", returnclass = "sf")

lima_city <- data.frame(
  name = "Lima",
  lon = -77.0428,
  lat = -12.0464
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(st_crs(peru_data))

canete_bbox <- st_bbox(canete_highlight) %>%
  st_as_sfc()

fig_peru <- ggplot() +
  geom_sf(data = peru_deps, fill = "gray90", color = "black", size = 0.3) +
  geom_sf(data = canete_highlight, fill = "#9D0208", color = "#9D0208", size = 0.3) +
  geom_sf(data = lima_city, color = "black", size = 2, shape = 19) +
  geom_sf_text(data = lima_city, aes(label = name), 
               nudge_x = -1.3, nudge_y = 0.18, size = 3.5, fontface = "bold") +
  geom_sf(data = canete_bbox, fill = NA, color = "black", size = 0.5) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotation_scale(location = "bl")

# 10. combine into final figure with insets
ggdraw() +
  draw_plot(fig_peru, x = -0.2) +
  draw_plot(fig_globe, x = 0.48, y = 0.05, width = 0.3, height = 0.3) +
  draw_plot(fig_lq, x = 0.4, y = 0.5, width = 0.4, height = 0.4) +
  draw_line(
    x = c(0.245, 0.5),
    y = c(0.345, 0.9),
    color = "black", 
    size = 0.5
  ) +
  draw_line(
    x = c(0.27, 0.5),
    y = c(0.3, 0.5),
    color = "black", 
    size = 0.5
  )

# 11. exported as pdf from plots pane