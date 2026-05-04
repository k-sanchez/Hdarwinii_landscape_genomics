##################
## Figure 4 map ##
##################

library(dplyr)
library(rnaturalearth)
library(sf)
library(ggplot2)
library(ggrepel)
library(stringr)
theme_set(theme_minimal(base_family = "CMU Sans Serif"))

data <- read.csv("~/proj/Hdar_postdoc1/data/coords_ddRAD.csv", header = TRUE) %>% filter(sp != "wil")
countries <- ne_countries(scale = "medium", returnclass = "sf")
prov_arg <- ne_states(country = "Argentina", returnclass = "sf")
bbox <- st_bbox(st_as_sf(data, coords = c("lon", "lat"), crs = 4326))
expanded_bbox <- bbox + c(-1, -1, 1, 1)

# province labels
provs <- c("Mendoza", "Neuquén", "Río Negro", "Chubut", "Santa Cruz")
prov_label_df <- prov_arg %>% filter(name %in% provs) %>% mutate(centroid = st_centroid(geometry)) %>% mutate(lon = st_coordinates(centroid)[, 1], lat = st_coordinates(centroid)[, 2])

# rivers
# rivers <- st_read("~/OneDrive/GIS/capas/IGN/lineas_de_aguas_continentales_perenne/lineas_de_aguas_continentales_perenneLine.shp")
# rivers_chu <- rivers %>% filter(str_detect(nam, regex("Chubut")), !str_detect(nam, regex("Arroyo"))) %>% st_crop(ymin = -46, ymax = -41, xmin = -65, xmax = -72)
# Patagonian steppe
# estepa <- st_read("~/OneDrive/GIS/capas/Roig_etal_2018_ZS_PatagonianSteppe/shape estepa por distritos.shp")
# Pa <- estepa[estepa$Eco.terr.0 == "Payunia", ] # Payunia

# map
p_map <- ggplot() +
  geom_sf(data = countries, fill = NA, linewidth = 1, color = "#636363") +
  geom_sf(data = prov_arg, fill = NA, linetype = "dashed", linewidth = 0.5, color = "#636363") +
#   geom_sf(data = rivers_chu, color = "#1f78b4", linewidth = 1) +
#   geom_sf(data = Pa, fill = NA, color = "black", linewidth = 1) +
  # plot groups
  geom_point(data = data, aes(x = lon, y = lat, fill = sp), shape = 21, size = 4, color = "#636363") +
  scale_fill_manual(values = c("#1b9e77","#d95f02","#7570b3"), guide = "none") +
  coord_sf(xlim = c(expanded_bbox["xmin"], expanded_bbox["xmax"]), ylim = c(expanded_bbox["ymin"], expanded_bbox["ymax"]), expand = FALSE) +
  scale_x_continuous(breaks = seq(floor(expanded_bbox["xmin"]), ceiling(expanded_bbox["xmax"]), by = 4)) +
  scale_y_continuous(breaks = seq(floor(expanded_bbox["ymin"]), ceiling(expanded_bbox["ymax"]), by = 4)) +
  # province labels
  geom_text_repel(data = prov_label_df, aes(x = lon, y = lat, label = postal), size = 5, family = "CMU Sans Serif") +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 11),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, linewidth = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.1)
    )
p_map # SVG 500 x 500
