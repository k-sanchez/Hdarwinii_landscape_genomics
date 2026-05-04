################################################
################################################
## Figure 1: cytb tree + maps
################################################
################################################

library(ape)
library(phytools)
library(ggtree)
library(dplyr)
library(ggplot2)
library(patchwork)
library(elevatr)
library(rnaturalearth)
library(sf)
library(terra)
library(ggnewscale)
library(grid)
library(ggrepel)
theme_set(theme_minimal(base_family = "CMU Sans Serif"))

##########
# Tree ----
##########

cytb_tree <- read.tree("~/proj/Hdar_postdoc1/analisis/cytb/raxml/out_part.raxml.support")
ddrad_data <- read.csv("~/proj/Hdar_postdoc1/analisis/popgenome/data.csv")
outgroup <- grep("^(Hand|Hche|Hwhi)", cytb_tree$tip.label, value = TRUE)
vp_tips <- c("Hdar_11731", "Hdar_11732") # Valdés Peninsula samples
tips_keep <- union(intersect(cytb_tree$tip.label, ddrad_data$voucher), c(outgroup, vp_tips)) # keep tips with ddrad data
cytb_pruned <- keep.tip(cytb_tree, tips_keep)
cytb_rooted <- root(cytb_pruned, outgroup = outgroup, resolve.root = TRUE)
# load metadata and ASAP groups
asap_data <- read.csv("~/proj/Hdar_postdoc1/analisis/cytb/data.csv")
asap_data$ASAP <- factor(asap_data$ASAP)
tip_df <- data.frame(label = cytb_rooted$tip.label) |> left_join(asap_data[, c("voucher", "ASAP")], by = c("label" = "voucher"))
hwill_tips <- grep("^Hwil", cytb_rooted$tip.label, value = TRUE)

# palette
asap_palette <- c(
  "6" = "#6a3d9a",
  "8" = "#e31a1c",
  "1" = "#33a02c",
  "3" = "#b15928",
  "4" = "#b2df8a",
  "2" = "#1f78b4",
  "7" = "#fdbf6f",
  "5" = "#ffff99",
  "10" = "#cab2d6",
  "9" = "#ff7f00"
  )

# plot tree
p_tree <- ggtree(cytb_rooted, right = TRUE, color = "#636363", size = 0.75) %<+% tip_df +
  geom_point2(aes(subset = grepl("^[0-9]+$", label) & as.numeric(label) >= 80), shape = 21, size = 1.5, fill = "white", color = "#636363") +
  geom_tippoint(aes(fill = ASAP), shape = 21, size = 2, show.legend = FALSE, na.rm = TRUE, color = "#636363") +
  scale_fill_manual(values = asap_palette) +
  geom_treescale(x = 0, y = 95, offset = 1, width = 0.1, family = "CMU Sans Serif")
# explore X/Y positions of tips to highlight
# tree_data <- ggtree(cytb_rooted, right = TRUE)$data
# H. williamsii tips
# tree_data %>% filter(label %in% hwill_tips) %>% arrange(x) %>% pull(x)
# tree_data %>% filter(label %in% hwill_tips) %>% arrange(y) %>% pull(y)
# outgtoups: H che = outgroup[1:2], Hand = outgroup[3:4], Hwhi = outgroup[5:6]
# tree_data %>% filter(label %in% outgroup[1:2]) %>% arrange(x) %>% pull(x)
# tree_data %>% filter(label %in% outgroup[1:2]) %>% arrange(y) %>% pull(y)
# H. darwiniii complex (Hdar + Hwil)
# tree_data %>% filter(node == findMRCA(cytb_rooted, tips = grep("^(Hdar|Hwil)", tree_data$label, value = TRUE), type = "node"))
# H. dar. type locality sample
# tree_data %>% filter(label %in% "Hdar_9813")
# highlight clades
p_tree2 <- p_tree + 
  # H. will
  annotate("segment", x = 0.42, xend = 0.42, y = 31, yend = 35, linewidth = 1) +
  annotate("text", x = 0.425, y = 33, label = expression(italic("H. williamsii")), hjust = 0, family = "CMU Sans Serif", size = 3.5) +
  # H. che
  annotate("text", x = 0.135, y = 122, label = expression(italic("H. chelemini")), hjust = 0, family = "CMU Sans Serif", size = 3.5) +
  # H. and
  annotate("text", x = 0.13, y = 125, label = expression(italic("H. andicola")), hjust = 0, family = "CMU Sans Serif", size = 3.5) +
  # H. whi
  annotate("text", x = 0.06, y = 126.5, label = expression(italic("H. whitii")), hjust = 0, family = "CMU Sans Serif", size = 3.5) +
  # H. dar complex
  annotate("text", x = 0.12, y = 105, label = expression(italic("H. darwinii") ~ "complex"), hjust = 0, family = "CMU Sans Serif", size = 3.5) +
  # H. dar TL
  annotate("text", x = 0.4, y = 71, label = "*", hjust = 0, family = "CMU Sans Serif", size = 6)
p_tree2
  
##########
# Main map ----
##########

countries <- ne_countries(scale = "medium", returnclass = "sf")
prov_arg <- ne_states(country = "Argentina", returnclass = "sf")
bbox <- st_bbox(st_as_sf(asap_data %>% filter(!is.na(lon) & !is.na(lat)), coords = c("lon", "lat"), crs = 4326))
expanded_bbox <- bbox + c(-1, -1, 1, 1)
bbox_polygon <- st_as_sf(st_sfc(st_polygon(list(matrix(c(
  expanded_bbox["xmin"], expanded_bbox["ymin"],
  expanded_bbox["xmin"], expanded_bbox["ymax"],
  expanded_bbox["xmax"], expanded_bbox["ymax"],
  expanded_bbox["xmax"], expanded_bbox["ymin"],
  expanded_bbox["xmin"], expanded_bbox["ymin"]
  ), ncol = 2, byrow = TRUE))), crs = 4326))
elev <- get_elev_raster(locations = bbox_polygon, z = 5, prj = st_crs(4326)$proj4string)
elev_df <- as.data.frame(mask(rast(elev), vect(st_crop(countries, expanded_bbox))), xy = TRUE, na.rm = TRUE)
colnames(elev_df) <- c("lon", "lat", "elevation")
elev_breaks <- c(0, 200, 400, 800, max(elev_df$elevation, na.rm = TRUE))
elev_df$elev_cat <- cut(elev_df$elevation, breaks = elev_breaks, include.lowest = TRUE, labels = c("0–200", "200–400", "400–800", "800+"))
elev_df <- elev_df[!is.na(elev_df$elev_cat), ]
elev_cols <- gray.colors(n = 4, start = 0.5, end = 0.95)

# province labels
provs <- c("Mendoza", "Neuquén", "Río Negro", "Chubut", "Santa Cruz", "Buenos Aires")
prov_label_df <- prov_arg %>% filter(name %in% provs) %>% mutate(centroid = st_centroid(geometry)) %>% mutate(lon = st_coordinates(centroid)[, 1], lat = st_coordinates(centroid)[, 2])
prov_label_df <- prov_label_df %>% mutate(label = if_else(postal == "DF", "BA", postal))
p_map <- ggplot() +
  geom_raster(data = elev_df, aes(x = lon, y = lat, fill = elev_cat)) +
  scale_fill_manual(values = elev_cols, name = "Elevation (m)") +
  new_scale_fill() +
  geom_sf(data = countries, fill = NA, linewidth = 1, color = "#636363") +
  geom_sf(data = prov_arg, fill = NA, linetype = "dashed", linewidth = 0.5, color = "#636363") +
  # plot all ASAP groups except VP
  geom_point(data = subset(asap_data, ASAP != "9"), aes(x = lon, y = lat, fill = ASAP), shape = 21, size = 4, color = "#636363") +
  # plot VP on top
  geom_point(data = subset(asap_data, ASAP == "9"), aes(x = lon, y = lat, fill = ASAP), shape = 21, size = 4, color = "#636363") +
  scale_fill_manual(values = asap_palette, guide = "none") +
  coord_sf(xlim = c(expanded_bbox["xmin"], expanded_bbox["xmax"]), ylim = c(expanded_bbox["ymin"], expanded_bbox["ymax"]), expand = FALSE) +
  scale_x_continuous(breaks = seq(floor(expanded_bbox["xmin"]), ceiling(expanded_bbox["xmax"]), by = 4)) +
  scale_y_continuous(breaks = seq(floor(expanded_bbox["ymin"]), ceiling(expanded_bbox["ymax"]), by = 4)) +
  # type locality H. darwinii
  geom_segment(aes(x = -65.5, y = -46.3, xend = -65.90, yend = -47.75), arrow = arrow(length = unit(0.5, "cm")), linewidth = 1) +
  # type locality H. williamsii
  geom_segment(aes(x = -63.5, y = -38.5, xend = -61.99, yend = -38.14), arrow = arrow(length = unit(0.5, "cm")), linewidth = 1) +
  # province labels
  geom_text_repel(data = prov_label_df, aes(x = lon, y = lat, label = label), size = 5, family = "CMU Sans Serif") +
  # highlight distrib. of H. williamsii
  annotate("point", x = -61.99, y = -38.14, shape = 21, size = 12, fill = NA, stroke = 1.5) +
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
p_map

##########
# Argentina with SA inset ----
########## 

argentina <- ne_countries(scale = "medium", country = "Argentina", returnclass = "sf")
arg_box <- st_bbox(argentina)
# shapefiles
hdar_shp <- st_read("~/OneDrive/GIS/capas/IUCN_shapefiles/Homonota/darwinii/data_0.shp", quiet = TRUE)
hwil_shp <- st_read("~/OneDrive/GIS/capas/IUCN_shapefiles/Homonota/williamsii/data_0.shp", quiet = TRUE)
hand_shp <- st_read("~/OneDrive/GIS/capas/IUCN_shapefiles/Homonota/andicola/data_0.shp", quiet = TRUE)
hwhi_shp <- st_read("~/OneDrive/GIS/capas/IUCN_shapefiles/Homonota/whitii/data_0.shp", quiet = TRUE)
hche_shp <- st_read("~/OneDrive/GIS/capas/capas_KS/Homonota_chelemini/distribution_smooth.shp", quiet = TRUE)
# Argentina map
p_arg <- ggplot() +
  geom_sf(data = argentina, fill = NA, color = "#636363", linewidth = 0.8) +
  geom_sf(data = st_crop(hdar_shp, arg_box), aes(fill = "H. darwinii"), alpha = 0.6, color = NA) +
  geom_sf(data = st_crop(hwil_shp, arg_box), aes(fill = "H. williamsii"), alpha = 0.6, color = NA) +
  geom_sf(data = st_crop(hand_shp, arg_box), aes(fill = "H. andicola"), alpha = 0.6, color = NA) +
  geom_sf(data = st_crop(hwhi_shp, arg_box), aes(fill = "H. whitii"), alpha = 0.6, color = NA) +
  geom_sf(data = st_crop(hche_shp, arg_box), aes(fill = "H. chelemini"), alpha = 0.6, color = NA) +
  annotate("text", x = -64.02, y = -25.33, label = "?", size = 6, fontface = "bold", family = "CMU Sans Serif") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set2"), labels = c(
    expression(italic("H. darwinii")),
    expression(italic("H. williamsii")),
    expression(italic("H. andicola")),
    expression(italic("H. whitii")),
    expression(italic("H. chelemini"))),
    name = NULL
    ) +
  theme_void() +
  theme(legend.text = element_text(family = "CMU Sans Serif", size = 11), plot.background = element_rect(fill = NA, colour = NA))
# South America inset
sa <- st_union(ne_countries(scale = "small", continent = "South America", returnclass = "sf"))
arg_small <- ne_countries(scale = "small", country = "Argentina", returnclass = "sf")
p_sa <- ggplot() +
  geom_sf(data = sa, fill = NA, color = "#636363", linewidth = 0.8) +
  geom_sf(data = arg_small, fill = "#636363", color = NA) +
  theme_void() +
  theme(plot.background  = element_rect(fill = NA, colour = NA))

##########
# Combine ----
##########

main_plot <- p_tree2 + plot_spacer() + p_map + plot_layout(widths = c(1.4, 0.1, 1))
# insert arg and SA inside the first panel (tree)
main_plot +
  inset_element(p_arg, left = -2.7, bottom = -0.03, right = 0.2, top = 0.6, on_top = TRUE) +
  inset_element(p_sa, left = -2.7, bottom = -0.02, right = 0.2, top = 0.2, on_top = TRUE) +
  theme(plot.background  = element_rect(fill = NA, colour = NA))
# PDF = 800 x 550

# Edits:
# - move province labels
# - pop H. wil and H. whi labels out of group