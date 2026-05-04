###############################
##############################
## Figure 2: population genomics
##############################
##############################

library(ape)
library(dplyr)
library(elevatr)
library(ggnewscale)
library(ggplot2)
library(ggrepel)
library(ggtree)
library(grDevices)
library(grid)
library(gridExtra)
library(LEA)
library(patchwork)
library(raster)
library(reemsplots2)
library(reshape)
library(rnaturalearth)
library(RColorBrewer)
library(sf)
library(terra)
library(tidyr)
library(tidyverse)
theme_set(theme_minimal(base_family = "CMU Sans Serif"))

#######################
# Tree + barplots ----
#######################

ddrad_data <- read.csv("~/proj/Hdar_postdoc1/analisis/popgenome/data.csv")
asap_data <- read.csv("~/proj/Hdar_postdoc1/analisis/cytb/data.csv")
# add missing sample (not included in ASAP analysis but is co-distributed with samples assigned to Santa Cruz + Chubut east group)
asap_data <- asap_data %>% add_row(voucher = "Hdar_12187", lat = -43.48325, lon = -67.743639, ASAP = 10)
asap_data$ASAP <- factor(asap_data$ASAP)

# nuclear ML tree
tree <- read.tree("~/proj/Hdar_postdoc1/analisis/phylogenetics/ML_whigroup_loci_unphased/out_part.raxml.support")
outgroup <- grep("Hwhi|Hand", tree$tip.label, value = TRUE)
root_node <- getMRCA(tree, outgroup)
tree <- root(tree, node = root_node, resolve.root = TRUE)
# keep only vouchers that match data
matches <- intersect(tree$tip.label, ddrad_data$voucher) # this removes JC55
tree <- keep.tip(tree, matches)
# create ggtree object
tree_gg <- ggtree(tree, color = "#636363", size = 0.75)
# merge ASAP data into tree object; this allows plotting tip data (ASAP groups) directly
tree_gg$data <- tree_gg$data %>% left_join(asap_data, by = c("label" = "voucher"))
# ASAP groups color palette
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
# expanded cluster palette to match population genetic structure results
cluster_palette <- asap_palette
names(cluster_palette)[names(cluster_palette) == "2"] <- "Pa1"
cluster_palette["Pa2"] <- "#a6cee3"
cluster_palette["pop9"] <- "#fb9a99"

# Plot
p_tree <- tree_gg +
  # node support
  geom_point2(aes(subset = grepl("^[0-9]+$", label) & as.numeric(label) >= 80), shape = 21, size = 2, fill = "white", color = "#636363") +
  # add dashed lines from tips to aligned point
  geom_segment(data = tree_gg$data %>% filter(isTip), aes(x = x, xend = max(tree_gg$data$x) + 2e-3, y = y, yend = y), linetype = "dashed", linewidth = 0.3, color = "#636363") +
  # aligned points with ASAP group fill
  geom_point(data = tree_gg$data %>% filter(isTip), aes(x = max(tree_gg$data$x) + 2e-3, y = y, fill = ASAP), shape = 21, size = 2.5, color = "#636363", show.legend = FALSE) +
  scale_fill_manual(values = asap_palette) +
  coord_flip() +
  scale_x_reverse()
p_tree

# LEA project
project <- load.snmfProject("~/proj/Hdar_postdoc1/analisis/popgenome/LEA/populations.snps.filt.snmfProject")
# barplot function
qbar <- function(K, group_names) {
  best_run <- which.min(cross.entropy(project, K = K))
  qfile <- sprintf("~/proj/Hdar_postdoc1/analisis/popgenome/LEA/populations.snps.filt.snmf/K%d/run%d/populations.snps.filt_r%d.%d.Q", K, best_run, best_run, K)
  qmatrix <- read.table(qfile, header = FALSE)
  qmatrix$voucher <- ddrad_data$voucher
  tip_order <- p_tree$data %>% filter(isTip) %>% arrange(desc(y)) %>% pull(label)
  # reverse to match top-to-bottom order of tips
  qmatrix <- qmatrix %>% mutate(voucher = factor(voucher, levels = rev(tip_order))) %>% arrange(voucher)
  colnames(qmatrix)[1:K] <- group_names
  q_long <- pivot_longer(qmatrix, cols = all_of(group_names), names_to = "group", values_to = "proportion")
  ggplot(q_long, aes(x = voucher, y = proportion, fill = group)) +
    geom_bar(stat = "identity", position = "stack", width = .7, color = "#636363", linewidth = 0.1) +
    scale_fill_manual(values = cluster_palette[group_names]) +
    theme_void() +
    theme(legend.position = "none")
    }

# cluster group names in correct order
group_names_k7 <- c("10", "8", "5", "1", "4", "6", "Pa1")
group_names_k8 <- c("6", "Pa1", "8", "1", "Pa2", "4", "10", "5")
group_names_k9 <- c("4", "pop9", "8", "5", "1", "Pa2", "6", "10", "Pa1")
# create barplots
p_k7 <- qbar(7, group_names_k7)
p_k8 <- qbar(8, group_names_k8)
p_k9 <- qbar(9, group_names_k9)
p_k9_labs <- p_k9 +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = "Pa1", hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0, bottom = -0.5, right = 0.05, top = 0) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = "Pa2", hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0, bottom = -0.5, right = 0.16, top = 0) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = "WRN", hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0, bottom = -0.5, right = 0.28, top = 0) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = "SNQ", hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0, bottom = -0.5, right = 0.42, top = 0) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = "MZ-NNQ", hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0, bottom = -0.5, right = 0.7, top = 0) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = "CRN-NEC", hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0, bottom = -0.5, right = 1.23, top = 0.03) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = "SEL", hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0, bottom = -0.5, right = 1.45, top = 0) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = "WCH", hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0, bottom = -0.5, right = 1.75, top = 0) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = "***", hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0, bottom = -0.5, right = 0.49, top = 0) # Hwil
p_k9_labs
# build full tree + bars
p_tree_tagged <- p_tree + labs(tag = "a") + theme(plot.tag = element_text(size = 18, face = "bold", family = "CMU Sans Serif"), plot.tag.position = c(0, 0.96)) 
tree_and_bars <- (p_tree_tagged / p_k7 / p_k8 / p_k9_labs) + plot_layout(heights = c(1, 0.2, 0.2, 0.2))
tree_and_bars



##########################################
# Map of pop structure ----
##########################################

# data for K = 8
traitsK8 <- read.table("~/proj/Hdar_postdoc1/analisis/popgenome/LEA/dar_wil_traits_K8.txt", header = TRUE)
locs <- read.csv("~/proj/Hdar_postdoc1/analisis/popgenome/data.csv")
dataK8 <- merge(locs, traitsK8, by.x = "voucher", by.y = "Traits")
# reshape the data into long format for ggplot
admix_long <- pivot_longer(dataK8, cols = K1:K8, names_to = "K", values_to = "proportion")

# basemap
countries <- ne_countries(scale = "medium", returnclass = "sf")
prov_arg <- ne_states(country = "Argentina", returnclass = "sf")
# these coordinates are based on the cytb map, I use the same bounding box to match both maps
bbox <- st_bbox(st_as_sf(asap_data, coords = c("lon", "lat"), crs = 4326))
expanded_bbox <- bbox + c(-1, -1, 1, 1)
bbox_polygon <- st_as_sf(st_sfc(st_polygon(list(matrix(c(
  expanded_bbox["xmin"], expanded_bbox["ymin"],
  expanded_bbox["xmin"], expanded_bbox["ymax"],
  expanded_bbox["xmax"], expanded_bbox["ymax"],
  expanded_bbox["xmax"], expanded_bbox["ymin"],
  expanded_bbox["xmin"], expanded_bbox["ymin"]
  ), ncol = 2, byrow = TRUE))), crs = 4326))
# elev <- get_elev_raster(locations = bbox_polygon, z = 5, prj = st_crs(4326)$proj4string)
# elev_df <- as.data.frame(mask(rast(elev), vect(st_crop(countries, expanded_bbox))), xy = TRUE, na.rm = TRUE)
# colnames(elev_df) <- c("lon", "lat", "elevation")
# elev_breaks <- c(0, 200, 400, 800, max(elev_df$elevation, na.rm = TRUE))
# elev_df$elev_cat <- cut(elev_df$elevation, breaks = elev_breaks, include.lowest = TRUE, labels = c("0–200", "200–400", "400–800", "800+"))
# elev_df <- elev_df[!is.na(elev_df$elev_cat), ]
# elev_cols <- gray.colors(n = 4, start = 0.5, end = 0.95)

# province labels
provs <- c("Mendoza", "Neuquén", "Río Negro", "Chubut", "Santa Cruz", "Buenos Aires")
prov_label_df <- prov_arg %>% filter(name %in% provs) %>% mutate(centroid = st_centroid(geometry)) %>% mutate(lon = st_coordinates(centroid)[, 1], lat = st_coordinates(centroid)[, 2])
prov_label_df <- prov_label_df %>% mutate(label = if_else(postal == "DF", "BA", postal))

# rivers
rivers <- st_read("~/OneDrive/GIS/capas/IGN/lineas_de_aguas_continentales_perenne/lineas_de_aguas_continentales_perenneLine.shp")
rivers_chu <- rivers %>% filter(str_detect(nam, regex("Chubut|Chico")), !str_detect(nam, regex("Arroyo"))) %>% st_crop(ymin = -46, ymax = -41, xmin = -65, xmax = -72)

# plot basemap
basemap <- ggplot() +
#  geom_raster(data = elev_df, aes(x = lon, y = lat, fill = elev_cat)) +
#  scale_fill_manual(values = elev_cols, name = "Elevation (m)") +
  geom_sf(data = countries, fill = "#f2f2f2", linewidth = 1, color = "#636363") +
  geom_sf(data = rivers_chu, color = "#1f78b4", linewidth = 1) +
  geom_sf(data = prov_arg, fill = NA, linetype = "dashed", linewidth = 0.5, color = "#636363") +
  coord_sf(xlim = c(expanded_bbox["xmin"], expanded_bbox["xmax"]), ylim = c(expanded_bbox["ymin"], expanded_bbox["ymax"]), expand = FALSE) +
  scale_x_continuous(breaks = seq(floor(expanded_bbox["xmin"]), ceiling(expanded_bbox["xmax"]), by = 4)) +
  scale_y_continuous(breaks = seq(floor(expanded_bbox["ymin"]), ceiling(expanded_bbox["ymax"]), by = 4)) +
  # province labels
  geom_text_repel(data = prov_label_df, aes(x = lon, y = lat, label = label), size = 3.5, family = "CMU Sans Serif") +
  # highlight distrib. of H. williamsii
  geom_segment(aes(x = -64, y = -39, xend = -62.5, yend = -38.14), arrow = arrow(length = unit(0.5, "cm")), linewidth = 1) +
  annotate("text", x = expanded_bbox["xmax"] - 1.9, y = expanded_bbox["ymax"] - 0.5, label = expression(bolditalic("K") * bold(" = 8")), family = "CMU Sans Serif", size = 5) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 9),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, linewidth = 0.5),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 9)
    )
basemap
# create pie charts
# define group names for K = 8
group_names_k8 <- c("6", "Pa1", "8", "1", "Pa2", "4", "10", "5")
# lookup vector: "K1" → "Cerro Alto", etc.
k_lookup <- setNames(group_names_k8, paste0("K", 1:8))
# convert K codes to actual group names for plotting
admix_long$K_name <- k_lookup[admix_long$K]
# function to create pie charts
pies <- function(admix_df, site, cols){
  site_data <- subset(admix_df, voucher == site)
  ggplot(site_data, aes(x = "", y = proportion, fill = K_name)) +
    geom_bar(stat = "identity", width = 5., color = "#636363", show.legend = FALSE) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
    }

# create pie charts for each sample
pies_list <- list()
for (i in unique(dataK8$voucher)) { pies_list[[i]] <- pies(admix_df = admix_long, site = i, cols = cluster_palette) }
# define coordinates for pie chart placement
radius = 0.35
coord_list <- list()
for (i in dataK8$voucher) { coord_list[[i]] <- c(subset(dataK8, voucher == i)$lon, subset(dataK8, voucher == i)$lat) }
# add pies to basemap
pie_map <- basemap
for (i in 1:nrow(dataK8)) {
  pie_map <- pie_map + annotation_custom(
    grob = ggplotGrob(pies_list[[i]]),
    xmin = coord_list[[i]][1] - radius,
    xmax = coord_list[[i]][1] + radius,
    ymin = coord_list[[i]][2] - radius,
    ymax = coord_list[[i]][2] + radius
    )
  }
pie_map


################
# EEMS ----
################

load("~/proj/Hdar_postdoc1/analisis/eems/plots.RData")
# eems_path <- "~/proj/Hdar_posdoc1/analisis/eems/"
# plots <- make_eems_plots(mcmcpath = paste0(eems_path, c("nDemes500_chain1/", "nDemes500_chain2/", "nDemes500_chain3/")), longlat = TRUE, add_grid = FALSE)
# save.image("~/proj/Hdar_posdoc1/analisis/eems/plots.RData")

# m rates
m <- plots$mrates02 +
  geom_sf(data = countries, fill = NA, linewidth = 1, color = "#636363", inherit.aes = FALSE) +
  geom_sf(data = prov_arg, fill = NA, linetype = "dashed", linewidth = 0.5, color = "#636363", inherit.aes = FALSE) +
  # geom_point(data = ddrad_data, aes(x = lon, y = lat), shape = 16, size = 1, show.legend = FALSE) +
  coord_sf(xlim = c(expanded_bbox["xmin"], expanded_bbox["xmax"]), ylim = c(expanded_bbox["ymin"], expanded_bbox["ymax"]), expand = FALSE) +
  # province labels
  geom_text_repel(data = prov_label_df, aes(x = lon, y = lat, label = label), size = 3.5, family = "CMU Sans Serif") +
  scale_x_continuous(breaks = seq(floor(expanded_bbox["xmin"]), ceiling(expanded_bbox["xmax"]), by = 4)) +
  scale_y_continuous(breaks = seq(floor(expanded_bbox["ymin"]), ceiling(expanded_bbox["ymax"]), by = 4)) +
  scale_fill_gradient2(name = expression(bold(log(bolditalic(m)))), low = "#fc8d59", mid = "white", high = "#91bfdb", midpoint = 0.5, breaks = c(0.1, 0.9), labels = c("-0.9", "0.9")) +
  theme(
    axis.text = element_text(size = 9),
    axis.text.y = element_blank(),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, linewidth = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.19)
    )
# q rates
q <- plots$qrates02 +
  geom_sf(data = countries, fill = NA, linewidth = 1, color = "#636363", inherit.aes = FALSE) +
  geom_sf(data = prov_arg, fill = NA, linetype = "dashed", linewidth = 0.5, color = "#636363", inherit.aes = FALSE) +
  # geom_point(data = ddrad_data, aes(x = lon, y = lat), shape = 16, size = 1, show.legend = FALSE) +
  # province labels
  geom_text_repel(data = prov_label_df, aes(x = lon, y = lat, label = label), size = 3.5, family = "CMU Sans Serif") +
  coord_sf(xlim = c(expanded_bbox["xmin"], expanded_bbox["xmax"]), ylim = c(expanded_bbox["ymin"], expanded_bbox["ymax"]), expand = FALSE) +
  scale_x_continuous(breaks = seq(floor(expanded_bbox["xmin"]), ceiling(expanded_bbox["xmax"]), by = 4)) +
  scale_y_continuous(breaks = seq(floor(expanded_bbox["ymin"]), ceiling(expanded_bbox["ymax"]), by = 4)) +
  scale_fill_gradient2(name = expression(bold(log(bolditalic(q)))), low = "#fc8d59", mid = "white", high = "#91bfdb", midpoint = 0.5, breaks = c(0.1, 0.9), labels = c("-0.9", "0.9")) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 9),
    axis.text.y = element_blank(),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, linewidth = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.19)
    )
  

##################
# Combine all ----
##################

# add tags inside each plot manually
pie_map_tag <- pie_map + labs(tag = "b") + theme(plot.tag = element_text(size = 18, face = "bold", family = "CMU Sans Serif"), plot.tag.position = c(0.00, 0.99))
m_tag <- m + labs(tag = "c") + theme(plot.tag = element_text(size = 18, face = "bold", family = "CMU Sans Serif"), plot.tag.position = c(-0.02, 0.99))
q_tag <- q + labs(tag = "d") + theme(plot.tag = element_text(size = 18, face = "bold", family = "CMU Sans Serif"), plot.tag.position = c(-0.02, 0.99))
# compose final layout
bottom_row <- pie_map_tag + plot_spacer() + m_tag + plot_spacer() + q_tag + plot_layout(widths = c(1, 0.05, 1, 0.05, 1))
final_plot_base <- tree_and_bars / plot_spacer() / bottom_row + plot_layout(heights = c(0.5, 0.1, 0.1, 0.1, 0.01, 1.5))
# add labels K = 7, 8, 9 as text grobs
final_plot <- final_plot_base +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = expression(italic("K =") ~ "7"), hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0.8, bottom = 0, right = 1.15, top = 2.45) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = expression(italic("K =") ~ "8"), hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0.8, bottom = 0, right = 1.15, top = 2.32) +
  inset_element(ggplot() + annotate("text", x = 0, y = 0, label = expression(italic("K =") ~ "9"), hjust = 0, size = 4, family = "CMU Sans Serif") + theme_void(), left = 0.8, bottom = 0, right = 1.15, top = 2.18)
final_plot # PDF 800 x 600

# Manual edits:
# - move group labels
# - move province labels
