#########################
# Code for Figure 5 ----
#########################

library(ggplot2)
library(raster)
library(cowplot)
library(viridis)
library(dplyr)
library(sf)
library(terra)
library(rnaturalearth)
library(stringr)
library(classInt)

# raster data
emberger <- raster("/Users/kevinsanchez/OneDrive/GIS/capas/ENVIREM/SAmerica_current_30arcsec_geotiff/current_30arcsec_embergerQ.tif")
cont <- raster("/Users/kevinsanchez/OneDrive/GIS/capas/ENVIREM/SAmerica_current_30arcsec_geotiff/current_30arcsec_continentality.tif")
therm <- raster("/Users/kevinsanchez/OneDrive/GIS/capas/ENVIREM/SAmerica_current_30arcsec_geotiff/current_30arcsec_thermicityIndex.tif")
# reproject rasters
crs_wgs84 <- CRS("+init=epsg:4326")
emberger <- projectRaster(emberger, crs = crs_wgs84)
cont <- projectRaster(cont, crs = crs_wgs84)
therm <- projectRaster(therm, crs = crs_wgs84)

data <- read.csv("~/proj/Hdar_postdoc1/data/coords_all.csv", header = TRUE)
data <- data[!data$sp == "Hwil", ]
bbox <- extent(min(data$lon), max(data$lon), min(data$lat), max(data$lat))
emberger_crop <- crop(emberger, bbox)
cont_crop <- crop(cont, bbox)
therm_crop <- crop(therm, bbox)

# political boundaries
countries <- ne_countries(scale = "medium", returnclass = "sf")
# convert rasters to dataframes for ggplot
raster_to_df <- function(r) { as.data.frame(rasterToPoints(r), xy = TRUE) %>% rename(value = 3) }
df_emberger <- raster_to_df(emberger_crop)
df_cont <- raster_to_df(cont_crop)
df_therm <- raster_to_df(therm_crop)



#### Maps ####
# PDF = 400 x 400

ci <- classIntervals(df_emberger$value, n = 6, style = "fisher", samp_prop = 1); breaks <- ci$brks
map_emberger <- ggplot() +
  geom_raster(data = df_emberger, aes(x = x, y = y, fill = value)) +
  geom_sf(data = countries, fill = NA, linewidth = 1) +
  geom_point(data = data, aes(x = lon, y = lat, shape = sp), size = 2, color = "black", stroke = 1) +
  scale_fill_stepsn(colors = rev(viridisLite::inferno(5)), breaks = breaks) +
  scale_x_continuous(breaks = seq(floor(bbox[1]), ceiling(bbox[2]), by = 2)) +
  scale_y_continuous(breaks = seq(floor(bbox[3]), ceiling(bbox[4]), by = 2)) +
  coord_sf(xlim = c(bbox[1], bbox[2]), ylim = c(bbox[3], bbox[4]), expand = FALSE) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(family = "CMU Sans Serif", colour = "black", size = 12),
    legend.text = element_text(family = "CMU Sans Serif", size = 12),
    plot.title = element_text(family = "CMU Sans Serif", size = 14)
  ) +
  guides(shape = "none") +
  ggtitle("Emberger's Q")

ci <- classIntervals(df_cont$value, n = 6, style = "fisher", samp_prop = 1); breaks <- ci$brks
map_cont <- ggplot() +
  geom_raster(data = df_cont, aes(x = x, y = y, fill = value)) +
  geom_sf(data = countries, fill = NA, linewidth = 1) +
  geom_point(data = data, aes(x = lon, y = lat, shape = sp), size = 2, color = "white", stroke = 1) +
  scale_fill_stepsn(colors = rev(viridisLite::viridis(5)), breaks = breaks) +
  scale_x_continuous(breaks = seq(floor(bbox[1]), ceiling(bbox[2]), by = 2)) +
  scale_y_continuous(breaks = seq(floor(bbox[3]), ceiling(bbox[4]), by = 2)) +
  coord_sf(xlim = c(bbox[1], bbox[2]), ylim = c(bbox[3], bbox[4]), expand = FALSE) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(family = "CMU Sans Serif", colour = "black", size = 12),
    legend.text = element_text(family = "CMU Sans Serif", size = 12),
    plot.title = element_text(family = "CMU Sans Serif", size = 14)
  ) +
  guides(shape = "none") +
  ggtitle("Continentality")

ci <- classIntervals(df_therm$value, n = 6, style = "fisher", samp_prop = 1); breaks <- ci$brks
map_therm <- ggplot() +
  geom_raster(data = df_therm, aes(x = x, y = y, fill = value)) +
  geom_sf(data = countries, fill = NA, linewidth = 1) +
  geom_point(data = data, aes(x = lon, y = lat, shape = sp), size = 2, color = "white", stroke = 1) +
  scale_shape_manual(values = 0:4) +
  scale_fill_stepsn(
    colors = rev(viridisLite::cividis(5)),
    breaks = 
  ) +
  scale_x_continuous(breaks = seq(floor(bbox[1]), ceiling(bbox[2]), by = 2)) +
  scale_y_continuous(breaks = seq(floor(bbox[3]), ceiling(bbox[4]), by = 2)) +
  coord_sf(xlim = c(bbox[1], bbox[2]), ylim = c(bbox[3], bbox[4]), expand = FALSE) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(family = "CMU Sans Serif", colour = "black", size = 12),
    legend.text = element_text(family = "CMU Sans Serif", size = 12),
    plot.title = element_text(family = "CMU Sans Serif", size = 14)
  ) +
  guides(shape = "none") +
  ggtitle("Compensated thermicity index")

#### SNP barplot ####

snp_data <- data.frame(
  Pair = c("CN–Pa", "CN–SO", "Pa-SO"),
  SNPs = c(25, 32, 21),
  Top_vars = c(
    "Emberger Q, compensated thermicity index, aridity index",
    "PET seasonality, PET wettest, continentality",
    "Continentality, PET driest quarter, climatic conture index"
    )
    )

snp_bar <- ggplot(snp_data, aes(y = reorder(Pair, SNPs), x = SNPs)) +
  geom_bar(stat = "identity", fill = "grey30", width = 0.6) +
  geom_text(aes(label = Top_vars), hjust = 0, size = 3.2, nudge_x = 1, family = "CMU Sans Serif") +
  labs(
    title = "Environmentally associated SNPs",
    x = "Number of candidate SNPs",
    y = "Population pair"
  ) +
  scale_x_continuous(breaks = seq(0, 60, by = 10)) +
  coord_cartesian(clip = "off", xlim = c(0, max(snp_data$SNPs) + 15)) +
theme_void() +
theme(
  plot.margin = margin(5.5, 0, 5.5, 5.5),
  axis.title.x = element_text(family = "CMU Sans Serif", colour = "black", size = 12),
  axis.title.y = element_text(family = "CMU Sans Serif", colour = "black", size = 12, angle = 90),
  axis.text = element_text(family = "CMU Sans Serif", colour = "black", size = 10),
  axis.line.x = element_line(colour = "black"),
  axis.line.y = element_blank(),
  axis.ticks.x = element_line(colour = "black"),
  axis.ticks.y = element_blank(),
  plot.title = element_text(family = "CMU Sans Serif", size = 14, hjust = 0)
  ) # PDF 600 x 180
