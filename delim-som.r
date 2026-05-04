################
# Delim-SOM ----
################

library(adegenet)
library(scales)
library(conStruct)
library(poppr)
library(kohonen)
library(lsr)
library(combinat)
library(viridis)
library(caret)
library(biomod2)
library(dplyr)
library(stringr)
library(plotly)
library(tidyr)
library(rnaturalearth)
library(sf)
library(terra)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
setwd("~/proj/Hdar_postdoc1/analisis/delim-som/")
source("~/OneDrive/code_data/species_delimitation/delim-SOM/R/kohonen_code.R")

#### Data

# Genetic data
vcf <- vcfR::read.vcfR("~/proj/Hdar_postdoc1/data/RADseq/stacks/denovo_map_whigroup/populations_eems/populations.snps.filt.vcf")
samples <- colnames(vcf@gt)[-1]
genotype <- read.structure("~/proj/Hdar_postdoc1/data/RADseq/stacks/denovo_map_whigroup/populations_eems/populations.snps.filt.stru",
                           n.ind = length(samples),
                           n.loc = nrow(vcf@fix),
                           onerowperind = FALSE,
                           col.lab = 1,
                           col.pop = 0,
                           col.others = 2,
                           row.marknames = 1,
                           NA.char = -9)
alleles <- makefreq(genotype)
colnames(alleles) <- paste0("allele", 1:dim(alleles)[2])
write.csv(alleles, "alleles.txt", row.names = samples, quote = FALSE)

# Geographic data
coords <- read.csv("~/proj/Hdar_postdoc1/data/coords_ddRAD.csv", header = TRUE)
# reorder coords to match genotype order
coords_ordered <- coords %>% filter(voucher %in% samples) %>% arrange(match(voucher, samples))
# verify order
identical(coords_ordered$voucher, samples) # TRUE
names(coords_ordered)[4:5] <- c("x", "y")
# Elevation
dem <- rast("~/OneDrive/GIS/capas/ENVIREM/elev_SAmerica_current_30arcsec_geotiff/current_30arcsec_tri.tif")
target_crs <- "+proj=longlat +datum=WGS84 +no_defs"
dat <- vect(coords_ordered, geom = c("x", "y"), crs = crs(target_crs))
elev <- terra::extract(dem, dat)
names(elev)[2] <- "z"
# integrate lon, lat, and elev
space <- data.frame(coords_ordered[, c("x", "y")], z = elev$z)
# normalize to 0-1
space <- space %>% mutate(across(c(x, y, z), ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))))
# write final spatial data
write.csv(space, "space.txt", row.names = samples, quote = FALSE)

# Climatic data
clim <- rast("~/proj/Hdar_postdoc1/analisis/lscape_genetics/algatr_K8/stacked_rast.tiff")
names(clim) <- sub(".*_", "", names(clim)) # rename
# reproject to same CRS
clim_proj <- project(clim, target_crs)
# check NAs
values <- terra::extract(clim_proj, dat)
sum(is.na(values)) # should return 0
# pairwise correlations
corr <- cor(values[, -1])
rmv_cor <- findCorrelation(corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact = TRUE)
sel <- values[, -rmv_cor]
# normalize
climate <- as.data.frame(lapply(sel, function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))))
write.csv(climate, "climate.txt", row.names = samples, quote = FALSE)


#### SOM-based species delimitation

# Load processed data

# allele frequencies
alleles <- as.matrix(read.csv("alleles.txt", header = TRUE, row.names = 1))
# spatial data
space <- as.matrix(read.csv("space.txt", header = TRUE, row.names = 1))
# climatic data
climate <- as.matrix(read.csv("climate.txt", header = TRUE, row.names = 1))

# size of grid
g <- round(sqrt(5*sqrt(length(rownames(alleles)))))
# grid of size sqrt(n)
som_grid <- kohonen::somgrid(xdim = g, ydim = g, topo = "hexagonal", neighbourhood.fct = "gaussian")
# number of replicates
n <- 100
# number of steps
m <- 100

#### Run

res <- Climate.SOM()
# save.image("res_ClimateSOM.RData")
# load
rm(list = ls())
load("res_ClimateSOM.RData")
# plots
labels <- match.labels(alleles)
plotLearning.Climate(res) # 500 x 400
plotLayers(res) # 400 x 400
plotK(res) # 400 x 400
plotModel(res) # 400 x 400
Climate.SOM.varImp(res) # 500 x 800


#### Map ####

library(rnaturalearth)
library(sf)
library(terra)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(grid)
library(gridExtra)
library(readr)

q <- data.frame(match.k(res, labels)) # get admixture coefficients
k <- ncol(q)
q$voucher <- rownames(q); rownames(q) <- NULL
data <- merge(coords_ordered, q, by = "voucher")

# basemap
countries <- ne_countries(scale = "medium", returnclass = "sf")
prov_arg <- ne_states(country = "Argentina", returnclass = "sf")
# create a bounding box to crop maps
bbox <- st_bbox(st_as_sf(data, coords = c("x", "y"), crs = 4326))
expanded_bbox <- bbox + c(-1, -1, 1, 1)
# plot basemap
basemap <- ggplot() +
  geom_sf(data = countries, fill = "#f0f0f0", linewidth = 1) +
  geom_sf(data = prov_arg, fill = NA, linetype = "dashed", linewidth = 0.5) +
  coord_sf(xlim = c(expanded_bbox["xmin"], expanded_bbox["xmax"]), ylim = c(expanded_bbox["ymin"], expanded_bbox["ymax"]), expand = FALSE) +
  scale_x_continuous(breaks = seq(floor(expanded_bbox["xmin"]), ceiling(expanded_bbox["xmax"]), by = 2)) +
  scale_y_continuous(breaks = seq(floor(expanded_bbox["ymin"]), ceiling(expanded_bbox["ymax"]), by = 3)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    panel.grid.minor = element_line(colour="grey90", linewidth = 0.5),
    panel.grid.major = element_line(colour="grey90", linewidth = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "right"
  )


# draw piecharts
pal <- RColorBrewer::brewer.pal(6, name = "Paired")
# reshape the data into long format for ggplot
admix_long <- pivot_longer(data, cols = X1:paste0("X", k), names_to = "K", values_to = "proportion")

# define function to create pie charts
pies <- function(admix_df, ind, cols){
  ind_data <- subset(admix_df, voucher == ind)
  ggplot(ind_data, aes(x = "", y = proportion, fill = K)) +
    geom_bar(stat = "identity", width = 5., color = "#525252", show.legend = FALSE) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
    }

pies_list <- list()
for (i in data$voucher) { pies_list[[i]] <- pies(admix_df = admix_long, ind = i, cols = pal) }
coord_list <- list()
for (i in data$voucher) { coord_list[[i]] <- c(subset(data, voucher == i)$x, subset(data, voucher == i)$y) }
# convert each pie chart to an annotation_custom grob
radius = 0.3
pies_ac <- list()
for (i in 1:nrow(data)) {
  pies_ac[[i]] <- annotation_custom(
    grob = ggplotGrob(pies_list[[i]]),
    xmin = coord_list[[i]][1] - radius,
    xmax = coord_list[[i]][1] + radius,
    ymin = coord_list[[i]][2] - radius,
    ymax = coord_list[[i]][2] + radius
    )
  }

# add pie charts to the basemap
pie_map <- basemap
for (i in 1:nrow(data)) {
  pie_map <- pie_map + annotation_custom(
    grob = ggplotGrob(pies_list[[i]]),
    xmin = coord_list[[i]][1] - radius,
    xmax = coord_list[[i]][1] + radius,
    ymin = coord_list[[i]][2] - radius,
    ymax = coord_list[[i]][2] + radius
    )
  }
pie_map # 500 x 500
