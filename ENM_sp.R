#############################################
# ENM for CN, Pa, and SO lineages – BIOMOD2 #
# Each lineage projected to 3 time periods  #
#############################################

library(terra)
library(rnaturalearth)
library(caret)
library(biomod2)
library(dplyr)
library(sf)
setwd("~/proj/Hdar_postdoc1/analisis/ENM_sp/")

#### Data

# occurrence
occur_all <- read.csv("~/proj/Hdar_postdoc1/data/coords_all.csv")
lineages <- c("cn", "so", "Pa")
bbox <- ext(min(occur_all$lon) - 2, max(occur_all$lon) + 2, min(occur_all$lat) - 2, max(occur_all$lat) + 2)

# environmental data folders
env_folders <- list(
    curr = "~/OneDrive/GIS/capas/ENVIREM/SAmerica_current_2.5arcmin_geotiff",
    lgm = "~/OneDrive/GIS/capas/ENVIREM/SAmerica_lgm_ccsm4_2.5arcmin_geotiff",
    hol = "~/OneDrive/GIS/capas/ENVIREM/SAmerica_holo_ccsm4_2.5arcmin_geotiff"
    )

#### Loop over each lineage

for (sp_i in lineages) {
  # filter for this lineage
  occur <- occur_all %>% filter(sp == sp_i) %>% distinct(lon, lat, .keep_all = TRUE)
  # load current climatic variables and crop
  clim_curr_files <- list.files(env_folders$curr, pattern = "\\.tif$", full.names = TRUE)
  clim_curr <- rast(clim_curr_files)
  names(clim_curr) <- sub(".*_", "", names(clim_curr))
  # crop raster to the extent of distribution
  clim_curr_crop <- crop(clim_curr, bbox)
  clim_curr_crop <- project(clim_curr_crop, "EPSG:4326")
  pres_vect <- vect(occur, geom = c("lon", "lat"), crs = "EPSG:4326")
  # extract environmental values at presence locations
  env_vals <- terra::extract(clim_curr_crop, pres_vect, ID = FALSE)
  corr <- cor(env_vals)
  rmv_cor <- findCorrelation(corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact = TRUE)
  if (length(rmv_cor) > 0) { clim_sel <- clim_curr_crop[[-rmv_cor]] } else { clim_sel <- clim_curr_crop }
  # load LGM and Holocene rasters
  clim_lgm_files <- list.files(env_folders$lgm, pattern = "\\.tif$", full.names = TRUE)
  clim_lgm <- rast(clim_lgm_files)
  names(clim_lgm) <- sub(".*_", "", names(clim_lgm))
  clim_lgm_crop <- crop(clim_lgm, bbox)
  clim_lgm_crop <- clim_lgm_crop[[intersect(names(clim_sel), names(clim_lgm_crop))]]
  clim_lgm_crop <- project(clim_lgm_crop, "EPSG:4326")
  clim_hol_files <- list.files(env_folders$hol, pattern = "\\.tif$", full.names = TRUE)
  clim_hol <- rast(clim_hol_files)
  names(clim_hol) <- sub(".*_", "", names(clim_hol))
  clim_hol_crop <- crop(clim_hol, bbox)
  clim_hol_crop <- clim_hol_crop[[intersect(names(clim_sel), names(clim_hol_crop))]]
  clim_hol_crop <- project(clim_hol_crop, "EPSG:4326")
  # prepare dataset for BIOMOD2
  pres <- occur[, c("lon", "lat")]
  dir.create(sp_i, showWarnings = FALSE, recursive = TRUE)
  myBiomodData <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(pres)),
    expl.var = clim_sel,
    resp.xy  = pres,
    resp.name = sp_i,
    dir.name = ".",
    PA.nb.rep = 10,
    PA.nb.absences = 10000,
    PA.strategy = "disk",
    PA.dist.max = 110000,
    filter.raster = TRUE
    )
  # model building
  myBiomodModelOut <- BIOMOD_Modeling(
    bm.format = myBiomodData,
    modeling.id = paste0("EcoMod_", sp_i),
    models = c("GAM", "GLM", "MAXNET", "XGBOOST"),
    CV.strategy = "kfold",
    CV.k = 5,
    CV.do.full.models = FALSE,
    OPT.strategy = "bigboss",
    var.import = 5,
    metric.eval = c("TSS", "ROC")
    )
  # ensemble model
  myBiomodEM <- BIOMOD_EnsembleModeling(
    bm.mod = myBiomodModelOut,
    models.chosen = "all",
    em.by = "all",
    em.algo = "EMmedian",
    metric.eval = c("TSS", "ROC"),
    var.import = 5
    )
  # projections for current, LGM, Holocene
  # current
  myBiomodEMProj_curr <- BIOMOD_EnsembleForecasting(
    bm.em = myBiomodEM,
    proj.name = "curr",
    new.env = clim_sel,
    models.chosen = "all",
    metric.binary = "TSS",
    dir.name = file.path(sp_i, "proj_curr")
    )
  # LGM
  myBiomodEMProj_lgm <- BIOMOD_EnsembleForecasting(
    bm.em = myBiomodEM,
    proj.name = "lgm",
    new.env = clim_lgm_crop,
    models.chosen = "all",
    metric.binary = "TSS",
    dir.name = file.path(sp_i, "proj_lgm")
    )
  # Holocene
  myBiomodEMProj_holo <- BIOMOD_EnsembleForecasting(
    bm.em = myBiomodEM,
    proj.name = "holo",
    new.env = clim_hol_crop,
    models.chosen = "all",
    metric.binary = "TSS",
    dir.name = file.path(sp_i, "proj_holo")
    )
  }

#### Plot models

library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(patchwork)
library(dplyr)
setwd("~/proj/Hdar_postdoc1/analisis/ENM_sp/")
occur_all <- read.csv("~/proj/Hdar_postdoc1/data/coords_ddRAD.csv") %>% select(-c(1:2))
lineages <- c("cn", "so", "Pa")

# basemap
countries <- ne_countries(scale = "medium", returnclass = "sf")
prov_arg <- ne_states(country = "Argentina", returnclass = "sf")
bbox <- st_bbox(st_as_sf(occur_all, coords = c("lon", "lat"), crs = 4326))
# time periods
times <- c("lgm", "holo", "curr")

# Find global min and max across all rasters
raster_paths <- c()
for (l in lineages) { for (t in times) { raster_paths <- c(raster_paths, paste0("./", l, "/proj_", t, "/proj_", t, "_", l, "_ensemble.tif")) } }
# extract global min/max (first layer of each raster)
vals_min <- c(); vals_max <- c()
for (path in raster_paths) {
  r <- rast(path)[[1]]
  vals_min <- c(vals_min, global(r, "min", na.rm = TRUE)[1, 1])
  vals_max <- c(vals_max, global(r, "max", na.rm = TRUE)[1, 1])
  }
global_min <- min(vals_min, na.rm = TRUE)
global_max <- max(vals_max, na.rm = TRUE)

# loop
plot_list <- list()
for (l in lineages) {
    for (t in times) {
        # load raster
        file_path <- paste0("./", l, "/proj_", t, "/proj_", t, "_", l, "_ensemble.tif")
        r <- rast(file_path)[[1]]
        df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
        names(df)[3] <- "suit"
        # normalize using global min/max
    df$suit_norm <- (df$suit - global_min) / (global_max - global_min)
    # keep only high suitability
    # thresh <- 0
    # df <- df %>% filter(suit_norm > thresh)
    # plot
    p <- ggplot() +
      geom_raster(data = df, aes(x = x, y = y, fill = suit_norm)) +
      geom_sf(data = countries, fill = NA, color = "#636363", linewidth = 0.8) +
      geom_sf(data = prov_arg, fill = NA, linetype = "dashed", color = "#636363", linewidth = 0.4) +
      coord_sf(xlim = c(bbox["xmin"] - 2, bbox["xmax"] + 2), ylim = c(bbox["ymin"] - 2, bbox["ymax"] + 2), expand = FALSE) +
      scale_fill_viridis_c(option = "D", limits = c(0, 1), na.value = "transparent") +
      theme_bw() +
      theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
      # store
      plot_name <- paste0(l, "_", t)
      plot_list[[plot_name]] <- p
      }
    }

panel_plot <- (plot_list[["Pa_lgm"]] | plot_list[["Pa_holo"]] | plot_list[["Pa_curr"]]) /
              (plot_list[["cn_lgm"]] | plot_list[["cn_holo"]] | plot_list[["cn_curr"]]) /
              (plot_list[["so_lgm"]] | plot_list[["so_holo"]] | plot_list[["so_curr"]])
panel_plot # SVG 800 x 800





df <- as.data.frame(rast("./so/proj_lgm/proj_lgm_so_ensemble.tif")[[1]], xy = TRUE, na.rm = TRUE)
range(df$so_EMmedianByTSS_mergedData_mergedRun_mergedAlgo)
