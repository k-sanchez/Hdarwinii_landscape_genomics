###########
# ENM -----
###########

setwd("~/proj/Hdar_postdoc1/analisis/lscape_genetics/ENM/")
library(terra)
library(rnaturalearth)
library(caret)
library(biomod2)

save.image("ENM.RData")
load("ENM.RData")

#### Data ####

# Occurrences
occur <- read.csv("data/occurrences.csv")[, -3]
# remove duplicated coordinates
occur_unique <- unique(occur)
# rows to modify (move "inland")
occur_unique[95, "lon"] <- occur_unique[95, "lon"] - 0.5
occur_unique[102, "lon"] <- occur_unique[102, "lon"] - 0.1

# Geo-climatic data
# map of Argentina
arg <- ne_countries(country = "argentina", returnclass = "sv")
clim_files <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
clim <- rast(clim_files)
names(clim) <- sub(".*_", "", names(clim)) # rename
# reproject to same CRS
target_crs <- "+proj=longlat +datum=WGS84 +no_defs"
clim <- project(clim, target_crs)
# crop to sampling extent
bbox <- ext(min(occur_unique$lon), max(occur_unique$lon), min(occur_unique$lat), max(occur_unique$lat))
extent <- bbox + c(-2, +2, -2, +2)
clim_crop <- crop(clim, extent)
# check NAs
dat <- vect(occur_unique, geom = c("lon", "lat"), crs = crs(target_crs))
values <- terra::extract(clim_crop, dat)
sum(is.na(values))

# EVI data
evi_files <- list.files("~/GIS/capas/OpenLandMap/EVI/", "\\.tif", full.names = TRUE)
evi <- rast(evi_files)
# crop and match resolution of climate data 
evi <- crop(evi, extent)
evi <- resample(evi, clim_crop, method = "average")
evi <- project(evi, target_crs)

# summarize the data on a pixel-by-pixel basis
evi <- app(evi, max, na.rm = TRUE)
names(evi) <- "evi" # rename

# set same No Data mask for all variables
na_clim <- app(is.na(clim_crop), sum)
plot(na_clim)
na_evi <- is.na(evi)
plot(na_evi)
plot(na_clim + na_evi)
na_mask <- (na_clim + na_evi) > 0
plot(na_mask)
# mask
clim_crop[na_mask] <- NA
evi[na_mask] <- NA

# export rasters
writeRaster(clim_crop, "clim.tiff", overwrite = TRUE)
writeRaster(evi, "evi.tiff", overwrite = TRUE)
# load
clim_crop <- rast("data/clim.tiff")
evi <- rast("data/evi.tiff")

# Ensure that no presence data is repeated in the same pixel

dt <- terra::extract(evi, occur_unique[, c("lon", "lat")], cells = TRUE)
# remove presences that fall in No Data
mask <- is.na(dt$evi)
sum(mask)
dt <- dt[!mask, ]
pres <- occur_unique[!mask, ]
# remove duplicates in each pixel
pres$duplicated <- NA
sp_dup <- duplicated(dt$cell)
pres$duplicated <- sp_dup
# how many duplicates?
sum(pres$duplicated)
# filter and export
final_pres <- pres[!pres$duplicated, 1:2]
write.table(final_pres, "data/final_pres.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
pres <- read.table("data/final_pres.tsv", sep = "\t", header = TRUE)[, c(2, 1)]

#### Variable selection ####

# stack all rasters
rst <- c(clim_crop, evi)
# define spatial vector
v <- vect(final_pres, geom = c("lon", "lat"))
# define training area with a buffer of 1 degree
bsize <- 1
buf <- buffer(v, bsize)
buf <- aggregate(buf)
# check
plot(evi); plot(buf, cex = 0.25, add = TRUE); plot(v, add = TRUE)
# extract environmental data
dt <- terra::extract(rst, buf, ID = FALSE)
# remove NAs coming from ocean pixels
dt <- dt[!is.na(dt[, 1]), ]

# Pairwise correlations
corr <- cor(dt)
rmv_cor <- findCorrelation(corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact = TRUE)
sel_rst <- rst[[-rmv_cor]]
# write final raster
writeRaster(sel_rst, "data/final_vars.tiff", overwrite = TRUE)
sel_rst <- rast("data/final_vars.tiff")

# Prepare dataset for projection
arg <- crop(arg, extent)
vars <- crop(sel_rst, arg, mask = TRUE)
writeRaster(vars, "data/proj_current.tiff", overwrite = TRUE)
curp <- rast("data/proj_current.tiff")

#### Model building ####

myBiomodData <- BIOMOD_FormatingData(resp.var = rep(1, nrow(pres)),
                                     expl.var = sel_rst,
                                     resp.xy = pres,
                                     resp.name = "Hdar",
                                     dir.name = "models",
                                     PA.nb.rep = 10, # recommended 10 replicates
                                     PA.nb.absences = 10000, # at least 3 * the number of presences
                                     PA.strategy = "disk",
                                     PA.dist.max = 110000,
                                     filter.raster = TRUE # filtra si aparecen varios puntos en la misma celda
                                     )

# modeling
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = "EcoMod",
                                    models = c("GAM", "GLM", "MAXNET", "XGBOOST"),
                                    CV.strategy = "kfold",
                                    CV.k = 5,
                                    CV.do.full.models = FALSE,
                                    OPT.strategy = "bigboss",
                                    var.import = 5,
                                    metric.eval = c("TSS", "ROC"))

plt <- bm_PlotEvalMean(myBiomodModelOut, dataset = "calibration")
plt <- bm_PlotEvalMean(myBiomodModelOut, dataset = "validation")
# compare the runs and algorithms
plt <- bm_PlotEvalBoxplot(myBiomodModelOut, group.by = c("algo", "run"))
# variable importance
plt <- bm_PlotVarImpBoxplot(myBiomodModelOut, group.by = c("expl.var", "run", "algo"))

# Model ensemble

# modelname <- load("models/Hdar/Hdar.EcoMod.models.out")
# myBiomodModelOut <- eval(str2lang(modelname))
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = "all",
                                      em.by = "all",
                                      em.algo = "EMmedian",
                                      metric.eval =  c("TSS", "ROC"),
                                      var.import = 5)
# check performance
plt <- bm_PlotEvalBoxplot(myBiomodEM, group.by = c("algo", "algo"))
# variable importance
plt <- bm_PlotVarImpBoxplot(myBiomodEM, group.by = c("expl.var", "algo", "algo"))
# response curves
plt <- bm_PlotResponseCurves(myBiomodEM, fixed.var = "median")

#### Projecting model ####

myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                             proj.name = "current",
                                             new.env = curp,
                                             models.chosen = "all",
                                             metric.binary = "TSS")
# continuous model: proj_current_Hdar_ensemble.tif
# binary model: proj_current_Hdar_ensemble_TSSbin.tif





####################
# algatr -----------
# Landscape genomics
####################

# devtools::install_github("TheWangLab/algatr")
library(algatr)
# install all of the packages for algatr:
# alazygatr_packages()
library(vcfR)
library(terra)
library(caret)
library(RStoolbox)
library(dplyr)
library(tess3r)
library(gdm)
library(vegan)

# repeat this for the 3 sp comparisons
setwd("~/proj/Hdar_postdoc1/analisis/lscape_genetics/algatr_sp/Pa-so/")
rm(list = ls())

#### Data ####

# Genotypes
vcf_file <- "~/proj/Hdar_postdoc1/data/RADseq/stacks/denovo_map_whigroup/populations_algatr/Pa-so_pops/populations.snps.filt.vcf"
vcf <- read.vcfR(vcf_file)
samples <- colnames(vcf@gt)[-1]

# Coordinates
coords <- read.csv("~/proj/Hdar_postdoc1/data/coords_ddRAD.csv")
# filter coordinates for samples in VCF and reorder
coords_sub <- coords %>% filter(voucher %in% samples) %>% arrange(match(voucher, samples))
# how many samples in each group?
summary(factor(coords_sub$cluster))
points <- coords_sub[, 4:5]
colnames(points) <- c("x", "y")
# define bounding box for cropping layers
bbox <- ext(min(points$x), max(points$x), min(points$y), max(points$y))
expanded_bbox <- bbox
expanded_bbox[1] <- bbox[1] - 2
expanded_bbox[2] <- bbox[2] + 2
expanded_bbox[3] <- bbox[3] - 2
expanded_bbox[4] <- bbox[4] + 2

# Environmental
stacked_rast <- rast("~/proj/Hdar_postdoc1/analisis/lscape_genetics/algatr_sp/stacked_rast.tiff")
stacked_rast_crop <- crop(stacked_rast, expanded_bbox)

# Resistance
# elevation
dem <- rast("~/OneDrive/GIS/capas/ENVIREM/elev_SAmerica_current_30arcsec_geotiff/current_30arcsec_tri.tif")
dem_crop <- crop(dem, expanded_bbox)
# ENM
enm <- rast("~/proj/Hdar_postdoc1/analisis/lscape_genetics/ENM/models/Hdar/proj_current/proj_current_Hdar_ensemble.tif")
enm_crop <- crop(enm[[1]], expanded_bbox)

# Process genetic data
# convert VCF to dosage matrix
dosage <- vcf_to_dosage(vcf)
# impute missing genotypes
str_dos <- str_impute(gen = dosage, K = 2, entropy = TRUE, repetitions = 100, quiet = FALSE, save_output = FALSE)
# explore the imputation
dosage[1:5, 1:10]; str_dos[1:5, 1:10]
freqs  <- str_dos / 2 # allele frequencies
# aggregate allele frequencies by locality
site_freqs <- aggregate(freqs, by = list(x = points$x, y = points$y), FUN = mean, na.rm = TRUE)
# genetic distances
gendist <- as.matrix(dist(site_freqs, method = "euclidean"))

# Process geographic data
loc_coords <- site_freqs[, c("x", "y")]
geodist <- geo_dist(loc_coords)

# Process environmental data
env_values <- raster::extract(stacked_rast_crop, loc_coords)[, -1]
colnames(env_values) <- gsub("current_30arcsec_", "", colnames(env_values))
# pairwise correlations:
corr <- cor(env_values)
rmv_cor <- findCorrelation(corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact = TRUE)
env_values_uncor <- env_values[, -rmv_cor]
# scale
env_scaled <- scale(env_values_uncor, center = TRUE, scale = TRUE)

#### Pop structure ####

# does this structure corresponds to the delimited entities?
tess3_obj <- tess3(str_dos, coord = as.matrix(points), K = 2, method = "projected.ls", ploidy = 2)
q <- qmatrix(tess3_obj, K = 2)
barplot(q, border = NA, space = 0, xlab = "Individuals", ylab = "Ancestry proportions",  main = "Ancestry matrix") -> bp
labels <- coords_sub$voucher[bp$order]
axis(1, at = 1:nrow(q), labels = labels, las = 3, cex.axis = .4) # PDF: 500/600 x 200


#### Isolation hypotheses: MMRR ####
# IBD/IBE/IBR

# Dependent variable (gen. distances):
Y <- gendist
# Independent variables:
# environmental distances
X <- env_dist(env_scaled)
# add geographic distance
X[["geodist"]] <- geodist
# add topographic distance
X[["topodist"]] <- geo_dist(loc_coords, type = "topographic", lyr = dem_crop)
# add resistance distance
X[["resistdist"]] <- geo_dist(loc_coords, type = "resistance", lyr = enm_crop)


# MMRR
mmrr_full <- mmrr_run(Y, X, nperm = 999, stdz = TRUE, model = "full")
mmrr_plot(Y, X, mod = mmrr_full$mod, plot_type = "all", stdz = TRUE) # PDF 400/500 x 400
mmrr_full$coeff_df %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% arrange(desc(abs(estimate)))

#### Isolation hypotheses: GDM ####

# gendist vs coords + env
gdm_full <- gdm_run(gendist = gendist, coords = loc_coords, env = env_scaled, model = "full", scale_gendist = TRUE)
gdm_plot_isplines(gdm_full$model, scales = "free_x") # SVG: 500/400 x 400
gdm_df(gdm_full) %>% filter(coefficient != 0) %>% arrange(desc(coefficient)) %>% bind_rows(data.frame(predictor = "%explained", coefficient = gdm_full$model$explained)) %>% print()
# gendist vs resist (topo) + env
gdm_full_topo <- gdm_run(gendist = gendist, coords = loc_coords, env = env_scaled, model = "full", scale_gendist = TRUE, geodist_type = "topographic", dist_lyr = dem_crop)
gdm_plot_isplines(gdm_full_topo$model, scales = "free_x")
gdm_df(gdm_full_topo) %>% filter(coefficient != 0) %>% arrange(desc(coefficient)) %>% bind_rows(data.frame(predictor = "%explained", coefficient = gdm_full_topo$model$explained)) %>% print()
# gendist vs resist (ENM) + env
gdm_full_res <- gdm_run(gendist = gendist, coords = loc_coords, env = env_scaled, model = "full", scale_gendist = TRUE, geodist_type = "resistance", dist_lyr = enm_crop)
gdm_plot_isplines(gdm_full_res$model, scales = "free_x")
gdm_df(gdm_full_res) %>% filter(coefficient != 0) %>% arrange(desc(coefficient)) %>% bind_rows(data.frame(predictor = "%explained", coefficient = gdm_full_res$model$explained)) %>% print()

# Variable importance:
gdmData <- gdm_format(gendist, loc_coords, env_scaled, scale_gendist = TRUE)
varimp <- gdm.varImp(gdmData, geo = TRUE, nPerm = 100, parallel = TRUE, cores = 8)
importance <- varimp$`Predictor Importance`
pvalues <- varimp$`Predictor p-values`
importance_pvalues <- merge(importance, pvalues, by = "row.names")
colnames(importance_pvalues) <- c("Predictor", "Importance", "P_Value")
importance_pvalues[order(-importance_pvalues$Importance), ]
# w/topo
gdmData_topo <- gdm_format(gendist, loc_coords, env_scaled, scale_gendist = TRUE, geodist_type = "topographic", dist_lyr = dem_crop)
varimp_topo <- gdm.varImp(gdmData_topo, geo = TRUE, nPerm = 100, parallel = TRUE, cores = 8)
importance <- varimp_topo$`Predictor Importance`
pvalues <- varimp_topo$`Predictor p-values`
importance_pvalues <- merge(importance, pvalues, by = "row.names")
colnames(importance_pvalues) <- c("Predictor", "Importance", "P_Value")
importance_pvalues[order(-importance_pvalues$Importance), ]
# w/ENM
gdmData_enm <- gdm_format(gendist, loc_coords, env_scaled, scale_gendist = TRUE, geodist_type = "resistance", dist_lyr = enm_crop)
varimp_enm <- gdm.varImp(gdmData_enm, geo = TRUE, nPerm = 100, parallel = TRUE, cores = 8)
importance <- varimp_enm$`Predictor Importance`
pvalues <- varimp_enm$`Predictor p-values`
importance_pvalues <- merge(importance, pvalues, by = "row.names")
colnames(importance_pvalues) <- c("Predictor", "Importance", "P_Value")
importance_pvalues[order(-importance_pvalues$Importance), ]

#### Genotype-environment associations: RDA ####

env_inds <- raster::extract(stacked_rast_crop, points)[, -1]
colnames(env_inds) <- gsub("current_30arcsec_", "", colnames(env_inds))
env_inds_uncor <- env_inds[, -rmv_cor]
env_inds_scaled <- scale(env_inds_uncor, center = TRUE, scale = TRUE)

# simple RDA
mod_full <- rda_run(str_dos, env_inds_scaled, model = "full")
RsquareAdj(mod_full)
# partial RDA with geography as covariable
mod_pRDA_geo <- rda_run(str_dos, env_inds_scaled, points, model = "full", correctGEO = TRUE, correctPC = FALSE)
RsquareAdj(mod_pRDA_geo)
# partial RDA with geography and pop structure as covariables
mod_pRDA_gs <- rda_run(str_dos, env_inds_scaled, points, model = "full", correctGEO = TRUE, correctPC = TRUE, nPC = 3)
RsquareAdj(mod_pRDA_gs)

# Variance partitioning with partial RDA
var_part <- rda_varpart(str_dos, env_inds_scaled, points, Pin = 0.05, R2permutations = 1000, R2scope = TRUE, nPC = 3)
options(width = 300)
var_part %>% select(-call) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))

# Candidate SNPs:
mod_pRDA <- rda_run(str_dos, env_inds_scaled, model = "full", correctPC = TRUE)
# Z-scores method
sig_z <- rda_getoutliers(mod_pRDA, naxes = "all", outlier_method = "z", z = 3, plot = FALSE)
length(sig_z$rda_snps)
# p-value method
sig_p <- rda_getoutliers(mod_pRDA, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.05, plot = FALSE)
length(sig_p$rda_snps)
# identify outliers that have q-values <= 0.1
sig_q <- sig_p$rdadapt %>%
  # make SNP names column from row names
  mutate(snp_names = row.names(.)) %>%  filter(q.values <= 0.1)
nrow(sig_q)
# SNPs detected by the three methods
rda_outliers <- Reduce(intersect, list(sig_z$rda_snps, sig_p$rda_snps, sig_q$snp_names))
length(rda_outliers)
 

#### Genotype-environment associations: LFMM ####

ridge_results <- lfmm_run(str_dos, env_inds_scaled, K = 2, lfmm_method = "ridge")
lasso_results <- lfmm_run(str_dos, env_inds_scaled, K = 2, lfmm_method = "lasso")

# extract significant SNPs based on both
lfmm_ridge_outliers <- ridge_results$df %>% filter(pvalue <= 0.05) %>% pull(snp)
length(lfmm_ridge_outliers)
lfmm_lasso_outliers <- lasso_results$df %>% filter(pvalue <= 0.05) %>% pull(snp)
length(lfmm_lasso_outliers)

lfmm_outliers <- intersect(lfmm_ridge_outliers, lfmm_lasso_outliers)
length(lfmm_outliers)

# Number of SNPs detected by RDA and LFMM
high_conf_outliers <- intersect(rda_outliers, lfmm_outliers)
length(high_conf_outliers)
# explore associations between SNPs and environment:
rda_gen <- str_dos[, high_conf_outliers]
cor_df <- rda_cor(rda_gen, env_inds_scaled)
cor_df %>%
  filter(snp %in% high_conf_outliers) %>%
  mutate(abs_r = abs(r)) %>% 
  group_by(snp) %>%
  slice_max(abs_r, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(abs_r)) %>%
  select(snp, var, r, p) %>%
  print(n = Inf)
