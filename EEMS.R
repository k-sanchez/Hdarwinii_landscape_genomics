#############################
# Create habitat polygon ----
#############################

library(ggplot2)
library(sf)
library(dplyr)
library(rnaturalearth)
setwd("~/proj/Hdar_posdoc1/analisis/eems/")

df <- read.csv("~/proj/Hdar_posdoc/analisis/lscape_genetics/coords_K8_fix.csv")
df_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
world <- ne_countries(scale = "medium", returnclass = "sf")
prov <- ne_states(country = "Argentina", returnclass = "sf")
# identify points outside land boundaries
out_points <- df_sf %>% filter(!apply(st_within(df_sf, world, sparse = FALSE), 1, any))
# shift westward only the longitude of the points that fall outside land boundaries
df_sf_fix <- df_sf
df_sf_fix$geometry[df_sf_fix$voucher %in% out_points$voucher] <- 
  st_geometry(df_sf_fix[df_sf_fix$voucher %in% out_points$voucher, ]) - c(0.1, 0)
# check again
df_sf_fix %>% filter(!apply(st_within(df_sf_fix, world, sparse = FALSE), 1, any))
# save
df_csv <- df_sf_fix %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()
write.csv(df_csv, "~/proj/Hdar_posdoc/analisis/lscape_genetics/coords_K8_fix.csv", row.names = FALSE)

bbox <- st_bbox(df_sf_fix)
expanded_bbox <- bbox
expanded_bbox["xmin"] <- bbox["xmin"] - 1
expanded_bbox["xmax"] <- bbox["xmax"] + 1
expanded_bbox["ymin"] <- bbox["ymin"] - 1
expanded_bbox["ymax"] <- bbox["ymax"] + 1

# create a buffer around points
buff <- df_sf %>%
    st_union() %>%  # merge points into a single geometry
    # st_convex_hull(),
    st_concave_hull(ratio = 0.6) %>%
    st_buffer(dist = 50000)
# crop hull to land boundaries
world_merge <- world %>% st_union()
buff_crop <- st_intersection(buff, world_merge)
# apply additional buffer
buff_crop <- st_buffer(buff_crop, dist = 20000)

ggplot() +
  geom_sf(data = world, fill = "lightblue", linewidth = 0.5) +
  geom_sf(data = prov, fill = NA, linetype = "dashed", linewidth = 0.2) +
  geom_sf(data = buff_crop, fill = "red", alpha = 0.2) +
  geom_sf(data = df_sf_fix, color = "black", size = 2) +
  coord_sf(xlim = c(expanded_bbox["xmin"], expanded_bbox["xmax"]),
           ylim = c(expanded_bbox["ymin"], expanded_bbox["ymax"]),
           expand = FALSE) +
  theme_bw()
  
# convert closed polygon to coordinates
coords <- as.data.frame(st_coordinates(buff)[, 1:2])
colnames(coords) <- c("lon", "lat")
write.table(coords, "samples.outer", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
# export coordinates df
# check same order between vcf and coords file
vcf <- vcfR::read.vcfR("~/proj/Hdar_posdoc/data/RADseq/stacks/denovo_map_whigroup/populations_eems/populations.snps.filt.vcf")
identical(colnames(vcf@gt)[-1], df_csv$voucher) # TRUE
write.table(df_csv[, 3:4], "samples.coord", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")