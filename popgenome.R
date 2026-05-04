# devtools::install_github("JLSteenwyk/ggpubfigs")
library(vcfR)
library(adegenet)
library(tidyverse)
library(rnaturalearth)
library(RColorBrewer)
library(ape)
library(plotly)

###########
# DAPC ----
###########

setwd("~/proj/Hdar_posdoc1/analisis/popgenome/DAPC/")
vcf <- read.vcfR("~/proj/Hdar_posdoc/data/RADseq/stacks/denovo_map_whigroup/populations_popstructure/populations.snps.filt.vcf")
gl <- vcfR2genlight(vcf)
nclust <- find.clusters(gl, n.pca = dim(gl)[1], max.n.clust = 15) # 6-8; SVG = 300 * 300
dapc <- dapc(gl, nclust$grp, n.pca = dim(gl)[1], n.da = length(levels(nclust$grp)))
pca <- glPca(gl, nf = dim(gl)[1])

# show sample membership
df <- data.frame(dapc$grp)
df <- tibble::rownames_to_column(df, var = "sample")
write.csv(df, "dapc_K8_grp.csv", row.names = FALSE, quote = FALSE)

# tidy the data for plotting
pca_df <- as_tibble(pca$scores, rownames = "individual") %>%
        mutate(group = dapc$grp)
groups <- unique(dapc$grp)
colors <- brewer.pal(length(groups), name = "Accent")
group_colors <- setNames(colors, groups)
# plot the data
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = group)) +
        geom_point(shape = 21, size = 3) +
        scale_fill_manual(values = group_colors) +
        theme_bw(base_size = 16)
pca_plot # SVG = 400 * 350

# 3D visualization
pca_plot_3d <- plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3,
                       color = ~group, colors = group_colors,
                       marker = list(size = 5),
                       text = ~individual, # show specimen 
                       hoverinfo = "text") %>%
add_markers() %>%
layout(scene = list(xaxis = list(title = "PC1"),
                    yaxis = list(title = "PC2"),
                    zaxis = list(title = "PC3")),
       title = "3D PCA Plot")
pca_plot_3d
# export
htmlwidgets::saveWidget(pca_plot_3d, "dapc_K8_3D.html")

# Map 

data <- read.csv("~/proj/Hdar_posdoc1/analisis/popgenome/data.csv")
# merge locality data with the genetic data
locs <- data %>% mutate(individual = voucher) %>% full_join(pca_df, by = "individual")
# read in a basemap for plotting
arg <- ne_countries(country = "argentina", returnclass = "sf")
# plot the map
map <- ggplot(data = arg) +
        geom_sf() +
        geom_point(data = locs, aes(x = lon, y = lat, fill = group, size = 1), shape = 21) +
        scale_fill_manual(values = group_colors) +
        guides(size = "none", fill = guide_legend(override.aes = list(size = 3))) +
        coord_sf(xlim = c(-75, -60), ylim = c(-55, -33)) +
        theme_bw(base_size = 16)
map # SVG: 400 * 500


##########
# LEA ----
##########

# BiocManager::install("LEA")
library(LEA)
library(phytools)
library(tidyverse)
library(ggpubfigs)
setwd("~/proj/Hdar_posdoc1/analisis/popgenome/LEA")
vcf2geno("path_to_vcf")

# sNMF
project <- NULL
project <- snmf("populations.snps.filt.geno", K = 1:15, ploidy = 2, entropy = TRUE, repetitions = 100, project = "new", alpha = 100, CPU = 6)
# to continue a previously interrupted: project = "continue"
# plot cross-entropy criterion for all runs in the SNMF project
project <- load.snmfProject("populations.snps.filt.snmfProject")
plot(project, col = "blue", pch = 19, cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, main = "Cross-entropy criterion")

#### Barplot ####

data <- read.csv("~/proj/Hdar_postdoc1/analisis/popgenome/data.csv")
# read the ML tree
# we will use it to set the order of the bars in the barplot
tree <- read.newick("~/proj/Hdar_postdoc1/analisis/phylogenetics/ML_whigroup_loci_unphased/out_part.raxml.support")
# get the tip labels for the samples that contain "Hwhi" or "Hand"
root_samples <- grep("Hwhi|Hand", tree$tip.label, value = TRUE)
# find their MRCA and root on this node
root_node <- getMRCA(tree, root_samples)
tree <- root(tree, node = root_node)
# find matches between the tree and the samples used in the pop analyzes
matches <- intersect(tree$tip.label, data$voucher)
prun_tree <- keep.tip(tree, matches)
prun_tree <- ladderize(prun_tree)

# select the best run for K
k <- 8
best <- which.min(cross.entropy(project, K = k))
file <- paste0("populations.snps.filt.snmf/K", k, "/run", best, "/populations.snps.filt_r", best, ".", k, ".Q")
qmatrix <- read.table(file, header = FALSE)
qmatrix$voucher <- data$voucher
# reorder qmatrix based on the tree labels
qmatrix <- qmatrix %>% mutate(voucher = factor(voucher, levels = prun_tree$tip.label)) %>% arrange(voucher)
# reshape qmatrix into long format
qmatrix_long <- qmatrix %>% pivot_longer(cols = -voucher, names_to = "cluster", values_to = "proportion")
# plot
palette <- friendly_pal("muted_nine", k)
ggplot(qmatrix_long, aes(x = voucher, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = palette) +
  theme_simple() +
  labs(x = "Samples", y = "Proportion", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
  
# export the qmatrix
traits <- qmatrix[, c(ncol(qmatrix), 1:(ncol(qmatrix) - 1))]
colnames(traits) <- c("Traits", paste0("K", 1:k))
write.table(traits, "dar_wil_traits_K8.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# extract the group with the highest value for each individual
qmatrix_dt <- data.table::as.data.table(qmatrix_long)
# find the group with the maximum value for each individual
result <- qmatrix_dt[, .SD[which.max(proportion)], by = voucher]
write.csv(result, "max_prop_samples_K8.csv", quote = FALSE, row.names = FALSE)

