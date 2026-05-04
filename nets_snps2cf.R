##################################################
# CF table for network estimation SNPs2CF v1.7----
##################################################

setwd("~/proj/Hdar_postdoc1/analisis/phylonetworks_K8")
# source the function
source("~/OneDrive/code_data/filogenia/SNPs2CF/functions_v1.7.R")
# VCF to phy
library(pegas)
vcfpath <- "/Users/kevinsanchez/Library/CloudStorage/OneDrive-Personal/PRODUCCION/ongoing/Hdar_postdoc1/data/RADseq/stacks/denovo_map_whigroup/populations_K8_net_allSNPs/populations.snps.filt.vcf" # absolute paths, otherwise it crashes
vcf <- read.vcf(vcfpath, to = 100000); vcf
vcf2phylip(vcf.name = vcfpath, total.SNPs = dim(vcf)[2], output.name = "populations.snps.filt.phy")
phy <- "populations.snps.filt.phy" # PUT .phy MATRIX IN WORKING DIRECTORY

# Generate CF table
SNPs2CF(seqMatrix = phy,
        multipleSNPs.per.locus = TRUE,
        kept.sites.file = "~/proj/Hdar_postdoc1/data/RADseq/stacks/denovo_map_whigroup/populations_K8_net_allSNPs/populations.snps.filt.kept.sites",
        ImapName = "popmap.txt",
        between.sp.only = TRUE,
        max.SNPs = NULL,
        bootstrap = TRUE,
        outputName = "CF.csv",
        save.progress = FALSE,
        n.quartets = 100,
        indels.as.fifth.state = FALSE,
        boots.rep = 100,
        starting.sp.quartet = 1,
        max.quartets = 100000,
        cores = 4)


#### CF histograms ####

# rows <- 1
# cols <- 5
# mt <- matrix(1:(rows * cols), nrow = rows, ncol = cols, byrow = TRUE)
# par(family = "CMU Sans Serif")
# plotCF(CF.table.name = "CF.csv", col = "gray", plot.stats = TRUE, asterisk.cex = 1.2, p.lines.cex = 1, y.line.adj = 0.02, xlabel.cex = 1)
