##########
# IBD ----
##########

# install.packages("prabclus")
# devtools::install_github("thierrygosselin/radiator")
library(prabclus)
library(radiator)
library(dplyr)

setwd("~/proj/Hdar_postdoc1/analisis/")
df <- read.csv("~/proj/Hdar_postdoc1/data/coords_ddRAD.csv")
# only keep clusters to test
df_filt <- df %>% filter(cluster %in% c("WCH", "WRN"))
# if there are grouped clusters also run this: 
df_filt <- df_filt %>% mutate(cluster = ifelse(cluster %in% c("MZ_NNQ", "CRN_SEC", "SNQ", "Pa1", "Pa2", "WRN"), "MZ_NNQ_CRN_SEC_SNQ_Pa12_WRN", cluster)
                              # cluster = ifelse(cluster %in% c("Pa1", "Pa2"), "Pa1_Pa2", cluster)
                              )

geodist <- log1p(coord2dist(coordmatrix = df_filt[, 5:4], file.format = "decimal2", output.dist = TRUE)) # lat, lon
vcf_file <- "~/proj/Hdar_postdoc1/data/RADseq/stacks/denovo_map_whigroup/populations_ibd/WCH-WRN_pops/populations.snps.filt.vcf"
# force same order between vcf and coords file
vcf <- vcfR::read.vcfR(vcf_file)
df_filt <- df_filt[order(factor(df_filt$voucher, levels = colnames(vcf@gt)[-1])), ]
identical(colnames(vcf@gt)[-1], df_filt$voucher) # TRUE?

genomic_converter(data = vcf_file, output = "genepop", filter.common.markers = FALSE, filter.monomorphic = FALSE, filename = "snps")
df_gp <- read.table(paste(list.dirs()[2], "/snps_genepop.gen", sep = ""), sep = "", skip = 3)[, -1]
unlink(grep("radiator", list.dirs(), value = TRUE), recursive = TRUE) # remove directory
snps_pc <- alleleconvert(strmatrix = df_gp, format.in = "genepop", format.out = "prabclus")
snps_ao <- alleleinit(allelematrix = snps_pc, neighborhood = "none")
gendist <- alleledist(unbuild.charmatrix(snps_ao$charmatrix, snps_ao$n.individuals, snps_ao$n.variables), snps_ao$n.individuals, snps_ao$n.variables)

#### Tests ####

{
       grouping <- df_filt$cluster
       groups <- levels(factor(df_filt$cluster))
       h01 <- regeqdist(dmx = geodist, dmy = gendist, grouping = grouping, groups = groups)
       if ((min(h01$pval) * 2) > 0.05){ # bonferroni = min(h01$pval) * 2
       print(paste("H01 not rejected, p-value =", signif(min(h01$pval * 2), 3), "→ testing H02"))
       h02 <- regdistbetween(dmx = geodist, dmy = gendist, grouping = grouping, groups = groups)
       ifelse(h02$pval < 0.05,
              paste("Barrier to gene flow, p-value =", signif(h02$pval, 3)),
              paste("No barrier to gene flow (IBD), p-value =", signif(h02$pval, 3))
              )
       } else {
              print(paste("H01 rejected, p-value =", signif(min(h01$pval * 2), 3), "→ testing H03"))
              h03_1 <- regdistbetweenone(dmx = geodist, dmy = geodist, grouping = grouping, groups = groups, rgroup = groups[1])
              h03_2 <- regdistbetweenone(dmx = geodist, dmy = geodist, grouping = grouping, groups = groups, rgroup = groups[2])
              ifelse(h03_1$pval < 0.05 | h03_2$pval < 0.05,
                     paste0("Barrier to gene flow, p-value =", signif(h03_1$pval, 3), "/", signif(h03_2$pval, 3)),
                     paste("No barrier to gene flow (IBD), p-value =", signif(h03_1$pval, 3), "/", signif(h03_2$pval, 3))
                     )
              }
}
