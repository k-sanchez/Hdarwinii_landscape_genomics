# QC of RADseq reads ----
## based on: https://github.com/DevonDeRaad/RADstackshelpR/blob/master/inst/extdata/fastqcr.Rmd

library(gridExtra)
library(ggplot2)
library(fastqcr)
library(knitr)
load("~/OneDrive/ongoing/Hdar_posdoc/data/RADseq/fastqcr/fastqcr.RData")


# specify full path to directory containing a .fastq.gz file for each sample
fq_dir <- "~/OneDrive/PRODUCCION/ongoing_papers/Hdar_posdoc/data/RADseq/fastq/"

# specify full path to the output directory where you want 
qc_dir <- "~/OneDrive/PRODUCCION/ongoing_papers/Hdar_posdoc/data/RADseq/fastqcr"

# Run fastqc on all .fastq.gz files, through r
# this only needs to be run once, if only tweaking downstream visualizations, you can comment out this step
fastqc(fastqc.path = "/usr/local/bin/fastqc",
       fq.dir = fq_dir, # FASTQ files directory
       qc.dir = qc_dir, # Results directory
       threads = 8      # Number of threads
       )

# List of files in the output directory to ensure fastqc worked
list.files(qc_dir)

# create a character vector where each value is the full path to the .zip created by fastqc() for a given sample
samps <- list.files(qc_dir, full.names = TRUE, pattern = "*.zip")

# plot qc test results for each sample
for (i in samps){
  # read info for given sample from the .zip file generated in the previous step
  samp_info <- qc_read(i)
  # open blank list to hold qc visualizations for the given sample
  plot <- list()
  # do qc for the given sample
  plot[[1]] <- qc_plot(samp_info, "Basic statistics")
  plot[[2]] <- qc_plot(samp_info, "Per sequence quality scores")
  plot[[3]] <- qc_plot(samp_info, "Sequence duplication levels")
  # visualize tables
  print(paste0("QC results for sample ", gsub(".*/", "", i)))

  cat('\n')

  print(kable(plot[[1]]))

  cat('\n')

  # visualize plots
  grid.arrange(plot[[2]], plot[[3]], ncol =2 )
  
  # clear plot to hold info for next sample
  rm(plot)
}

# aggregate the reports by pointing this function to the folder holding output of fastqc()
qc <- qc_aggregate(qc_dir, progressbar = TRUE)

# stats per sample
knitr::kable(qc_stats(qc))

### solid red line = median sample value
### dashed red line = 10% of median sample value

# save stats info as an object
stats_info <- qc_stats(qc)
# make tot.seq numeric
stats_info$tot.seq <- as.numeric(stats_info$tot.seq)

# make histogram of number of sequence reads for each sample
ggplot(stats_info, aes(x = tot.seq)) +
  geom_histogram(color = "black", fill = "white", bins = 20) +
  geom_vline(aes(xintercept = median(tot.seq)), color = "red") +
  geom_vline(aes(xintercept = median(tot.seq)*.1), color = "red", lty = 14) +
  theme_classic() +
  xlab("number") +
  ggtitle("Number of sequence reads per sample")

# solid red line = median sample value
# dashed red line = 10% of median sample value
ggplot(stats_info, aes(x = tot.seq)) +
  geom_histogram(color = "black", fill = "white", bins = 200) +
  geom_vline(aes(xintercept = median(tot.seq)), color = "red") +
  geom_vline(aes(xintercept = median(tot.seq)*.1), color = "red", lty = 14) +
  theme_classic() +
  xlab("Number of sequencing reads")

# show me the samples that have less than 10% of the number of reads as the median sample from this experiment (these should be dropped immediately)
print(paste("Median sample contains", median(stats_info$tot.seq), "reads. The following samples contain less than", median(stats_info$tot.seq)*.1, "reads (10% of the median), and should likely be dropped")); knitr::kable(stats_info[stats_info$tot.seq < median(stats_info$tot.seq)*.1, ])
# remove Hdar_Ca_15972 and Hdar_Cu_15243 from the dataset

save.image("~/OneDrive/PRODUCCION/ongoing_papers/Hdar_posdoc/data/RADseq/fastqcr/fastqcr.RData")


#######################################
# Parameter optimization in Stacks ----
#######################################

library(tidyverse)

setwd("~/OneDrive/PRODUCCION/ongoing/Hdar_posdoc/data/RADseq/stacks/param_opt_whigroup")
df <- read.table("result.tsv", header = TRUE)
df <- df %>% mutate(change = R80_loci - lag(R80_loci))

# line plot
ggplot(df, aes(x = M_n, y = change)) +
  scale_x_continuous(breaks = 2:6, labels = c("1/2", "2/3", "3/4", "4/5", "5/6")) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_point(color = "red", size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Change in ustacks M",y = "Change in R80 loci") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) # SVG: 250 * 250


####################
# SNP filtering ----
####################

# Filtros aplicados a las matrices:
# DAPC/LEA: bi/depth/genotype/inds./MAC/thin
# IBD/EEMS/landscape genomics: bi/depth/genotype/thin
# Snapper/networks: bi/depth/genotype/inds./thin
# networks (for SNPs2CF v1.7): bi/depth

library(vcfR)
library(SNPfiltR)
library(stringr)
library(dplyr)
setwd("~/proj/Hdar_postdoc1/data/RADseq/stacks/denovo_map_whigroup/populations_K8_wil_net_allSNPs")
vcf <- read.vcfR("populations.snps.vcf")

#### Biallelic and depth filters ####

# keep only bi-allelic SNPs
vcf_bi <- filter_biallelic(vcf)
# generate exploratory visualizations of depth for all called genotypes
hard_filter(vcfR = vcf_bi)
# hard filter to minimum depth of 10
vcf_bi_minD <- hard_filter(vcfR = vcf_bi, depth = 10)

# visualize and pick appropriate max depth cutoff
max_depth(vcf_bi_minD)
# a good rule of thumb is mean depth * 2
vcf_bi_minD_maxD <- max_depth(vcf_bi_minD, maxdepth = 200)

#### Missing data filters ####

# Round 1: 
# visualize proportion of missing data per genotype and sample
missing_by_snp(vcfR = vcf_bi_minD_maxD)
missing_by_sample(vcfR = vcf_bi_minD_maxD)
# genotype filtering: 50 % allowed missing data per SNP
vcf_bi_minD_maxD_gen1 <- missing_by_snp(vcf_bi_minD_maxD, cutoff = 0.5)
# sample filtering: 90 % allowed missing data per sample
vcf_bi_minD_maxD_gen1_ind1 <- missing_by_sample(vcf_bi_minD_maxD_gen1, cutoff = 0.9)

# Round 2:
# visualize proportion of missing data per sample and genotype
missing_by_snp(vcfR = vcf_bi_minD_maxD_gen1_ind1)
missing_by_sample(vcfR = vcf_bi_minD_maxD_gen1_ind1)
# genotype filtering: 60 % allowed missing data per SNP
vcf_bi_minD_maxD_gen2_ind1 <- missing_by_snp(vcf_bi_minD_maxD_gen1_ind1, cutoff = 0.6)
# sample filtering: 70 % allowed missing data per sample
vcf_bi_minD_maxD_gen2_ind2 <- missing_by_sample(vcf_bi_minD_maxD_gen2_ind1, cutoff = 0.7)

# Round 3:
# visualize proportion of missing data per sample and genotype
missing_by_snp(vcfR = vcf_bi_minD_maxD_gen2_ind2)
missing_by_sample(vcfR = vcf_bi_minD_maxD_gen2_ind2)
# genotype filtering: 70 % allowed missing data per SNP
vcf_bi_minD_maxD_gen3_ind2 <- missing_by_snp(vcf_bi_minD_maxD_gen2_ind2, cutoff = 0.7)
# sample filtering: 50 % allowed missing data per sample
vcf_bi_minD_maxD_gen3_ind3 <- missing_by_sample(vcf_bi_minD_maxD_gen3_ind2, cutoff = 0.5)

# Round 4:
# visualize proportion of missing data per sample and genotype
missing_by_snp(vcfR = vcf_bi_minD_maxD_gen3_ind3)
missing_by_sample(vcfR = vcf_bi_minD_maxD_gen3_ind3)
# genotype filtering: 80 % allowed missing data per SNP
vcf_bi_minD_maxD_gen4_ind3 <- missing_by_snp(vcf_bi_minD_maxD_gen3_ind3, cutoff = 0.8)
# sample filtering: 30 % allowed missing data per sample
vcf_bi_minD_maxD_gen4_ind4 <- missing_by_sample(vcf_bi_minD_maxD_gen4_ind3, cutoff = 0.3)


#### MAC and thinning ####

vcf_filt_mac <- min_mac(vcf_bi_minD_maxD_gen4_ind3, min.mac = 1)
vcf_filt_mac_unlink <- distance_thin(vcf_filt_mac, min.distance = 50)

#### Assess clustering ####

sample_names <- colnames(extract.gt(vcf_filt_mac_unlink))
pop_names <- str_extract(sample_names, ".*(?=_[^_]+$)")
popmap <- data.frame(id = sample_names, pop = pop_names)
# pca
assess_missing_data_pca(vcf_filt_mac_unlink, popmap, clustering = FALSE)
# t-sne
assess_missing_data_tsne(vcf_filt_mac_unlink, popmap, clustering = FALSE)

#### Final matrices ####

# plot depth per SNP and per sample
dp <- extract.gt(vcf_filt_mac_unlink, element = "DP", as.numeric = TRUE)
heatmap.bp(dp, rlabels = FALSE)

# compare SNPs matrices
get_vcf_info <- function(vcf) {
  num_individuals <- ncol(vcf@gt) - 1  # exclude "FORMAT" column
  num_snps <- nrow(vcf@fix)
  return(c(num_individuals, num_snps))
}
# list VCF objects
vcf_list <- list(raw = vcf,
                 round1 = vcf_bi_minD_maxD_gen1_ind1, 
                 round2 = vcf_bi_minD_maxD_gen2_ind2,
                 round3 = vcf_bi_minD_maxD_gen3_ind3,
                 round4 = vcf_bi_minD_maxD_gen4_ind4,
                 final = vcf_filt_mac_unlink)
# extract information and create data frame
vcf_info <- do.call(rbind, lapply(vcf_list, get_vcf_info))
df <- data.frame(n_inds = vcf_info[,1], 
                 n_SNPs = vcf_info[,2],
                 row.names = names(vcf_list))
df

# export VCF
# write.vcf(vcf_bi_minD_maxD_gen4_ind4, file = "populations.snps.filt.vcf.gz")
write.vcf(vcf_filt_mac_unlink, file = "populations.snps.filt.vcf.gz")


#################################
# Sample PHYLIP MSAs for BPP ----
#################################

library(ape)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)

msa_dir <- "~/proj/Hdar_postdoc1/data/RADseq/stacks/denovo_map_whigroup/populations_inds_dar_wil_loci/loci/" # folder with .phy files
seq_table_file <- "~/proj/Hdar_postdoc1/analisis/bpp/seq_assignment.txt" # table: voucher cluster
out_dir <- "~/proj/Hdar_postdoc1/analisis/hhsd/tmp"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(2442)
# mapping table
seq_table <- read.table(seq_table_file, sep = "\t", header = TRUE)

# function to process one PHYLIP file
process_phylip <- function(file, seq_table, n_per_lineage = 2) {
  aln <- read.dna(file, format = "sequential")
  seq_names <- rownames(aln)  
  # join with lineage assignment
  df <- tibble(voucher = seq_names) %>% left_join(seq_table, by = c("voucher"))
  # sample up to n_per_lineage sequences per lineage
  sampled <- df %>%
    group_by(cluster) %>%
    group_modify(~{
      n_rows <- nrow(.x)
      # if lineage has fewer than n_per_lineage, take all
      if (n_rows <= n_per_lineage) {
        .x
      } else {
        slice_sample(.x, n = n_per_lineage)
      }
    }) %>%
    mutate(new_name = paste0(cluster, "_", row_number())) %>%
    ungroup()  
  # subset alignment to sampled sequences
  aln_sampled <- aln[sampled$voucher, , drop = FALSE]
  rownames(aln_sampled) <- sampled$new_name  
  # write output
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(file)), "_samp.phy"))
  write.dna(aln_sampled, file = out_file, format = "sequential", nbcol = -1, colsep = "")
  return(out_file)
  }

# apply to all MSAs
phy_files <- list.files(msa_dir, pattern = "\\.phylip$", full.names = TRUE)
results <- list()
for (i in seq_along(phy_files)) {
  file <- phy_files[i]
  out_file <- process_phylip(file, seq_table = seq_table, n_per_lineage = 2)
  results[[file]] <- out_file
  }
