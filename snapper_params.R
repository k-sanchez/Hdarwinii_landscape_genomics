#################
# Snapper params ----
################# 

# Estimate the parameters (mean and sd) of a log-normal distribution
# from given mean and 95 % HPD

# this is based on a previous time-calibrated phylogeny
mean_ln <- 6.8
HPD_lower <- 5.1
HPD_upper <- 8.3
# approximate z-score for a 95% confidence interval
z <- 1.96
# estimate the median in real space
median_ln <- sqrt(HPD_lower * HPD_upper)
# solve for mu (log of the median)
mu <- log(median_ln)
# solve for sigma (s.d. in log space)
sigma <- (log(HPD_upper) - mu) / z  
# compute mean in real space
mean_real <- exp(mu + (sigma^2) / 2)

cat("Mean (real space):", mean_real, "\n")
cat("Standard deviation (log-transformed distribution):", sigma, "\n")


#### Densitree ####

library(ggtree)
library(deeptime)
mcc <- treeio::read.beast("~/proj/Hdar_postdoc1/analisis/phylogenetics/snapper_K8/snapper_whigroup_msc1_mcc_median.tre")
posterior <- treeio::read.newick("~/Downloads/A01nucOGadmR1.mcmc.txt")
posterior_subset <- posterior[sample(seq_along(posterior), 100)]
ggdensitree(posterior_subset, alpha = .3, colour = 'steelblue') +
    geom_tiplab(size = 3) +
    #coord_geo(pos = list("bottom", "bottom"), xlim = c(-8, 0), ylim = c(0, 10), neg = TRUE, dat = list("periods", "epochs"), abbrv = list(TRUE, FALSE), expand = TRUE, size = "auto") +
    #labs(x = "Time (Ma)", y = NULL) +
    #scale_x_continuous(breaks = seq(-8, 0, 2), labels = abs(seq(-8, 0, 2))) +
    theme_tree2()

