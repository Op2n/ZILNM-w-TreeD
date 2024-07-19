###############
## Libraries ##
library(ggplot2)
library(dplyr)

###################
## load the data ##

## Scenarios 3: small T
load("running_time_smallT.RData")
load("Frobenius_norm_error_smallT.RData")
load("Simpson_index_mse_smallT.RData")
load("Wasserstein_distance_error_smallT.RData")
load("Frobenius_norm_error_smallT_2.RData")
load("Simpson_index_mse_smallT_2.RData")
load("Wasserstein_distance_error_smallT_2.RData")

######################################
## mean of the Frobenius norm error ##
Frob1_mean_smallT <- array(round(apply(Frobenius_norm_error_smallT, 2, mean), digits = 8), dim = c(1, 5))
colnames(Frob1_mean_smallT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Frob1_mean_smallT_noGBM <- array(Frob1_mean_smallT[, -4], dim = c(1, 4))
colnames(Frob1_mean_smallT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

Frob2_mean_smallT <- array(round(apply(Frobenius_norm_error_smallT_2, 2, mean), digits = 8), dim = c(1, 5))
colnames(Frob2_mean_smallT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Frob2_mean_smallT_noGBM <- array(Frob2_mean_smallT[, -4], dim = c(1, 4))
colnames(Frob2_mean_smallT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

######################################
## mean of the mse of Simpson index ##
Simpson1_mean_smallT <- array(round(apply(Simpson_index_mse_smallT, 2, mean), digits = 8), dim = c(1, 5))
colnames(Simpson1_mean_smallT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Simpson1_mean_smallT_noGBM <- array(Simpson1_mean_smallT[, -4], dim = c(1, 4))
colnames(Simpson1_mean_smallT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

Simpson2_mean_smallT <- array(round(apply(Simpson_index_mse_smallT_2, 2, mean), digits = 8), dim = c(1, 5))
colnames(Simpson2_mean_smallT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Simpson2_mean_smallT_noGBM <- array(Simpson2_mean_smallT[, -4], dim = c(1, 4))
colnames(Simpson2_mean_smallT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

############################################
## mean of the Wasserstein distance error ##
WassD1_mean_smallT <- array(round(apply(Wasserstein_distance_error_smallT, 2, mean), digits = 8), dim = c(1, 5))
colnames(WassD1_mean_smallT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
WassD1_mean_smallT_noGBM <- array(WassD1_mean_smallT[, -4], dim = c(1, 4))
colnames(WassD1_mean_smallT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

WassD2_mean_smallT <- array(round(apply(Wasserstein_distance_error_smallT_2, 2, mean), digits = 8), dim = c(1, 5))
colnames(WassD2_mean_smallT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
WassD2_mean_smallT_noGBM <- array(WassD2_mean_smallT[, -4], dim = c(1, 4))
colnames(WassD2_mean_smallT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

#############################
## running time in minutes ##
running_time_smallT_min <- round((running_time_smallT / 60), digits = 1) 
running_time_smallT_min_noGBM <- array(running_time_smallT_min[, -4], dim = c(1, 4))
colnames(running_time_smallT_min_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")