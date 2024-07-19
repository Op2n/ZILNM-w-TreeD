###############
## Libraries ##
library(ggplot2)
library(dplyr)

###################
## load the data ##

## Scenarions 2: 
load("running_time_moderateT.RData")
load("Frobenius_norm_error_moderateT.RData")
load("Simpson_index_mse_moderateT.RData")
load("Wasserstein_distance_error_moderateT.RData")
load("Frobenius_norm_error_moderateT_2.RData")
load("Simpson_index_mse_moderateT_2.RData")
load("Wasserstein_distance_error_moderateT_2.RData")

######################################
## mean of the Frobenius norm error ##
Frob1_mean_moderateT <- array(round(apply(Frobenius_norm_error_moderateT, 2, mean), digits = 8), dim = c(1, 5))
colnames(Frob1_mean_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Frob1_mean_moderateT_noGBM <- array(Frob1_mean_moderateT[, -4], dim = c(1, 4))
colnames(Frob1_mean_moderateT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

Frob2_mean_moderateT <- array(round(apply(Frobenius_norm_error_moderateT_2, 2, mean), digits = 8), dim = c(1, 5))
colnames(Frob2_mean_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Frob2_mean_moderateT_noGBM <- array(Frob2_mean_moderateT[, -4], dim = c(1, 4))
colnames(Frob2_mean_moderateT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

######################################
## mean of the mse of Simpson index ##
Simpson1_mean_moderateT <- array(round(apply(Simpson_index_mse_moderateT, 2, mean), digits = 8), dim = c(1, 5))
colnames(Simpson1_mean_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Simpson1_mean_moderateT_noGBM <- array(Simpson1_mean_moderateT[, -4], dim = c(1, 4))
colnames(Simpson1_mean_moderateT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

Simpson2_mean_moderateT <- array(round(apply(Simpson_index_mse_moderateT_2, 2, mean), digits = 8), dim = c(1, 5))
colnames(Simpson2_mean_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Simpson2_mean_moderateT_noGBM <- array(Simpson2_mean_moderateT[, -4], dim = c(1, 4))
colnames(Simpson2_mean_moderateT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

############################################
## mean of the Wasserstein distance error ##
WassD1_mean_moderateT <- array(round(apply(Wasserstein_distance_error_moderateT, 2, mean), digits = 8), dim = c(1, 5))
colnames(WassD1_mean_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
WassD1_mean_moderateT_noGBM <- array(WassD1_mean_moderateT[, -4], dim = c(1, 4))
colnames(WassD1_mean_moderateT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

WassD2_mean_moderateT <- array(round(apply(Wasserstein_distance_error_moderateT_2, 2, mean), digits = 8), dim = c(1, 5))
colnames(WassD2_mean_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
WassD2_mean_moderateT_noGBM <- array(WassD2_mean_moderateT[, -4], dim = c(1, 4))
colnames(WassD2_mean_moderateT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

#############################
## running time in minutes ##
running_time_moderateT_min <- round((running_time_moderateT / 60), digits = 1) 
running_time_moderateT_min_noGBM <- array(running_time_moderateT_min[, -4], dim = c(1, 4))
colnames(running_time_moderateT_min_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")


##########################
## Frobenius Norm Error ##
df_frb <- data.frame(
  Method = rep(c('zComp-SQ', 'mbImpute', 'phyloMDA-MIX'), each = 2),
  Tree = rep(c('Tree1', 'Tree2'), times = 5),
  Mean_time = c(time_zComp_SQ_tr1, time_zComp_SQ_tr2,
                time_zComp_GBM_tr1, time_zComp_GBM_tr2,
                time_mbImpute_tr1, time_mbImpute_tr2,
                time_phyloMDA_05_tr1, time_phyloMDA_05_tr2,
                time_phyloMDA_MIX_tr1, time_phyloMDA_MIX_tr2)
)









