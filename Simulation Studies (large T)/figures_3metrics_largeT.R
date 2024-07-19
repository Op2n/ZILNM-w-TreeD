###############
## Libraries ##
library(ggplot2)
library(dplyr)

###################
## load the data ##

## Scenarios 1: large T
load("running_time_largeT.RData")
load("Frobenius_norm_error_largeT.RData")
load("Simpson_index_mse_largeT.RData")
load("Wasserstein_distance_error_largeT.RData")
load("Frobenius_norm_error_largeT_2.RData")
load("Simpson_index_mse_largeT_2.RData")
load("Wasserstein_distance_error_largeT_2.RData")

######################################
## mean of the Frobenius norm error ##
Frob1_mean_largeT <- array(round(apply(Frobenius_norm_error_largeT, 2, mean), digits = 8), dim = c(1, 5))
colnames(Frob1_mean_largeT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Frob1_mean_largeT_noGBM <- array(Frob1_mean_largeT[, -4], dim = c(1, 4))
colnames(Frob1_mean_largeT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

Frob2_mean_largeT <- array(round(apply(Frobenius_norm_error_largeT_2, 2, mean), digits = 8), dim = c(1, 5))
colnames(Frob2_mean_largeT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Frob2_mean_largeT_noGBM <- array(Frob2_mean_largeT[, -4], dim = c(1, 4))
colnames(Frob2_mean_largeT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

######################################
## mean of the mse of Simpson index ##
Simpson1_mean_largeT <- array(round(apply(Simpson_index_mse_largeT, 2, mean), digits = 8), dim = c(1, 5))
colnames(Simpson1_mean_largeT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Simpson1_mean_largeT_noGBM <- array(Simpson1_mean_largeT[, -4], dim = c(1, 4))
colnames(Simpson1_mean_largeT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

Simpson2_mean_largeT <- array(round(apply(Simpson_index_mse_largeT_2, 2, mean), digits = 8), dim = c(1, 5))
colnames(Simpson2_mean_largeT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Simpson2_mean_largeT_noGBM <- array(Simpson2_mean_largeT[, -4], dim = c(1, 4))
colnames(Simpson2_mean_largeT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

############################################
## mean of the Wasserstein distance error ##
WassD1_mean_largeT <- array(round(apply(Wasserstein_distance_error_largeT, 2, mean), digits = 8), dim = c(1, 5))
colnames(WassD1_mean_largeT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
WassD1_mean_largeT_noGBM <- array(WassD1_mean_largeT[, -4], dim = c(1, 4))
colnames(WassD1_mean_largeT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

WassD2_mean_largeT <- array(round(apply(Wasserstein_distance_error_largeT_2, 2, mean), digits = 8), dim = c(1, 5))
colnames(WassD2_mean_largeT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
WassD2_mean_largeT_noGBM <- array(WassD2_mean_largeT[, -4], dim = c(1, 4))
colnames(WassD2_mean_largeT_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")

#############################
## running time in minutes ##
running_time_largeT_min <- round((running_time_largeT / 60), digits = 1) 
running_time_largeT_min_noGBM <- array(running_time_largeT_min[, -4], dim = c(1, 4))
colnames(running_time_largeT_min_noGBM) <- c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm")


###########
## Graph ##

##########################
## Frobenius Norm Error ##

# Create the data frame
df_frb_noGBM <- data.frame(
  Method = rep(c("mbImpute", "phyloMDA", "zComp-SQ", "ZIPPCAlnm"), each = 3),
  Eval_metrics = rep(c("Frobenius Norm Error", "MSE of Simpson index", "Wasserstein distance error"), times = 4),
  Frobenius1_mean_largeT_noGBM = c(Frob1_mean_largeT_noGBM[1, 1], Simpson1_mean_largeT_noGBM[1, 1], WassD1_mean_largeT_noGBM[1, 1],
                                   Frob1_mean_largeT_noGBM[1, 2], Simpson1_mean_largeT_noGBM[1, 2], WassD1_mean_largeT_noGBM[1, 2],
                                   Frob1_mean_largeT_noGBM[1, 3], Simpson1_mean_largeT_noGBM[1, 3], WassD1_mean_largeT_noGBM[1, 3],
                                   Frob1_mean_largeT_noGBM[1, 4], Simpson1_mean_largeT_noGBM[1, 4], WassD1_mean_largeT_noGBM[1, 4])
)

# Plotting with ggplot2
ggplot(df_frb_noGBM, aes(x=Method, y=Frobenius1_mean_largeT_noGBM, fill=Eval_metrics)) + 
  geom_bar(stat='identity', position=position_dodge()) + 
  geom_text(aes(label=Frobenius1_mean_largeT_noGBM), position=position_dodge(width=1), vjust=-0.5) + 
  theme_minimal() +
  labs(title="Large T scenario", x="Methods", y="", fill="Eval_metrics") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(hjust = 0.5), 
        axis.title.y = element_text(hjust = 0.5))


























