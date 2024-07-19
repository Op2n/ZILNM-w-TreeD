###############
## Libraries ##
library(ape)
library(mbImpute)
library(phyloMDA)
library(phyloseq)
library(zCompositions)
library(ZIPPCAlnm)
library(mvtnorm)
library(eBay)

###########################################
## load the True Abundance and OTU table ##
load("truAb_OTU_moderateT_1000times.RData")

#########################################################
## ZIPPCAlnm                                           ##
## Need to remove the columns and rows with all zeros  ##
ZIPPCAlnm_result_moderateT <- list()
time7 <- Sys.time() ## start time
for (i in 1: 1000){
  ZIPPCAlnm_result_moderateT[[i]] <- ZIPPCAlnm(X = y_noColRowzero_moderateT[[i]])
}
time8 <- Sys.time() ## End time
ZIPPCAlnm_time_moderateT <- difftime(time8, time7, units = "secs")#; ZIPPCAlnm_time_moderateT
# save the result
save(ZIPPCAlnm_result_moderateT, file = "ZIPPCAlnm_result_moderateT.RData")
save(ZIPPCAlnm_time_moderateT, file = "ZIPPCAlnm_time_moderateT.RData")