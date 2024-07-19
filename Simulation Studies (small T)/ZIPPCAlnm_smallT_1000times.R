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
load("truAb_OTU_smallT_1000times.RData")

#########################################################
## ZIPPCAlnm                                           ##
## Need to remove the columns and rows with all zeros  ##
ZIPPCAlnm_result_smallT <- list()
time7 <- Sys.time() ## start time
for (i in 1: 1000){
  ZIPPCAlnm_result_smallT[[i]] <- ZIPPCAlnm(X = y_noColRowzero_smallT[[i]])
}
time8 <- Sys.time() ## End time
ZIPPCAlnm_time_smallT <- difftime(time8, time7, units = "secs")#; ZIPPCAlnm_time_smallT
# save the result
save(ZIPPCAlnm_result_smallT, file = "ZIPPCAlnm_result_smallT.RData")
save(ZIPPCAlnm_time_smallT, file = "ZIPPCAlnm_time_smallT.RData")