##############################
## Working Directory Set Up ##
current_wd<- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_wd)

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
load("truAb_OTU_largeT_1000times.RData")

#########################################################
## ZIPPCAlnm                                           ##
## Need to remove the columns and rows with all zeros  ##
ZIPPCAlnm_result_largeT <- list()
time7 <- Sys.time() ## start time
for (i in 1: 1000){
  ZIPPCAlnm_result_largeT[[i]] <- ZIPPCAlnm(X = y_noColRowzero_largeT[[i]])
}
time8 <- Sys.time() ## End time
ZIPPCAlnm_time_largeT <- difftime(time8, time7, units = "secs")#; ZIPPCAlnm_time_largeT
# save the result
save(ZIPPCAlnm_result_largeT, file = "ZIPPCAlnm_result_largeT.RData")
save(ZIPPCAlnm_time_largeT, file = "ZIPPCAlnm_time_largeT.RData")