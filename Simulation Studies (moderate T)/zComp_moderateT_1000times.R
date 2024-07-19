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


#################################################################
## zComposition                                                ##
## Need to remove the columns with less than 2 positive value  ##

# SQ method
zComposition_SQ_result_moderateT <- list()
time1 <- Sys.time() ## start time
for (i in 1: 1000){
  zComposition_SQ_result_moderateT[[i]] <- as.matrix(cmultRepl(y_positive2value_moderateT[[i]], method= "SQ", z.warning = 1))
}
time2 <- Sys.time() ## End time
zComp_SQ_time_moderateT <- difftime(time2, time1, units = "secs")#; zComp_SQ_time_moderateT
# save the result
save(zComposition_SQ_result_moderateT, file = "zComp_SQ_result_moderateT.RData")
save(zComp_SQ_time_moderateT, file = "zComp_SQ_time_moderateT.RData")

# GBM method
zComposition_GBM_result_moderateT <- list()
time3 <- Sys.time() ## start time
for (i in 1: 1000){
  zComposition_GBM_result_moderateT[[i]] <- as.matrix(cmultRepl(y_positive2value_moderateT[[i]], method= "GBM", z.warning = 1))
}
time4 <- Sys.time() ## End time
zComp_GBM_time_moderateT <- difftime(time4, time3, units = "secs")#; zComp_GBM_time_moderateT
# save the result
save(zComposition_GBM_result_moderateT, file = "zComp_GBM_result_moderateT.RData")
save(zComp_GBM_time_moderateT, file = "zComp_GBM_time_moderateT.RData")