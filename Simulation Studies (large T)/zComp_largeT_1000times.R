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


#################################################################
## zComposition                                                ##
## Need to remove the columns with less than 2 positive value  ##

# SQ method
zComposition_SQ_result_largeT <- list()
time1 <- Sys.time() ## start time
for (i in 1: 1000){
  zComposition_SQ_result_largeT[[i]] <- as.matrix(cmultRepl(y_positive2value_largeT[[i]], method= "SQ", z.warning = 1))
}
time2 <- Sys.time() ## End time
zComp_SQ_time_largeT <- difftime(time2, time1, units = "secs")#; zComp_SQ_time_largeT
# save the result
save(zComposition_SQ_result_largeT, file = "zComp_SQ_result_largeT.RData")
save(zComp_SQ_time_largeT, file = "zComp_SQ_time_largeT.RData")

# GBM method
zComposition_GBM_result_largeT <- list()
time3 <- Sys.time() ## start time
for (i in 1: 1000){
  zComposition_GBM_result_largeT[[i]] <- as.matrix(cmultRepl(y_positive2value_largeT[[i]], method= "GBM", z.warning = 1))
}
time4 <- Sys.time() ## End time
zComp_GBM_time_largeT <- difftime(time4, time3, units = "secs")#; zComp_GBM_time_largeT
# save the result
save(zComposition_GBM_result_largeT, file = "zComp_GBM_result_largeT.RData")
save(zComp_GBM_time_largeT, file = "zComp_GBM_time_largeT.RData")

