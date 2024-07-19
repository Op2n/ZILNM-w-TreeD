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

##############
## mbImpute ##
mbImpute_result_largeT <- list()
time9 <- Sys.time() ## start time
for (i in 1: 1000){
  mbImpute_result_largeT[[i]] <- mbImpute(otu = y_largeT_list[[i]], D = D_phylo_combo2, 
                                          parallel = T, ncores = 8)
}
time10 <- Sys.time() ## End time
mbImpute_time_largeT <- difftime(time10, time9, units = "secs")#; mbImpute_time_largeT
# save the result
save(mbImpute_result_largeT, file = "mbImpute_result_largeT.RData")
save(mbImpute_time_largeT, file = "mbImpute_time_largeT.RData")

#################################
## Convert into the proportion ##
mbImpute_prop_result_largeT <- list()
for (l in 1: 1000){
  P <- mbImpute_result_largeT[[l]]$imp_count_mat_origlibsize
  for (i in 1:98){
    temp <- sum(P[i, ])
    for (j in 1:62){
      P[i, j] <- P[i, j] / temp
    }
  }
  mbImpute_prop_result_largeT[[l]] <- P
}
# save the result
save(mbImpute_prop_result_largeT, file = "mbImpute_prop_result_largeT.RData")