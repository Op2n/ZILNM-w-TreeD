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

##############
## mbImpute ##
mbImpute_result_smallT <- list()
time9 <- Sys.time() ## start time
for (i in 1: 1000){
  mbImpute_result_smallT[[i]] <- mbImpute(otu = y_smallT_list[[i]], D = D_phylo_combo2, 
                                          parallel = T, ncores = 4)
}
time10 <- Sys.time() ## End time
mbImpute_time_smallT <- difftime(time10, time9, units = "secs")#; mbImpute_time_smallT
# save the result
save(mbImpute_result_smallT, file = "mbImpute_result_smallT.RData")
save(mbImpute_time_smallT, file = "mbImpute_time_smallT.RData")

#################################
## Convert into the proportion ##
mbImpute_prop_result_smallT <- list()
for (l in 1: 1000){
  P <- mbImpute_result_smallT[[l]]$imp_count_mat_origlibsize
  for (i in 1:98){
    temp <- sum(P[i, ])
    for (j in 1:62){
      P[i, j] <- P[i, j] / temp
    }
  }
  mbImpute_prop_result_smallT[[l]] <- P
}
# save the result
save(mbImpute_prop_result_smallT, file = "mbImpute_prop_result_smallT.RData")