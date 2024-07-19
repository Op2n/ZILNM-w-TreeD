##################
## Loading data ##
load("truAb_OTU_moderateT_1000times.RData")
load("mbImpute_prop_result_moderateT.RData")
load("mbImpute_time_moderateT.RData")
load("phyloMDA_result_moderateT.RData")
load("phyloMDA_time_moderateT.RData")
load("zComp_GBM_result_moderateT.RData")
load("zComp_GBM_time_moderateT.RData")
load("zComp_SQ_result_moderateT.RData")
load("zComp_SQ_time_moderateT.RData")
load("ZIPPCAlnm_result_moderateT.RData")
load("ZIPPCAlnm_time_moderateT.RData")


##################
## Running time ##
mbImpute_time_moderateT; phyloMDA_time_moderateT; zComp_SQ_time_moderateT; zComp_GBM_time_moderateT; ZIPPCAlnm_time_moderateT
## time (list)
running_time_moderateT <- array(c(17799.16, 1604.231, 37.77921, 38.81938, 9854.954), dim = c(1, 5)) # unit: seconds
colnames(running_time_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")


#########################
## Evaluation function ##
## Note: the final step of the Wasserstein distance calculation is dividing the number of taxa
## For mbImpute: The columns of all zeros are including in the final results, so the number of taxa will be varied in the Wasserstein distance comparison
## For phyloMDA: The columns of all zeros are imputed, so the number of taxa will be varied in the Wasserstein distance comparison

## Function 1: 
evaluate_function <- function(est_value, true_ab){
  
  ## number of samples
  n <- dim(est_value)[1]
  ## number of taxa
  k<- dim(est_value)[2]
  
  ## Frobenius norm error
  Frobenius_norm_error <- sqrt(sum((est_value - true_ab)^2))
  
  ## Mean square error of the Simpson's index
  Simpson_index_mse <- (sum((rowSums(true_ab^2)- rowSums((est_value)^2))^2))/ n
  
  ## Wasserstein Distance error
  ## mean of each taxon
  true_ab_mean<- (apply(true_ab, 2, sum))/n
  est_value_mean<- (apply(est_value, 2, sum))/n
  ## standard deviation of each taxon
  true_ab_sd<- sqrt(apply((true_ab- true_ab_mean)^2, 2, sum)/ (n-1))
  est_value_sd<- sqrt(apply((est_value- est_value_mean)^2, 2, sum)/ (n-1))
  ## mean/ sd in order statstics
  true_ab_mean_sd<- sort(true_ab_mean/ true_ab_sd)
  est_value_mean_sd<- sort(est_value_mean/ est_value_sd)
  ## The Wasserstein distance error between the distribution of mean/sd
  Wasserstein_distance_error<- sum(abs(true_ab_mean_sd- est_value_mean_sd)) / k
  
  ## Results
  evaluation_results<- list(Frobenius_norm_error= Frobenius_norm_error, 
                            Simpson_index_mse= Simpson_index_mse, 
                            Wasserstein_distance_error= Wasserstein_distance_error)
  
}

## Function 2: only for mbImpute and phyloMDA to remove the columns with all zeros
evaluate_mbImpute_phyloMDA <- function(est_value, true_ab){
  
  ## index: no column with all zeros in the true pi
  index1 <- which(colSums(true_ab) != 0)
  est_value <- est_value[, index1]
  true_ab <- true_ab[, index1]
  
  ## number of samples
  n <- dim(est_value)[1]
  ## number of taxa
  k<- dim(est_value)[2]
  
  ## Frobenius norm error
  Frobenius_norm_error <- sqrt(sum((est_value - true_ab)^2))
  
  ## Mean square error of the Simpson's index
  Simpson_index_mse <- (sum((rowSums(true_ab^2)- rowSums((est_value)^2))^2))/ n
  
  ## Wasserstein Distance error
  ## mean of each taxon
  true_ab_mean<- (apply(true_ab, 2, sum))/n
  est_value_mean<- (apply(est_value, 2, sum))/n
  ## standard deviation of each taxon
  true_ab_sd<- sqrt(apply((true_ab- true_ab_mean)^2, 2, sum)/ (n-1))
  est_value_sd<- sqrt(apply((est_value- est_value_mean)^2, 2, sum)/ (n-1))
  ## mean/ sd in order statstics
  true_ab_mean_sd<- sort(true_ab_mean/ true_ab_sd)
  est_value_mean_sd<- sort(est_value_mean/ est_value_sd)
  ## The Wasserstein distance error between the distribution of mean/sd
  Wasserstein_distance_error<- sum(abs(true_ab_mean_sd - est_value_mean_sd)) / k
  
  ## Results
  evaluation_results<- list(Frobenius_norm_error= Frobenius_norm_error, 
                            Simpson_index_mse= Simpson_index_mse, 
                            Wasserstein_distance_error= Wasserstein_distance_error)
}


#######################
## Evaluation result ##

## Function 1: 
Frobenius_norm_error_moderateT <- array(NA, dim = c(1000, 5))
colnames(Frobenius_norm_error_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Simpson_index_mse_moderateT <- array(NA, dim = c(1000, 5))
colnames(Simpson_index_mse_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Wasserstein_distance_error_moderateT <- array(NA, dim = c(1000, 5))
colnames(Wasserstein_distance_error_moderateT) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
for (l in 1:1000){
  mbImpute_compare <- evaluate_function(est_value = mbImpute_prop_result_moderateT[[l]], 
                                        true_ab = true_pi_moderateT)
  phyloMDA_compare <- evaluate_function(est_value = phyloMDA_result_moderateT[[l]], 
                                        true_ab = true_pi_moderateT)
  zComp_SQ_compare <- evaluate_function(est_value = zComposition_SQ_result_moderateT[[l]], 
                                        true_ab = true_pi_positive2value_moderateT[[l]])
  zComp_GBM_compare <- evaluate_function(est_value = zComposition_GBM_result_moderateT[[l]], 
                                         true_ab = true_pi_positive2value_moderateT[[l]])
  ZIPPCAlnm_compare <- evaluate_function(est_value = (ZIPPCAlnm_result_moderateT[[l]])$Q, 
                                         true_ab = true_pi_noColRowzero_moderateT[[l]])
  Frobenius_norm_error_moderateT[l, 1] <- mbImpute_compare$Frobenius_norm_error
  Frobenius_norm_error_moderateT[l, 2] <- phyloMDA_compare$Frobenius_norm_error
  Frobenius_norm_error_moderateT[l, 3] <- zComp_SQ_compare$Frobenius_norm_error
  Frobenius_norm_error_moderateT[l, 4] <- zComp_GBM_compare$Frobenius_norm_error
  Frobenius_norm_error_moderateT[l, 5] <- ZIPPCAlnm_compare$Frobenius_norm_error
  Simpson_index_mse_moderateT[l, 1] <- mbImpute_compare$Simpson_index_mse
  Simpson_index_mse_moderateT[l, 2] <- phyloMDA_compare$Simpson_index_mse
  Simpson_index_mse_moderateT[l, 3] <- zComp_SQ_compare$Simpson_index_mse
  Simpson_index_mse_moderateT[l, 4] <- zComp_GBM_compare$Simpson_index_mse
  Simpson_index_mse_moderateT[l, 5] <- ZIPPCAlnm_compare$Simpson_index_mse
  Wasserstein_distance_error_moderateT[l, 1] <- mbImpute_compare$Wasserstein_distance_error
  Wasserstein_distance_error_moderateT[l, 2] <- phyloMDA_compare$Wasserstein_distance_error
  Wasserstein_distance_error_moderateT[l, 3] <- zComp_SQ_compare$Wasserstein_distance_error
  Wasserstein_distance_error_moderateT[l, 4] <- zComp_GBM_compare$Wasserstein_distance_error
  Wasserstein_distance_error_moderateT[l, 5] <- ZIPPCAlnm_compare$Wasserstein_distance_error
}


## Function 1 and Function 2: 
Frobenius_norm_error_moderateT_2 <- array(NA, dim = c(1000, 5))
colnames(Frobenius_norm_error_moderateT_2) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Simpson_index_mse_moderateT_2 <- array(NA, dim = c(1000, 5))
colnames(Simpson_index_mse_moderateT_2) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
Wasserstein_distance_error_moderateT_2 <- array(NA, dim = c(1000, 5))
colnames(Wasserstein_distance_error_moderateT_2) <- c("mbImpute", "phyloMDA", "zComp-SQ", "zComp-GBM", "ZIPPCAlnm")
for (l in 1:1000){
  mbImpute_compare_2 <- evaluate_mbImpute_phyloMDA(est_value = mbImpute_prop_result_moderateT[[l]], 
                                                   true_ab = true_pi_moderateT)
  phyloMDA_compare_2 <- evaluate_mbImpute_phyloMDA(est_value = phyloMDA_result_moderateT[[l]], 
                                                   true_ab = true_pi_moderateT)
  zComp_SQ_compare_2 <- evaluate_function(est_value = zComposition_SQ_result_moderateT[[l]], 
                                          true_ab = true_pi_positive2value_moderateT[[l]])
  zComp_GBM_compare_2 <- evaluate_function(est_value = zComposition_GBM_result_moderateT[[l]], 
                                           true_ab = true_pi_positive2value_moderateT[[l]])
  ZIPPCAlnm_compare_2 <- evaluate_function(est_value = (ZIPPCAlnm_result_moderateT[[l]])$Q, 
                                           true_ab = true_pi_noColRowzero_moderateT[[l]])
  Frobenius_norm_error_moderateT_2[l, 1] <- mbImpute_compare_2$Frobenius_norm_error
  Frobenius_norm_error_moderateT_2[l, 2] <- phyloMDA_compare_2$Frobenius_norm_error
  Frobenius_norm_error_moderateT_2[l, 3] <- zComp_SQ_compare_2$Frobenius_norm_error
  Frobenius_norm_error_moderateT_2[l, 4] <- zComp_GBM_compare_2$Frobenius_norm_error
  Frobenius_norm_error_moderateT_2[l, 5] <- ZIPPCAlnm_compare_2$Frobenius_norm_error
  Simpson_index_mse_moderateT_2[l, 1] <- mbImpute_compare_2$Simpson_index_mse
  Simpson_index_mse_moderateT_2[l, 2] <- phyloMDA_compare_2$Simpson_index_mse
  Simpson_index_mse_moderateT_2[l, 3] <- zComp_SQ_compare_2$Simpson_index_mse
  Simpson_index_mse_moderateT_2[l, 4] <- zComp_GBM_compare_2$Simpson_index_mse
  Simpson_index_mse_moderateT_2[l, 5] <- ZIPPCAlnm_compare_2$Simpson_index_mse
  Wasserstein_distance_error_moderateT_2[l, 1] <- mbImpute_compare_2$Wasserstein_distance_error
  Wasserstein_distance_error_moderateT_2[l, 2] <- phyloMDA_compare_2$Wasserstein_distance_error
  Wasserstein_distance_error_moderateT_2[l, 3] <- zComp_SQ_compare_2$Wasserstein_distance_error
  Wasserstein_distance_error_moderateT_2[l, 4] <- zComp_GBM_compare_2$Wasserstein_distance_error
  Wasserstein_distance_error_moderateT_2[l, 5] <- ZIPPCAlnm_compare_2$Wasserstein_distance_error
}


#####################
## save the result ##
save(running_time_moderateT, file = "running_time_moderateT.RData")
save(Frobenius_norm_error_moderateT, file = "Frobenius_norm_error_moderateT.RData")
save(Simpson_index_mse_moderateT, file = "Simpson_index_mse_moderateT.RData")
save(Wasserstein_distance_error_moderateT, file = "Wasserstein_distance_error_moderateT.RData")
save(Frobenius_norm_error_moderateT_2, file = "Frobenius_norm_error_moderateT_2.RData")
save(Simpson_index_mse_moderateT_2, file = "Simpson_index_mse_moderateT_2.RData")
save(Wasserstein_distance_error_moderateT_2, file = "Wasserstein_distance_error_moderateT_2.RData")
