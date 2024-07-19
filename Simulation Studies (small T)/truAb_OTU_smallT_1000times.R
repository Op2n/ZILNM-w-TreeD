# updated on May 28, 2024
# using Distance matrix to set up the P(structural zero)
# three choices of Temp, T= 5, 12, 30

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

########################
## Parameters Setting ##
# Real dataset: COMBO from phyloMDA package
COMBO_OTU<- as.matrix(t(combo.phyloseq.obj@otu_table@.Data))
COMBO_tree1<- phy_tree(combo.phyloseq.obj) ## Original edge length
COMBO_tree2<- COMBO_tree1
COMBO_tree2$edge.length<- rep(1, 122) ## Edge length equal to 1
# D: phylogenetic distance matrix
D_phylo_combo1<- cophenetic.phylo(COMBO_tree1) ## sum of all true edge lengths
D_phylo_combo2<- cophenetic.phylo(COMBO_tree2) ## sum of the number of edges

# number of taxa
K <- dim(D_phylo_combo2)[1]

# total count for each of n samples
nraw <- dim(COMBO_OTU)[1]
ni_combo_raw <- array(apply(COMBO_OTU, 1, sum), c(nraw, 1))

REP <- 1
n <- REP * nraw
ni_combo <- rep(ni_combo_raw, REP)

###############################
##  OTU count data Generator ##

###############################################################################
## Logit-normal-distribution function based on log-ratio transformation ##
###############################################################################
mu_prob<- function(x){ # input are L - 1 ratios and output is L \pi's please noted the sum of L ratios equal to 1
  pl<- length(x) 
  tmp<- sum(exp(x))
  ans<- c(exp(x[1: pl])/ (1+ tmp), 1/ (1+ tmp))
}

true_abundance<- function(sigma, rho, D_phylo, n, K, Ni, Temp, seednum = 123, num_replicate){
  
  set.seed(seednum)
  
  # larger T, smaller P(biological zero)
  # PIMS VXML choices: T = 5, 15, 50
  p_delta1<- rep(0, K) # using D_phylo to model the structural zero probability
  # K is reserved for the taxon that must exist
  delta_param<- rep(NA, K)
  for (j in 1:K ){
    p_delta1[j]<- 1 - exp(-(D_phylo[j, K]) /Temp ) # d=0 P(biological 0) = 1, larger d, larger P(biological 0)
    delta_param[j]<- rbinom(1, 1, p_delta1[j])
  }
  
  
  ##############################################################################
  ## Define var-cov matrix with the phylogenentic tree information            ##
  ## rho: evolution rate; D_phy: phylogenentic tree distance matrix           ##
  ##############################################################################
  var_cov_matrix<- sigma^2* exp(-2* rho* D_phylo)
  
  ### sample theta
  
  theta <- rmvnorm(1, mean= rep(0, K) , sigma = 10 * var_cov_matrix)
  
  ######################################
  ## Generate the true abundance 
  ######################################
  true_pi<- array(0, c(n, K))
  otu <- array(0, c(n, K))
  colnames(true_pi)<- colnames(otu)<- colnames(D_phylo)
  ind<- which(delta_param[-K] == 0) # 0 for existing species including the reference taxon (labeled by K)
  
  ## testing purpose: check whether the true abundance remains the same
  ## true_abundance_list<- list() 
  set.seed(seednum)
  for (i in 1:n) {
    ## ONLY FOR EXISTING TAXA
    U_ind<- rmvnorm(1, mean= theta[ind] , sigma= var_cov_matrix[ind, ind])
    true_pi[i, c(ind, K)]<- mu_prob(U_ind)
  }
  
  otu_list<- list()  ## This step is using for simulation studies, the 
  for (l in 1:num_replicate){
    set.seed(l+seednum)
    for (i in 1:n){
      ## log(True abundance ratio) ~ multivariate normal
      otu[i, c(ind, K)] <- rmultinom (1, Ni[i], true_pi[i, c(ind, K)])
    }
    otu_list[[l]]<- otu
    ## testing purpose: check whether the true abundance remains the same
    ## true_abundance_list[[l]]<- true_pi
  }
  
  #############################
  ## output                  ##
  #############################
  results<- list(p_delta1= p_delta1, 
                 delta_param= delta_param, 
                 time= time, theta = theta,
                 var_cov_matrix= var_cov_matrix, 
                 true_pi= true_pi, true_otu = otu_list)
  return(results)
}


#######################
## Scenarios 3: small T
simulated_combo_data2_smallT <- true_abundance(sigma = 1, rho = 1, D_phylo = D_phylo_combo2, 
                                               n= n, K = 62, Ni = ni_combo, Temp = 5, num_replicate = 1000)
# Save the simulated result
# save(simulated_combo_data2_smallT, file = "simulated_combo_data2_smallT.RData")

# True abundance
true_pi_smallT <- simulated_combo_data2_smallT$true_pi
# save the true abundance table
# save(true_pi_smallT, file = "true_pi_smallT.RData")

# OTU counts list
y_smallT_list <- simulated_combo_data2_smallT$true_otu
# save the OTU counts table
# save(y_smallT_list, file = "y_smallT_list.RData")

# Proportion of zeros
prop_zero_smallT <- c()
for (i in 1: 1000){
  prop_zero_smallT[i] <- signif(mean(y_smallT_list[[i]] == 0), digits = 6)
}
biological_smallT <- simulated_combo_data2_smallT$delta_param
prop_bio_zero_smallT<- signif(mean(biological_smallT == 1), digits = 6)
# save the proportion of zeros
zero_smallT <- list(prop_zero_smallT, prop_bio_zero_smallT)
# save(zero, file = "zero.RData")

# Remove columns of the OTU tables with all zeros
y_noColzero_smallT <- list()
# Remove columns of the true abundance with all zeros
true_pi_noColzero_smallT <- list()
for (i in 1: 1000){
  X <- y_smallT_list[[i]]
  Z1 <- true_pi_smallT
  # remove the columns with all zero
  zerocol <- which(colSums(X)==0)
  if(length(zerocol) >0 ){
    X <- X[, -zerocol]
    Z1 <- Z1[, -zerocol]
  }
  y_noColzero_smallT[[i]] <- X
  true_pi_noColzero_smallT[[i]] <- Z1
}
# save the OTU counts table with no columns with all zeros
# save(y_noColzero_smallT, file = "y_noColzero_smallT.RData")
# save the true abundance table with no columns with all zeros
# save(true_pi_noColzero_smallT, file = "true_pi_noColzero_smallT.RData")


# Remove columns and rows of the OTU tables with all zeros
y_noColRowzero_smallT <- list()
true_pi_noColRowzero_smallT <- list()
for (i in 1: 1000){
  Y <- y_smallT_list[[i]]
  Z2 <- true_pi_smallT
  # remove the rows with all zero
  zerorow <- which(rowSums(Y)==0)
  if(length(zerorow) >0 ){
    Y <- Y[-zerorow, ]
    Z2 <- Z2[-zerorow, ]
  }
  # remove the columns with all zero
  zerocol2 <- which(colSums(Y)==0)
  if(length(zerocol2) >0 ){
    Y <- Y[, -zerocol2]
    Z2 <- Z2[, -zerocol2]
  }
  y_noColRowzero_smallT[[i]] <- Y
  true_pi_noColRowzero_smallT[[i]] <- Z2
}

# save the OTU counts table with no columns or rows with all zeros
# save(y_noColRowzero_smallT, file = "y_noColRowzero_smallT.RData")
# save the true abundance table with no columns or rows with all zeros
# save(true_pi_noColRowzero_smallT, file = "true_pi_noColRowzero_smallT.RData")

# Remove columns with less than two positive values (for zComposition package)
y_positive2value_smallT <- list()
true_pi_positive2value_smallT <- list()
for (i in 1:1000){
  X2 <- y_smallT_list[[i]]
  R <- X2 > 0
  positive_counts <- colSums(R)
  y_positive2value_smallT[[i]] <- X2[, positive_counts >= 2]
  true_pi_positive2value_smallT[[i]] <- true_pi_smallT[, positive_counts >= 2]
}

# Remove those unnecessary information
rm(X, X2, Y, Z1, Z2, R, positive_counts, current_wd)

# save the result for the whole file
save.image(file = "truAb_OTU_smallT_1000times.RData")