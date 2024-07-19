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
## phyloMDA ##

###########################################################################
## Error Debug for phyloMDA: zidm methods 
## change the while loop condition from && to |
## remove the the if-else condition and remove case for number of taxa=2
## Error message will comes up when the number of taxa is large
## Modified function: 
est.zidm.EM_modified <- function(X, init.a=NULL, init.pi=NULL, iter=100, conv=1e-6) {
  print("modified function used")
  # for only binary
  Y <- X
  Y <- Y[rowSums(Y)>0, ]
  p <- ncol(Y)
  if(p > 2) {
    est.zidm <- function(Y, convergence=1e-4, iter = 1000) {
      f_t2 <- function(x,child,delta_e,A_e,B_e) {
        J=length(child)-1
        b<-c()
        for(j in 1:(J+1)){
          b<-c(b,x[j])
        }
        fx=0
        I=nrow(delta_e)
        for(i in 1:I){
          for(j in 1:J){
            fx=fx+(1-delta_e[i,child[j]])*(-log(beta(b[j],sum(b[(j+1):(J+1)])))+(b[j]-1)*A_e[i,child[j]]+(sum(b[(j+1):(J+1)])-1)*B_e[i,child[j]])
            #print(fx)
          }
        }
        -fx
      }
      
      g_t2 <- function(x,child,delta_e,A_e,B_e) {
        J=length(child)-1
        b<-c()
        for(j in 1:(J+1)){
          b<-c(b,x[j])
        }
        gx.all<-c()
        
        ### 1
        j=1
        gx=0
        I=nrow(delta_e)
        for(i in 1:I){
          gx=gx+(1-delta_e[i,child[j]])*(digamma(sum(b[j:(J+1)]))-digamma(b[j])+A_e[i,child[j]])
        }
        gx.all<-c(gx.all,gx)
        
        ### 2-J
        if(J>1){
          for(k in 2:J){
            gx=0
            for(j in 1:(k-1)){
              for(i in 1:I){
                gx=gx+(1-delta_e[i,child[j]])*(digamma(sum(b[j:(J+1)]))-digamma(sum(b[(j+1):(J+1)]))+B_e[i,child[j]])
              }
            }
            j=k
            for(i in 1:I){
              gx=gx+(1-delta_e[i,child[j]])*(digamma(sum(b[j:(J+1)]))-digamma(b[j])+A_e[i,child[j]])
            }
            gx.all<-c(gx.all,gx)
          }
        }
        ### J+1
        gx=0
        for(j in 1:J){
          for(i in 1:I){
            gx=gx+(1-delta_e[i,child[j]])*(digamma(sum(b[j:(J+1)]))-digamma(sum(b[(j+1):(J+1)]))+B_e[i,child[j]])
          }
        }
        gx.all<-c(gx.all,gx)
        -gx.all
      }
      tree <- ape::rtree(p)
      tree$edge <- cbind(p+1, 1:p)
      tree$Nnode <- 1
      tree$tip.label <- as.character(1:p)
      
      I=nrow(Y)
      M=nrow(tree$edge)
      leaf=length(tree$tip.label)
      
      pai_e=array(rep(0.5,I*M),c(I,M))
      A_e=array(rep(NA,I*M),c(I,M))
      B_e=array(rep(NA,I*M),c(I,M))
      delta_e=array(rep(NA,I*M),c(I,M))
      #delta_e=array(rep(0,I*M),c(I,M))
      a_e=array(rep(5,I*M),c(I,M))
      colnames(pai_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      colnames(A_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      colnames(B_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      colnames(delta_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      colnames(a_e)=paste(tree$edge[,1],tree$edge[,2],sep="->")
      
      Q_a=100
      Q_ab=1000
      num=0
      
      ###########################################################################################
      ## Error found here, change && to or !!!!
      while( (abs(Q_a-Q_ab)/abs(Q_ab))>conv | num<iter){
        
        # print((abs(Q_a-Q_ab)/abs(Q_ab)))
        # num=num+1
        # print(num)
        
        for(i in 1:I){
          for(j in (leaf + 1) : (leaf + tree$Nnode)){
            child=which(tree$edge[, 1] == j)
            for(c in child[1:(length(child)-1)]){
              
              if(Y[i,c]>0){
                delta_e[i,c]=0
              }else{
                mid=beta(a_e[i,c]+Y[i,c],sum(a_e[i,child[(which(child==c)+1):length(child)]]+
                                               Y[i,child[(which(child==c)+1):length(child)]]))/
                  beta(a_e[i,c],sum(a_e[i,child[(which(child==c)+1):length(child)]]))
                delta_e[i,c]=pai_e[i,c]/(pai_e[i,c]+(1-pai_e[i,c])*mid)
              }
            }
          }
        }
        
        for(i in 1:I){
          for(j in (leaf + 1) : (leaf + tree$Nnode)){
            child=which(tree$edge[, 1] == j)
            for(c in child[1:(length(child)-1)]){
              
              A_e[i,c]=digamma(a_e[i,c]+Y[i,c])-
                digamma(sum(a_e[i,child[which(child==c):length(child)]]+Y[i,child[which(child==c):length(child)]]))
              B_e[i,c]=digamma(sum(a_e[i,child[(which(child==c)+1):length(child)]]+Y[i,child[(which(child==c)+1):length(child)]]))-
                digamma(sum(a_e[i,child[which(child==c):length(child)]]+Y[i,child[which(child==c):length(child)]]))
            }
          }
        }
        
        # pai_e[1,]=constrOptim(rep(0.2,M), f_t1, g_t1, ui =rbind(diag(1,M,M),diag(-1,M,M)), ci=c(rep(0.000001,M),rep(-1,M)),
        #                       method="BFGS",tree=tree,delta_e=delta_e,outer.iterations = 100)$par
        
        pai_e[1,]<- c(colMeans(as.matrix(delta_e[,1:(M-1)])),0)
        pai_e[1,]=pai_e[1,]+1e-6
        for(i in 1:I){
          pai_e[i,]=pai_e[1,]
        }
        
        
        for(j in (leaf + 1) : (leaf + tree$Nnode)){
          child=which(tree$edge[, 1] == j)
          H=length(child)
          a_e[1,child]=constrOptim(rep(2,H),f_t2,grad=g_t2,ui = diag(1,H,H),ci=c(rep(0,H)),method="BFGS",
                                   child=child,delta_e=delta_e,A_e=A_e,B_e=B_e)$par
        }
        
        for(i in 2:I){
          a_e[i,]=a_e[1,]
        }
        
        Q_ab=Q_a
        Q_a=0
        
        for(i in 1:I){
          for(j in (leaf + 1) : (leaf + tree$Nnode)){
            child=which(tree$edge[, 1] == j)
            for(c in child[1:(length(child)-1)]){
              Q_a=Q_a+delta_e[i,c]*log(pai_e[i,c])+(1-delta_e[i,c])*log(1-pai_e[i,c])
              +(1-delta_e[i,c])*(-log(beta(a_e[i,c],sum(a_e[i,child[(which(child==c)+1):length(child)]])))+
                                   (a_e[i,c]-1)*A_e[i,c]+
                                   (sum(a_e[i,child[(which(child==c)+1):length(child)]])-1)*B_e[i,c])
            }
          }
        }
      }
      #print("done")
      return(list(alpha=a_e[1,], pi=pai_e[1,], loglik=Q_a))
    }
    res <- est.zidm(Y)
    return(res)
  }else{
    N  <- nrow(Y)
    yl <- Y[, 1]
    yr <- Y[, 2]
    zr <- which(yl==0)
    Ql <- function(a,b,A,B,d) {
      q1 <- sum(d*log(b)+(1-d)*log(1-b))
      q2 <- sum((1-d)*(-lbeta(a[1],a[2]) + (a[1]-1)*A + (a[2]-1)*B))
      q1+q2
    }
    
    # init
    b <- init.pi
    if(is.null(b))  b <-  length(zr)/N
    a <- init.a
    if(is.null(a)) {
      if(sum(yl>0)<2) {
        a <- c(1, 1)
      }else{
        a <-  mom.dm(Y)
      }
    }
    
    if(any(is.infinite(a))) a <- c(1, 1)
    lq <- loglik0 <- loglik_zidm(Y, a, b)
    #lq <- loglik0 <- -Inf
    conv0 <- 1
    delta<- rep(0, N)
    if(b<1e-5) {
      conv0 <- conv/10
      ta <- est.dm.NR(Y)
      a <- ta[[1]]; b <- 0; loglik <- ta[[2]]
    }
    if(any(is.na(a)) | b>0.95) { # can not estimate the value by EM
      conv0 <- conv/10
      ta <- est.dm.NR(Y)
      a <- c(0.5, 0.5); b <- 0; loglik <- -Inf
      if(!any(is.na(ta[[1]]))){
        a<- ta[[1]]; loglik <- ta[[2]]
      }
    }
    Iter <- 0
    ## Error found here, change and to or !!!!
    while(conv0 > conv | Iter < iter) {
      Iter <- Iter+1
      # e-step
      A <- digamma(a[1]+yl) - digamma(sum(a)+rowSums(Y))
      B <- digamma(a[2]+yr) - digamma(sum(a)+rowSums(Y))
      delta[zr] <- 1/(1+(1/b-1)*exp(lbeta(a[1], yr[zr]+a[2])-lbeta(a[1],a[2])))
      
      # m-step
      b <- mean(delta)
      da<- digamma(sum(a))-digamma(a)
      g <- c(N*(1-b)*da[1] + sum((1-delta)*A), N*(1-b)*da[2] + sum((1-delta)*B))
      ta<- trigamma(sum(a))
      H <- N*(1-b)*matrix(c(ta - trigamma(a[1]), ta, ta, ta - trigamma(a[2])),2,2)
      t1 <- try(H1 <- solve(H), silent = T)
      if(class(t1)[[1]]=='try-error') H1 <- matrix(c(1,0,0,1), 2)
      a <- a - as.numeric(H1%*%t(t(g)))
      a[a<0] <- 0.5
      loglik1 <- loglik_zidm(Y, a, b)
      #loglik1 <- Ql(a, b, A, B, delta)
      lq <- c(lq, loglik1)
      conv0 <- abs(loglik1-loglik0)/abs(loglik1)
      loglik0 <-loglik1
    }
    res <- list(alpha=a, pi=b, loglik=loglik0)
    return(res)
  }
}
environment(est.zidm.EM_modified) <- asNamespace("adaANCOM")
assignInNamespace("est.zidm.EM",est.zidm.EM_modified,ns="adaANCOM")

#############
## Fitting ##
phyloMDA_result_smallT <- list()
time5 <- Sys.time() ## start time
for (i in 1: 1000){
  phyloMDA_result_smallT[[i]] <- eBay_comps(otu_tab= y_smallT_list[[i]], 
                                            prior = "zero-inflated-Dirichlet-tree",
                                            tree= COMBO_tree2, 
                                            model= "MIX")
}
time6 <- Sys.time() ## End time
phyloMDA_time_smallT <- difftime(time6, time5, units = "secs")#; phyloMDA_time_smallT
# save the result
save(phyloMDA_result_smallT, file = "phyloMDA_result_smallT.RData")
save(phyloMDA_time_smallT, file = "phyloMDA_time_smallT.RData")
