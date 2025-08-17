###--------------------------------------------------------###
###        Code for Bayesian Functional Principal          ###
###            Components Analysis with AOP                ###
###--------------------------------------------------------###
library(MASS)
library(MCMCpack)
library(splines)
library(progress)

## Input
# Y: (n, TT)-matrix of observed functions
# tt: grid of observed points (same for all the functional observations)
# K: number of PCA functions
# L: number of basis functions for modeling PCA functions
# constraint: "fix", "IG-global" or "IG-local"


FPCA <- function(Y, tt, K=5, L=7, draw=2000, bn=1000, constraint="IG-global", tau=0.1){
  # preliminary 
  mc <- draw - bn    # number of posterior samples
  n <- dim(Y)[1]
  TT <- dim(Y)[2]
  Phi <- bs(tt, df=L, degree=3, intercept=T)
  Om <- t(Phi)%*%Phi    # correlation (inner product)
  H <- list()
  for(k in 1:K){
    H[[k]] <- cbind(matrix(0, L-k+1, k-1), diag(L-k+1))
  }
  
  # parameter for orthogonality constraint
  tau <- rep(tau, K)
  tau[1] <- NA
  
  # initial values
  gam <- 10     # prior variance
  B0 <- gam
  lam <- seq(1, 0, length=K+1)[1:K]
  Beta <- matrix(0, L, K)
  xi <- matrix(0, n, K)
  for(k in 1:K){
    xi[,k] <- 1-(k-1)/K
  }
  PC <- Phi%*%Beta
  Fn <- xi%*%t(PC)
  sig <- 1
  
  # objects for posterior samples
  FPC_pos <- array(NA, c(mc, TT, K))   
  Xi_pos <- array(NA, c(mc, n, K))
  Sig_pos <- c()
  Tau_pos <- matrix(NA, mc, K)
  
  # MCMC iterations 
  pb <- progress_bar$new(total=draw)
  for(itr in 1:draw){
    # Beta & xi
    for(k in 1:K){
      sY <- Y - xi[,-k]%*%t(PC[,-k])
      vec_sY <- c(t(sY))
      
      # Beta
      xi_Phi <- kronecker(xi[,k], Phi)
      beta_B <- t(xi_Phi)%*%vec_sY/sig^2
      if(k==1){
        inv_beta_A <- t(xi_Phi)%*%xi_Phi/sig^2 + (1/B0)*diag(L) 
      }else{
        H_k <- H[[k]]
        inv_beta_A <- t(xi_Phi)%*%xi_Phi/sig^2 + t(H_k)%*%H_k/B0
      }
      for(j in 1:K){   # soft orthogonal constraint
        if(j<k){
          inv_beta_A <- inv_beta_A + Om%*%Beta[,j]%*%t(Beta[,j])%*%Om/tau[k]^2
        }
        if(j>k){
          inv_beta_A <- inv_beta_A + Om%*%Beta[,j]%*%t(Beta[,j])%*%Om/tau[j]^2
        }
      }
      beta_A <- solve( inv_beta_A )
      Beta[,k] <- mvrnorm(1, beta_A%*%beta_B, beta_A)
      pc <- c(Phi%*%Beta[,k])
      PC[,k] <- pc
      
      # xi
      A <- sum(pc^2)/sig^2 + rep(1/lam[k], n)
      B <- apply(t(sY)*pc, 2, sum)/sig^2
      xi[,k] <- rnorm(n, B/A, sqrt(1/A))
    }
    
    # tau (local parameter for orthogonality)
    if(constraint=="IG-global"){
      ss <- 0
      for(k in 2:K){
        for(j in 1:(k-1)){
          ss <- ss + (t(Beta[,j])%*%Om%*%Beta[,k])^2
        }
      }
      a_tau <- 3
      b_tau <- 2*(1/K)^2
      tau_new <- sqrt( rinvgamma(1, a_tau+K*(K-1)/4, b_tau+ss/2) )
      tau[2:K] <- tau_new
    }
    
    if(constraint=="IG-local"){
      for(k in 2:K){
        ss <- 0
        for(j in 1:(k-1)){
          ss <- ss + (t(Beta[,j])%*%Om%*%Beta[,k])^2
        }
        a_tau <- 3
        b_tau <- 2*(1/K)^2
        tau[k] <- sqrt( rinvgamma(1, a_tau+(k-1)/2, b_tau+ss/2) )
      }
    }
    
    # gamma 
    ss_gam <- 0
    for(k in 1:K){
      ss_gam <- ss_gam + sum((H[[k]]%*%Beta[,k])^2)
    }
    gam <- rinvgamma(1, 1+(2*L-K+1)*K/2, 1+ss_gam/2)
    B0 <- gam
    
    # lam 
    for(k in 2:K){
      lam[k] <- rinvgamma(1, 1+n/2, 1+sum(xi[,k]^2)/2)
      if(lam[k]>lam[k-1]){ lam[k] <- lam[k-1] }
    }
    
    # sigma (error variance)
    resid <- Y - xi%*%t(PC)
    sig <- sqrt( rinvgamma(1, 1+n*TT/2, 1+sum(resid^2)/2) )
    
    # save samples
    if(itr>bn){
      FPC_pos[itr-bn,,] <- PC
      Xi_pos[itr-bn,,] <- xi
      Sig_pos[itr-bn] <- sig
      Tau_pos[itr-bn,] <- tau
    }
    
    pb$tick()
  }
  
  # output 
  Result <- list(FPC=FPC_pos, Xi=Xi_pos, Sig=Sig_pos, Tau=Tau_pos)
  return(Result)
}

