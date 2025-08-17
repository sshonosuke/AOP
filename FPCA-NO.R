###---------------------------------------------------------###
###         Code for Bayesian Functional Principle          ###
###    Components Analysis without Orthogonal Constraint    ###
###---------------------------------------------------------###
library(MASS)
library(MCMCpack)
library(splines)

## Input
# Y: (n, TT)-matrix of observed functions
# tt: grid of observed points (same for all the functional observations)
# K: number of PCA functions
# L: number of basis functions for modeling PCA functions


FPCA_no <- function(Y, tt, K=5, L=7, shrink=F, draw=2000, bn=1000){
  # preliminary 
  mc <- draw - bn    # number of posterior samples
  n <- dim(Y)[1]
  TT <- dim(Y)[2]
  Phi <- bs(tt, df=L, degree=3, intercept=T)
  Om <- t(Phi)%*%Phi    # correlation (inner product)
  
  # initial values 
  gam <- 10
  U <- matrix(1, L, K)    # local shrinkage parameter
  lam <- seq(1, 0, length=K+1)[1:K]
  Beta <- matrix(0, L, K)
  xi <- matrix(0, n, K)
  for(k in 1:K){
    xi[,k] <- 1-(k-1)/K
  }
  PC <- Phi%*%Beta
  Fn <- xi%*%t(PC)
  sig <- 1
  
  # objets for posterior samples
  FPC_pos <- array(NA, c(mc, TT, K))   
  Xi_pos <- array(NA, c(mc, n, K))
  
  # MCMC iterations 
  for(itr in 1:draw){
    # Beta & xi
    for(k in 1:K){
      sY <- Y - xi[,-k]%*%t(PC[,-k])
      vec_sY <- c(t(sY))
      
      # Beta
      xi_Phi <- kronecker(xi[,k], Phi)
      beta_B <- t(xi_Phi)%*%vec_sY/sig^2
      inv_beta_A <- t(xi_Phi)%*%xi_Phi/sig^2 + (1/gam)*diag(1/U[,k]) 
      beta_A <- solve( inv_beta_A )
      Beta[,k] <- mvrnorm(1, beta_A%*%beta_B, beta_A)
      pc <- c(Phi%*%Beta[,k])
      PC[,k] <- pc
      
      # xi
      A <- sum(pc^2)/sig^2 + rep(1/lam[k], n)
      B <- apply(t(sY)*pc, 2, sum)/sig^2
      xi[,k] <- rnorm(n, B/A, sqrt(1/A))
    }
    
    # gam & u
    if(shrink){
      gam <- rinvgamma(1, 1+L*K/2, 1+sum(Beta^2/U)/2)
      for(k in 1:K){
        V <- rinvgamma(L, 1, 1+1/U[,k])
        U[,k] <- rinvgamma(L, 1, 1/V+0.5*Beta[,k]^2/gam) 
      }
    }
    
    # lam 
    for(k in 2:K){
      lam[k] <- rinvgamma(1, 1+n/2, 1+sum(xi[,k]^2)/2)
      if(lam[k]>lam[k-1]){ lam[k] <- lam[k-1] }
    }
    
    # sigma 
    resid <- Y - xi%*%t(PC)
    sig <- sqrt( rinvgamma(1, 1+n*TT/2, 1+sum(resid^2)/2) )
    
    # save samples
    if(itr>bn){
      FPC_pos[itr-bn,,] <- PC
      Xi_pos[itr-bn,,] <- xi
    }
  }
  
  # output 
  Result <- list(FPC=FPC_pos, Xi=Xi_pos)
  return(Result)
}

