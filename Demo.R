###-------------------------------------------------------------###
###               Code for oneshot simulation of                ###
###      Bayesian functional principal component analysis       ###
###-------------------------------------------------------------###
rm(list=ls())

## load packages
library(Rcpp)
library(tidyverse)
library(reshape2)
library(GGally)
library(rootSolve)

source("FPCA-NO.R")
source("FPCA-AOP.R")
set.seed(123) 

scenario <- 1   # 1 or 2

## Generate time points
TT <- 30
tt <- seq(0, 1, length=TT)

## Define basis functions 
if(scenario==1){
  # Legendre polynomials
  phi_1 <- sqrt(3)*(2*tt-1)
  phi_2 <- sqrt(5)*(6*tt^2-6*tt+1)
  phi_3 <- sqrt(7)*(20*tt^3-30*tt^2+12*tt-1)
}
if(scenario==2){
  # Haar wavelet 
  phi_1 <- ifelse(tt<0.5, 1, -1)                  
  phi_2 <- sqrt(2)*ifelse(tt<0.25, 1, ifelse(tt<0.50, -1, 0))           
  phi_3 <- sqrt(2)*ifelse(tt<0.50, 0, ifelse(tt<0.75, 1, -1))  
}


## Number of samples
n <- 100

## Generate random scores
xi_1 <- rnorm(n, 0, 1)
xi_2 <- rnorm(n, 0, 0.7)
xi_3 <- rnorm(n, 0, 0.5)
sigma <- 1

# Generate functional data
Fn_true <- matrix(0, nrow=n, ncol=TT)
Y <- matrix(0, nrow=n, ncol=TT)
for(i in 1:n){
  Fn_true[i,] <- xi_1[i]*phi_1 + xi_2[i]*phi_2 + xi_3[i]*phi_3
  Y[i,] <- Fn_true[i,] + rnorm(TT, 0, sigma)
}


## Plot sample functions
matplot(tt, t(Y), type="l", lty=1, col=rainbow(n), xlab="t", 
        ylab="Y(t)", main="Simulated Functional Data")



## MCMC 
K <- 10
L <- 12
draw <- 2000     # length of MCMC
bn <- 1000       # length of burn^in

Fit <- list()
Fit[[1]] <- FPCA_no(Y=Y, tt=tt, K=K, L=L, shrink=F, draw=draw, bn=bn)
Fit[[2]] <- FPCA_no(Y=Y, tt=tt, K=K, L=L, shrink=T, draw=draw, bn=bn)
Fit[[3]] <- FPCA(Y=Y, tt=tt, K=K, L=L, constraint="IG-global", draw=draw, bn=bn)
Fit[[4]] <- FPCA(Y=Y, tt=tt, K=K, L=L, constraint="IG-local", draw=draw, bn=bn)
M <- length(Fit)

## Principal component functions
lower_bound <- 0.05
NC <- OG <- c()
hOm <- list()
for(m in 1:M){
  FPC <- apply(Fit[[m]]$FPC, c(2,3), mean)
  hOm[[m]] <- abs(t(FPC)%*%FPC)
  NC[m] <- sum( diag(hOm[[m]])>lower_bound )          # Number of components
  OG[m] <- sum( abs(hOm[[m]][upper.tri(hOm[[m]])]) )  # Orthogonality
}


## Heatmap of orthogonality
library(ggplot2)
library(gridExtra)

hOm_max <- max(hOm[[1]], hOm[[2]], hOm[[3]], hOm[[4]])
plot <- list()
for(j in 1:M){
  data <- data.frame(x=rep(1:K, each=K), y=rep(1:K, times=K), value=as.vector(hOm[[j]]))
  plot[[j]] <- ggplot(data, aes(x=x, y=y, fill=value)) +
    geom_tile() + scale_fill_gradient(low="#FAFAFA", high="red", limits=c(0,hOm_max)) +
    theme_minimal() + labs(title = "", x="Principle component", y="Principle component") +
    scale_x_continuous(breaks=1:K, labels=1:K) +  
    scale_y_continuous(breaks=1:K, labels=1:K) +  
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title=element_text(hjust=0.5, size=25),  
          axis.title=element_text(size=14),             
          legend.title=element_text(size=13),           
          legend.text=element_text(size=14)) 
}

grid.arrange(plot[[1]], plot[[3]], ncol=2)



## Plot of principal functions
par(mfcol=c(1,2))
# Standard (non-orthogonal) prior
ind <- order(diag(hOm[[1]]), decreasing=T)[1:3]
col <- lty <- lwd <- rep(1, K)
col[ind] <- 2:4
lty[-ind] <- 2
lwd[ind] <- 2
matplot(tt, apply(Fit[[1]]$FPC, c(2,3), mean), type="l", lty=lty, col=col, lwd=lwd, 
        xlab="", ylab="", main="NO (non-orthogonal prior)", ylim=c(-4, 4))
legend("bottomright", legend=c("1st principal function", "2nd principal function", "3rd principal function"), 
       col=2:4, lty=1, lwd=2)
# AOP prior
ind <- order(diag(hOm[[3]]), decreasing=T)[1:3]
col <- lty <- lwd <- rep(1, K)
col[ind] <- 2:4
lty[-ind] <- 2
lwd[ind] <- 2
matplot(tt, apply(Fit[[3]]$FPC, c(2,3), mean), type="l", lty=lty, col=col, lwd=lwd, 
        xlab="", ylab="", main="AOP (adaptive orthogonal prior)", ylim=c(-4, 4))
legend("bottomright", legend=c("1st principal function", "2nd principal function", "3rd principal function"), 
       col=2:4, lty=1, lwd=2)


