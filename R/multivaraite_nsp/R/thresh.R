
#-------------------------------------------------------------------------------
#
# Functions for obtaining thresholds for panel NSP via Monte Carlo simulations
#
#--------------------------------------------------------------------------------

library(docstring)
library(foreach)
library(parallel)
library(doParallel)
library(mvtnorm)
library(tcltk)



multi.sup.norms <- function(E) {
  
  #' Multi-Sup Norms
  #' 
  #' Calculates $\ell_2$ and $\ell_\infty$ aggregated multi-sup norms for panel E. 
  #'
  #'@param E matrix (n x T) of noise terms
  
  
  tt <- ncol(E)
  
  n <- nrow(E)
  
  
  E.cumsum <- cbind(rep(0,n), t(apply(E,1,cumsum)))
  
  run.L2 <- rep(0,n)
  
  run.L.infty <- 0 
  
  
  for (i in 0:(tt-1)) for (j in (i+1):tt) {
    
    
    region <- abs(E.cumsum[,j+1] - E.cumsum[,i+1]) / sqrt(j - i)
    
    
    local.L.infty <- max(region)
    
    if (local.L.infty > run.L.infty) run.L.infty <- local.L.infty
    
    
    run.L2 <- run.L2 * (run.L2 > region) + region * (run.L2 <= region) 
    
  }
  
  c(run.L.infty, sum(run.L2 ** 2))
  
} 



simulate.panel.max <- function(n, tt, seed = 42, K = 1000, do.par = TRUE) {
  
  #' Simulate Panel Maxima
  #'
  #' Simulate $\ell_2$ and $\ell_\infty$ maxima of panel with standard Gaussian entries. Called by the function `MC.thresh`. 
  #' 
  #' @param n int, number of channels
  #' @param tt int, number of time points
  #' @param seed int
  #' @param K int, number of Monte Carlo draws
  #' @param do.par bool, whether calculations are carried out in paralell
  
  
  set.seed(seed)
  
  
  
  if (do.par) {
    
    
    maxima <- foreach(k=1:K, .inorder = FALSE, .export = c("multi.sup.norms", "rmvnorm","tkProgressBar","setTkProgressBar")) %dopar% {
      
      if(!exists("pb")) pb <- tkProgressBar("para pb", min=1, max=K)

      setTkProgressBar(pb, k)
      
      E <- rmvnorm(n, mean = rep(0,tt), sigma = diag(tt))
      
      multi.sup.norms(E)
      
    }
    
    
  } else {
    
    pb <- txtProgressBar(min = 1, max = K)
    
    maxima <- foreach(k=1:K) %do% {
      
      E <- rmvnorm(n, mean = rep(0,tt), sigma = diag(tt))
      
      setTxtProgressBar(pb, k)
      
      multi.sup.norms(E)
      
    }
    
  }
  
  
  max.sim <- matrix(unlist(maxima), ncol = 2, byrow = TRUE)
  
  if (!dir.exists("../temp")) dir.create("../temp")
  
  write.csv(max.sim, file.path("../temp",paste("n",n,"T",tt,"K",K,".csv", sep = "_")))
  
  
  return(max.sim)
  
}




MC.thresh <- function(n, tt, alpha, seed = 42, K = 1000, do.par = FALSE) {
  
  #' Monte Carlo Threshold
  #'
  #' Generates thresholds which controll the global type I error 
  #'
  #' @param n int, number of channels
  #' @param tt int, number of time points
  #' @param alpha float, type I error prob. 
  #' @param seed int
  #' @param K int, number of Monte Carlo draws
  #' @param do.par bool, whether calculations are carried out in paralell
  

  sim.path <- file.path("../temp",paste("n",n,"T",tt,"K",K,".csv", sep = "_"))
  
  
  if (file.exists(sim.path)) {
    
    max.sim <- read.csv(sim.path) 
    
    return(c(quantile(max.sim[,2], 1-alpha), quantile(max.sim[,3], 1-alpha)))
    
    
  } else {
    
    max.sim <- simulate.panel.max(n, tt, seed, K, do.par)
    
    return(c(quantile(max.sim[,1], 1-alpha), quantile(max.sim[,2], 1-alpha)))
    
  }
  
}

