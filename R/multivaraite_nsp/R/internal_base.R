
#--------------------------------------------------
#
# Internal functions: lowest level of abstraction
#
#--------------------------------------------------

library(docstring)
library(foreach)

build_poly_constraint <- function(tt, deg) {
  
  #' Build matrix of contrainats based on linear model
  #'
  #'@param tt int, number of time points in the panel
  #'@param deg int, degree of polynomial in linear model 
  
  Z <- matrix(,tt, deg+1)
  
  for (i in 1:(deg+1)) Z[,i] <- seq(from = 0, to =  1, length = tt) ** (i-1)
  
  return(Z)
  
}



build_lpSolve_LHS <- function(Z, shifts, tt) {
  
  #'Constraints, `const.mat` in the function `lp`
  #'
  #'@param Z array, dim(tt x deg+1 x shifts)
  #'@param shifts vector
  #'@param tt int
  
  
  scales <- length(shifts)
  
  
  if (dim(Z)[2] > 1) {
    
    const.mat <- matrix(0,0,dim(Z)[2])
    
    for (i in 1:scales) const.mat <- rbind(const.mat,Z[1:(tt-shifts[i]),,i])
    
    const.mat <- rbind(const.mat, -const.mat)
    
    q <- nrow(const.mat)
    
    
  } else {
    
    const.mat <- c()
    
    for (i in 1:scales) const.mat <- c(const.mat,Z[1:(tt-shifts[i]),,i])
    
    const.mat <- c(const.mat, -const.mat)
    
    q <- length(const.mat)
    
  }
  
  
  return(cbind(rep(1,q), const.mat, -const.mat))
  
}



build_lpSolve_RHS <- function(X, shifts, tt) {
  
  #' Constraints, `const.rhs` in  the function `lp`
  #'
  #'@param X array, dim (1 x tt x scales)
  #'@param shifts vector
  #'@param tt int
  
  
  if (!is.array(X)) return(c(X, -X))
  
  
  scales <- length(shifts)
  
  tt <- dim(X)[1]
  
  
  const.rhs <- c()
  
  for (i in 1:scales) const.rhs <- c(const.rhs, X[1:(tt-shifts[i]),i])
  
  return(c(const.rhs,-const.rhs))
  
}



all_intervals_flat <- function(n) {
  
  #' Piotr's code: https://github.com/pfryz/nsp/blob/master/NSP_for_Github_v4.R#L322 

  if (n == 2) ind <- matrix(1:2, 2, 1) else {
    M <- (n-1)*n/2
    ind <- matrix(0, 2, M)
    ind[1,] <- rep(1:(n-1), (n-1):1)
    ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)
  }
  ind

}


all_intervals_sorted <- function(n) {
  
  #' Piotr's code: https://github.com/pfryz/nsp/blob/master/NSP_for_Github_v4.R#L322

  d <- all_intervals_flat(n)
  d.ord <- order(d[2,] - d[1,])
  d[,d.ord, drop=FALSE]

}


grid_intervals_sorted <- function(tt, M) {
  
  #' grid of dyadic intervals sorted by length
  #'
  #'@param M int, max number of intervals to draw
  #'@param tt int, intervals will span c(1,tt)
  
  
  if ((tt==2) || (M == 0)) {
    
    ind <- matrix(c(1, tt), 2, 1)
    
    
  } else if (M >= (tt-1)*tt/2) {
    
    ind <- all_intervals_sorted(tt)
    
    
  } else {
    
    k <- 1
    
    while (k*(k-1)/2 < M) k <- k+1
    
    ind2 <- all_intervals_sorted(k)
    
    ind2.mx <- max(ind2)
    
    ind <- round((ind2 - 1) * ((tt-1) / (ind2.mx-1)) + 1)
  
  }

  return(ind)
  
}
