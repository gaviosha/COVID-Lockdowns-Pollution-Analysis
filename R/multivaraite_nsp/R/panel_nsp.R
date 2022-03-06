
#-----------------------------------------------------------------
#
# Search for narrowest regions of structural change in panel data 
#
# Motivation: http://stats.lse.ac.uk/fryzlewicz/nsp/nsp.pdf
#
#-----------------------------------------------------------------

# cl <- makeCluster(detectCores()-1)
#   
# registerDoParallel(cl)


panel_nsp <- function(X, deg = 0, thresh = NA, M = c(100,1000), do.par = FALSE, K = 1000) {
  
  #' Narrowest Significance Pursuit for Panel Data
  #'
  #'@param X matrix, panel data with dimension (n x tt)
  #'@param deg int, degree of piecewise polynomila signal
  #'@param M int, number of intervals to draw
  #'@param do.par bool, determines whether calculations are carried out in paralell 
  
  X <- t(apply(X,1, function(ii) ii / mad(ii)))

  Z <- build_poly_constraint(ncol(X),deg)
  
  ads.array <- all_dyadic_scans(X,Z)
  
  nsp.out <- iter_scan_array(c(1,ncol(X)), ads.array, M, thresh, do.par)
  
  return(list(thresh.val = thresh, nsp.out = nsp.out))
  
}