
#-----------------------------------------------
#
# Internal functions: for scanning sub-intervals
#
#-----------------------------------------------



library(lpSolveAPI)



all_dyadic_scans <- function(X,Z) {
  
  #' Build all dyadic scans based on panel and constrain matrix
  #'
  #'@param X matrix, dim (n x tt) panel data 
  #'@param Z matrix, dim (tt x deg+1) linear model constraints
  
  
  tt <- dim(X)[2]
  
  if (tt > 0) {
    
    add.scales <- floor(logb(tt,2))
    
    shifts <- rep(0, add.scales+1)
    
    res.X <- array(X, c(dim(X), add.scales+1))
    
    res.Z <- array(Z, c(dim(Z), add.scales+1))
    
    
    if (add.scales) for (i in 1:add.scales) {
      
      
      res.X[,1:(tt-2^i+1),(i+1)] <- 2^(-1/2) * (res.X[,1:(tt-2^i+1),i] + res.X[,(2^(i-1)+1):(tt-2^i+1+2^(i-1)),i])
      
      res.X[,(tt-2^i+2):(tt),(i+1)] <- 0
      
      
      res.Z[1:(tt-2^i+1),,(i+1)] <- 2^(-1/2) * (res.Z[1:(tt-2^i+1),,i] + res.Z[(2^(i-1)+1):(tt-2^i+1+2^(i-1)),,i])
      
      res.Z[(tt-2^i+2):(tt),,(i+1)] <- 0
      
      
      shifts[i+1] <- 2^i-1
      
      
    }
    
  }
  
  
  return(list(res.X = res.X, res.Z = res.Z, shifts = shifts))
  
  
}



iter_scan_array <- function(ind, ads.array, M, thresh, do.par) {
  
  #' What it does...
  #'
  #'@param ind vector, start and end of interval
  #'@param ads.array list, as returned by `all_dyadic_scans`
  #'@param M int, number of intervals to draw
  #'@param thresh vector, thresholds forcomparing $\ell_2$ and $\ell_\infty$ aggregation 
  #'@param do.par bool, determines whether calculations are carried out in paralell 
  
  
  s <- ind[1]
  
  e <- ind[2]
  
  tt <- e - s + 1
  
  
  indices <- ((ads.array$shifts+1) <= (tt/2))
  
  ads.array$res.X <- ads.array$res.X[,s:e,indices,drop=F]
  
  ads.array$res.Z <- ads.array$res.Z[s:e,,indices,drop=F]
  
  ads.array$shifts <- ads.array$shifts[indices]
  
  
  if (tt > 1) {
    
    next.int <- scan_2_stage_array(c(1,tt), ads.array, M, thresh, do.par)
    
    
    if (!is.na(next.int$selected.ind))  {
      
      
      if (next.int$selected.val[1,1] >= 2) {
        
        left <- iter_scan_array(c(1, next.int$selected.val[1,1]), ads.array, M, thresh, do.par)
        
      } else {
        
        left <- matrix(NA, 4, 0)
        
      }
      
      
      if (tt - next.int$selected.val[1,2] >= 1) {
        
        right <- iter_scan_array(c(next.int$selected.val[1,2], tt), ads.array, M, thresh, do.par)
        
        if (dim(right)[2]) right <- right + c(rep(next.int$selected.val[1,2]-1, 2), 0,0)
        
      } else {
        
        right <- matrix(NA, 4, 0)
        
      } 
      
      
      return(cbind(t(next.int$selected.val), left, right))
      
      
    } else {
      
      return(matrix(NA, 4, 0))
      
    }
    
    
  } else {
    
    return(matrix(NA, 4, 0))
    
  }
}




scan_2_stage_array <- function(ind, ads.array, M, thresh, do.par) {
  
  #' Scan aray over interval, if deviation is detected perform secondary check 
  #'
  #'@param ind vector, start and end of interval
  #'@param ads.array list, as returned by `all_dyadic_scans`
  #'@param M int, number of intervals to draw
  #'@param thresh vector, thresholds forcomparing $\ell_2$ and $\ell_\infty$ aggregation 
  #'@param do.par bool, determines whether calculations are carried out in paralell 
  
  
  s1 <- scan_array_1_by_1(ind, ads.array, M[1], thresh, do.par)
  
  
  if (!is.na(s1$selected.ind)) {
    
    s <- s1$selected.val[1,1] + ind[1] - 1
    
    e <- s1$selected.val[1,2] + ind[1] - 1
    
    s2 <- scan_array_1_by_1(c(s,e), ads.array, M[2], thresh, do.par)
    
    
    if (!is.na(s2$selected.ind)) {
      
      replacement <- s2$selected.val
      
      replacement[1,1:2] <- replacement[1,1:2] + s - ind[1]
      
      s1$selected.val <- replacement
      
    }
    
  }
  
  s1	
  
}



scan_array_1_by_1 <- function(ind, ads.array, M, thresh, do.par) {
  
  #' Scan sub-intervals 1-by-1 looking for narrowest region of deviation
  #' 
  #'@param ind vector, start and end of interval
  #'@param ads.array list, as returned by `all_dyadic_scans`
  #'@param M int, number of intervals to draw
  #'@param thresh vector, thresholds forcomparing $\ell_2$ and $\ell_\infty$ aggregation 
  #'@param do.par bool, determines whether calculations are carried out in paralell 
  
  
  s <- ind[1]
  
  e <- ind[2]
  
  tt <- e - s + 1
  
  
  if (tt > 1) {
    
    
    indices <- ((ads.array$shifts+1) <= (tt/2))
    
    ads.array$res.X <- ads.array$res.X[,s:e,indices,drop=F]
    
    ads.array$res.Z <- ads.array$res.Z[s:e,,indices,drop=F]
    
    ads.array$shifts <- ads.array$shifts[indices]
    
    
    M <- min(M, (tt-1)*tt/2)
    
    grid_intervals <- grid_intervals_sorted(tt, M)
    
    res <- matrix(0, dim(grid_intervals)[2], 4)
    
    res[,1:2] <- t(grid_intervals)
    
    
    M <- nrow(res)
    
    thresh.check <- TRUE
    
    j <- 1
    

    while (thresh.check && (j <= M)) {
      
      res[j,3:4] <- check_interval_array(res[j,1:2], ads.array, do.par, thresh)
      
      thresh.check <- all(res[j,3:4] <= thresh)
      
      j <- j + 1 
      
    }
    
    
    if (thresh.check) {
      
      selected.ind <- NA
      
      selected.val <- matrix(0, 0, 4)
      
      
    } else {
      
      selected.ind <- j-1
      
      selected.val <- res[selected.ind,,drop=FALSE]
      
    }
    
    
  } else {
    
    selected.val <- matrix(0, 0, 4)
    
    selected.ind <- NA
    
    M <- 0
    
  }
  
  
  return(list(selected.ind = selected.ind, selected.val = selected.val, M.eff=max(1, M)))
  
}



check_interval_array <- function(ind, ads.array, do.par, thresh) {
  
  #' Check interval for deviation from linearity
  #'
  #'@param ind vector, start and end of interval
  #'@param ads.array list, as returned by `all_dyadic_scans`
  #'@param do.par bool, determines whether calculations are carried out in paralell 
  #'@param thresh vector, thresholds forcomparing $\ell_2$ and $\ell_\infty$ aggregation 
  
  
  s <- ind[1]
  
  e <- ind[2]
  
  tt <- e - s + 1
  
  
  indices <- ((ads.array$shifts+1) <= (tt/2))
  
  ads.array$res.X <- ads.array$res.X[,s:e,indices,drop=F]
  
  ads.array$res.Z <- ads.array$res.Z[s:e,,indices,drop=F]
  
  ads.array$shifts <- ads.array$shifts[indices]
  
  
  d.X <- dim(ads.array$res.X)
  
  d.Z <- dim(ads.array$res.Z)
  
  
  objective.in <- c(1,rep(0,2*d.Z[2]))
  
  const.mat <- build_lpSolve_LHS(ads.array$res.Z, ads.array$shifts, tt)
  
  const.dir <- rep(">=", nrow(const.mat))
  
  
  if (do.par) {
    
    
    max.sum.comb <- function(a,b) return(c(max(a[1],b[1]), sum(a[2],b[2])))
    
    scan.vals <- foreach(i=1:d.X[1], .combine = "max.sum.comb") %do% {
      
      i.const.rhs <- build_lpSolve_RHS(ads.array$res.X[i,,], ads.array$shifts, tt)
      
      i.MR.norm <- lp("min", objective.in, const.mat, const.dir, i.const.rhs)$solution[1]
      
      c(i.MR.norm, i.MR.norm**2)
      
    }
   
     
  } else {
    
    
    scan.max <- 0 
    
    scan.sum <- 0 
    
    for (i in 1:d.X[1]) {
      
      i.const.rhs <- build_lpSolve_RHS(ads.array$res.X[i,,], ads.array$shifts, tt)
      
      lps.model <- make.lp(nrow = nrow(const.mat), ncol = ncol(const.mat))
      
      for (j in 1:ncol(const.mat)) set.column(lps.model, j, const.mat[,j])
      
      set.rhs(lps.model, i.const.rhs)
      
      set.objfn(lps.model,objective.in)
      
      set.constr.type(lps.model,const.dir)
      
      solve(lps.model)
      
      i.MR.norm <- get.objective(lps.model)
      
      
      if (i.MR.norm > scan.max) scan.max <- i.MR.norm
      
      scan.sum <- scan.sum + (i.MR.norm ** 2)
      
      if (any(c(scan.max, scan.sum) >= thresh)) break()
      
    }
    
    scan.vals <- c(scan.max, scan.sum)
    
  }
  
  
  return(scan.vals)
  
}



