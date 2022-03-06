
min_width <- as.numeric(
  readline(prompt = "enter minimum segment width: ")
)

random_checks_scan_array_1by1 <- function(ind, ads.array, M, thresh) {
  
  s <- ind[1]
  e <- ind[2]
  n <- e - s + 1
  
  if (n >= min_width) {
    
    indices <- ((ads.array$shifts+1) <= (n/2))
    
    ads.array$res <- ads.array$res[s:e,,indices,drop=F]
    
    ads.array$shifts <- ads.array$shifts[indices]
    
    M <- min(M, (n-1)*n/2)
    
    ind <- grid_intervals_sorted(n, M)
    
    M <- dim(ind)[2]
    
    res <- matrix(0, M, 3)
    
    res[,1:2] <- t(ind)
    
    zero.check <- TRUE
    j <- 1
    
    while (zero.check && (j <= M)) {
      
      res[j,3] <- check_interval_array(res[j,1:2], ads.array, thresh)
      zero.check <- (res[j,3] == 0)
      j <- j + 1
      
    }
    
    if (zero.check) {
      
      selected.ind <- NA
      selected.val <- matrix(0, 0, 3)
      
    }
    
    else {
      
      selected.ind <- j-1
      selected.val <- res[selected.ind,,drop=FALSE]
      
    }
    
    
  }
  
  else {
    
    selected.val <- matrix(0, 0, 3)
    selected.ind <- NA
    M <- 0
    
  }
  
  
  list(selected.ind = selected.ind, selected.val = selected.val, M.eff=max(1, M))
  
}