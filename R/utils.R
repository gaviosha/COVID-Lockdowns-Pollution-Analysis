
plot.with.dates.axis <- function(x, dates, xlab = "", ylab = "", main = "", col = "black", len.out = 10, ylim = NULL)
{
  #'Plot time series with correct dates on x-axis
  #'
  #'@param x vector, the time series
  #'@param dates vector, the dates
  
  plot(x, 
       type = "l",
       xaxt = "n",
       xlab = xlab, 
       ylab = ylab, 
       main = main, 
       col = col,
       ylim = ylim
       )
  
  axis(1, 
       at = seq(1, length(x), length.out = len.out), 
       label = format.Date(dates, "%Y-%m")[seq(1, length(x), length.out = len.out)]
       )
  
}


df.plot.with.dates.axis <- function(df, xlab = "", ylab = "", main = "", col = "black", len.out = 10, ylim = NULL)
{
  #'Plot all time series in data frame with date on x-axis
  #'
  #'@param df data frame 
  
  max.sd <- max(sapply(1:ncol(df), function(ii) sd(df[,ii])))
  
  ylim <- c(min(df) - max.sd, max(df) + max.sd)
  
  plot.with.dates.axis(df[,1], rownames(df), xlab, ylab, main, col, len.out, ylim)
  
  if (ncol(df) > 1) for (ii in 2:ncol(df)) lines(df[,ii], col = col)
  
}


fit.pcwsLin <- function(x)
{
  #'Fit pieicewise linear function to signal using competing methods
  #'
  #' Methods used: 
  #'    * Narrowest Over Threshold
  #'    * Isolate Detect 
  #'    * Trend Filtering
  #'    * Least Squares Estimation
  #'
  #'@param x vector

  res <- list()
  
  not.obj <- not(x, contrast = "pcwsLinContMean")
  
  res[["NOT"]] <- list(
   fit = predict(not.obj), 
   cpt = features(not.obj)[["cpt"]]
  )
  
  id.obj <- ID(x, contrast = "slope")
  
  res[["ID"]] <- list(
    fit = id.obj$fit, 
    cpt = id.obj$cpt
  )
  
  
  tf.obj <- trendfilter(x, ord = 1)

  res[["trend filtering"]] <- list(
    fit = tf.obj$fit[,100], 
    cpt = which(abs(diff(tf.obj$fit[,100], differences=2)) > sqrt(.Machine$double.eps))
  )
  
  tt <- 1:length(x)
  struc.obj <- breakpoints(x ~ tt, h = 0.08)
  
  res[["struc-change"]] <- list(
    fit = as.vector(fitted(struc.obj)), 
    cpt = struc.obj$breakpoints
  )
  
  return(res)
}




fit.multivaraute.changepoints <- function(x)
{
  #'Find changepoints in panel data using the following methods
  #'    * Sparse Projection
  #'    * 
  #'    * Double CUSUM 
  #'    * Sparsified BinSeg
  #'
  #'@param x matrix,l each time series is a 
  
  res <- list()
  
  inspect.obj <- inspect(x)
  res[["inspect"]] <- inspect.obj$changepoints[,1]
  
  geom.obj <- geomcp(t(x))
  res[["geometric mapping"]] <- geom.obj@ang.cpts
  
  dcbs.obj <- dcbs.alg(x)
  res[["double cusum"]] <- dcbs.obj$ecp 
  
  sbs.obj <- sbs.alg(x)
  res[["sparsified binary segmentation"]] <- sbs.obj$ecp
  
  return(res)
}



draw_mnsp_rects <- function(mnsp.obj, yrange, density = 10, col = "red")
{
  #'
  #'
  #'@param mnps.obj
  #'@param yrange
  #'@param density
  #'@param col 
  
  for (ii in 1:ncol(mnsp.obj))
  {
    rect(mnsp.obj[,ii][1], 
         yrange[1], 
         mnsp.obj[,ii][2],
         yrange[2], 
         density = density, 
         col = col
    )
  }
}


























