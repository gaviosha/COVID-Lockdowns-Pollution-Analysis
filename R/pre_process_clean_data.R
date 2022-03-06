
stabilize.varaince.detrend <- function(path.to.data, 
                                   cities, 
                                   target.year = "2020", 
                                   trend.fit.years = c("2019", "2018"),
                                   max.na.run = 31, 
                                   max.na.prop = 0.1, 
                                   order.max = 3, 
                                   show.prog.bar = TRUE
                                   )
{
  #' Square-root transform + trend removal by regressing on Fourier frequencies
  #'
  #'@param path.to.data string, path to roo of pollutant data dir 
  #'@param cities vector, which cities to pre-process
  #'@param target.year string, which year to whiten and detrend 
  #'@param trend.fit.years vector, which years to use for model fitting
  #'@param order.max int, max order for AR model used in pre-whitening
  #'@param max.na.run int, maximum run of missing values permitted before time series is dropped
  #'@param max.na.prop float, maximum proportion of missing values permitted before time series is dropped
  #'@param show.prog.bar bool, whether to display pb  
  
  all.datasets <- list()
  selected.cities <- c()
  all.countries <- list.files(path.to.data)
  
  all.datasets.detrended <- list()
  all.datasets.whitened <- list()
  
  if (show.prog.bar && length(cities) > 1)
  {
    pb <- txtProgressBar(min = 1, max = length(cities))
    ii <- 1
  }

  for (country in all.countries)
  {
    country.dir <- file.path(path.to.data, country)
    all.cities <- list.files(country.dir)
    
    for (city in all.cities)
    {
      
      match.city <- agrep(city, cities)
      
      if (identical(integer(0), match.city)) next()
      city.data <- load.test.train.data(file.path(country.dir,city), target.year, trend.fit.years, max.na.run, max.na.prop)
      if (is.null(city.data)) next()

      selected.cities <- c(selected.cities, city)
      temp <- detrend.by.fourier.frequencies(city.data)
      all.datasets.detrended[[length(all.datasets.detrended)+1]] <- temp$detrended
      all.datasets.whitened[[length(all.datasets.whitened)+1]] <- whiten.detreneded.data(temp$detrended, temp$training.residuals, order.max)
      
      if (show.prog.bar && length(cities) > 1)
      {
        setTxtProgressBar(pb, ii)
        ii <- ii + 1 
      }
    }
  }
  
  all.datasets.detrended <- tidy.up(selected.cities, all.datasets.detrended, max.na.run, max.na.prop)
  all.datasets.whitened <- tidy.up(selected.cities, all.datasets.whitened, max.na.run, max.na.prop)
  
  return(
    list(detrended = all.datasets.detrended, whitened = all.datasets.whitened)
  )
}

load.test.train.data <- function(path.to.files,
                                 target.year,
                                 trend.fit.years,
                                 max.na.run, 
                                 max.na.prop)
{
  #' Load target data and some training data to fit trend component
  #'
  #'@param path.to.files string 
  #'@param target.year string 
  #'@param trend.fit.years string 
  
  all.files <- list.files(path.to.files)
  city.data <- list()
  trend.fit.data <- list()
  
  ind <- grep(target.year, all.files)
  if (identical(integer(0), ind)) return(NULL)
  if (length(ind) > 1) ind <- ind[1]
  city.data[["target.year"]] <- read.csv(file.path(path.to.files, all.files[ind]), row.names = 1)
  city.data[["target.year"]]$Concentration <- special.na.approx(city.data[["target.year"]]$Concentration)
  
  for (year in trend.fit.years)
  {
    ind <- grep(year, all.files)
    if (identical(integer(0),ind)) next()
    if (length(ind) > 1) ind <- ind[1]
    
    temp <- read.csv(file.path(path.to.files, all.files[ind]), row.names = 1)
    if (dim(temp)[1] < 180 ) next()
    if (longest.na.run(temp[["Concentration"]]) > max.na.run) next()
    if (na.prop(temp[["Concentration"]]) > max.na.prop) next()

    trend.fit.data[[length(trend.fit.data)+1]] <- temp
    trend.fit.data[[length(trend.fit.data)]]$Concentration <- special.na.approx(trend.fit.data[[length(trend.fit.data)]]$Concentration)
  }
  
  if (length(trend.fit.data) == 0) return(NULL)
  if (length(trend.fit.data) == 1) {
    city.data[["trend.fit.data"]] <- trend.fit.data[[1]]
    return(city.data)
  }
  
  suppressWarnings(
    city.data[["trend.fit.data"]] <- Reduce(function(df1, df2) rbind(df1, df2), trend.fit.data)
  )
  
  return(city.data)
  
}


detrend.by.fourier.frequencies <- function(city.data)
{
  #' Remove seasonal component from target data
  #'
  #'@param city.data list, output from `load.test.train.data`
  
  n.fit <- nrow(city.data$trend.fit.data)
  n.target <- nrow(city.data$target.year)
  
  days <- as.Date(city.data$trend.fit.data$Date, format = "%Y-%m-%d")
  sat <- ifelse(weekdays(days) == "Saturday",1, 0)
  sun <- ifelse(weekdays(days) == "Sunday", 1, 0)
  
  tt <- 1:n.fit
  
  s.qt <- sin(2*pi*tt/365/4)
  c.qt <- cos(2*pi*tt/365/4)
  s.yr <- sin(2*pi*tt/365)
  c.yr <- cos(2*pi*tt/365)
  
  lm.seasonal <- lm(
    city.data$trend.fit.data$Concentration ~ sat + sun + s.qt + c.qt + s.yr + c.yr
    )
  
  
  days <- as.Date(city.data$target.year$Date, format = "%Y-%m-%d")
  sat <- ifelse(weekdays(days) == "Saturday",1, 0)
  sun <- ifelse(weekdays(days) == "Sunday", 1, 0)
  
  tt <- 1:n.target
  
  s.qt <- sin(2*pi*tt/365/4)
  c.qt <- cos(2*pi*tt/365/4)
  s.yr <- sin(2*pi*tt/365)
  c.yr <- cos(2*pi*tt/365)
  
  e <- city.data$target.year$Concentration - predict.lm(
    lm.seasonal, 
    newdata = data.frame(sat, sun, s.qt, c.qt, s.yr, c.yr)
    )
  
  city.data$target.year$Concentration <- e
  
  return(
    list(detrended = city.data$target.year, training.residuals = lm.seasonal$residuals)
  )
  
}


special.na.approx <- function(x)
{
  #' Linear interpolation for missing values including NAs at start and end
  #'
  #'@param x vector 

  x <- suppressWarnings(
    na.approx(x, na.rm = FALSE)
  )
  
  if (is.na(x[1])) x <- interpolate.LHS(x)
   
  if (is.na(tail(x,1))) x <- interpolate.RHS(x)
  
  return (x)
   
}


interpolate.LHS <- function(x)
{
  #'Linear interpolation for run of NAs at start of vector only
  #'
  #'@param x vector
  
  zz <- cumsum(is.na(x))
  n <- length(x)
  
  for (ii in 2:length(zz)) if(zz[ii] == zz[ii-1]) ind <- zz[ii]
  
  if (ind == 1)
  {
    x[ind] <- x[ind+1]
    return(x)
  }
  
  tt <- (ind+1):min(n,(2*ind))
  y <- x[tt]
  
  b.1 <- sum((y-mean(y))*(tt-mean(tt))) / sum((tt - mean(tt))**2)
  b.0 <- mean(y) - b.1 * mean(tt)
  
  x[1:ind] <- b.0 + 1:ind * b.1
  
  return(x)
}


interpolate.RHS <- function(x)
{
  #'Linear interpolation for run of NAs at end of vector only
  #'
  #'@param x vector
  
  zz <- cumsum(is.na(rev(x)))
  n <- length(x)
  
  for (ii in 2:length(zz)) if(zz[ii] == zz[ii-1]) ind <- zz[ii]
  
  if (ind == 1)
  {
    x[n] <- x[n-ind]
    return(x)
  }
  
  n <- length(x)
  tt <- max(1,(n - 2*ind+1)):(n-ind)
  y <- x[tt]
  
  b.1 <- sum((y-mean(y))*(tt-mean(tt))) / sum((tt - mean(tt))**2)
  b.0 <- mean(y) - b.1 * mean(tt)
  
  x[(n-ind+1):n] <- b.0 + (n-ind+1):n * b.1
  
  return(x)
  
}


whiten.detreneded.data <- function(detrended, training.residuals, order.max)
{
  #'White data by saving resids from AR fit
  #'
  #'@param detrended data fram
  #'@param training.residuals vector
  #'@param order.max int, max order for model fit
  
  ar.fit <- ar(training.residuals, order.max = order.max)
  ar.coefs <- ar.fit$ar
  ar.order <- ar.fit$order
  
  z <- c(rep(0,ar.order), detrended$Concentration)
  n <- nrow(detrended)
  
  for (ii in 1:n) detrended$Concentration[ii] <- detrended$Concentration[ii] - sum(ar.coefs * z[(ii+ar.order):(ii+1)])
  
  return(detrended)
  
}


tidy.up <- function(selected.cities, df.list, max.na.run, max.na.prop)
{
  #'Re-format to datafram with cities as row names and dates as col names
  #'
  #'@param selected.cities vector, name of selected cities 
  #'@param df.list list, data frame of clean time series
  #'@param max.na.run int
  #'@param max.na.prop float

  if (length(selected.cities) == 1)
  {
    dates <- df.list[[1]][,1]
    df.list <- data.frame(df.list[[1]][,2])
    names(df.list) <- selected.cities
    rownames(df.list) <- dates
    return(df.list)
  }
  
  suppressWarnings(
    df.list <- Reduce(function(df1, df2) merge(df1, df2, by = "Date", all.x = TRUE, all.y = TRUE), df.list)
  )
  
  rownames(df.list) <- df.list[,1]
  df.list <- df.list[,-1]
  
  drop.ind <- union(
    which(apply(df.list,2, longest.na.run) > max.na.run),
    which(apply(df.list, 2, na.prop) > max.na.prop)
  )
  
  if (!identical(drop.ind, integer(0)))
  {
    df.list <- df.list[,-drop.ind]
    selected.cities <- selected.cities[-drop.ind]
  }
  
  names(df.list) <- selected.cities
  for (city in selected.cities) df.list[[city]] <- special.na.approx(df.list[[city]])
  
  return(df.list)
  
}


longest.na.run <- function(x)
{
  #'Finds the longest run of NAs
  #'
  #'@param x vector
  
  max.run <- 0
  cur.run <- 0 
  n <- length(x)
  
  for (ii in 1:n)
  {
    if (is.na(x[ii]))
    {
      cur.run <- cur.run + 1
    }
    else 
    {
      if (cur.run > max.run) max.run <- cur.run
      cur.run <- 0 
    }
  }
  
  max.run
}


na.prop <- function(x)
{
  #'Finds proportion of missing vals 
  #'
  #'@param x vector

  sum(is.na(x)) / length(x)
    
}
