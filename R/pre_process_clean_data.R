

stabilize.varaince.detrend <- function(path.to.data, 
                                   largest.cities, 
                                   target.year = "2020", 
                                   trend.fit.years = c("2019", "2018")
                                   )
{
  #' Square-root transform + trend removal by regressing on Fourier frequencies
  #'
  #'@param path.to.data string
  #'@param largest.cities vector
  #'@param target.year string
  #'@param trend.fit.years vector
  
  all.datasets <- list()
  all.countries <- list.files(path.to.data)

  for (country in all.countries)
  {
    country.dir <- file.path(path.to.data, country)
    all.cities <- list.files(country.dir)
    
    for (city in all.cities)
    {
      
      match.city <- agrep(city, largest.cities)
      if (identical(integer(0), match.city)) next()
      
      city.data <- load.test.train.data(file.path(country.dir,city), target.year, trend.fit.years)
      if (is.null(city.data)) next()
      
      all.datasets[[length(all.datasets)+1]] <- detrend.by.fourier.frequencies(city.data)
    }
  }
  
  suppressWarnings(
    all.datasets <- Reduce(function(df1, df2) merge(df1, df2, by = "Date", all.x = TRUE, all.y = TRUE), all.datasets)
  )
  
  return(all.datasets)
  
}

load.test.train.data <- function(path.to.files,
                                 target.year,
                                 trend.fit.years)
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
    if (length(ind) > 1) ind <- ind[1]
    trend.fit.data[[length(trend.fit.data)+1]] <- read.csv(file.path(path.to.files, all.files[ind]), row.names = 1)
    trend.fit.data[[length(trend.fit.data)]]$Concentration <- special.na.approx(trend.fit.data[[length(trend.fit.data)]]$Concentration)
  }
  
  if (length(trend.fit.data) == 0) return(NULL)
  if (length(trend.fit.data) == 1) {
    city.data[["trend.fit.data"]] <- trend.fit.data[[1]]
    return(city.data)
  }
  
  city.data[["trend.fit.data"]] <- Reduce(function(df1, df2) rbind(df1, df2), trend.fit.data)
  return(city.data)
  
}


detrend.by.fourier.frequencies <- function(city.data)
{
  #' Remove seasonal component from target data
  #'
  #'@param city.data list, output from `load.test.train.data`
  
  
  n.fit <- nrow(city.data$trend.fit.data)
  n.target <- nrow(city.data$target.year)
  tt <- 1:n.fit
  
  c1 <- cos(2*pi*tt/7); s1 <- sin(2*pi*tt/7)
  c2 <- cos(2*pi*tt/30); s2 <- sin(2*pi*tt/30)
  c3 <- cos(2*pi*tt/365); s3 <- sin(2*pi*tt/365)
  
  lm.seasonal <- lm(city.data$trend.fit.data$Concentration ~ c1 + s1 + c2 + s2 + c3 + s3)
  
  e <- city.data$target.year$Concentration - predict.lm(lm.seasonal)[1:n.target]
  city.data$target.year$Concentration <- e
  
  return(city.data$target.year)
  
}


special.na.approx <- function(x)
{
  #' Linear interpolation for missing values including NAs at start and end
  #'
  #'@param x vector 

  x <- suppressWarnings(
    na.approx(x, na.rm = FALSE)
  )
  
  ind <- which(is.na(x))
  if (identical(integer(0), ind)) return(x)
  
  
  nn <- length(x)
  for (ii in ind)
  {
    if (ii == 1) x[ii] <- x[ii+1]
    if (ii == nn) x[ii] <- x[ii-1]
  }
   
  return (x)
   
}


