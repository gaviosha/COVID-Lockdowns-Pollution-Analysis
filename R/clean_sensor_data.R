
aggregate.sensor.data <- function(path.to.data)
{
  #'Clean and aggregate mid-frequency sensor measurments
  #'
  #' * Average over hourly measurements at sensor locations to obtain daily measirements
  #' * Average over sensor measure measurements for to obtain single concentration measure for city
  #'
  #'@param path.to.data string 
  
  list.datasets <- list.files(path.to.data)
  
  if (identical(list.datasets, character(0))) return(character(0))

  all.datasets <- list()
  
  
  for (dataset in list.datasets)
  {
    
    dataset.path <- file.path(path.to.data, dataset)
    
    temp <- read.csv(dataset.path)
    if (!all(c("DatetimeBegin", "Concentration") %in% names(temp)))
    {
      next()
    } 
    else
    {
      temp <- temp[c("DatetimeBegin", "Concentration")]
    }
    
    temp$DatetimeBegin <- sapply(temp$DatetimeBegin, function(ii) substr(toString(ii),1,10))
    temp <- aggregate(temp$Concentration, by = list(temp$DatetimeBegin), function(ii) mean(ii, na.rm = TRUE))
    names(temp) <- c("Date", "Concentration")

    all.datasets[[length(all.datasets)+1]] <- temp
  }
  
  if (length(all.datasets) == 1) return(all.datasets[[1]])
  
  suppressWarnings(
    all.datasets <- Reduce(function(df1, df2) merge(df1, df2, by = "Date", all.x = TRUE, all.y = TRUE), all.datasets)
  )
  
  all.datasets[,2] <- sapply(
    1:nrow(all.datasets),
    function(ii) ifelse(all(is.na(all.datasets[,-1][ii,])), NA, rowMeans(all.datasets[,-1][ii,], na.rm = TRUE))
  )

  all.datasets <- all.datasets[,1:2]
  names(all.datasets) <- c("Date", "Concentration")

  return(all.datasets)
  
}


clean.pollutant.data <- function(pollutant, data.dir = NULL)
{
  
 #' Clean pollution data
 #' 
 #' * Make dir structure to store clean data
 #' * Produce daily data from mid-frequency data from many sensors 
 #' 
 #'@param pollutant string
 #'@param data.dir string

  
  if (is.null(data.dir)) data.dir <- getwd()
  
  pollutant.dir <- file.path(data.dir, "raw", pollutant)
  if (!file.exists(pollutant.dir)) stop(paste("No directory found for pollutant: ",pollutant))
  
  
  if (!file.exists("clean")){
    dir.create("clean")
  }

  
  clean.dir <- file.path(data.dir, "clean", pollutant)
  if (!file.exists(clean.dir)){
    dir.create(clean.dir)
  }
  
  
  countries.list <- list.files(pollutant.dir)
  
  n <- length(countries.list)
  pb <- txtProgressBar(1, n, style = 3)
  ii <- 0
  
  
  for (country in countries.list)
  {
    
    clean.country.dir <- file.path(clean.dir, country)
    if (!file.exists(clean.country.dir)) dir.create(clean.country.dir)

    country.dir <- file.path(pollutant.dir, country)
    cities.list <- list.files(country.dir)
    
    for (city in cities.list)
    {
      
      clean.city.dir <- file.path(clean.country.dir, city)
      if (!file.exists(clean.city.dir)) dir.create(clean.city.dir)
      
      city.dir <- file.path(country.dir, city)
      years.list <- list.files(city.dir)
      
      for (year in years.list)
      {
        
        path.to.data <- file.path(city.dir, year)
        cleaned.data <- aggregate.sensor.data(path.to.data)
        
        if (identical(character(0), cleaned.data)) next()
        
        file.name <- file.path(
          clean.city.dir,
          paste(
            paste(pollutant, country, city, year, sep = "_"),
            ".csv",
            sep = ""
          )
        )
        
        write.csv(cleaned.data, file = file.name)
        
      }
    }
    
    ii <- ii + 1
    setTxtProgressBar(pb, ii)
  }
}


























