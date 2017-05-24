# Library handle ----------------------------------------------------------

set_lib_path_local <- function(my.lib = "Packages"){
  .libPaths(c(.libPaths(), my.lib))
}

install_lib_local <- function(package, my.lib = "Packages"){
  set_lib_path_local(my.lib)
  install.packages(package, lib = my.lib)
}

set_library <- function(){
  if(!require(MSGARCH) || !require(data.table)){
    set_lib_path_local()
    library(MSGARCH)
    library(data.table)
  }
}

# Raw data processing -----------------------------------------------------

read.data.closing <- function(in_file, length_data, head_data = F){
  data <- read.csv(paste0("Input/", in_file, ".csv"), sep = ",", stringsAsFactors = F)
  data <- data[order(as.Date(data$Date)), ]
  
  if(!missing(length_data)){
    if(nrow(data)<=length_data) stop("Input file has smaller length than specified by length_data.")
    
    if(head_data){
      data <- data[1:(length_data+1),]
    } else {
      data <- data[(nrow(data)-length_data):nrow(data),]
    }
  }
  
  return(list(date = as.Date(data$Date),
              closing = data$Close,
              log.return = 100*diff(log(data$Close))))
}

plot.data <- function(data, name, save_plot = F){
  data$date <- as.Date(data$date)
  
  if(save_plot) jpeg(paste0("Output/", name, "_historical_plots.jpg"), width = 800, height = 500, quality = 100)
  
  par(mfrow = c(2,1))
  plot(data$date, data$closing, type = "l", xlab = "date", ylab = "closing value")
  plot(data$date[-1], data$log.return, type = "l", xlab = "date", ylab = "log return value")
  title(main = paste0(name, "\n" , min(data$date), " - ", max(data$date)), outer = T, line = -3)
  
  if(save_plot) dev.off()
}

summary.data <- function(data, name, write_table = F){
  dat <- data$log.return
  df <- data.frame(Mean = mean(dat), Variance = var(dat), N = length(dat), Min = min(dat), 
                   Max = max(dat), Skewness = skewness(dat), Kurtosis = kurtosis(dat))
  if(write_table) write.csv(df, file = paste0("Output/", name, "_stats.csv"), row.names = F)
  return(df)
}

skewness <- function(data){
  mean( ( (data-mean(data))/sd(data) )^3 )
}

kurtosis <- function(data){
  mean( ( (data-mean(data))/sd(data) )^4 )
}

# Utility functions -------------------------------------------------------

simahead_exclude_inf <- function(object, n, m, theta, y){
  y <- as.matrix(y)
  theta <- as.vector(theta)
  draws <- state <- matrix(data = NA, nrow = m, ncol = n)
  
  for(j in 1:n){
    for(i in 1:m){
      rand = object$rcpp.func$rnd_Rcpp(1, theta, c(y, draws[i, 0:(j-1)]))
      while(length(which(is.infinite(rand$draws))) > 0) rand = object$rcpp.func$rnd_Rcpp(1, theta, c(y, draws[i, 0:(j-1)]))
      draws[i,j] = rand$draws
      state[i,j] = rand$state
    }
  }
  
  return(list(draws = draws,
              state = state))
}

check.alpha <- function(alpha){
  if((!is.numeric(alpha)) || (alpha < 0) || (alpha > 1)) stop("Input argument (alpha) must be numeric between 0 and 1.")
  return(alpha)
}

check.hit <- function(hit){
  if(prod(hit %in% c(T,F))) stop("Input argument (hit) should be a vector consists of only 0 or FALSE or 1 or TRUE.")
  return(hit)
}