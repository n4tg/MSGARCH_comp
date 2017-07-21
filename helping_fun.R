# Library handle ----------------------------------------------------------

## Function: install_lib_local
# Description: Install R package to the specified folder.
# Input:  package -  A string represents an R package name.
#         my.lib - A string for folder path. Default as "Packages".
# Output: None.
install_lib_local <- function(package, my.lib = "Packages"){
  set_lib_path_local(my.lib)
  install.packages(package, lib = my.lib)
}

## Function: set_library
# Description: Call require packages for the project.
# Input:  user_lib - A logical input indicates whether the packages are 
#                    called from specified folder, default as True. If 
#                    false, packages are called from default R library.
#         my.lib - A string for the folder path, not required when user_lib 
#                  is False. Default as "Packages".
# Output: None.
set_library <- function(user_lib = T, my.lib = "Packages"){
  if(user_lib) .libPaths(c(.libPaths(), my.lib))
  library(MSGARCH)
  library(data.table)
  library(MCS)
  library(MASS)
  library(plyr)
  library(expm)
  library(dplyr)
}

# Raw data processing -----------------------------------------------------

## Function: read.data.log.return
# Description: Read input data from folder "Input" and provide date, adjusted closing prices, 
#              and log return, on a daily basis, as an output.
# Input:  input - A string indicates input data, i.e. "DAX", "SP500", or "Nikkei".
#         start_date - (Optional) A date format of start date for the output data.
#         end_date - (Optional) A date format of end date for the output data.
#         length_lr - (Optional) A numeric input specifies the fixed log-return output length.
# Output: A list of
#         + Date - A date object of the format: yyyy-mm-dd.
#         + Adj.Close - A numeric vector of daily adjusted closing price at the specific date.
#         + LogReturn - A numeric vector of daily log returns of the adjusted closing price.
read.data.log.return <- function(input, start_date, end_date, length_lr = NULL){
  data = read.csv(paste0("Input/", input, ".csv"), sep = ",", stringsAsFactors = F)
  data$Date = as.Date(data$Date)
  data = data[order(data$Date), ]
  
  if(!missing(start_date)) {
    data = data[which(data$Date >= as.Date(start_date)),]
    
    if(!missing(end_date)) {
      data = data[which(data$Date <= as.Date(end_date)),]
      if(!is.null(length_lr)) {
        warning("Parameters: start_date and end_date were given, so neglecting parameter: length_lr")
        length_lr = NULL
      }
    }
    
    if(!is.null(length_lr)) {
      if(nrow(data) <= length_lr) {
        stop(paste0("Number of input data from time ", start_date, " is not larger than the given length_lr (", length_lr, ")."))
      } else data = data[1:(length_lr+1),] # nrow(data) = length_lr+1
    }
  }
  
  if(!missing(end_date)) {
    data = data[which(data$Date <= as.Date(end_date)),]
    if(!is.null(length_lr)){
      if(nrow(data) <= length_lr) {
        stop(paste0("Number of input data up to time ", end_date, " is not larger than the given length_lr (", length_lr, ")."))
      } else data = data[(nrow(data) - length_lr):nrow(data),] # nrow(data) = length_lr+1
    }
  }
  
  return(list(Date = data$Date,
              Adj.Close = data$Adj.Close,
              LogReturn = 100*diff(log(data$Adj.Close))))
}


## Function: plot.data
# Description: Plot raw data and log-return data.
# Input:  data - A list similar to output from function "read.data.log.return".
#         input - A string indicates input data, i.e. "DAX", "SP500", or "Nikkei".
#         save_plot - A logical input states whether the graphical output will be saved or not. Default as False.
#         filename - (Optional) A path and name to which the graphical output is saved.
# Output: A figure consists of two plots, raw data and log return data at the top and the low panels, respectively.
plot.data <- function(data, input, save_plot = F, ...){
  opt_arg = list(...)
  
  data$Date = as.Date(data$Date)
  
  if(save_plot){
    if(!is.null(opt_arg$filename)){
      filename = opt_arg$filename
    } else filename = paste0("Output/", input, "_data_plots.jpg")
    jpeg(filename, width = 800, height = 500, quality = 100)
  }
  
  par(mfrow = c(2,1))
  if(!is.null(opt_arg$N.trunc)){
    N.trunc = opt_arg$N.trunc
    data$Date = data$Date[1:(N.trunc+1)]
    data$Adj.Close = data$Adj.Close[1:(N.trunc+1)]
    data$LogReturn = data$LogReturn[1:N.trunc]
  }
  plot(data$Date, data$Adj.Close, type = "l", xlab = "date", ylab = "Adjusted closing value")
  plot(data$Date[-1], data$LogReturn, type = "l", xlab = "date", ylab = "Log return value")
  title(main = paste0(input, "\n" , min(data$Date), " - ", max(data$Date)), outer = T, line = -3)
  
  if(save_plot) dev.off()
}

## Function: summary.data
# Description: Display data statistics as a table.
# Input:  data - A list similar to output from function "read.data.log.return".
#         input - A string indicates input data, i.e. "DAX", "SP500", or "Nikkei".
#         stat - A string determines whether the data statistics are calculated for the whole period ("all"),
#                or based on rolling windows ("roll_win"). Default as stat = "all".
#         write_table - A logical input states whether the output table will be saved or not. Default as False.
#         ... - optional argument:
#           + filename - A path and name to which the output table is saved.
#           + roll_win - A numeric input indicates window length, required when stat = "roll_win".
# Output: A dataframe consists of input string and its following statistics: mean, variance, length, 
#         minimum, maximum, skewness, and kurtosis.
summary.data <- function(data, input, stat = "all", write_table = F, ...){
  opt_arg = list(...)
  
  if(stat=="all"){
    df = data.frame(Data = input, Mean = mean(data), Variance = var(data), N = length(data), 
                    Min = min(data), Max = max(data), Skewness = skewness(data), Kurtosis = kurtosis(data))
  } else if(stat=="roll_win"){
    if(is.null(opt_arg$roll_win)) stop("Please specify the rolling window length (roll_win) for statistics calculation; otherwise, use stat = 'all'.")
    N.data = length(data)
    roll_win = opt_arg$roll_win
    df = lapply(1:(N.data-roll_win), 
                function(i){
                  dat = data[i:(i+roll_win-1)]
                  data.frame(Mean = mean(dat), Variance = var(dat), Skewness = skewness(dat), Kurtosis = kurtosis(dat))}
                )
    df = rbindlist(df, use.names = T, fill = T)
    df = cbind(Data = input, data.frame(t(colMeans(df))))
  }
  
  if(write_table) {
    if(!is.null(opt_arg$filename)){
      filename = opt_arg$filename
    } else filename = paste0("Output/", input, "_stats_", stat, ".csv")
    write.table(df, file = filename, sep = ";", row.names = F)
  }
  return(df)
}


## Function: skewness
# Description: Calculate skewness of data.
# Input:  data - A numeric vector from which skewness is calculated.
# Output: A number represents skewness.
skewness <- function(data){
  mean( ( (data-mean(data))/sd(data) )^3 )
}

## Function: kurtosis
# Description: Calculate kurtosis of data.
# Input:  data - A numeric vector from which kurtosis is calculated.
# Output: A number represents kurtosis.
kurtosis <- function(data){
  mean( ( (data-mean(data))/sd(data) )^4 )
}

# Utility functions -------------------------------------------------------

## Function: simahead_exclude_inf
# Description: Re-simulate the draws that are non-numerics (NA, NaN) or (-/+)Inf.
# Input:  object, n, m, theta, y as in MSGARCH::simahead. See MSGARCH package manual.
# Output: A list of draws and states, in which all draws are valid numeric.
simahead_exclude_inf <- function(object, n, m, theta, y){
  y = as.matrix(y)
  theta = as.vector(theta)
  draws <- state <- matrix(data = NA, nrow = m, ncol = n)
  
  for(j in 1:n){
    print(paste0("Random draws for time ", j, " out of ", n))
    for(i in 1:m){
      rand = object$rcpp.func$rnd_Rcpp(1, theta, c(y, draws[i, 0:(j-1)]))
      tmp.it = 0
      while(length(which(is.infinite(rand$draws)))+length((which(is.nan(rand$draws)))) > 0) {
        rand = object$rcpp.func$rnd_Rcpp(1, theta, c(y, draws[i, 0:(j-1)]))
        tmp.it = tmp.it + 1
        if(tmp.it > 10000){
          rand$draws = NA
          break
        }
      }
      draws[i,j] = rand$draws
      state[i,j] = rand$state
    }
  }
  
  return(list(draws = draws,
              state = state))
}

## Function: check.VaR
# Description: Check whether the VaR is of the correct form. If not, 
#              e.g. contains NA, NaN, (+/-)Inf, correct it.
# Input:  VaR - A vector or matrix of VaR.
#         tau - A number indicates forecasting steps for VaR prediction.
#         save.out - A logical input indicates whether the corrected VaR is saved. Default as True.
#         input - A string indicates input data, i.e. "DAX", "SP500", or "Nikkei".
#         interval - A numeric vector of time index.
# Output: A numeric vector or matrix of corrected VaR.
check.VaR <- function(VaR, tau, save.out = T, input, interval){
  tau = is.numeric(tau)
  
  N.period = nrow(VaR)
  n.fix = find.n.fix(VaR)
  n.err = which(which(abs(VaR) > 100) > N.period)
  n.fix = c(n.fix, n.err)
  
  for(i in n.fix){
    print(paste0("Start fix ", input, ", tau = ", tau, ", n.fix = ", i, " (", which(n.fix == i), " of ", length(n.fix), ")"))
    n.col.fix = find.n.col.fix(i, N.period)
    spec.fix = spec[[colnames(VaR[n.col.fix])]]
    
    n.row.fix = find.n.row.fix(i, N.period)
    
    y = input.data$LogReturn[(n.row.fix-tau+1):(n.row.fix+window_size-tau)]
    
    print(paste0("Fit ", input, ", tau = ", tau, ", n.fix = ", i, " (", which(n.fix == i), " of ", length(n.fix), ")"))
    fit = fit.bayes(spec = spec.fix, y = y, ctr = ctr.bayes)
    
    print(paste0("Start draw ", input, ", tau = ", tau, ", n.fix = ", i, " (", which(n.fix == i), " of ", length(n.fix), ")"))
    draws = simahead_exclude_inf(spec.fix, n = tau, m = N.sim, theta = colMeans(fit$theta), y = y)$draws
    print(paste0("Finish draw ", input, ", tau = ", tau, ", n.fix = ", i, " (", which(n.fix == i), " of ", length(n.fix), ")"))
    VaR[n.row.fix, n.col.fix] = quantile(draws[, tau], na.rm = T, probs = alpha, names = F)
    if(save.out) write.table(VaR, file = paste0("Output/", input, "_VaR_", tau, "_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F)
  }
  
  VaR$iteration = NULL
  return(VaR)
}

## Function: find.n.fix
# Description: Find NA, NaN, (+/-)Inf, in the given vector/matrix. 
# Input:  VaR - A vector or matrix of interest.
# Output: A numeric vector or a number indicates the location of NA, NaN, (+/-)Inf.
find.n.fix <- function(VaR){
  n.fix.na = which(is.na(VaR))
  n.fix.inf = which(is.infinite(as.matrix(VaR)))
  n.fix.nan = which(is.nan(as.matrix(VaR)))
  n.fix = c(n.fix.na, n.fix.inf, n.fix.nan)
  return(n.fix)
}

## Function: find.n.col.fix
# Description: Find column of the given location and the total row size.
# Input:  n.fix - A number of location (e.g. output from function "find.n.fix").
#         N.period - A number of total row size of the considered matrix.
# Output: A number indicates the column index, regarding the given location and total row size.
find.n.col.fix <- function(n.fix, N.period){
  n.fix = as.numeric(n.fix)
  N.period = as.numeric(N.period)
  
  return(ceiling(n.fix/N.period))
}

## Function: find.n.row.fix
# Description: Find row of the given location and the total row size.
# Input:  n.fix - A number of location (e.g. output from function "find.n.fix").
#         N.period - A number of total row size of the considered matrix.
# Output: A number indicates the row index, regarding the given location and total row size.
find.n.row.fix <- function(n.fix, N.period){
  n.fix = as.numeric(n.fix)
  N.period = as.numeric(N.period)
  
  n.row.fix = n.fix %% N.period
  if(n.row.fix == 0)  n.row.fix = N.period
  return(n.row.fix)
}

## Function: check.Mat
# Description: Check if the input matrix contains NA, NaN, or (+/-)Inf. 
#              If yes, change NA, NaN to 0 and (+/-)Inf to (+/-)1e6.
# Input:  mat - A matrix.
# Output: A matrix without NA, NaN, (+/-)Inf.
check.Mat <- function(mat){
  n.period = nrow(mat)
  n.fix = find.n.fix(mat)
  if(length(n.fix) > 0){
    n.col.fix = find.n.col.fix(n.fix, n.period)
    n.row.fix = find.n.row.fix(n.fix, n.period)
    
    if(is.na(mat[n.row.fix,n.col.fix]) || is.nan(mat[n.row.fix,n.col.fix])) mat[n.row.fix,n.col.fix] = 0
    if(is.infinite(mat[n.row.fix,n.col.fix])){
      if(mat[n.row.fix,n.col.fix] < 0){mat[n.row.fix,n.col.fix] = -1e6} else{mat[n.row.fix,n.col.fix] = 1e6}
    }
  }
  
  return(mat)
}

## Function: check.alpha
# Description: Check if the input (alpha) is the valid probability: 0<=alpha<=1. 
# Input:  alpha - A number or character indicates (shortfall) probability.
# Output: A number of (shortfall) probability.
check.alpha <- function(alpha){
  alpha = as.numeric(alpha)
  if((alpha < 0) || (alpha > 1)) stop("Input argument (alpha) must be numeric between 0 and 1.")
  return(alpha)
}

## Function: check.hit
# Description: Check if the input vector (hit) contains only 0 or 1 (False or True). 
# Input:  hit - A considered vector.
# Output: A numeric vector consists of 0 or 1.
check.hit <- function(hit){
  hit = as.numeric(hit)
  if(!prod(hit %in% c(0,1))) stop("Input argument (hit) should be a vector consists of only 0 or FALSE or 1 or TRUE.")
  return(hit)
}

## Function: create.outfile
# Description: Create an output file for the criteria tests and VaR forecasting.
# Input:  input - A string indicates input data, i.e. "DAX", "SP500", or "Nikkei".
#         specs - A vector of strings consist of all competitor models.
#         taus - A numeric vector or a number of the forecasting period.
# Output: None.
create.outfile <- function(input, specs, taus){
  taus = as.numeric(taus)
  write.table(matrix(c("test", "iteration", specs), nrow = 1), file = paste0("Output/", input, "_AIC_BIC_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F)
  l_ply(taus, function(t) write.table(matrix(c("iteration", specs), nrow = 1), file = paste0("Output/", input, "_VaR_", t, "_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F))
}
