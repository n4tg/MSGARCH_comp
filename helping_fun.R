# Library handle ----------------------------------------------------------

set_lib_path_local <- function(my.lib = "Packages"){
  .libPaths(c(.libPaths(), my.lib))
}

install_lib_local <- function(package, my.lib = "Packages"){
  set_lib_path_local(my.lib)
  install.packages(package, lib = my.lib)
}

set_library <- function(){
  set_lib_path_local()
  library(MSGARCH)
  library(data.table)
  library(MCS)
  library(MASS)
  library(plyr)
  library(expm)
}

# Raw data processing -----------------------------------------------------

read.data.log.return <- function(in_file, start_date, end_date, length_lr = NULL){
  data = read.csv(paste0("Input/", in_file, ".csv"), sep = ",", stringsAsFactors = F)
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

plot.data <- function(data, name, save_plot = F){
  data$Date = as.Date(data$Date)
  
  if(save_plot) jpeg(paste0("Output/", name, "_historical_plots.jpg"), width = 800, height = 500, quality = 100)
  
  par(mfrow = c(2,1))
  plot(data$Date, data$closing, type = "l", xlab = "date", ylab = "closing value")
  plot(data$Date[-1], data$log.return, type = "l", xlab = "date", ylab = "log return value")
  title(main = paste0(name, "\n" , min(data$Date), " - ", max(data$Date)), outer = T, line = -3)
  
  if(save_plot) dev.off()
}

summary.data <- function(data, input, stat = 'all', write_table = F, ...){
  #dat = data$log.return
  opt_arg = list(...)
  
  if(stat=='all'){
    df = data.frame(Data = input, Mean = mean(data), Variance = var(data), N = length(data), 
                    Min = min(data), Max = max(data), Skewness = skewness(data), Kurtosis = kurtosis(data))
  } else if(stat=='roll_win'){
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

skewness <- function(data){
  mean( ( (data-mean(data))/sd(data) )^3 )
}

kurtosis <- function(data){
  mean( ( (data-mean(data))/sd(data) )^4 )
}

# Utility functions -------------------------------------------------------

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

check.VaR <- function(VaR, tau, save.out = T){
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

find.n.fix <- function(VaR){
  n.fix.na = which(is.na(VaR))
  n.fix.inf = which(is.infinite(as.matrix(VaR)))
  n.fix.nan = which(is.nan(as.matrix(VaR)))
  n.fix = c(n.fix.na, n.fix.inf, n.fix.nan)
  return(n.fix)
}

find.n.col.fix <- function(n.fix, N.period){
  n.fix = as.numeric(n.fix)
  N.period = as.numeric(N.period)
  
  return(ceiling(n.fix/N.period))
}

find.n.row.fix <- function(n.fix, N.period){
  n.fix = as.numeric(n.fix)
  N.period = as.numeric(N.period)
  
  n.row.fix = n.fix %% N.period
  if(n.row.fix == 0)  n.row.fix = N.period
  return(n.row.fix)
}

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

check.alpha <- function(alpha){
  alpha = as.numeric(alpha)
  if((alpha < 0) || (alpha > 1)) stop("Input argument (alpha) must be numeric between 0 and 1.")
  return(alpha)
}

check.hit <- function(hit){
  hit = as.numeric(hit)
  if(!prod(hit %in% c(0,1))) stop("Input argument (hit) should be a vector consists of only 0 or FALSE or 1 or TRUE.")
  return(hit)
}

create.outfile <- function(input, specs, taus){
  taus = as.numeric(taus)
  write.table(matrix(c("test", "iteration", specs), nrow = 1), file = paste0("Output/", input, "_AIC_BIC_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F)
  l_ply(taus, function(t) write.table(matrix(c("iteration", specs), nrow = 1), file = paste0("Output/", input, "_VaR_", t, "_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F))
}