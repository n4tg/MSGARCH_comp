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

Wald.test <- function(theta_hat,
                      Vn,
                      R,
                      r = 0){
  W <- t(R %*% theta_hat - r) %*% solve(R %*% Vn %*% t(R)) %*% (R %*% theta_hat - r)
  df <- nrow(R)
  pval <- 1-pchisq(as.numeric(W), df)
  return(data.frame(W = W, df = df, p = pval))
}