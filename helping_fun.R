read.data.closing <- function(in_file, length_data, head_data = T){
  data <- read.csv(paste0("Input/", in_file, ".csv"), sep = ",", stringsAsFactors = F)
  if(!missing(length_data)){
    if(nrow(data)<=length_data) stop("Input file has smaller length than specified by length_data.")
    
    if(head_data){
      data <- data[1:length_data+1,]
    } else {
      data <- data[(nrow(data)-length_data+1):nrow(data),]
    }
  }
  
  return(list(date = as.Date(data$Date),
              closing = data$Close,
              log.return = diff(log(data$Close))))
}

plot.data <- function(data, name, save_plot = F){
  if(save_plot) jpeg(paste0("Output/", name, "_historical_plots.jpg"), width = 800, height = 500, quality = 100)
  
  par(mfrow = c(1,2))
  plot(data$date, data$closing, type = "l", xlab = "date", ylab = "closing value")
  plot(data$date, data$log.return, type = "l", xlab = "date", ylab = "log return value")
  title(main = paste0(name, "\n" , min(data$date), " - ", max(data$date)), outer = T, line = -3)
  
  if(save_plot) dev.off()
}

summary.data <- function(data, name, write_table = F){
  dat <- data$log.return
  df <- data.frame(Mean = mean(dat), Variance = var(dat), N = length(dat), Min = min(dat), 
                   Max = max(dat), Skewness = skewness(dat), Kurtosis = kurtosis(dat))
  if(write_table) write.csv(df, file = paste0("Output/", name, "_stats.csv"), row.names = F)
}