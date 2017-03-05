rm(list = ls())

library(MSGARCH)
library(moments)
# library(ggplot2)
# library(gridExtra)

source("helping_fun.R")

input <- "DAX"
# input <- "S&P500"
# input <- "Nikkei"

data <- read.data.closing(input)

# statistical properties
plot.data(data, name = input, save_plot = T)
summary.data(data, name = input, write_table = T)
