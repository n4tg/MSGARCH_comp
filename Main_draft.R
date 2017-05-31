# Set up ------------------------------------------------------------------
install.packages("devtools")
library(devtools)
install_github("keblu/MSGARCH/Package")
source("helping_fun.R")
source("testing_fun.R")
set_library()

# Input -------------------------------------------------------------------

input = "DAX"
interval = 1:500

# VaR forecasting parameter
alpha = 0.01
tau = 22 # the tau period ahead of VaR prediction
window_size = 5000

# Bayesian parameters
N.sim = 10000 # number of draws
N.thin = 50 # every N.thin'th draws are kept
N.burn = 5000 # number of burn in draws

# Main --------------------------------------------------------------------
N.data <- 

# Preprocessing
input.data <- read.data.closing(input, )

for(i in interval){
  
}