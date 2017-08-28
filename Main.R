####################################################
# Start MAIN.R of MSGARCH_comp
####################################################

rm(list = ls()) # clear the environment before the script starts

# Set up ------------------------------------------------------------------
source("helping_fun.R") # call functions in the file helping_fun.R
source("testing_fun.R") # call functions in the file testing_fun.R
set_library() 
############## END OF SECTION "Set up" #############

# Input -------------------------------------------------------------------
input = "DAX" # choice of "DAX", "SP500", and "Nikkei"
roll_win_stat = 250 # window size for computing data statistics
window_size = 3000 # window length
interval = 1:2000 # rolling windows index

# Model specifications
models <- c("sGARCH", "eGARCH", "gjrGARCH")
distributions <- c("norm", "std")

# Bayesian parameters
N.sim = 50000 # number of draws
N.thin = 50 # draws are kept at every N.thin'th
N.burn = 50000 # number of burn-in draws

# VaR forecasting parameter
alpha = 0.01 # shortfall probability
taus <- c(1, 3, 10, 22) # forecasting horizon: VaR_t+taus
############## END OF SECTION "Input" #############

# Preprocessing -----------------------------------------------------------
input.data <- read.data.log.return(input)

plot.data(input.data, input, N.trunc = (window_size+max(interval)), save_plot = T)
summary.data(input.data$LogReturn[(window_size+1):(window_size+max(interval))], input, write_table = T)
summary.data(input.data$LogReturn[(window_size+1):(window_size+max(interval))], input, stat = 'roll_win', roll_win = roll_win_stat, write_table = T)

spec <- list() # create an empty list to keep specification objects for each model
for(model in models){
  for(distr in distributions){ 
    
    # create models specification objects
    spec[[paste0(model, ".", distr)]] <- MSGARCH::create.spec(model = c(model, model),
                                                              distribution = c(distr, distr), do.skew = c(F, F), 
                                                              do.mix = F, do.shape.ind = F)
  }
}
ctr.bayes <- list(N.burn = N.burn, N.mcmc = N.sim, N.thin = N.thin, do.enhance.theta0 = T, adapt = T, acc.rate = 0.4) # a list of control parameters for Bayesian estimation
specs <- names(spec) # get all specifications name
bayes <- list() # create an empty list to keep model estimation and forecasting results
create.outfile(input, specs, taus)
############ END OF SECTION "Preprocessing" ###########

# Fit data and forecasting VaR --------------------------------------------
for(i in interval){
  # show the computing status on the console panel
  print(paste0("start computation for time ", i, " out of ", max(interval))) 
  
  # select the input log return data in according to relevant window
  y <- input.data$LogReturn[i:(window_size+i-1)]   
  
  set.seed(123) # set seed to allow reproduction of results
  # estimated models by Bayesian for each model spec regarding to the input log return data
  bayes$fit <- lapply(spec, function(s) MSGARCH::fit.bayes(spec = s, y = y, ctr = ctr.bayes))
  
  # save the AIC, BIC of each models fitting
  bayes$ic <- data.frame(i, sapply(bayes$fit, function(f) list(AIC = MSGARCH::AIC(f), BIC = MSGARCH::BIC(f)) ))
  write.table(as.matrix(bayes$ic), file = paste0("Output/", input, "_AIC_BIC_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = T, col.names = F, append = T)
  
  # show the computing status on the console panel
  print(paste0("Computing 1-period-ahead VaR at time ", i, " out of ", max(interval)))
  
  # compute and save the 1-period-ahead VaR with the alpha shortfall probability
  bayes$VaR_1 <- sapply(bayes$fit, function(f) MSGARCH::risk(object = f, level = 1-alpha, ES = F, do.its = F)$VaR)
  write.table(matrix(c(i,bayes$VaR_1), nrow = 1), file = paste0("Output/", input, "_VaR_1_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F, append = T)
  
  # show the computing status on the console panel
  print(paste0("Computing multi-period-ahead VaR at time ", i, " out of ", max(interval)))
  
  set.seed(123) # set seed to allow reproduction of result
  # simulate the data for multi-period VaR computation
  bayes$sim.data <- lapply(bayes$fit, function(f) MSGARCH::simahead(f, n = max(taus)))
  
  # compute and save the 3-, 10-, and 22-period-ahead VaR (alpha-quantile)
  bayes$VaR_multi <- lapply(taus[-1], function(tau) matrix(c(i, tau, sapply(bayes$sim.data, function(sim) quantile(sim$draws[,tau], prob = alpha, na.rm = T))), nrow = 1))
  l_ply(bayes$VaR_multi, function(v) write.table(matrix(v[-2], nrow = 1), file = paste0("Output/", input, "_VaR_", v[2], "_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F, append = T))
}
print("END OF LOOP") # show the computing status on the console panel
####### END OF SECTION "Fit data and forecasting VaR" ######

# AIC and BIC -------------------------------------------------------------
# read in the AIC an BIC results
filename <- paste0("Output/", input, "_AIC_BIC_", min(interval), "_", max(interval), ".csv")
df <- read.csv(filename, sep = ";", stringsAsFactors = F)
# select AIC values only and keep them in a data frame called aic
aic <- df %>% filter(test == "AIC") %>% select(-c(test,iteration))
# find the average AIC for each model
colMeans(aic) 
# find the model that has minimum AIC per each rolling window
aic$min <- sapply(interval, function(i) which.min(aic[i,]))
table(aic$min)

# select BC values only and keep them in a data frame called aic
bic <- df %>% filter(test == "BIC") %>% select(-c(test,iteration))
# find the average BIC for each model
colMeans(bic)
# find the model that has minimum BIC per each rolling window
bic$min <- sapply(interval, function(i) which.min(bic[i,]))
table(bic$min)
############ END OF SECTION "AIC and BIC" ###########

# Forecasting evaluation --------------------------------------------------
VaR <- list() # create an empty list to keep VaR results
VaR_t <- paste0("VaR_", taus) # a string vector that keeps VaR forecasting names
# create a string vector that keeps the output filenames for each VaR forecasting result
filename <- paste0(input, "_VaR_", taus, "_", min(interval), "_", max(interval), ".csv")
# read the forecasted VaR in as a list of data frame
VaR <- lapply(filename, function(fn) read.csv(paste0("Output/", fn), header = T, stringsAsFactors = F, sep = ";"))
names(VaR) <- VaR_t # change the list name
# check the forecasted VaR at each forecasting horizon and fix it if error e.g. NA, NaN
VaR <- lapply(taus, function(tau) check.VaR(VaR[[paste0("VaR_", tau)]], tau = tau, input = input, interval = interval))
names(VaR) <- VaR_t # change the list name

# Backtesting and MCS
l_ply(taus, function(tau){
  
  VaR_tau = VaR[[paste0("VaR_", tau)]] # Select the relevant VaR
  N.period = nrow(VaR_tau) # find the total period considered in VaR forecasting
  
  # Select the relevant log return input data
  r = input.data$LogReturn[(window_size+tau):(window_size+tau+N.period-1)]
  
  # Evaluation the forecasting performances by backtesting and MCS
  Backtesting(alpha = alpha, r = r, VaR = VaR_tau, filename = paste0("Output/", input, "_Backtest_VaR_", tau, ".csv"))
  test.MCS(r = r, VaR = VaR_tau, LossFn = MCS::LossVaR, tau = alpha, type = 'differentiable', filename = paste0("Output/", input, "_MCS_VaR_", tau, ".txt"))})
####### END OF SECTION "Forecasting evaluation" ######

print("END OF PROGRAM") # show the computing status on the console panel

####################################################
# End MAIN.R of MSGARCH_comp
####################################################
