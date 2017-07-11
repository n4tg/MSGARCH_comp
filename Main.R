rm(list = ls())

# Set up ------------------------------------------------------------------
# install.packages("devtools")
# library(devtools)
# install_github("keblu/MSGARCH/Package")
source("helping_fun.R")
source("testing_fun.R")
set_library()

# Input -------------------------------------------------------------------

input = "DAX"
roll_win_stat = 250 # window size for computing data statistics

# Models specification
models <- c("sGARCH", "eGARCH", "gjrGARCH")
distributions <- c("norm", "std")

# Bayesian parameters
N.sim = 50000 # number of draws
N.thin = 50 # every N.thin'th draws are kept
N.burn = 50000 # number of burn in draws

# VaR forecasting parameter
alpha = 0.01
taus <- c(1, 3, 10, 22) # VaR_t+taus
window_size <- 3000
interval = 1:2000 # rolling windows

# Preprocessing -----------------------------------------------------------

input.data <- read.data.log.return(input)

plot.data(input.data, input, N.trunc = (window_size+max(interval)), save_plot = T)
summary.data(input.data$LogReturn[(window_size+1):(window_size+max(interval))], input, write_table = T)
summary.data(input.data$LogReturn[(window_size+1):(window_size+max(interval))], input, stat = 'roll_win', roll_win = roll_win_stat, write_table = T)

spec <- list()
for(model in models){
  for(distr in distributions){
    spec[[paste0(model, ".", distr)]] <- MSGARCH::create.spec(model = c(model, model),
                                                              distribution = c(distr, distr),
                                                              do.skew = c(F, F),
                                                              do.mix = F,
                                                              do.shape.ind = F)
  }
}

ctr.bayes <- list(N.burn = N.burn, N.mcmc = N.sim, N.thin = N.thin, do.enhance.theta0 = T, adapt = T, acc.rate = 0.4)
specs <- names(spec)
bayes <- list()
create.outfile(input, specs, taus)

# Fit data and forecasting VaR --------------------------------------------

for(i in interval){
  print(paste0("start computation for time ", i, " out of ", max(interval)))
  y <- input.data$LogReturn[i:(window_size+i-1)]
  
  set.seed(123)
  bayes$fit <- lapply(spec, function(s) MSGARCH::fit.bayes(spec = s, y = y, ctr = ctr.bayes))
  bayes$ic <- data.frame(i, sapply(bayes$fit, function(f) list(AIC = MSGARCH::AIC(f), BIC = MSGARCH::BIC(f)) ))
  write.table(as.matrix(bayes$ic), file = paste0("Output/", input, "_AIC_BIC_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = T, col.names = F, append = T)

  #ptm = proc.time()
  print(paste0("Computing 1-period-ahead VaR at time ", i, " out of ", max(interval)))
  bayes$VaR_1 <- sapply(bayes$fit, function(f) MSGARCH::risk(object = f, level = 1-alpha, ES = F, do.its = F)$VaR)
  write.table(matrix(c(i,bayes$VaR_1), nrow = 1), file = paste0("Output/", input, "_VaR_1_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F, append = T)
  
  print(paste0("Computing multi-period-ahead VaR at time ", i, " out of ", max(interval)))
  set.seed(123)
  bayes$sim.data <- lapply(bayes$fit, function(f) MSGARCH::simahead(f, n = max(taus)))
  #proc.time()-ptm
  
  bayes$VaR_multi <- lapply(taus[-1], function(tau) matrix(c(i, tau, sapply(bayes$sim.data, function(sim) quantile(sim$draws[,tau], prob = alpha, na.rm = T))), nrow = 1))
  l_ply(bayes$VaR_multi, function(v) write.table(matrix(v[-2], nrow = 1), file = paste0("Output/", input, "_VaR_", v[2], "_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F, append = T))
}

print("END OF LOOP")

# AIC and BIC -------------------------------------------------------------

filename <- paste0("Output/", input, "_AIC_BIC_", min(interval), "_", max(interval), ".csv")
df <- read.csv(filename, sep = ";", stringsAsFactors = F)
aic <- df %>% filter(test == "AIC") %>% select(-c(test,iteration))
colMeans(aic)
aic$min <- sapply(interval, function(i) which.min(aic[i,]))
table(aic$min)

bic <- df %>% filter(test == "BIC") %>% select(-c(test,iteration))
colMeans(bic)
bic$min <- sapply(interval, function(i) which.min(bic[i,]))
table(bic$min)

# Forecasting evaluation --------------------------------------------------

VaR <- list()
VaR_t <- paste0("VaR_", taus)
filename <- paste0(input, "_VaR_", taus, "_", min(interval), "_", max(interval), ".csv")
VaR <- lapply(filename, function(fn) read.csv(paste0("Output/", fn), header = T, stringsAsFactors = F, sep = ";"))
names(VaR) <- VaR_t
VaR <- lapply(taus, function(tau) check.VaR(VaR[[paste0("VaR_", tau)]], tau = tau, input = input, interval = interval))
names(VaR) <- VaR_t

# Backtesting and MCS
l_ply(taus, function(tau){
  VaR_tau = VaR[[paste0("VaR_", tau)]]
  N.period = nrow(VaR_tau)
  r = input.data$LogReturn[(window_size+tau):(window_size+tau+N.period-1)]
  Backtesting(alpha = alpha, r = r, VaR = VaR_tau, filename = paste0("Output/", input, "_Backtest_VaR_", tau, ".csv"))
  test.MCS(r = r, VaR = VaR_tau, LossFn = MCS::LossVaR, tau = alpha, type = 'differentiable', filename = paste0("Output/", input, "_MCS_VaR_", tau, ".txt"))
  })

print("END OF PROGRAM")
