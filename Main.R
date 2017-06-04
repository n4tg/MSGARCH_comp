rm(list = ls())

input = "DAX"
interval = 1:200

# NO CHANGES BELOW THIS LINE ----------------------------------------------

# Set up ------------------------------------------------------------------
# install.packages("devtools")
#library(devtools)
#install_github("keblu/MSGARCH/Package")
source("helping_fun.R")
source("testing_fun.R")
set_library()

# Input -------------------------------------------------------------------

# input = "DAX"
# interval = 1:200

# Models specification
models <- c("sGARCH", "eGARCH", "gjrGARCH")
# models <- c("sGARCH", "gjrGARCH")
distributions <- c("norm", "std")

# MLE parameters
# do.init = T
# itermax = 5000

# Bayesian parameters
N.sim = 50000 # number of draws
N.thin = 50 # every N.thin'th draws are kept
N.burn = 5000 # number of burn in draws

# VaR forecasting parameter
alpha = 0.01
tau <- 22 # the tau period ahead of VaR prediction
taus <- c(3, 10 ,22)
window_size = 3000

# Main --------------------------------------------------------------------

input.data <- read.data.log.return(input)

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

fit <- list()
# ctr.mle <- list(do.init = do.init, itermax = itermax)
ctr.bayes <- list(N.burn = N.burn, N.mcmc = N.sim, N.thin = N.thin, do.enhance.theta0 = T, adapt = T, acc.rate = 0.4)
specs <- names(spec)

for(i in interval){
  print(paste0("start computation for time ", i, " out of ", max(interval)))
  y <- input.data$LogReturn[i:(window_size+i-1)]
  # fit.mle <- lapply(spec, function(s) fit.mle(spec = s, y = y))
  # ctr.bayes.mle <- lapply(fit.mle, function(f) list(N.burn = N.burn, N.mcmc = N.sim, N.thin = N.thin, do.enhance.theta0 = T, adapt = TRUE, acc.rate = 0.4, theta0 = f$theta))
  # fit.bayes.mle <- lapply(specs, function(s) fit.bayes(spec = spec[[s]], y = y, ctr = ctr.bayes.mle[[s]]))
  # names(fit.bayes.mle) <- names(spec)
  
  set.seed(123)
  fit <- lapply(spec, function(s) MSGARCH::fit.bayes(spec = s, y = y, ctr = ctr.bayes))
  ic <- sapply(fit, function(f) list(AIC = MSGARCH::AIC(f), BIC = MSGARCH::BIC(f)) )

  write.table(ic, file = paste0("Output/", input, "_AIC_BIC_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = T, col.names = F, append = T)
  
  wald <- sapply(specs, function(t) test.Wald(theta_hat = matrix(colMeans(fit[[t]]$theta)), 
                                             Vn = var(fit[[t]]$theta), 
                                             R = matrix(c(rep(0 , length(fit[[t]]$theta)/N.sim -2), 1, -1), nrow = 1) ))
  
  write.table(as.matrix(wald), file = paste0("Output/", input, "_Wald_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = T, col.names = F, append = T)
  
  #ptm = proc.time()
  print(paste0("Computing 1-period-ahead VaR at time ", i, " out of ", max(interval)))
  VaR_1 <- sapply(fit, function(f) MSGARCH::risk(object = f, level = 1-alpha, ES = F, do.its = F)$VaR)
  write.table(matrix(c(i,VaR_1), nrow = 1), file = paste0("Output/", input, "_VaR_1_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F, append = T)
  
  print(paste0("Computing multi-period-ahead VaR at time ", i, " out of ", max(interval)))
  set.seed(123)
  sim.data <- lapply(fit, function(f) MSGARCH::simahead(f, n=tau))
  #proc.time()-ptm
  
  VaR_multi <- lapply(taus, function(tau) matrix(c(i, tau, sapply(sim.data, function(sim) quantile(sim$draws[,tau], prob = alpha, na.rm = T))), nrow = 1))
  l_ply(VaR_multi, function(v) write.table(matrix(v[-2], nrow = 1), file = paste0("Output/tmp_", input, "_VaR_", v[2], "_", min(interval), "_", max(interval), ".csv"), sep = ";", row.names = F, col.names = F, append = T))
} 
print("END OF LOOP")
