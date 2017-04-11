rm(list = ls())

library(MSGARCH)
#library(moments)
#library(numDeriv)
#library(Matrix)
# library(ggplot2)
# library(gridExtra)

source("helping_fun.R")

input <- "DAX"
#input <- "SP500"
#input <- "Nikkei"

window_size <- 3000
N.data <- 5000
N.sim <- 5000
N.mcmc <- 10000
N.thin <- 50
N.burn <- 10000

alpha <- 0.01

#data <- read.data.closing(input)

dax <- read.data.closing("DAX", N.data)
sp <- read.data.closing("SP500", N.data)
nikkei <- read.data.closing("Nikkei", N.data)

data <- dax

# statistical properties
plot.data(data, name = input, save_plot = F)
summary.data(data, name = input, write_table = F)

# MSGARCH

spec.sGARCH.norm <- create.spec(model = c("sGARCH", "sGARCH"), distribution = c("norm", "norm"))
spec.sGARCH.skewnorm <- create.spec(model = c("sGARCH", "sGARCH"), distribution = c("norm", "norm"), do.skew = c(T,T))
spec.sGARCH.t <- create.spec(model = c("sGARCH", "sGARCH"), distribution = c("std", "std"))
spec.eGARCH.norm <- create.spec(model = c("eGARCH", "eGARCH"), distribution = c("norm", "norm"))
spec.eGARCH.skewnorm <- create.spec(model = c("eGARCH", "eGARCH"), distribution = c("norm", "norm"), do.skew = c(T,T))
spec.eGARCH.t <- create.spec(model = c("eGARCH", "eGARCH"), distribution = c("std", "std"))
spec.gjrGARCH.norm <- create.spec(model = c("gjrGARCH", "gjrGARCH"), distribution = c("norm", "norm"))
spec.gjrGARCH.skewnorm <- create.spec(model = c("gjrGARCH", "gjrGARCH"), distribution = c("norm", "norm"), do.skew = c(T,T))
spec.gjrGARCH.t <- create.spec(model = c("gjrGARCH", "gjrGARCH"), distribution = c("std", "std"))

fit_test = data.frame()
VaR <- data.frame()
violation <- data.frame()
VaR_3 <- data.frame()
violation_3 <- data.frame()
VaR_10 <- data.frame()
violation_10 <- data.frame()
VaR_22 <- data.frame()
violation_22 <- data.frame()

i = 1

## start loop

for(i in 1:(N.data-1)){
  print(i)
  set.seed(i)
  
  fit.mle.sGARCH.norm <- fit.mle(spec = spec.sGARCH.norm, data$log.return[i:(window_size+i-1)])
  fit.mle.sGARCH.skewnorm <- fit.mle(spec = spec.sGARCH.skewnorm, data$log.return[i:(window_size+i-1)])
  fit.mle.sGARCH.t <- fit.mle(spec = spec.sGARCH.t, data$log.return[i:(window_size+i-1)])
  fit.mle.eGARCH.norm <- fit.mle(spec = spec.eGARCH.norm, data$log.return[i:(window_size+i-1)])
  fit.mle.eGARCH.skewnorm <- fit.mle(spec = spec.eGARCH.skewnorm, data$log.return[i:(window_size+i-1)])
  fit.mle.eGARCH.t <- fit.mle(spec = spec.eGARCH.t, data$log.return[i:(window_size+i-1)])
  fit.mle.gjrGARCH.norm <- fit.mle(spec = spec.gjrGARCH.norm, data$log.return[i:(window_size+i-1)])
  fit.mle.gjrGARCH.skewnorm <- fit.mle(spec = spec.gjrGARCH.skewnorm, data$log.return[i:(window_size+i-1)])
  fit.mle.gjrGARCH.t <- fit.mle(spec = spec.gjrGARCH.t, data$log.return[i:(window_size+i-1)])
  
  fit.bayes.sGARCH.norm <- fit.bayes(spec = spec.sGARCH.norm, data$log.return[i:(window_size+i-1)], ctr = list(N.burn = N.burn, N.mcmc = N.mcmc, N.thin = N.thin))
  fit.bayes.sGARCH.skewnorm <- fit.bayes(spec = spec.sGARCH.skewnorm, data$log.return[i:(window_size+i-1)], ctr = list(N.burn = N.burn, N.mcmc = N.mcmc, N.thin = N.thin))
  fit.bayes.sGARCH.t <- fit.bayes(spec = spec.sGARCH.t, data$log.return[i:(window_size+i-1)], ctr = list(N.burn = N.burn, N.mcmc = N.mcmc, N.thin = N.thin))
  fit.bayes.eGARCH.norm <- fit.bayes(spec = spec.eGARCH.norm, data$log.return[i:(window_size+i-1)], ctr = list(N.burn = N.burn, N.mcmc = N.mcmc, N.thin = N.thin))
  fit.bayes.eGARCH.skewnorm <- fit.bayes(spec = spec.eGARCH.skewnorm, data$log.return[i:(window_size+i-1)], ctr = list(N.burn = N.burn, N.mcmc = N.mcmc, N.thin = N.thin))
  fit.bayes.eGARCH.t <- fit.bayes(spec = spec.eGARCH.t, data$log.return[i:(window_size+i-1)], ctr = list(N.burn = N.burn, N.mcmc = N.mcmc, N.thin = N.thin))
  fit.bayes.gjrGARCH.norm <- fit.bayes(spec = spec.gjrGARCH.norm, data$log.return[i:(window_size+i-1)], ctr = list(N.burn = N.burn, N.mcmc = N.mcmc, N.thin = N.thin))
  fit.bayes.gjrGARCH.skewnorm <- fit.bayes(spec = spec.gjrGARCH.skewnorm, data$log.return[i:(window_size+i-1)], ctr = list(N.burn = N.burn, N.mcmc = N.mcmc, N.thin = N.thin))
  fit.bayes.gjrGARCH.t <- fit.bayes(spec = spec.gjrGARCH.t, data$log.return[i:(window_size+i-1)], ctr = list(N.burn = N.burn, N.mcmc = N.mcmc, N.thin = N.thin))
  
  fit_test <- rbind(fit_test,
                    c(AIC(fit.mle.sGARCH.norm), AIC(fit.mle.sGARCH.skewnorm), AIC(fit.mle.sGARCH.t),
                      AIC(fit.mle.eGARCH.norm), AIC(fit.mle.eGARCH.skewnorm), AIC(fit.mle.eGARCH.t),
                      AIC(fit.mle.gjrGARCH.norm), AIC(fit.mle.gjrGARCH.skewnorm), AIC(fit.mle.gjrGARCH.t),
                      BIC(fit.mle.sGARCH.norm), BIC(fit.mle.sGARCH.skewnorm), BIC(fit.mle.sGARCH.t),
                      BIC(fit.mle.eGARCH.norm), BIC(fit.mle.eGARCH.skewnorm), BIC(fit.mle.eGARCH.t),
                      BIC(fit.mle.gjrGARCH.norm), BIC(fit.mle.gjrGARCH.skewnorm), BIC(fit.mle.gjrGARCH.t),
                      AIC(fit.bayes.sGARCH.norm), AIC(fit.bayes.sGARCH.skewnorm), AIC(fit.bayes.sGARCH.t),
                      AIC(fit.bayes.eGARCH.norm), AIC(fit.bayes.eGARCH.skewnorm), AIC(fit.bayes.eGARCH.t),
                      AIC(fit.bayes.gjrGARCH.norm), AIC(fit.bayes.gjrGARCH.skewnorm), AIC(fit.bayes.gjrGARCH.t),
                      BIC(fit.bayes.sGARCH.norm), BIC(fit.bayes.sGARCH.skewnorm), BIC(fit.bayes.sGARCH.t),
                      BIC(fit.bayes.eGARCH.norm), BIC(fit.bayes.eGARCH.skewnorm), BIC(fit.bayes.eGARCH.t),
                      BIC(fit.bayes.gjrGARCH.norm), BIC(fit.bayes.gjrGARCH.skewnorm), BIC(fit.bayes.gjrGARCH.t),
                      DIC(fit.bayes.sGARCH.norm)$DIC, DIC(fit.bayes.sGARCH.skewnorm)$DIC, DIC(fit.bayes.sGARCH.t)$DIC,
                      DIC(fit.bayes.eGARCH.norm)$DIC, DIC(fit.bayes.eGARCH.skewnorm)$DIC, DIC(fit.bayes.eGARCH.t)$DIC,
                      DIC(fit.bayes.gjrGARCH.norm)$DIC, DIC(fit.bayes.gjrGARCH.skewnorm)$DIC, DIC(fit.bayes.gjrGARCH.t)$DIC))
  
  VaR <- rbind(VaR,
               c(risk(fit.mle.sGARCH.norm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.mle.sGARCH.skewnorm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.mle.sGARCH.t, level = 1-alpha, ES = F)$VaR,
                 risk(fit.mle.eGARCH.norm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.mle.eGARCH.skewnorm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.mle.eGARCH.t, level = 1-alpha, ES = F)$VaR,
                 risk(fit.mle.gjrGARCH.norm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.mle.gjrGARCH.skewnorm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.mle.gjrGARCH.t, level = 1-alpha, ES = F)$VaR,
                 risk(fit.bayes.sGARCH.norm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.bayes.sGARCH.skewnorm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.bayes.sGARCH.t, level = 1-alpha, ES = F)$VaR,
                 risk(fit.bayes.eGARCH.norm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.bayes.eGARCH.skewnorm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.bayes.eGARCH.t, level = 1-alpha, ES = F)$VaR,
                 risk(fit.bayes.gjrGARCH.norm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.bayes.gjrGARCH.skewnorm, level = 1-alpha, ES = F)$VaR,
                 risk(fit.bayes.gjrGARCH.t, level = 1-alpha, ES = F)$VaR))
  
  violation <- rbind(violation, (data$log.return[window_size+i] < VaR))
  
  if(i <= N.data-3){
    VaR_3 <- rbind(VaR_3,
                   c(quantile(rowSums(simahead(fit.mle.sGARCH.norm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.mle.sGARCH.skewnorm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.mle.sGARCH.t, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.mle.eGARCH.norm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.mle.eGARCH.skewnorm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.mle.eGARCH.t, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.mle.gjrGARCH.norm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.mle.gjrGARCH.skewnorm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.mle.gjrGARCH.t, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.bayes.sGARCH.norm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.bayes.sGARCH.skewnorm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.bayes.sGARCH.t, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.bayes.eGARCH.norm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.bayes.eGARCH.skewnorm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.bayes.eGARCH.t, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.bayes.gjrGARCH.norm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.bayes.gjrGARCH.skewnorm, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                     quantile(rowSums(simahead(fit.bayes.gjrGARCH.t, n = 3, m = N.sim)$draws, na.rm = T), alpha, na.rm = T)))
    
    violation_3 <- rbind(violation_3, (data$log.return[window_size+i+2] < VaR_3))
  }
  
  if(i <= N.data-10){
    VaR_10 <- rbind(VaR_10,
                    c(quantile(rowSums(simahead(fit.mle.sGARCH.norm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.sGARCH.skewnorm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.sGARCH.t, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.eGARCH.norm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.eGARCH.skewnorm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.eGARCH.t, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.gjrGARCH.norm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.gjrGARCH.skewnorm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.gjrGARCH.t, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.sGARCH.norm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.sGARCH.skewnorm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.sGARCH.t, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.eGARCH.norm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.eGARCH.skewnorm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.eGARCH.t, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.gjrGARCH.norm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.gjrGARCH.skewnorm, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.gjrGARCH.t, n = 10, m = N.sim)$draws, na.rm = T), alpha, na.rm = T)))
    
    violation_10 <- rbind(violation_10, (data$log.return[window_size+i+9] < VaR_10))
  }
  
  if(i <= N.data-22){
    VaR_22 <- rbind(VaR_22,
                    c(quantile(rowSums(simahead(fit.mle.sGARCH.norm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.sGARCH.skewnorm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.sGARCH.t, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.eGARCH.norm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.eGARCH.skewnorm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.eGARCH.t, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.gjrGARCH.norm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.gjrGARCH.skewnorm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.mle.gjrGARCH.t, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.sGARCH.norm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.sGARCH.skewnorm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.sGARCH.t, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.eGARCH.norm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.eGARCH.skewnorm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.eGARCH.t, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.gjrGARCH.norm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.gjrGARCH.skewnorm, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T),
                      quantile(rowSums(simahead(fit.bayes.gjrGARCH.t, n = 22, m = N.sim)$draws, na.rm = T), alpha, na.rm = T)))
    
    violation_22 <- rbind(violation_22, (data$log.return[window_size+i+21] < VaR_22))
  }
}

## end loop
names(violation) <- names(VaR) <- names(violation_3) <- names(VaR_3) <-
  names(violation_10) <- names(VaR_10) <- names(violation_22) <- names(VaR_22) <-
  c("fit.mle.sGARCH.norm", "fit.mle.sGARCH.skewnorm", "fit.mle.sGARCH.t",
    "fit.mle.eGARCH.norm", "fit.mle.eGARCH.skewnorm", "fit.mle.eGARCH.t",
    "fit.mle.gjrGARCH.norm", "fit.mle.gjrGARCH.skewnorm", "fit.mle.gjrGARCH.t",
    "fit.bayes.sGARCH.norm", "fit.bayes.sGARCH.skewnorm", "fit.bayes.sGARCH.t",
    "fit.bayes.eGARCH.norm", "fit.bayes.eGARCH.skewnorm", "fit.bayes.eGARCH.t",
    "fit.bayes.gjrGARCH.norm", "fit.bayes.gjrGARCH.skewnorm", "fit.bayes.gjrGARCH.t")

names(fit_test) <- c("AIC.mle.sGARCH.norm", "AIC.mle.sGARCH.skewnorm", "AIC.mle.sGARCH.t",
                     "AIC.mle.eGARCH.norm", "AIC.mle.eGARCH.skewnorm", "AIC.mle.eGARCH.t",
                     "AIC.mle.gjrGARCH.norm", "AIC.mle.gjrGARCH.skewnorm", "AIC.mle.gjrGARCH.t",
                     "BIC.mle.sGARCH.norm", "BIC.mle.sGARCH.skewnorm", "BIC.mle.sGARCH.t",
                     "BIC.mle.eGARCH.norm", "BIC.mle.eGARCH.skewnorm", "BIC.mle.eGARCH.t",
                     "BIC.mle.gjrGARCH.norm", "BIC.mle.gjrGARCH.skewnorm", "BIC.mle.gjrGARCH.t",
                     "AIC.bayes.sGARCH.norm", "AIC.bayes.sGARCH.skewnorm", "AIC.bayes.sGARCH.t",
                     "AIC.bayes.eGARCH.norm", "AIC.bayes.eGARCH.skewnorm", "AIC.bayes.eGARCH.t",
                     "AIC.bayes.gjrGARCH.norm", "AIC.bayes.gjrGARCH.skewnorm", "AIC.bayes.gjrGARCH.t",
                     "BIC.bayes.sGARCH.norm", "BIC.bayes.sGARCH.skewnorm", "BIC.bayes.sGARCH.t",
                     "BIC.bayes.eGARCH.norm", "BIC.bayes.eGARCH.skewnorm", "BIC.bayes.eGARCH.t",
                     "BIC.bayes.gjrGARCH.norm", "BIC.bayes.gjrGARCH.skewnorm", "BIC.bayes.gjrGARCH.t",
                     "DIC.bayes.sGARCH.norm", "DIC.bayes.sGARCH.skewnorm", "DIC.bayes.sGARCH.t",
                     "DIC.bayes.eGARCH.norm", "DIC.bayes.eGARCH.skewnorm", "DIC.bayes.eGARCH.t",
                     "DIC.bayes.gjrGARCH.norm", "DIC.bayes.gjrGARCH.skewnorm", "DIC.bayes.gjrGARCH.t")

## Automation

# models <- c("sGARCH", "eGARCH", "gjrGARCH")
# distr <- c("norm", "std")

# in.sample.stat <- data.frame(model = character(1), dist = character(1),
#                              alpha0_1 = numeric(1), alpha1_1 = numeric(1),
#                              alpha2_1 = numeric(1), beta_1 = numeric(1), 
#                              nu_1 = numeric(1), alpha0_2 = numeric(1), 
#                              alpha1_2 = numeric(1), alpha2_2 = numeric(1),
#                              beta_2 = numeric(1), nu_2 = numeric(1),
#                              P = numeric(1), P = numeric(1))

#in.sample.stats <- list(case = list())
# i = 1

#theta_hat <- data.frame()
# fit_data <- list()
# fit_data_bayes <- list()
# test <- list()
# 
# for(model in models){
#   for(dis in distr){
#     spec = create.spec(model = c(model, model), distribution = c(dis, dis))
#     fit = fit.mle(spec = spec, y = data$log.return)
#     fit_bayes = fit.bayes(spec = spec, y = data$log.return)
#     #theta_hat <- rbind(theta_hat, fit$theta)
#     spec[[i]] <- spec 
#     fit_data[[i]] <- fit
#     #fit_data_bayes[[i]] <- fit_bayes
#     post_cov = var(fit_bayes$theta)
#     
#     hess = hessian(function(theta){-MSGARCH::kernel(spec, theta, data$log.return)}, fit$theta)
#     hess = hessian( func = pdf_fun, x = fit$theta)
#     hess_pd = nearPD(hess)
#     cov_mat = solve(hess_pd$mat)
#     #test[[i]] = Wald.test(theta_hat = as.matrix(colMeans(fit_bayes$theta)), Vn = var(fit_bayes$theta), R = matrix(c(1, rep(0,11)), nrow = 1))
# 
#     #in.sample.stats$case[[i]] <- list(model = model, dist = dis, fit = fit, test = test)
#     i = i+1
#   }
# }
# 
# pdf_fun <- function(theta){
#   pdf = MSGARCH::pdf(object = spec, theta = theta, y = data$log.return, log = T, do.its = T)
#   -sum(pdf$pdf[which(!is.na(pdf$pdf))])
# }