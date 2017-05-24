## SPA

loss_fun <- function(alpha, r, VaR, delta){
  (alpha - 1/(1+exp(delta*(r-VaR))))*(r-VaR)
}

SPA <- function(alpha, r, VaR, benchmark = 1, delta = 25, N.boot = 50000, q = 1){
  N.period = length(r)
  N.model = ncol(VaR)-1
  
  L = loss_fun(alpha, r, VaR, delta)
  d = L[,benchmark]-L
  d = d[,-benchmark]
  d_bar = colMeans(d, na.rm = T)
  #mat.d_bar = matrix(rep(d_bar, each = N.boot), nrow = N.boot)
  
  tmp.g = matrix(NA, N.period, N.model)
  for(k in 1:N.model){
    for(i in 0:(N.period-1)){
      tmp.g[i+1,k] <- 1/N.period * sum( (d[1:(N.period-i),k]-d_bar[k])*(d[(1+i):N.period,k]-d_bar[k]) )
    }
  }
  
  tmp.k = (N.period - 1:(N.period-1))/N.period * (1-q)^(1:(N.period-1)) + (1:(N.period-1))/N.period * (1-q)^(N.period-(1:(N.period-1)))
  var_hat = tmp.g[1,] + colSums(2*tmp.k*tmp.g[2:N.period,], na.rm = T)
  #mat.var_hat = matrix(rep(var_hat, each = N.boot), nrow = N.boot)
  
  # d_boot_bar = rep(0, N.model)
  # for(i in 1:N.period){
  #   d_boot_bar = d_boot_bar + colMeans(stationaryBoot(d, N.boot, q))
  # }
  # d_boot_bar = d_boot_bar/N.period
  
  d_boot_bar = colMeans(stationaryBoot(d, N.boot, q))
  
  df <- data.frame(d_boot_bar = d_boot_bar, d_bar, const = -sqrt(2*var_hat/N.period *log(log(N.period))), cond = (d_bar>=-sqrt(2*var_hat/N.period *log(log(N.period)))))
  
  SPA.c = sqrt(N.period/var_hat) * (d_boot_bar - (d_bar>=-sqrt(2*var_hat/N.period *log(log(N.period))))*d_bar)
  SPA.u = sqrt(N.period/var_hat) * (d_boot_bar - d_bar)
  SPA.l = sqrt(N.period/var_hat) * (d_boot_bar - pmax(0, d_bar))

  df.SPA.out <- data.frame(Model = colnames(d), SPA.c = pmax(0, SPA.c), SPA.l = pmax(0, SPA.l), SPA.u = pmax(0, SPA.u), stringsAsFactors = F)
  View(df.SPA.out)
  
  tmp.df = df
  tmp.SPA.c = sqrt(N.period/var_hat) * (tmp.df$d_boot_bar - (d_bar>=-sqrt(2*var_hat/N.period *log(log(N.period))))*d_bar)
  tmp.SPA.u = sqrt(N.period/var_hat) * (tmp.df$d_boot_bar - d_bar)
  tmp.SPA.l = sqrt(N.period/var_hat) * (tmp.df$d_boot_bar - pmax(0, d_bar))
  tmp.df.SPA.out <- data.frame(Model = colnames(d), SPA.c = pmax(0, tmp.SPA.c), SPA.l = pmax(0, tmp.SPA.l), SPA.u = pmax(0, tmp.SPA.u), stringsAsFactors = F)
  View(tmp.df.SPA.out)
  
  # t.SPA = matrix(rep(pmax(0, sqrt(N.period)*d_bar/var_hat), each = N.boot), nrow = N.boot)
  # pval.SPA.c = colMeans(pmax(0, SPA.c)>t.SPA)
  # pval.SPA.u = colMeans(pmax(0, SPA.u)>t.SPA)
  # pval.SPA.l = colMeans(pmax(0, SPA.l)>t.SPA)
  # 
  # df.out <- data.frame(Model = colnames(d), SPA.c = colMeans(SPA.c), pval.SPA.c = pval.SPA.c, SPA.l = colMeans(SPA.l),  
  #                      pval.SPA.l = pval.SPA.l, SPA.u = colMeans(SPA.u), pval.SPA.u = pval.SPA.u, stringsAsFactors = F)
  
  return(rbind(c(colnames(VaR[benchmark]), rep("benchmark", ncol(df.SPA.out))), df.SPA.out))
}

## Stationary bootstrap

stationaryBoot <- function(d, 
                           N.boot = 50000, #B
                           q = 1){
  N.model = ncol(d) #m
  N.period = nrow(d) #n
  
  U = replicate(N.period, runif(N.boot))
  V = replicate(N.period, runif(N.boot))
  
  while(length(which(U == 0))) U[which(U == 0)] = runif(1)
  while(length(which(V == 0))) V[which(V == 0)] = runif(1)
  
  tau = matrix(NA, N.boot, N.period)
  tau[,1] = ceiling(N.period*U[,1])
  
  for(t in 2:N.period){
    tau[,t] = (V[,t]<q)*ceiling(N.period*U[,t]) + (V[,t]>=q)*((tau[,t-1]<N.period)*tau[,t-1]+1)
  }
  
  d_boot = matrix(NA, N.boot, N.model)
  
  for(b in 1:N.boot){
    for(k in 1:N.model){
      d_boot[b,k] = mean(d[tau[b,],k])
    }
  }
  
  return(d_boot)
}
