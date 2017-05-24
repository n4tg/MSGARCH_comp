# Wrapper function --------------------------------------------------------

Backtesting <- function(alpha, r, VaR, lags = 4){
  hit = (r < VaR)
  
  tests = c("test.UC", "test.Ind", "test.CC", "test.DQ")
  res_table = data.table(Model = colnames(VaR),
                         Violation.percentage = 100*colMeans(hit))
  
  for(test in tests){
    if(test == "test.DQ"){
      res = apply(X = VaR, MARGIN = 2, FUN = test, alpha = alpha, r = r, lags = lags)
    } else res = apply(X = hit, MARGIN = 2, FUN = test, alpha = alpha)
    res_table = cbind(res_table, rbindlist(res, use.names = T, fill = T))
  }
  
  return(res_table)
}

# Conditional coverage test -----------------------------------------------

## Unconditional coverage test

test.UC <- function(alpha, hit){
  # alpha <- check.alpha(alpha)
  # hit <- check.hit(hit)
  
  n0 = length(which(hit == 0))
  n1 = length(which(hit == 1))
  alpha_hat = n1/(n0+n1)
  
  LR.UC = -2*log( ((alpha^n1)*(1-alpha)^n0)/((alpha_hat^n1)*(1-alpha_hat)^n0) )
  pval.UC = 1 - pchisq(LR.UC, df = 1)
  
  return(list(LR.UC = LR.UC, pval.UC = pval.UC))
}

## Independence test

test.Ind <- function(alpha, hit){
  # alpha <- check.alpha(alpha)
  # hit <- check.hit(hit)
  
  n00 <- n01 <- n10 <- n11 <- 0
  
  for(i in 1:(length(hit)-1)){
    if(hit[i] == 0){
      if(hit[i+1] == 0){
        n00 = n00+1
      } else{
        n01 = n01+1
      }
    } else{
      if(hit[i+1] == 0){
        n10 = n10+1
      } else{
        n11 = n11+1
      }
    }
  }
  
  p01 = n01/(n00+n01)
  p11 = n11/(n10+n11)
  p1 = (n01+n11)/(n00+n01+n10+n11)
  
  LR.Ind = -2*log( ( (1-p1)^(n00+n10) * p1^(n01+n11) )/( (1-p01)^n00 * p01^n01 * (1-p11)^n10 * p11^n11 ) )
  pval.Ind = 1 - pchisq(LR.Ind, df = 1)
  
  return(list(LR.Ind = LR.Ind, pval.Ind = pval.Ind))
}

## Conditional coverage test
test.CC <- function(alpha, hit){
  # alpha <- check.alpha(alpha)
  # hit <- check.hit(hit)
  
  LR.CC = test.UC(alpha, hit)$LR.UC + test.Ind(alpha, hit)$LR.Ind
  pval.CC = 1- pchisq(LR.CC, df = 2)
  
  return(list(LR.CC = LR.CC, pval.CC = pval.CC))
}

# Dynamic quantile test ---------------------------------------------------

test.DQ <- function(alpha, r, VaR, lags = 4){
  N.period = length(r)
  H = (r<VaR)-alpha
  Hit = H[(lags+1):N.period]
  VaR_hat = VaR[(lags+1):N.period]
  
  const = rep(1, times = N.period-lags)
  tmp_Mat = matrix(0, nrow = N.period-lags, ncol = lags)
  
  for(i in 1:lags){
    tmp_Mat[,i] = H[i:(N.period-lags+i-1)]
  }
  
  X = cbind(const, VaR_hat, tmp_Mat)
  
  DQ = (t(Hit) %*% X %*% solve(t(X) %*% X) %*% t(X) %*% Hit)/(alpha*(1-alpha))
  pval.DQ = 1 - pchisq(DQ, df = ncol(X))
  
  return(list(DQ = DQ, pval.DQ = pval.DQ))
}

# Superior predictive ability test ---------------------------------------
# 
# loss_fun <- function(alpha, r, VaR, delta = 25){
#   mean( ( alpha - 1/(1+exp( delta*(r-VaR) ) ) ) * (r-VaR) )
# }

# test.rc <- function()


# Wald test ---------------------------------------------------------------

# test.Wald <- function(theta_hat,
#                       Vn,
#                       R,
#                       r = 0){
#   W <- t(R %*% theta_hat - r) %*% solve(R %*% Vn %*% t(R)) %*% (R %*% theta_hat - r)
#   df <- nrow(R)
#   pval <- 1-pchisq(as.numeric(W), df)
#   return(data.frame(W = W, df = df, p = pval))
# }
