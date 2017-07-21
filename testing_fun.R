# MCS ---------------------------------------------------------------------

## Function: test.MCS
# Description: A wrapper function for doing MCS.
# Input:  r - A numeric vector of the real data.
#         VaR - A numeric matrix of the forecasted data.
#         LossFn - A function indicates loss function used for MCS evaluation. 
#                  See MCS package manual for possible choices.
#         save.out - A logical input indicates whether the output MCS is saved.
#         filename - A string for path and filename to which the output is saved.
#         ... - Optional argument for the LossFn. See MCS package manual.
# Output: An SSM object. See MCS package manual.
test.MCS <- function(r, VaR, LossFn, 
                     save.out = T, filename = paste0("Output/MCS.txt"),
                     ...){
  r = as.numeric(r)

  Loss = apply(X = VaR, MARGIN = 2, FUN = function(x) LossFn(realized = r, evaluated = x, ...))
  SSM = MCSprocedure(Loss = Loss, alpha = alpha, verbose = F)
  if(save.out){
    sink(filename)
    print(SSM)
    sink()
  }
  return(SSM)
}

# Bactesting --------------------------------------------------------------

## Function: Backtesting
# Description: A wrapper function for backtesting.
# Input:  alpha - A number indicates shortfall probability.
#         r - A numeric vector of the real data.
#         VaR - A numeric matrix of the forecasted data.
#         lags - A number indicates lags of Hit in DQ test. See Engle and Mangenelli (2004).
#         save.out - A logical input indicates whether the backtesting output table is saved or not.
#         filename - A string for path and filename to which the output is saved.
# Output: A data table consists of model names, %violations, test statistics (UC, Ind, CC, DQ) and their p-values.
Backtesting <- function(alpha, r, VaR, lags = 4,
                        save.out = T, filename = paste0("Output/Backtest.csv")){
  alpha = check.alpha(alpha)
  r = as.numeric(r)

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
  
  if(save.out) write.table(res_table, file = filename, sep = ";", col.names = T, row.names = F)
  return(res_table)
}

## Function: test.UC
# Description: Perform unconditional coverage test of Kupiec (1995).
# Input:  alpha - A number indicates shortfall probability.
#         hit - A numeric vector of hit sequence.
# Output: A list consists of likelihood ratio of UC test (LR.UC) and its p-value (pval.UC).
test.UC <- function(alpha, hit){
  alpha = check.alpha(alpha)
  hit = check.hit(hit)
  
  n0 = length(which(hit == 0))
  n1 = length(which(hit == 1))
  alpha_hat = n1/(n0+n1)
  
  LR.UC = -2*log( ((alpha^n1)*(1-alpha)^n0)/((alpha_hat^n1)*(1-alpha_hat)^n0) )
  pval.UC = 1 - pchisq(LR.UC, df = 1)
  
  return(list(LR.UC = LR.UC, pval.UC = pval.UC))
}

## Function: test.Ind
# Description: Perform independent test as in Christoffersen (1998).
# Input:  alpha - A number indicates shortfall probability.
#         hit - A numeric vector of hit sequence.
# Output: A list consists of likelihood ratio of Ind test (LR.Ind) and its p-value (pval.Ind).
test.Ind <- function(alpha, hit){
  alpha = check.alpha(alpha)
  hit = check.hit(hit)
  
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


## Function: test.CC
# Description: Perform conditional coverage test of Christoffersen (1998).
# Input:  alpha - A number indicates shortfall probability.
#         hit - A numeric vector of hit sequence.
# Output: A list consists of likelihood ratio of CC test (LR.CC) and its p-value (pval.CC).
test.CC <- function(alpha, hit){
  alpha = check.alpha(alpha)
  hit = check.hit(hit)
  
  LR.CC = test.UC(alpha, hit)$LR.UC + test.Ind(alpha, hit)$LR.Ind
  pval.CC = 1- pchisq(LR.CC, df = 2)
  
  return(list(LR.CC = LR.CC, pval.CC = pval.CC))
}

## Function: test.DQ
# Description: Perform dynamic quantile test as in Engle and Mangenelli (2004).
# Input:  alpha - A number indicates shortfall probability.
#         r - A numeric vector of the real data.
#         VaR - A numeric matrix of the forecasted data.
#         lags - A number indicates lags of Hit. Default as 4.
# Output: A list consists of DQ test statistics (DQ) and its p-value (pval.DQ).
test.DQ <- function(alpha, r, VaR, lags = 4){
  alpha = check.alpha(alpha)
  r = as.numeric(r)

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
  mat = t(X) %*% X
  mat = check.Mat(mat)
  DQ = (t(Hit) %*% X %*% ginv(mat) %*% t(X) %*% Hit)/(alpha*(1-alpha))
  pval.DQ = 1 - pchisq(DQ, df = ncol(X))
  
  return(list(DQ = DQ, pval.DQ = pval.DQ))
}

# Wald test ---------------------------------------------------------------

## Function: test.Wald
# Description: Perform Wald test.
# Input:  theta_hat - A numeric vector of the estimated parameters.
#         Vn - A numeric matrix represents parameters covariance matrix.
#         R - A numeric vector or matrix of hypothesis coefficients.
#         r - A numeric vector or number in accordance to the hypothesis.
# Output: A list consists of Wald test statistics (W), its degree of freedom (df), and its p-value (pval).
test.Wald <- function(theta_hat,
                      Vn,
                      R,
                      r = 0){
  W = t(R %*% theta_hat - r) %*% ginv(R %*% Vn %*% t(R)) %*% (R %*% theta_hat - r)
  df = nrow(R)
  pval = 1-pchisq(as.numeric(W), df)
  return(list(W = W, df = df, p = pval))
}
