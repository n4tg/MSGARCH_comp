backtest_res_pack = matrix(NA, nrow = 18, ncol = 7)
for(i in 1:18){
  backtest_tmp = GAS::BacktestVaR(r, VaR[,i], alpha)
  backtest_res_pack[i,] = c(colnames(VaR[i]), backtest_tmp$LRuc[1], backtest_tmp$LRuc[2], 
                            backtest_tmp$LRcc[1], backtest_tmp$LRcc[2], backtest_tmp$DQ$stat, backtest_tmp$DQ$pvalue)
}
colnames(backtest_res_pack) <- c("Model", "LR.UC", "pval.UC", "LR.CC", "pval.CC", "DQ", "pval.DQ")
