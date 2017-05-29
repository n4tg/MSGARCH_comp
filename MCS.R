install.packages("MCS")
library(MCS)

Loss <- apply(X = VaR, MARGIN = 2, FUN = function(x) MCS::LossVaR(realized = r, evaluated = x, tau = alpha, type = 'differentiable'))
SSM <- MCSprocedure(Loss = Loss, alpha = alpha, verbose = F)
sink("Output/DAX_MCS.txt")
print(SSM)
sink()


