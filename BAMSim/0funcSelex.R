# =================================================================================
# Selectivity functions.
# Current options include logistic and double logistic
# KWS, MDD Feb 2025
#---------------------------------------------------------------------------------

#####################################################################################
logistic <- function(x, a = 100, b = 0) {
  ## INPUT PARAMETERS:
  ##    x = vector or scalar of independent variable
  ##    a = slope parameter, default is 100 ("knife edge" for most applications)
  ##    b = location parameter: x giving 50% probability, default is 0
  ## OUTPUT PARAMETERS:
  ##    P = values of fcn at x
  P <- 1 / (1 + exp(-a * (x - b)))
  return(P)
}

#####################################################################################
double.logistic <- function(x, a1 = 100, b1 = 0, a2 = 100, b2 = 0) {
  ## INPUT PARAMETERS:
  ##    x = vector or scalar of independent variable
  ##    a = slope parameter, default is 100 ("knife edge" for most applications)
  ##    b = location parameter: x giving 50% probability, default is 0
  ## OUTPUT PARAMETERS:
  ##    P = values of fcn at x

  P.temp1 <- 1 / (1 + exp(-a1 * (x - b1)))
  P.temp2 <- 1 - (1 / (1 + exp(-a2 * (x - (b1 + b2)))))
  P <- P.temp1 * P.temp2
  P <- P / max(P)
  return(P)
}
#####################################################################################
