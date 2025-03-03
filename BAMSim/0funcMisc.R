# =================================================================================
# Miscellaneous functions.
# Currently includes functions for geometric mean, power function, baranov catch eqn, and Lorenzen M
# KWS Feb 2025
#---------------------------------------------------------------------------------


#####################################################################################
geomean <- function(x) {
  ## compute geometric mean of a time series x
  nx <- length(x)
  y <- prod(x)^(1 / nx)
  return(y)
}

#####################################################################################
pow <- function(x, a, b) {
  ## INPUT PARAMETERS:
  ##    x = vector or scalar of independent variable
  ##    a = coefficient
  ##    b = exponent
  ## OUTPUT PARAMETERS:
  ##    y = values of fcn at x

  y <- a * x^b
  return(y)
}

#####################################################################################
baranov <- function(F.age, Z.age, N.age, wgt.age = NULL) {
  ## INPUT PARAMETERS:
  # F.age = vector of fishing mortality rate at age
  # Z.age = vector of total mortality rate at age including dead discards
  # N.age = vector of abundance at age
  # wgt.age = optional vector of weight at age
  ## OUTPUT PARAMETERS:
  # A list with total landings in numbers and optionally in weight
  # This will be "removals" with the addition of dead discards - MDD

  L.age.n <- L.age.w <- rep(NA, times = length(N.age))

  L.age.n <- F.age * N.age * (1.0 - exp(-Z.age)) / Z.age
  L.n <- sum(L.age.n)

  if (!is.null(wgt.age)) {
    L.age.w <- wgt.age * L.age.n
    L.w <- sum(L.age.w)
  } else {
    L.w <- NA
  }

  return(list(L.N = L.n, L.W = L.w))
} # end fcn baranov

#####################################################################################
M.lorenzen <- function(length, a = 0, b = -1) {
  ## FCN based on Lorenzen et al. 2022. Fisheries Research 106327
  ## INPUT PARAMETERS:
  # length = vector of length at age in cm
  # a = intercept parameter
  # b = slope parameter
  ## OUTPUT
  # M at age, which likely needs to be scaled for stock specific values
  Mage <- exp(a + b * log(length))
  return(Mage)
}
#####################################################################################
