# =================================================================================
# Spawner recruit functions.
# Options include Bev-Holt, Ricker, and Null models
# KWS Feb 2025
#---------------------------------------------------------------------------------

#####################################################################################
SR.fcn <- function(S, SR.switch, h = 0.75, R0 = 100, Phi.0 = 1) {
  # Beverton-Holt (B-H) Stock-Recruitment (S-R) Function
  ## INPUT PARAMETERS:
  ##    S = vector or scalar of independent variable (spawners)
  ##    SR.switch (1=BevHolt, 2=Ricker, 3=Null)
  ##    h = steepness, default is 0.75
  ##    R0 = recruitment of unfished popn or equilibrium recruitment for the Null model, default is 100
  ##    Phi.0 = spawners per recruit of unfished popn, default is 1
  ## OUTPUT PARAMETERS:
  ##    recruits = values of fcn at S

  recruits <- switch(SR.switch,
    ((0.8 * R0 * h * S) / (0.2 * R0 * Phi.0 * (1 - h) + S * (h - 0.2))),
    ((S / Phi.0) * exp(h * (1 - S / (R0 * Phi.0)))),
    (R0)
  )
  return(recruits)
}

#####################################################################################
SR.eq.fcn <- function(SR.switch, BC = 1, h = 0.75, R0 = 100, Phi.0 = 1, Phi.F = 1) {
  # Equilibrium spawning, conditional on F, with option for bias correction
  ## INPUT PARAMETERS:
  ##    S = vector or scalar of independent variable (spawners)
  ##    SR.switch (1=BevHolt, 2=Ricker, 3=Null)
  ##    BC = bias correction, default is none (=1)
  ##    h = steepness, default is 0.75
  ##    R0 = recruitment of unfished popn or equilibrium recruitment for the Null model, default is 100
  ##    Phi.0 = spawners per recruit of unfished popn, default is 1
  ##    Phi.F = spawners per recruit at fishing rate F
  ## OUTPUT PARAMETERS:
  ##    recruits = values of fcn at S

  recruits <- switch(SR.switch,
    ((R0 / ((5.0 * h - 1.0) * Phi.F)) * (BC * 4.0 * h * Phi.F - Phi.0 * (1.0 - h))),
    ((R0 / (Phi.F / Phi.0) * (1.0 + log(BC * Phi.F / Phi.0) / h))),
    (BC * R0)
  )
  return(recruits)
}
#####################################################################################
