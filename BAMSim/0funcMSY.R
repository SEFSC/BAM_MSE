
# =================================================================================
# Compute Maximum Sustainable Yield (MSY) reference points
# Output based on Bev-Holt, Ricker, or Null recruitment models
# KWS March 2006
# Last update: Feb 2025
#---------------------------------------------------------------------------------

MSY.func <- function(steep, R0, M, wgt, sp.frac = 0.0, reprod.age, selL, selD, selZ,
                     SR.switch = 1, SR.bc = 1, maxF = 2.0, step = 0.001, verbose = TRUE) {

  ## INPUT:
  ##    steep = steepness parameter
  ##    R0    = virgin recruitment parameter
  ##    M     = natural mortality, may be constant or vector
  ##    wgt   = vector of weight at age, supply vector of 1's if values in numbers desired
  ##    sp.frac = fraction of year when peak spawning occurs, default is 0.0 (i.e., Jan 1)
  ##    reprod.age = vector of reproductive output at age
  ##    selL  = selectivity at age for landings
  ##    selD  = selectivity at age for dead discards; set to vector of 0's if no discards
  ##    selZ  = selectivity at age for dead fish
  ##    SR.switch = switch for SR function (1=BH (default), 2=Ricker, 3=Null)
  ##    OPTIONAL INPUT
  ##    SR.bc = lognormal bias correction -- e.g., exp(sigma^2/2)
  ##    maxF  = maximum F examined, default is 2.0
  ##    step  = accuracy in MSY calculations (default is 0.001)
  ##
  ## OUTPUT:
  ##    MSY   = maximum sustainable yield in weight, units same as wgt input
  ##    Lmsy_num = landings in numbers at MSY
  ##    Fmsy  = fishing rate at MSY
  ##    Dmsy_wgt  = dead discards at MSY, units same as wgt input
  ##    Dmsy_num  = dead discards in numbers at MSY
  ##    spr_msy= spawners per recruit at MSY
  ##    SPRmsy= spawning potential ratio at MSY (spr_msy/spr_virgin)
  ##    SSBmsy= spawning output at MSY, in units of fecundity input
  ##    Rmsy  = equilibrium recruitment at MSY
  ##    Bmsy  = total biomass (male and female) at MSY
  ##    Emsy  = exploitation rate at MSY (total kills in number / abundance of fish in number)
  ##    F = vector of F's used to compute MSY
  ##    Vectors of landings and discards in numbers and weight and vector of SSB as functions of F

  #####################################################################################
  # Dependencies: 0funcSR.R
  #####################################################################################

  #  if (! is.numeric(amin)) stop ("Non-numeric value for minimum age!")
  #  if (! is.numeric(amax)) stop ("Non-numeric value for maximum age!")
  if (!is.numeric(steep)) stop("Non-numeric value for steepnes!")
  if (!is.numeric(R0)) stop("Non-numeric value for R0!")
  if (!is.numeric(M)) stop("Non-numeric value for M!")
  if (!is.numeric(wgt)) stop("Non-numeric value for weight at age!")
  if (!is.numeric(selL)) stop("Non-numeric value for L selectivity at age!")
  if (!is.numeric(selD)) stop("Non-numeric value for D selectivity at age!")
  if (!is.numeric(selZ)) stop("Non-numeric value for Z selectivity at age!")

  if (verbose) {
    if (BC == 1) {
      cat("*** MSY NOTE: Estimates contain no bias correction.\n")
    } else {
      cat("*** NOTE: Estimates contain bias correction (mean unbiased).\n")
    }
  }

  ## INITIALIZATION
  if (length(M) > 1) {
    M_age <- M
  } # natural mortality at age (may be constant)
  else {
    M_age <- rep(M, nages)
  }

  F <- seq(0.0, maxF, by = step)
  spr <- rep(0.0, length(F)) # equilibrium spr at F
  S_eq <- rep(0.0, length(F)) # equilibrium reproductive output (e.g., SSB) at F
  R_eq <- rep(0.0, length(F)) # equilibrium recruitment at F
  B_eq <- rep(0.0, length(F)) # equilibrium biomass at F
  L_num_eq <- rep(0.0, length(F)) # equilibrium landings at F
  D_num_eq <- rep(0.0, length(F)) # equilibrium dead discards at F
  L_wgt_eq <- rep(0.0, length(F)) # equilibrium landings at F
  D_wgt_eq <- rep(0.0, length(F)) # equilibrium dead discards at F
  E_eq <- rep(0.0, length(F)) # equilibrium exploitation rate at F (landings only)

  L_age <- rep(0.0, nages) # landings at age
  D_age <- rep(0.0, nages) # dead discards at age
  F_age <- rep(0.0, nages) # F at age
  Z_age <- rep(0.0, nages) # Z at age

  ## Compute virgin spr
  N0_sp <- rep(1.0, times = nages)
  N0_sp[1] <- N0_sp[1] * exp(-M_age[1] * sp.frac)
  for (iage in 2:nages) {
    N0_sp[iage] <- N0_sp[iage - 1] * exp(-1.0 * ((M_age[iage - 1] * (1.0 - sp.frac)) + (M_age[iage] * sp.frac)))
  }
  N0_sp[nages] <- N0_sp[nages] / (1. - exp(-1.0 * M_age[nages]))

  spr_F0 <- sum(N0_sp * reprod.age)

  for (i in 1:length(F)) {
    FL_age <- F[i] * selL
    FD_age <- F[i] * selD
    Z_age <- M_age + F[i] * selZ

    N_age <- N_age_sp <- N_age_msy <- rep(1.0, nages) # N at age
    for (iage in 2:nages) {
      N_age[iage] <- N_age[iage - 1] * exp(-1.0 * Z_age[iage - 1])
    }
    # last age is pooled
    N_age[nages] <- N_age[nages - 1] * exp(-1. * Z_age[nages - 1]) /
      (1. - exp(-1.0 * Z_age[nages]))
    N_age_sp[1:(nages - 1)] <- N_age[1:(nages - 1)] * exp(-sp.frac * Z_age[1:(nages - 1)])
    N_age_sp[nages] <- N_age_sp[(nages - 1)] * (exp(-(1 - sp.frac) * Z_age[(nages - 1)]) + sp.frac * Z_age[nages]) / (1 - exp(-Z_age[nages]))

    spr[i] <- sum(N_age_sp * reprod.age)

    R_eq[i] <- SR.eq.fcn(
      SR.switch = SR.switch, BC = SR.bc, h = steep, R0 = R0,
      Phi.0 = spr_F0, Phi.F = spr[i]
    )
    if (R_eq[i] < 0.0000001) R_eq[i] <- 0.0000001

    N_age_msy <- R_eq[i] * N_age
    N_age_msy_sp <- R_eq[i] * N_age_sp

    for (iage in 1:nages) {
      L_age[iage] <- N_age_msy[iage] *
        (FL_age[iage] / Z_age[iage]) * (1. - exp(-1.0 * Z_age[iage]))
      D_age[iage] <- N_age_msy[iage] *
        (FD_age[iage] / Z_age[iage]) * (1. - exp(-1.0 * Z_age[iage]))
    }
    S_eq[i] <- sum(N_age_msy_sp * reprod.age)
    B_eq[i] <- sum(N_age_msy * wgt)
    L_num_eq[i] <- sum(L_age)
    D_num_eq[i] <- sum(D_age)
    L_wgt_eq[i] <- sum(L_age * wgt)
    D_wgt_eq[i] <- sum(D_age * wgt)
    E_eq[i] <- (sum(L_age) + sum(D_age)) / sum(N_age_msy)
  } # END F loop

  msy_out <- max(L_wgt_eq)
  F_msy_out <- F[L_wgt_eq == msy_out]
  spr_msy_out <- spr[L_wgt_eq == msy_out]
  SR_msy_out <- spr_msy_out / spr_F0
  Lnum_msy_out <- L_num_eq[L_wgt_eq == msy_out]
  Dnum_msy_out <- D_num_eq[L_wgt_eq == msy_out]
  Dwgt_msy_out <- D_wgt_eq[L_wgt_eq == msy_out]
  R_msy_out <- R_eq[L_wgt_eq == msy_out]
  S_msy_out <- S_eq[L_wgt_eq == msy_out]
  B_msy_out <- B_eq[L_wgt_eq == msy_out]
  E_msy_out <- E_eq[L_wgt_eq == msy_out]

  if (F_msy_out == maxF) {
    cat("*** Fmsy reached a bound.\n")
  }

  return(list(
    msy = msy_out, Lmsy_num = Lnum_msy_out, Fmsy = F_msy_out, Dmsy_wgt = Dwgt_msy_out, Dmsy_num = Dnum_msy_out,
    spr_msy = spr_msy_out, SPRmsy = SR_msy_out, SSBmsy = S_msy_out, Rmsy = R_msy_out,
    Bmsy = B_msy_out, Emsy = E_msy_out,
    F = F, L_wgt_eq = L_wgt_eq, L_num_eq = L_num_eq, D_wgt_eq = D_wgt_eq, D_num_eq = D_num_eq, SSB_eq = S_eq
  ))
}
