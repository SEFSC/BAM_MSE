#####################################################################################
# Population simulation code
# Written: 18 Sept 2019 by KWS
# Last updated: Feb 2025 by MDD and KWS


#####################################################################################
popsim <- function(x) {
  ## INPUT PARAMETERS:
  # x = a list containing the following
  # nyr = number of yrs to simulate,
  # F.init are initialization F's for each fleet, F = vector of F's each year,
  # F time series for each fleet
  # ages = vector of ages,
  # nages = number of ages modeled,
  # SR.switch = switch for SR function (1=BH, 2=Ricker, 3=Null)
  # SR.bc = bias correction in equilibrium SR
  # SPR.proxy = SPR proxy as a fraction, computed even if not needed
  # SPR.Req = equilibrium recruitment assumed for computing spawn output put at FX%
  # sp.frac = fraction of year when spawning occurs. Must be on the interval [0,1)
  # R0 = unfished equilibrium recruitment, h = steepnes, Phi.0 = spr unfished,
  # M.age = natural mortality ate age,
  # W.mt = weight at age in mt
  # reprod = vector of reproductive output at age, per individual
  # selex1,2,3 = landings selectivity at age of fleets 1,2,3
  # selexL.wgt, selexD.wgt = weighted selex of landings and discards for computing MSY quantities
  # logR.resid = recruitment residuals (lognormal)
  # len.prob = probability of length at age, only necessary if fitting length comps

  #####################################################################################
  # Dependencies: 0funcSR.R, 0funcMSY.R, 0funcSPR.R
  #####################################################################################

  nyr <- x$nyr
  F1.init <- x$F1.init
  F2.init <- x$F2.init
  F3.init <- x$F3.init
  F1 <- x$F1
  F2 <- x$F2
  F3 <- x$F3
  ages <- x$ages
  nages <- x$nages
  R0 <- x$R0
  h <- x$h
  Phi.0 <- x$Phi.0
  SR.switch <- x$SR.switch
  SR.bc <- x$SR.bc
  SPR.proxy <- x$SPR.proxy
  SPR.Req <- x$SPR.Req
  sp.frac <- x$sp.frac
  M.age <- x$M.age
  W.mt <- x$W.mt
  reprod <- x$reprod
  selex1 <- x$selex1
  selex2 <- x$selex2
  selex3 <- x$selex3
  selexL.wgt <- x$selexL.wgt
  selexD.wgt <- x$selexD.wgt
  logR.resid <- x$logR.resid
  len.prob <- x$len.prob

  Ltot.mt <- Ltot.knum <- Dtot.mt <- Dtot.knum <- L1.mt <- L1.knum <- L2.mt <- L2.knum <- D.mt <- D.knum <- SSB <- biomass.mt <- abundance <- rep(0, nyr) # quantities of interest (annual)
  F.age <- N.age <- L1.age <- L2.age <- D.age <- matrix(0, nrow = nyr, ncol = nages)
  nlenbins <- ncol(len.prob)
  D.len <- matrix(0, nrow = nyr, ncol = nlenbins)
  D.dum <- matrix(0, nrow = nages, ncol = nlenbins)

  # Initial conditions assumes equilibrium age structure given initial F
  N.pr1 <- N.pr1.sp <- rep(1, nages) # Number of spawners per recruit at age at start of yr and at peak spawn
  F.age[1, ] <- F1.init * selex1 + F2.init * selex2 + F3.init * selex3
  Z <- F.age[1, ] + M.age
  for (a in 1:(nages - 1))
  {
    N.pr1[a + 1] <- N.pr1[a] * exp(-Z[a])
  }
  N.pr1[nages] <- N.pr1[nages] / (1 - exp(-Z[nages])) # Plus group
  N.pr1.sp[1:(nages - 1)] <- N.pr1[1:(nages - 1)] * exp(-sp.frac * Z[1:(nages - 1)])
  N.pr1.sp[nages] <- N.pr1.sp[(nages - 1)] * (exp(-(1 - sp.frac) * Z[(nages - 1)]) + sp.frac * Z[nages]) / (1 - exp(-Z[nages]))

  Phi.F <- sum(N.pr1.sp * reprod) # Spawners per recruit based on mature female biomass

  R.eq <- SR.eq.fcn(SR.switch = SR.switch, BC = SR.bc, h = h, R0 = R0, Phi.0 = Phi.0, Phi.F = Phi.F)
  if (R.eq < 1) {
    R.eq <- 1
  } # Catch numerical possibility that equilibrium R is negative
  N.age[1, ] <- R.eq * N.pr1.sp
  N.age[1, 1] <- N.age[1, 1] * exp(logR.resid[1])

  for (i in 1:(nyr - 1)) {
    F.age[i, ] <- F1[i] * selex1 + F2[i] * selex2 + F3[i] * selex3
    Z <- F.age[i, ] + M.age

    SSB[i] <- sum(N.age[i, ] * reprod)
    N.age[(i + 1), 1] <- SR.fcn(S = SSB[i], SR.switch = SR.switch, h = h, R0 = R0, Phi.0 = Phi.0) * exp(logR.resid[i + 1])
    for (a in 1:(nages - 1))
    {
      N.age[(i + 1), (a + 1)] <- N.age[i, a] * exp(-Z[a])
    } # Abundance at age in each year
    N.age[(i + 1), nages] <- N.age[(i + 1), nages] + N.age[i, nages] * exp(-Z[nages]) # Plus group correction

    L1.age[i, ] <- F1[i] * selex1 / (Z) * N.age[i, ] * (1 - exp(-Z))
    L1.knum[i] <- sum(L1.age[i, ]) / 1000
    L1.mt[i] <- sum(L1.age[i, ] * W.mt)

    L2.age[i, ] <- F2[i] * selex2 / (Z) * N.age[i, ] * (1 - exp(-Z))
    L2.knum[i] <- sum(L2.age[i, ]) / 1000
    L2.mt[i] <- sum(L2.age[i, ] * W.mt)

    D.age[i, ] <- F3[i] * selex3 / (Z) * N.age[i, ] * (1 - exp(-Z))
    D.knum[i] <- sum(D.age[i, ]) / 1000
    D.mt[i] <- sum(D.age[i, ] * W.mt)
    for (a in 1:nages)
    {
      D.dum[a, ] <- D.age[i, a] * len.prob[a, ]
    }
    D.len[i, ] <- colSums(D.dum) / sum(D.age[i, ])

    Ltot.knum[i] <- L1.knum[i] + L2.knum[i] # could be an issue not including discards here
    Ltot.mt[i] <- L1.mt[i] + L2.mt[i]
    Dtot.knum[i] <- D.knum[i]
    Dtot.mt[i] <- D.mt[i]
    abundance[i] <- sum(N.age[i, ])
    biomass.mt[i] <- sum(N.age[i, ] * W.mt)
  }
  F.age[nyr, ] <- F1[nyr] * selex1 + F2[nyr] * selex2 + F3[nyr] * selex3
  Z <- F.age[nyr, ] + M.age
  SSB[nyr] <- sum(N.age[nyr, ] * reprod)
  L1.age[nyr, ] <- F1[nyr] * selex1 / (Z) * N.age[nyr, ] * (1 - exp(-Z))
  L1.knum[nyr] <- sum(L1.age[nyr, ]) / 1000
  L1.mt[nyr] <- sum(L1.age[nyr, ] * W.mt)
  L2.age[nyr, ] <- F2[nyr] * selex2 / (Z) * N.age[nyr, ] * (1 - exp(-Z))
  L2.knum[nyr] <- sum(L2.age[nyr, ]) / 1000
  L2.mt[nyr] <- sum(L2.age[nyr, ] * W.mt)
  D.age[nyr, ] <- F3[nyr] * selex3 / (Z) * N.age[nyr, ] * (1 - exp(-Z))
  for (a in 1:nages)
  {
    D.dum[a, ] <- D.age[nyr, a] * len.prob[a, ]
  }
  D.len[nyr, ] <- colSums(D.dum) / sum(D.age[nyr, ])
  D.knum[nyr] <- sum(D.age[nyr, ]) / 1000
  D.mt[nyr] <- sum(D.age[nyr, ] * W.mt)
  Ltot.knum[nyr] <- L1.knum[nyr] + L2.knum[nyr]
  Ltot.mt[nyr] <- L1.mt[nyr] + L2.mt[nyr]
  Dtot.knum[nyr] <- D.knum[nyr]
  Dtot.mt[nyr] <- D.mt[nyr]
  abundance[nyr] <- sum(N.age[nyr, ])
  biomass.mt[nyr] <- sum(N.age[nyr, ] * W.mt)

  # Compute reference points
  selex.Z <- selexL.wgt + selexD.wgt
  msy <- MSY.func(
    steep = h, R0 = R0, M = M.age, wgt = W.mt, sp.frac = sp.frac, reprod.age = reprod, selL = selexL.wgt, selD = selexD.wgt, selZ = selex.Z,
    SR.switch = SR.switch, SR.bc = SR.bc, maxF = 1.0, step = 0.001, verbose = FALSE
  )

  spr <- SPR.func(
    nages = nages, max.F = 1.0, R.eq = SPR.Req, SPR.input = SPR.proxy, sp.frac = sp.frac,
    reprod.age = reprod, M.age = M.age, selex.age = selex.Z
  )


  return(list(
    yr = 1:nyr, SSB = SSB, abundance = abundance, biomass.mt = biomass.mt,
    N.age = N.age, F.age = F.age,
    L1.age = L1.age, L2.age = L2.age, D.age = D.age, D.len = D.len,
    F.fleet1 = F1, F.fleet2 = F2, F.fleet3 = F3,
    L1.knum = L1.knum, L1.mt = L1.mt, L2.knum = L2.knum, L2.mt = L2.mt,
    D.knum = D.knum, D.mt = D.mt, msy = msy, spr = spr, selex.Z = selex.Z, reprod.age = reprod
  ))
} # end popsim

#####################################################################################
