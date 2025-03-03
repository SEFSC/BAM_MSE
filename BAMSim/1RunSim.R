
#####################################################################################
# Run the simulation code, write data input files
# Written: 11 Oct 2019 by KWS
# Last updated: 10 Feb 2025 by MDD and KWS

# Simulation details: three fleets (fleets 1&2 are landings, fleet3 is discards),
#                     single survey, logistic selectivity for fleets 1 &2,
#                     double logistic for fleet 3
#                     some details based on red snapper (e.g. fecundity), but easily modified
#####################################################################################
# Clear any existing code
graphics.off()
rm(list = ls(all = TRUE))

source("1PopSim.R") # Population simulator
source("1ObservationModel.R") # Add observation error to the simulation output (operating model)
source("1BAM.write.dat.R") # Create dat file for BAM
source("0funcMSY.R") # Compute MSY-related benchmarks
source("0funcSPR.R") # Compute SPR proxies
source("0funcLenProb.R") # Compute probability of length at age (only required if fitting length comps)
source("0funcSR.R") # Spawner recruit functions (Available options are Bev-Holt, Ricker, Null)
source("0funcSelex.R") # Selectivity functions (Available options are logistic and double-logistic)
source("0funcMisc.R") # Miscellaneous functions (geometric mean, power fcn, baranov catch eqn, Lorenzen M)


# Compile BAM?
compile <- as.logical(F)
# admb file to run, without tpl extension
filename <- "BAM-Sim"
# switches applied during admb execution, list with space in between
admb.switch <- "-nox"

# set the RN seed for repeatability if generating stochastic results
# set.seed(4886)

#####################################################################################
# Define start years of assessment and data sources
# Landings/discards assumed to begin in styr of the assessment
styr <- 1 # first year of the assessment
styr.Rdevs <- 1 # first year in the assessment to estimate rec devs
styr.survey <- 21 # year to start survey (for index and age comps)
styr.ac1 <- 27 # year to start age comps of the fleet1 landings
styr.ac2 <- 31 # year to start age comps of the fleet2 landings
styr.d <- 61 # year to start length comps for discards (fleet 3)

if (any(c(styr.Rdevs, styr.survey, styr.ac1, styr.ac2, styr.d) < styr)) {
  stop("Error: All data styrs must be greater than or equal to assessment styr.")
}
#####################################################################################
# Define number of years and time series of mean F values for popn simulation
nyr <- 70

# Define a mean pattern of general recreational F
F.breakyr <- 50
F1.ts1 <- logistic(1:F.breakyr, 0.2, 20) * 0.1 # similar pattern to Siegfried et al. but without the descent
F1.ts2 <- 0.1 - (logistic((F.breakyr + 1):nyr, 0.2, (F.breakyr + 10)) * 0.1)
F1.ts <- c(F1.ts1, F1.ts2)
# plot(F1.ts) # check
F1.mu <- F1.ts

# Define a mean pattern of general commercial F
F2.ts1 <- logistic(1:F.breakyr, 0.1, 20) * 0.05
F2.ts2 <- 0.05 - (logistic((F.breakyr + 1):nyr, 0.1, (F.breakyr + 10)) * 0.05)
F2.ts <- c(F2.ts1, F2.ts2)
F2.mu <- F2.ts

# Define a mean pattern of discard F
F3.ts <- logistic(1:nyr, 0.2, 50) * 0.1
F3.mu <- F3.ts

Ftotal <- F1.mu + F2.mu + F3.mu

F1.prop <- F1.mu / Ftotal
F2.prop <- F2.mu / Ftotal
F3.prop <- F3.mu / Ftotal


sd.F1 <- 0.2 # process error: SD of F1 in log space
sd.F2 <- 0.2 # process error: SD of F2 in log space
sd.F3 <- 0.2 # process error: SD of F3 in log space
F1.devs <- rnorm(nyr, mean = 0, sd = sd.F1)
F1.devs <- F1.devs - mean(F1.devs)
F2.devs <- rnorm(nyr, mean = 0, sd = sd.F2)
F2.devs <- F2.devs - mean(F2.devs)
F3.devs <- rnorm(nyr, mean = 0, sd = sd.F3)
F3.devs <- F3.devs - mean(F3.devs)
F1 <- F1.mu * exp(F1.devs) # sequence of annual F values for fleet one (Rec)
F2 <- F2.mu * exp(F2.devs) # sequence of annual F values for fleet two (Comm)
F3 <- F3.mu * exp(F3.devs) # dead discard fleet
## NOTE, the F devs are different from those in BAM. These surround time-varying true vals, BAM surrounds an overall mean

F1.init <- F2.init <- F3.init <- 0.0 # initialize the popn at start of yr1 at unfished condition


# Noise levels in observation model and recruitment
cv.L1 <- 0.2 # CV of landings from fleet 1
cv.L2 <- 0.2 # CV of landings from fleet 2
cv.D <- 0.2 # CV of dead discards (fleet 3)
cv.survey <- 0.3 # CV of survey
n.ac.L1 <- 100 # annual sample size (nfish) of recreational age comps
n.ac.L2 <- 100 # annual sample size of commercial age comps
n.ac.survey <- 100 # annual sample size of survey age comps
n.lc.D <- 100 # sample size of length comps of discard fleet (fleet 3)

# CVs used for fitting, need not be same as CVs in observation error
cv.L1.fit <- 0.05
cv.L2.fit <- 0.05
cv.D.fit <- 0.05
cv.survey.fit <- 0.1

# Ages and length comp setup
nages <- 20 # Max Age, Number ages modeled
ages <- 1:nages # Ages modeled in the popn

binw <- 30 # bin width: 30mm
lenbins <- seq(210, 990, binw) # must be in same units as length variable len (defined below)
nlenbins <- length(lenbins)

#####################################################################################
# Define selectivity
A50.sel1 <- 2.7 # Age at 50% selection for landings fleet 1
slope.sel1 <- 5.0 # Slope of selectivity for landings fleet 1
A50.sel2 <- 3.3 # Age at 50% selectivity for landings fleet 2
slope.sel2 <- 2.9 # Slope of selectivity for landings fleet 2
A50a.sel3 <- 2.7 # Age at 50% selectivity for dead discards (fleet 3) ascending limb
slopea.sel3 <- 5.0 # Slope of selectivity for dead discards (fleet 3) ascending limb
A50d.sel3 <- 5.0 # A50a + A50d = age at 50% selectivity for dead discards (fleet 3) descending limb
sloped.sel3 <- 1.0 # Slope of selectivity for dead discards (fleet 3) descending limb
A50.sel.survey <- 1.95 # Age at 50% selectivity for survey
slope.sel.survey <- 2.5 # Slope of selectivity for survey

selex1 <- logistic(ages, slope.sel1, A50.sel1) # selectivity of fleet 1
selex2 <- logistic(ages, slope.sel2, A50.sel2) # selectivity of fleet 2
selex3 <- double.logistic(ages, slopea.sel3, A50a.sel3, sloped.sel3, A50d.sel3) # selectivity of fleet 3
selex.survey <- logistic(ages, slope.sel.survey, A50.sel.survey)

F1.mu.end <- geomean(F1[(nyr - 2):nyr])
F2.mu.end <- geomean(F2[(nyr - 2):nyr])
F3.mu.end <- geomean(F3[(nyr - 2):nyr])
tot.mu.end <- F1.mu.end + F2.mu.end + F3.mu.end
wgt.1 <- F1.mu.end / tot.mu.end
wgt.2 <- F2.mu.end / tot.mu.end
wgt.3 <- F3.mu.end / tot.mu.end
selexL.wgt <- (wgt.1 * selex1 + wgt.2 * selex2) / max((wgt.1 * selex1 + wgt.2 * selex2 + wgt.3 * selex3)) # wgted landings selex, used for MSY
selexD.wgt <- (wgt.3 * selex3) / max((wgt.1 * selex1 + wgt.2 * selex2 + wgt.3 * selex3)) # wgted discard selex, used for MSY calcs

#####################################################################################
# Life-history parms
SR.switch <- 3 # SR function switch (1=BH, 2=Ricker, 3=Null)
SR.bc.switch <- FALSE # TRUE or FALSE to turn on/off bias correction.
spawn.frac <- 0.0 # Fraction of year when spawning occurs. Must be on the interval [0,1)
R0 <- 500000 # Average annual unfished recruitment (scales the popn)
h <- 0.7 # Steepness of the spawner-recruit relationship; not used in Null model
logR.sd <- 0.3 # Standard deviation of log recruitment

if (SR.bc.switch) {
  SR.bc <- exp(logR.sd^2 / 2.0)
} else {
  SR.bc <- 1.0
}
SPR.proxy <- 0.4 # SPR proxy as a fraction. MSY values are also computed, conditional on BH SR fcn.
SPR.Req <- R0 # recruitment level used to compute spawning output at FX% (set ot 0 if not needed)

# Multiplier on age-dependent natural mortality; should be in the range (0,2)
M.mult <- 1
# NOTE: M.mult is specific to a study on M estimation, may not be needed in general.


Linf <- 800 # Asymptotic average length
K <- 0.15 # Growth coefficient
a0 <- -0.33 # Theoretical age at size 0
cv.laa <- 0.15 # CV of length at age; only needed if fitting length comps
a.lw <- 0.000000025 # Length-weight coefficient
b.lw <- 3.0 # Length-weight exponent
A50.mat <- 2.25 # Age at 50% maturity
slope.mat <- 3 # Slope of maturity ogive

# fecundity parameters; NOTE these values are from SEDAR 73
fec.a <- -271137.0
fec.b <- 0.0235
fec.c <- 2.7404
fec.thresh <- 400.0 # mm
fec.min <- 55523.0 # eggs at minimum threshold

#####################################################################################
# Create and fill vectors to be used in the population model
len <- Linf * (1 - exp(-K * (ages - a0))) # von Bertalanffy growth
len.cm <- len / 10
W.kg <- a.lw * len^b.lw # Weight-length relationship, assumed to be in kg
W.mt <- W.kg / 1000 # Weight in mt
# age-specific M from SEDAR73
Mlorenzen <- M.lorenzen(length = len.cm, a = 0.0, b = -1.0) # parameter values from Lorenzen et al. 2022 paper

M.age.unscaled <- 0.287 * Mlorenzen / max(Mlorenzen) # 0.287 scales similarly to red snapper

M.age <- M.mult * M.age.unscaled # natural mortality at age
mat.age <- logistic(ages, slope.mat, A50.mat) # maturity at age
proportion.female <- rep(0.5, nages) # proportion female at age
# fecundity
fecundity <- rep(NA, length = nages)
for (i in 1:nages) {
  if (len[i] < fec.thresh) {
    fecundity[i] <- fec.min
  } else {
    fecundity[i] <- fec.a + fec.b * len[i]^fec.c # sets fecundity for any fish <400mm to min. fecundity
  }
}

#####################################################################################
# Compute age-based vector of spawning output, which multiplies with N.age to compute total spawning output
reprod <- proportion.female * mat.age * fecundity

#####################################################################################
# Compute the number of spawners per recruit of an unfished population (Phi.0)
N0.sp <- rep(1.0, times = nages)
N0.sp[1] <- N0.sp[1] * exp(-M.age[1] * spawn.frac)
for (iage in 2:nages) {
  N0.sp[iage] <- N0.sp[iage - 1] * exp(-1.0 * ((M.age[iage - 1] * (1.0 - spawn.frac)) + (M.age[iage] * spawn.frac)))
}
N0.sp[nages] <- N0.sp[nages] / (1. - exp(-1.0 * M.age[nages]))

Phi.0 <- sum(N0.sp * reprod)

logR.resid <- rnorm(nyr, mean = 0, sd = logR.sd) # generate recruitment residuals to be passed to sim module
logR.resid <- logR.resid - mean(logR.resid) # center on zero, because that's a constraint in BAM (admb dev vectors)
# brings residuals downward to center on zero - this is the bias correction to lognormal error

#####################################################################################
# Compute probability of length at age
lenprob <- lenprob.fcn(len.mu = len, len.cv = cv.laa, len.bins = lenbins, len.binw = binw)

# Simulate the stock dynamics
input1 <- list(
  nyr = nyr, F1.init = F1.init, F2.init = F2.init, F3.init = F3.init, F1 = F1, F2 = F2, F3 = F3,
  ages = ages, nages = nages, SR.switch = SR.switch, SR.bc = SR.bc, SPR.proxy = SPR.proxy, SPR.Req = SPR.Req,
  sp.frac = spawn.frac, R0 = R0, h = h, Phi.0 = Phi.0, M.age = M.age, W.mt = W.mt,
  reprod = reprod,
  selex1 = selex1, selex2 = selex2, selex3 = selex3, selexL.wgt = selexL.wgt,
  selexD.wgt = selexD.wgt, logR.resid = logR.resid,
  len.prob = lenprob
)
sim <- popsim(x = input1)
true.vals <- list(
  nyr = nyr, Ftotal = Ftotal, styr.survey = styr.survey,
  sd.F1 = sd.F1, sd.F2 = sd.F2, sd.F3 = sd.F3, cv.L1 = cv.L1, cv.L2 = cv.L2, cv.D = cv.D, cv.survey = cv.survey,
  n.ac.L1 = n.ac.L1, n.ac.L2 = n.ac.L2, n.lc.D = n.lc.D, n.ac.survey = n.ac.survey, cv.L1.fit = cv.L1.fit, cv.L2.fit = cv.L2.fit, cv.D.fit = cv.D.fit,
  cv.survey.fit = cv.survey.fit, logR.sd = logR.sd, SR.switch = SR.switch, sp.frac = spawn.frac, R0 = R0, h = h, M.mult = M.mult, Linf = Linf, K = K, a0 = a0, len.cv = cv.laa, a.lw = a.lw, b.lw = b.lw,
  A50.mat = A50.mat, slope.mat = slope.mat, fec.a = fec.a, fec.b = fec.b, fec.c = fec.c, fec.thresh = fec.thresh, fec.min = fec.min,
  A50.sel1 = A50.sel1, slope.sel1 = slope.sel1, A50.sel2 = A50.sel2, slope.sel2 = slope.sel2, A50a.sel3 = A50a.sel3,
  slopea.sel3 = slopea.sel3, A50d.sel3 = A50d.sel3, sloped.sel3 = sloped.sel3,
  A50.sel.survey = A50.sel.survey, slope.sel.survey = slope.sel.survey,
  selex1 = selex1, selex2 = selex2, selex3 = selex3, selex.survey = selex.survey, selexL.wgt = selexL.wgt, selexD.wgt = selexD.wgt,
  M.age = M.age, reprod = reprod, F1.devs = F1.devs, F2.devs = F2.devs, F3.devs = F3.devs, F1 = F1, F2 = F2, F3 = F3, logR.resid = logR.resid
)
sim$true.vals <- true.vals
dput(sim, file = "simvals.rdat")

#####################################################################################
# Observation model; add noise to "true" data
survey.sim.age <- sim$N.age[styr.survey:nyr, ] %*% diag(selex.survey)
survey.sim.raw <- rowSums(survey.sim.age)
survey.sim <- survey.sim.raw / mean(survey.sim.raw)
sim.q <- 1 / mean(survey.sim.raw)
dat.sim <- ObsModel(
  L1 = sim$L1.mt, L2 = sim$L2.mt, D = sim$D.mt, survey = survey.sim,
  L1.age = sim$L1.age, L2.age = sim$L2.age, D.len = sim$D.len,
  survey.age = survey.sim.age,
  cv.L1 = cv.L1, cv.L2 = cv.L2, cv.D = cv.D, cv.survey = cv.survey,
  n.L1 = n.ac.L1, n.L2 = n.ac.L2, n.D = n.lc.D, n.survey = n.ac.survey,
  a1.yr1 = styr.ac1, a2.yr1 = styr.ac2, d.yr1 = styr.d
)

#####################################################################################
# Write dat files for BAM
survey.yrs <- styr.survey:nyr
survey.nyr <- length(survey.yrs)
survey.cv.yr <- rep(cv.survey.fit, survey.nyr)
nsamp.survey <- nfish.survey <- rep(n.ac.survey, survey.nyr)
dat.survey <- list(
  nyr = survey.nyr, yrs = survey.yrs,
  vals.obs = dat.sim$survey.obs, cv = survey.cv.yr,
  nyr.ages = length(survey.yrs), yrs.age = survey.yrs,
  nsamp = nsamp.survey, nfish = nfish.survey, acomp = dat.sim$survey.age.obs
)

L1.yrs <- styr:nyr
L1.nyr <- length(L1.yrs)
L1.cv.yr <- rep(cv.L1.fit, L1.nyr) # apply low CVs to fit the landings closely
L.ac1.yrs <- styr.ac1:nyr
L.ac1.nyr <- length(L.ac1.yrs)
nsamp.L1 <- nfish.L1 <- rep(n.ac.L1, L.ac1.nyr)
dat.L1 <- list(
  styr = styr, endyr = nyr, vals.obs = dat.sim$L1.obs[L1.yrs], cv = L1.cv.yr,
  nyr.ages = L.ac1.nyr, yrs.age = L.ac1.yrs,
  nsamp = nsamp.L1, nfish = nfish.L1, acomp = dat.sim$L1.age.obs[styr.ac1:nyr, ]
)

L2.yrs <- styr:nyr
L2.nyr <- length(L2.yrs)
L2.cv.yr <- rep(cv.L2.fit, L2.nyr) # apply low CVs to fit the landings closely
L.ac2.yrs <- styr.ac2:nyr
L.ac2.nyr <- length(L.ac2.yrs)
nsamp.L2 <- nfish.L2 <- rep(n.ac.L2, L.ac2.nyr)
dat.L2 <- list(
  styr = styr, endyr = nyr, vals.obs = dat.sim$L2.obs[L2.yrs], cv = L2.cv.yr,
  nyr.ages = L.ac2.nyr, yrs.age = L.ac2.yrs,
  nsamp = nsamp.L2, nfish = nfish.L2, acomp = dat.sim$L2.age.obs[styr.ac2:nyr, ]
)

D.yrs <- styr:nyr
D.nyr <- length(D.yrs)
D.cv.yr <- rep(cv.D.fit, D.nyr) # apply low CVs to fit the landings closely
D.lc.yrs <- styr.d:nyr
D.lc.nyr <- length(D.lc.yrs)
nsamp.D <- nfish.D <- rep(n.lc.D, D.lc.nyr)

dat.D <- list(
  styr = styr, endyr = nyr, vals.obs = dat.sim$D.obs[D.yrs], cv = D.cv.yr,
  D.lc.nyr = D.lc.nyr, yrs.len = D.lc.yrs,
  nsamp = nsamp.D, nfish = nfish.D, lcomp = dat.sim$D.len.obs[styr.d:nyr, ]
)

# Create parameter block for BAM
dat.parms <- data.frame(matrix(NA, nrow = 1, ncol = 7))
colnames(dat.parms) <- c("Guess", "LB", "UB", "Phase", "Prior.mu", "Prior.var", "pdf")
parms.Linf <- c(Linf, 500.0, 1200.0, -4, Linf, -0.25, 1)
dat.parms[1, ] <- parms.Linf
parms.K <- c(K, .01, 1.0, -4, K, -0.25, 1)
dat.parms <- rbind(dat.parms, parms.K)
parms.a0 <- c(a0, -5.0, -0.001, -4, a0, -0.25, 1)
dat.parms <- rbind(dat.parms, parms.a0)
parms.cvlen <- c(cv.laa, 0.05, 0.3, -4, cv.laa, -0.25, 1)
dat.parms <- rbind(dat.parms, parms.cvlen)
parms.M.scale <- c(0.0, -10.0, 10.0, -4, 0.0, 5.0, 1) # logit transform
dat.parms <- rbind(dat.parms, parms.M.scale)

parms.h <- c(h, 0.21, 0.99, -3, h, -0.25, 1)
dat.parms <- rbind(dat.parms, parms.h)
parms.logR0 <- c(log(R0), 10.0, 16.0, -1, log(R0), -0.25, 1)
dat.parms <- rbind(dat.parms, parms.logR0)
parms.Rrho <- c(0.0, -1.0, 1.0, -4, 0.0, -0.25, 1)
dat.parms <- rbind(dat.parms, parms.Rrho)
parms.logRsd <- c(logR.sd, 0.001, 1.0, -4, logR.sd, -0.25, 1)
dat.parms <- rbind(dat.parms, parms.logRsd)
if (SR.bc.switch) {
  bam.bc.switch <- -1.0
} else {
  bam.bc.switch <- 1.0
}


# dirichlet-multinomial pars by fleet or survey
parms.DMf1 <- c(2.0, -6.0, 6.0, 6, 0.0, 36.0, 3)
dat.parms <- rbind(dat.parms, parms.DMf1)
parms.DMf2 <- c(2.0, -6.0, 6.0, 6, 0.0, 36.0, 3)
dat.parms <- rbind(dat.parms, parms.DMf2)
parms.DMf3 <- c(2.0, -6.0, 6.0, 6, 0.0, 36.0, 3)
dat.parms <- rbind(dat.parms, parms.DMf3)
parms.DMsurvey <- c(2.0, -6.0, 6.0, 6, 0.0, 36.0, 3)
dat.parms <- rbind(dat.parms, parms.DMsurvey)

# selecitivity parameters
parms.f1a50 <- c(A50.sel1, 0.1, 10.0, 5, A50.sel1, -0.5, 3)
dat.parms <- rbind(dat.parms, parms.f1a50)
parms.f1slope <- c(slope.sel1, 0.1, 10.0, 5, slope.sel1, -0.5, 3)
dat.parms <- rbind(dat.parms, parms.f1slope)
parms.f2a50 <- c(A50.sel2, 0.1, 10.0, 5, A50.sel2, -0.5, 3)
dat.parms <- rbind(dat.parms, parms.f2a50)
parms.f2slope <- c(slope.sel2, 0.1, 10.0, 5, slope.sel2, -0.5, 3)
dat.parms <- rbind(dat.parms, parms.f2slope)
parms.f3a50a <- c(A50a.sel3, 0.1, 10.0, 5, A50a.sel3, -0.25, 3) # fix dome-shaped sels for now
dat.parms <- rbind(dat.parms, parms.f3a50a)
parms.f3slopea <- c(slopea.sel3, 0.1, 10.0, 5, slopea.sel3, -0.25, 3)
dat.parms <- rbind(dat.parms, parms.f3slopea)
parms.f3a50d <- c(A50d.sel3, 0.1, 10.0, 6, A50d.sel3, -0.25, 3)
dat.parms <- rbind(dat.parms, parms.f3a50d)
parms.f3sloped <- c(sloped.sel3, 0.1, 10.0, 6, sloped.sel3, -0.25, 3)
dat.parms <- rbind(dat.parms, parms.f3sloped)
parms.surveya50 <- c(A50.sel.survey, 0.1, 10.0, 5, A50.sel.survey, -0.5, 3)
dat.parms <- rbind(dat.parms, parms.surveya50)
parms.surveyslope <- c(slope.sel.survey, 0.1, 10.0, 5, slope.sel.survey, -0.5, 3)
dat.parms <- rbind(dat.parms, parms.surveyslope)

parms.logq <- c(log(sim.q), -18, -10, 2, log(sim.q), -1, 1)
dat.parms <- rbind(dat.parms, parms.logq)
if (styr == 1) {
  parms.Finitmult <- c(0.0, 0.0, 2.0, -1, 1.0, -1, 1)
} # Initialize at virgin if styr=1, to match popn model
if (styr > 1) {
  parms.Finitmult <- c(1.0, 0.0, 2.0, -1, 1.0, -1, 1)
} # Initialize using initial F

dat.parms <- rbind(dat.parms, parms.Finitmult)
parms.f1mu <- c(log(mean(F1)), -10, 1, 1, log(mean(F1)), -1, 1)
dat.parms <- rbind(dat.parms, parms.f1mu)
parms.f2mu <- c(log(mean(F2)), -10, 1, 1, log(mean(F2)), -1, 1)
dat.parms <- rbind(dat.parms, parms.f2mu)
parms.f3mu <- c(log(mean(F3)), -10, 1, 1, log(mean(F3)), -1, 1)
dat.parms <- rbind(dat.parms, parms.f3mu)

# Create initial dev vectors for recruits and F's. These can be true values or zeros.
logF1.resid <- c((log(F1[1:nyr]) - log(mean(F1[1:nyr]))))
logF2.resid <- c((log(F2[1:nyr]) - log(mean(F2[1:nyr]))))
logF3.resid <- c((log(F3[1:nyr]) - log(mean(F3[1:nyr]))))

F1.phase <- 3 # estimation phase for F1 devs
F2.phase <- 3 # estimation phase for F2 devs
F3.phase <- 3 # estimation phase for F3 devs
R.phase <- 2 # estimation phase for rec devs

dev.tseries <- list(
  logR.dev = logR.resid[styr.Rdevs:nyr],
  logF1.dev = logF1.resid[styr:nyr],
  logF2.dev = logF2.resid[styr:nyr],
  logF3.dev = logF3.resid[styr:nyr],
  F1.phase = F1.phase, F2.phase = F2.phase, F3.phase = F3.phase, R.phase = R.phase
)

# Miscellanous inputs
spawn.parms <- list(SR.switch = SR.switch, bam.bc.switch = bam.bc.switch, SPR.proxy = SPR.proxy, sp.frac = spawn.frac)
fec.parms <- list(fec.a = fec.a, fec.b = fec.b, fec.c = fec.c, fec.thresh = fec.thresh, fec.min = fec.min)
wgt.parms <- list(a.lw = a.lw, b.lw = b.lw)
lcomp.config <- list(nlenbins = nlenbins, binw = binw, lenbins = lenbins)

# create the BAM input file
BAM.write.dat(
  fname = "BAM-Sim.dat", styr = styr, styr.Rdevs = styr.Rdevs, nyr = nyr, nages = nages,
  spawn.parms <- spawn.parms,
  lcomp.config <- lcomp.config,
  dat.survey = dat.survey,
  dat.L1 = dat.L1, dat.L2 = dat.L2, dat.D = dat.D, parms = dat.parms, dev.tseries = dev.tseries,
  prop.f = proportion.female, mat.age = mat.age,
  wgt.parms = wgt.parms, fec.parms = fec.parms,
  M.age = M.age.unscaled
)

#####################################################################################
### Run BAM on simulated data
if (compile) {
  shell(paste("admb ", filename, sep = ""))
}
filename.dat <- paste(filename, ".dat", sep = "")
filename.rdat <- paste(filename, ".rdat", sep = "")
bamrun <- paste(basename(filename), ".exe", sep = "")
run.command <- paste(bamrun, admb.switch, "-ind", filename.dat, sep = " ")
shell(run.command)
shell(paste("01cleanup.bat"))

bam <- dget(filename.rdat)

#####################################################################################
# Some data plots

#####################################################################################
bam.nyr <- length(bam$t.series$year) - 1 # remove the projection year

windows(height = 10, width = 8, record = T)
mat <- matrix(1:8, ncol = 2, nrow = 4)
layout(mat = mat, widths = rep.int(1, ncol(mat)), heights = rep.int(1, nrow(mat)))
par(las = 1, mar = c(4.5, 4.5, 1, 0.5))

plot(sim$yr, sim$N.age[, 1] / 1000,
  type = "o", lwd = 2, ylab = "Recruits (1000 fish)",
  xlab = "", panel.first = grid(lty = 1), col = "blue",
  ylim = c(0, max(c(sim$N.age[, 1], bam$t.series$recruits)) / 1000)
)
points(bam$t.series$year[1:bam.nyr], bam$t.series$recruits[1:bam.nyr] / 1000, pch = 3, col = "black", cex = 1.2)
legend("top", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))

plot(sim$yr, logR.resid,
  type = "o", lwd = 2, ylab = "Recruitment residuals",
  xlab = "", panel.first = grid(lty = 1), col = "blue",
  ylim = c(min(c(logR.resid, bam$t.series$logR.dev)), max(c(logR.resid, bam$t.series$logR.dev)))
)
points(bam$t.series$year[1:bam.nyr], bam$t.series$logR.dev[1:bam.nyr], pch = 3, col = "black", cex = 1.2)
legend("top", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))

plot(sim$yr, sim$abundance / 1000,
  type = "o", lwd = 2, ylab = "Abundance (1000 fish)",
  xlab = "", panel.first = grid(lty = 1), col = "blue",
  ylim = c(0, max(c(sim$abundance, bam$t.series$N[1:bam.nyr])) / 1000)
)
points(bam$t.series$year[1:bam.nyr], bam$t.series$N[1:bam.nyr] / 1000, pch = 3, col = "black", cex = 1.2)
legend("top", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))


plot(sim$yr, sim$SSB,
  type = "o", lwd = 2, ylab = "Spawning biomass (mt)",
  xlab = "", panel.first = grid(lty = 1), col = "blue",
  ylim = c(0, max(c(sim$SSB, bam$t.series$SSB[1:bam.nyr])))
)
abline(h = sim$msy$SSBmsy, lty = 2)
text(sim$yr[4], sim$msy$SSBmsy, pos = 3, "SSBmsy")
points(bam$t.series$year[1:bam.nyr], bam$t.series$SSB[1:bam.nyr], pch = 3, col = "black", cex = 1.2)
legend("top", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))


plot(sim$yr, F1,
  type = "o", lwd = 2, ylab = "Fleet 1 fishing rate (/yr)", xlab = "",
  ylim = c(0, max(c(F1, F2, F3, bam$t.series$F.fleet1[1:bam.nyr], bam$t.series$F.fleet2[1:bam.nyr]))),
  panel.first = grid(lty = 1), col = "blue"
)
points(bam$t.series$year[1:bam.nyr], bam$t.series$F.fleet1[1:bam.nyr], pch = 3, col = "black", cex = 1.2)
legend("top", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))


plot(sim$yr, F2,
  type = "o", lwd = 2, ylab = "Fleet 2 fishing rate (/yr)", xlab = "",
  ylim = c(0, max(c(F1, F2, F3, bam$t.series$F.fleet1[1:bam.nyr], bam$t.series$F.fleet2[1:bam.nyr]))),
  panel.first = grid(lty = 1), col = "blue"
)
points(bam$t.series$year[1:bam.nyr], bam$t.series$F.fleet2[1:bam.nyr], pch = 3, col = "black", cex = 1.2)
legend("top", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))

plot(sim$yr, F3,
  type = "o", lwd = 2, ylab = "Fleet 3 fishing rate (/yr)", xlab = "",
  ylim = c(0, max(c(F1, F2, F3, bam$t.series$F.fleet3[1:bam.nyr], bam$t.series$F.fleet3[1:bam.nyr]))),
  panel.first = grid(lty = 1), col = "blue"
)
points(bam$t.series$year[1:bam.nyr], bam$t.series$F.fleet3[1:bam.nyr], pch = 3, col = "black", cex = 1.2)
legend("top", legend = c("True", "Predicted"), lty = c(1, -1), pch = c(-1, 3), col = c("blue", "black"))

#####################################################################################
library(FishGraph)
ptype <- NULL # (no quotes) for no plots saved; other options: "pdf", "wmf", "eps"
dtype <- TRUE # draft type
ctype <- TRUE # color type
windows(width = 8, height = 8, record = TRUE)
Comp.plots(bam, draft = dtype, use.color = ctype, graphics.type = ptype, c.min = 0.25, corr = FALSE)
windows(width = 8, height = 10, record = TRUE)
Comp.yearly.plots(bam,
  draft = dtype, use.color = ctype, graphics.type = ptype, plot.neff = FALSE, print.neff = TRUE,
  compact = TRUE, print.n = TRUE
)

#####################################################################################
RE.h <- switch(SR.switch,
  ((bam$parms$BH.steep - h) / h),
  ((bam$parms$Ricker.steep - h) / h),
  NA
)
RE.R0 <- switch(SR.switch,
  ((bam$parms$BH.R0 - R0) / R0),
  ((bam$parms$Ricker.R0 - R0) / R0),
  ((bam$parms$Mean.R0 - R0) / R0)
)
RE.M <- (bam$parms$M.scale - M.mult) / M.mult
RE.Fproxy <- (bam$parms$Fproxy - sim$spr$Fproxy) / sim$spr$Fproxy
RE.SSBFproxy <- (bam$parms$SSB.Fproxy - sim$spr$SSB.Fproxy) / sim$spr$SSB.Fproxy
RE.MSY <- (bam$parms$msy.mt - sim$msy$msy) / sim$msy$msy
RE.SSBmsy <- (bam$parms$SSBmsy - sim$msy$SSBmsy) / sim$msy$SSBmsy
RE.Fmsy <- (bam$parms$Fmsy - sim$msy$Fmsy) / sim$msy$Fmsy

if (SR.switch == 1 || SR.switch == 2) {
  RE <- c(RE.h, RE.R0, RE.M, RE.MSY, RE.SSBmsy, RE.Fmsy)
  RE.names <- c("h", "R0", "M", "MSY", "SSBmsy", "Fmsy")
} else if (SR.switch == 3) {
  RE <- c(RE.R0, RE.M, RE.Fproxy, RE.SSBFproxy)
  RE.names <- c("R0", "M", "Fproxy", "SSBFproxy")
}
windows(height = 5, width = 6, record = T)
plot(1:length(RE), RE,
  type = "p", pch = 15, col = "blue", ylim = c(-max(c(0.1, abs(RE))), max(c(0.1, abs(RE)))), xaxt = "n",
  ylab = "Relative Error", xlab = ""
)
axis(1, at = 1:length(RE), labels = RE.names)
abline(h = 0)

# RE <- c(RE.h, RE.R0, RE.M, RE.MSY, RE.SSBmsy, RE.Fmsy)
# RE.names <- c("h", "R0", "M", "MSY", "SSBmsy", "Fmsy")
# windows(height = 5, width = 6, record = T)
# plot(1:length(RE), RE,
#      type = "p", pch = 15, col = "blue", ylim = c(-max(abs(RE)), max(abs(RE))), xaxt = "n",
#      ylab = "Relative Error", xlab = ""
# )
# axis(1, at = 1:length(RE), labels = RE.names)
# abline(h = 0)


####################################################################################
# Some more plots
library(colorspace)
cols <- sequential_hcl(5, "Viridis")
windows(height = 8, width = 6, record = T)
# tiff(filename="Fig-F.tiff", width=5.5, height=4, units="in", compression="none", res=600)
mat <- matrix(1:2, ncol = 1, nrow = 2)
layout(mat = mat, widths = rep.int(1, ncol(mat)), heights = rep.int(1, nrow(mat)))
par(las = 1, mar = c(4.1, 4.25, 0.5, 0.5), cex = 0.75)

plot(1:nyr, Ftotal, type = "l", lwd = 2, xlab = "Year", ylab = "Fishing mortality rate", ylim = c(0, max(Ftotal, F1, F2)), col = cols[1])
lines(1:nyr, F1.mu, lwd = 2, col = cols[2])
lines(1:nyr, F2.mu, lwd = 2, col = cols[3])
lines(1:nyr, F3.mu, lwd = 2, col = cols[4])
lines(1:nyr, F1, lty = 2, col = cols[2])
lines(1:nyr, F2, lty = 2, col = cols[3])
lines(1:nyr, F3, lty = 2, col = cols[4])
legend("left",
  legend = c("Total mean", "Fleet1 mean", "Fleet2 mean", "Fleet3 mean", "Fleet1 realized", "Fleet2 realized", "Fleet3 realized"),
  lty = c(1, 1, 1, 1, 2, 2, 2), lwd = c(2, 2, 2, 2, 1, 1, 1), col = c(cols[1], cols[2], cols[3], cols[4], cols[2], cols[3], cols[4])
)
text(1, 0.98 * max(Ftotal, F1, F2, F3), "(A)", cex = 1.1)

plot(1:nages, selex1, type = "l", lwd = 2, xlab = "Age", ylab = "Selectivity", ylim = c(0, 1), col = cols[2], xlim = c(1, 7))
lines(1:nages, selex2, lty = 1, lwd = 2, col = cols[3])
lines(1:nages, selex3, lty = 1, lwd = 2, col = cols[4])
lines(1:nages, selex.survey, lty = 1, lwd = 2, col = cols[5])
points(1:nages, bam$sel.age$sel.m.fleet1[bam.nyr, ], pch = 3, col = cols[2], cex = 1.3, lwd = 3)
points(1:nages, bam$sel.age$sel.m.fleet2[bam.nyr, ], pch = 3, col = cols[3], cex = 1.3, lwd = 3)
points(1:nages, bam$sel.age$sel.m.fleet3[bam.nyr, ], pch = 3, col = cols[4], cex = 1.3, lwd = 3)
points(1:nages, bam$sel.age$sel.m.survey1[bam.nyr, ], pch = 3, col = cols[5], cex = 1.3, lwd = 3)

legend("right",
  legend = c("True Fleet1", "True Fleet2", "True Fleet3", "True Survey", "Predicted"),
  lty = c(1, 1, 1, 1, -1), lwd = 2, col = c(cols[2], cols[3], cols[4], cols[5], "black"), pch = c(-1, -1, -1, -1, 3)
)
text(1, 0.98, "(B)", cex = 1.1)

# savePlot(filename=paste("Fig-Fselex.tiff", sep=""), type="tiff")

cols <- sequential_hcl(6, "Viridis")
windows(height = 8, width = 10, record = T)
mat <- matrix(1:4, ncol = 2, nrow = 2)
layout(mat = mat, widths = rep.int(1, ncol(mat)), heights = rep.int(1, nrow(mat)))
par(las = 1, mar = c(4.5, 4.5, 1, 0.5))

plot(sim$yr, sim$L1.mt,
  type = "l", lwd = 2, ylab = "Fleet 1 landings (mt)", xlab = "",
  ylim = c(0, max(sim$L1.mt, dat.sim$L1.obs, bam$t.series$L.fleet1.pr, na.rm = T)), panel.first = grid(lty = 1), col = cols[1]
)
points(sim$yr, dat.sim$L1.obs, pch = 16, col = cols[3], cex = 1.3)
points(bam$t.series$year[1:bam.nyr], bam$t.series$L.fleet1.pr[1:bam.nyr], pch = 3, col = cols[5], cex = 1.3)
legend("topright", legend = c("True", "Observed", "Predicted"), lty = c(1, -1, -1), pch = c(-1, 16, 3), col = c(cols[1], cols[3], cols[5]))
text(1, 4100, "(A)", cex = 1.1)

plot(sim$yr, sim$L2.mt,
  type = "l", lwd = 2, ylab = "Fleet 2 landings (mt)", xlab = "",
  ylim = c(0, max(sim$L2.mt, dat.sim$L2.obs, bam$t.series$L.fleet2.pr, na.rm = T)), panel.first = grid(lty = 1), col = cols[1]
)
points(sim$yr, dat.sim$L2.obs, pch = 16, col = cols[3], cex = 1.3)
points(bam$t.series$year[1:bam.nyr], bam$t.series$L.fleet2.pr[1:bam.nyr], pch = 3, col = cols[5], cex = 1.3)
# legend("top", legend = c("True", "Observed", "Predicted"), lty = c(1, -1, -1), pch = c(-1, 16, 3), col = c(cols[1], cols[3], cols[5]))
text(1, 1300, "(B)", cex = 1.1)

plot(sim$yr, sim$D.mt,
  type = "l", lwd = 2, ylab = "Fleet 3 discards (mt)", xlab = "",
  ylim = c(0, max(sim$D.mt, dat.sim$D.obs, bam$t.series$L.fleet3.pr, na.rm = T)), panel.first = grid(lty = 1), col = cols[1]
)
points(sim$yr, dat.sim$D.obs, pch = 16, col = cols[3], cex = 1.3)
points(bam$t.series$year[1:bam.nyr], bam$t.series$L.fleet3.pr[1:bam.nyr], pch = 3, col = cols[5], cex = 1.3)
# legend("top", legend = c("True", "Observed", "Predicted"), lty = c(1, -1, -1), pch = c(-1, 16, 3), col = c(cols[1], cols[3], cols[5]))
text(1, 1300, "(B)", cex = 1.1)

plot(sim$yr, rep(0, nyr),
  type = "n", lwd = 2, ylab = "Survey (scaled)", xlab = "",
  ylim = c(0, max(survey.sim, dat.sim$survey.obs, bam$t.series$U.survey1.pr, na.rm = T)), panel.first = grid(lty = 1), col = cols[1]
)
lines(survey.yrs, survey.sim, lwd = 2, col = cols[1])
points(survey.yrs, dat.sim$survey.obs, pch = 16, col = cols[3], cex = 1.3)
points(survey.yrs, bam$t.series$U.survey1.pr[(styr.survey - styr + 1):bam.nyr], pch = 3, col = cols[5], cex = 1.3)
# legend("top", legend = c("True", "Observed", "Predicted"), lty = c(1, -1, -1), pch = c(-1, 16, 3), col = c(cols[1], cols[3], cols[5]))
text(1, 3, "(C)", cex = 1.1)

# savePlot(filename=paste("Fig-dynamics.tiff", sep=""), type="tiff")

# dev.off()
