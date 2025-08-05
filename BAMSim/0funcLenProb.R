# =================================================================================
# Compute matrix of proportion of length at age
# Matrix has nages rows and nlenbins columns, where columns sum to one
# KWS and MDD  Jan 2025
#---------------------------------------------------------------------------------

lenprob.fcn <- function(len.mu, len.cv, len.bins, len.binw) {
  ## INPUT:
  # len.mu = mean length at age
  # len.cv = cv of len at age
  # len.bins = length bins
  # len.binw = bin width

  ## OUTPUT:
  #    lenprob = matrix (nages by nlenbins) of proportion at length by age

  #####################################################################################
  nages <- length(len.mu)
  nlenbins <- length(len.bins)
  len_sd <- array(NA, dim = c(nages))
  cprob_lenvec <- array(NA, dim = c(nlenbins))
  lenprob <- array(NA, dim = c(nages, nlenbins))
  for (i in 1:nages) {
    len_sd[i] <- len.cv * len.mu[i]
    zscore_lzero <- (0.0 - len.mu[i]) / len_sd[i]
    cprob_lzero <- pnorm(zscore_lzero)
    # first length bin
    zscore_len <- ((len.bins[1] + 0.5 * len.binw) - len.mu[i]) / len_sd[i]
    cprob_lenvec[1] <- pnorm(zscore_len)
    lenprob[i, 1] <- cprob_lenvec[1] - cprob_lzero
    # most other length bins
    for (l in 2:(nlenbins - 1)) {
      zscore_len <- ((len.bins[l] + 0.5 * len.binw) - len.mu[i]) / len_sd[i]
      cprob_lenvec[l] <- pnorm(zscore_len)
      lenprob[i, l] <- cprob_lenvec[l] - cprob_lenvec[l - 1]
    }
    # last length bin is a plus group
    zscore_len <- ((len.bins[nlenbins] - 0.5 * len.binw) - len.mu[i]) / len_sd[i]
    lenprob[i, nlenbins] <- 1.0 - pnorm(zscore_len)
    lenprob[i, ] <- lenprob[i, ] / (1.0 - cprob_lzero) # renormalize to account for any prob mass below 0
  }
  return(lenprob)
}
