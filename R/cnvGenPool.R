##----------------------------------------------------------------------------##
## cnvGenPool: generate a simulated pool of samples
##----------------------------------------------------------------------------##

#' @name cnvGenPool
#' @title Generate a simulated pool of samples
#' 
#' @param ns integer of length 1, the number of samples
#' @param wndw integer of length 2, the bounds for the number of molecules per
#' sample
#' @param ne integer of length 1, the number of exons
#' @param cs numeric, the possible copy states
#' @param pc numeric, the prior probability for the copy state
#' @param cw integer of length 1, the width (number of exons) the cnv spans
#' @param seed numeric of length 1, the starting seed for the random number
#' generator, can be left NULL
#' 
#' @details 
#' The number of molecules per sample is drawn from a uniform distribution with
#' the bounds given by `wndw`.
#' 
#' For each pool a psuedo genome is created by modeling the probability of 
#' an exon receving a molecule during the capture process as a multinomial 
#' distribution, where each exon has its own probability. The multinomial 
#' distribution is given by drawing from a log-normal 
#' distribution with mean -12.36 and sd 0.7393. See 
#' \code{vignette("lognormDerivation", "cnvR")} for more information.
#' 
#' `cs`, `pc`, `cw` are all passed to \code{\link{cnvGenSmpl}}. Leaving `cs` &
#' `pc` NULL will use the default values. 
#' 
#' `seed` sets the seed for defining the multinomial distribution of 
#' probabilites that each exon gets a molecule. 
#' 
#' The number of reads per sample is passed as the seed to
#' \code{cnvGenSmpl}. This allows the simulated datasets to be completely
#' reproducible. 
#' 
#' @importFrom stats rlnorm runif
#' @import data.table
#' @export 

cnvGenPool <- function(ns, ne, cw, wndw, pc = NULL, cs = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  nr <- runif(ns, min = min(wndw), max = max(wndw))
  exonProbs <- rlnorm(ne, -12.36, 0.7393)
  exonProbs <- exonProbs/sum(exonProbs)
  p <- list(pe = exonProbs, cw = cw)
  if(!is.null(pc)) p[c("pc", "cs")] <- list(pc = pc, cs = cs)
  smpls <- mapply(cnvGenSmpl, nr = nr, seed = nr, 
                  SIMPLIFY = FALSE, MoreArgs = p)
  smpls <- rbindlist(smpls)
  smpls[ , sbj := rep(sprintf("sbj%0.2d", 1:ns), each = ne)]
  smpls[ , ref := sprintf("ref%0.6d", ref)]
  smpls[]
}
