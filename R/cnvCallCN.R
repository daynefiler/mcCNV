##----------------------------------------------------------------------------##
## cnvCallCN
##----------------------------------------------------------------------------##

#' @name cnvCallCN
#' @title Call copy number-state for the given data
#' 
#' @param counts data.table [counts object][validObjects]
#' @param width integer, window of contiguous intervals to calculate CN on; can 
#' provide a vector with multiple widths (defaults to 1)
#' @param prior numeric of length 1, the prior probability used to estimate
#' copy number variants; see details
#' @param delta integer of length 1, minimum copy number state changes to stop
#' the algorithm; see details
#' @param iterations integer of length 1, the maximum number of iterations 
#' before stopping algorithm; see details
#' @param verbose TRUE/FALSE
#' 
#' @details 
#' Need to add
#' 
#' @import data.table
#' @export

cnvCallCN <- function(counts, 
                      prior, 
                      width = 1L, 
                      delta = 20L, 
                      iterations = 30L, 
                      verbose = TRUE) {
  
  stopifnot(cnvValidCounts(counts))
  stopifnot(is.numeric(prior) && length(prior) == 1)
  stopifnot(is.integer(width))
  stopifnot(is.integer(delta) && length(delta) == 1)
  stopifnot(is.integer(iterations) && length(iterations) == 1)
  stopifnot(is.logical(verbose) && length(verbose) == 1)
  
  ## Create interval names as seqnames:start-end
  counts[ , intName := sprintf("%s:%d-%d", seqnames, start, end)]
  
  if (verbose) cat("Collapsing intervals...")
  counts <- rbindlist(lapply(width, .clpsExon, dat = counts))
  if (verbose) cat("done.\n")
  
  if (verbose) cat("Calling CNVs...")
  calls <- .callCN(cnts = counts, 
                   delta = delta, 
                   iterations = iterations, 
                   prior = prior,
                   shrink = TRUE)
  if (verbose) cat("done.\n")
  
  ## Overwrite placeholder CN == 0.001; negate calls where a mean was not 
  ## established
  calls[CN == 0.001, CN := 0.0]
  calls[is.na(mn), CN := NA]
  
  calls[ , c("adjN", "use", "geomn", "seqnames") := NULL]
  n0 <- c("phi", "mn", "vr", "width", "sf", "llk", "lp", "lp1")
  n1 <- c("intPhi", "intMean", "intSD", "intWidth", "sbjSizeFactor", 
          "cnLogLik", "cnLogP", "cn1LogP")
  setnames(calls, n0, n1)
  
  calls[]
  
}