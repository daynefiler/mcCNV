##----------------------------------------------------------------------------##
## cnvSimPool: generate a simulated pool of samples
##----------------------------------------------------------------------------##

#' @name cnvSimPool
#' @title Generate a simulated pool of samples
#' @description \code{cnvSimPool} is a convenience wrapper to generate a pool
#' of samples from the same interval object to mimic a multiplexed capture.
#' 
#' @param nSubjects integer of length 1, the number of samples
#' @param countRange integer of length 2, the bounds for the number of 
#' molecules per sample, see details
#' @inheritParams cnvSimCounts
#' @param ... passed to [cnvSimCounts()]
#' 
#' @details 
#' The number of molecules per sample is drawn from a uniform distribution with
#' the bounds given by 'countRange'.
#' 
#' The number of reads per sample is passed as the seed to [cnvSimCounts]. 
#' This allows the simulated pools to be completely reproducible.
#' 
#' @seealso cnvSimCounts cnvSimPool
#' 
#' @importFrom stats runif
#' @import data.table
#' @export 

cnvSimPool <- function(nSubjects = 16L,
                       countRange = c(8e6L, 14e6L),
                       interval = cnvSimInterval(),
                       seed = NULL,
                       ...) {
  
  stopifnot(is.integer(nSubjects) && length(nSubjects) == 1)
  stopifnot(is.integer(countRange))
  stopifnot(cnvValidInterval(interval))
  stopifnot(is.numeric(interval$captureProb))
  
  if (!is.null(seed)) set.seed(seed)
  nr <- runif(nSubjects, min = min(countRange), max = max(countRange))
  nr <- as.integer(nr)
  sn <- sprintf("simSbj-%03d", seq(nSubjects))
  counts <- mapply(cnvSimCounts, 
                   totalMolecules = nr, 
                   subject = sn,
                   seed = nr,
                   SIMPLIFY = FALSE, 
                   MoreArgs = list(interval = interval, ...))
  
  counts <- rbindlist(counts)
  counts[]
  
}
