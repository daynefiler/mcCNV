##----------------------------------------------------------------------------##
## cnvSimInterval
##----------------------------------------------------------------------------##

#' @name cnvSimInterval
#' @title Generate a simulated interval
#' 
#' @param n integer of length 1, the number of intervals
#' @param seqnames character of length 1, 'seqnames' for the simulated intervals
#' @param meanlog numeric of length 1, passed to [rlnorm][stats::rlnorm]
#' to define distribution of interval capture probabilities
#' @param sdlog numeric of length 1, passed to [rlnorm][stats::rlnorm]
#' to define distribution of interval capture probabilities
#' @param seed integer, passed to [set.seed()][base::set.seed()]
#' 
#' @details 
#' \code{cnvSimInterval} generates an [interval object][validObjects] with  
#' randomly-generated capture probabilities required for [cnvSimCounts()]. 
#' Interval capture probabilities are drawn initially from a lognormal 
#' distribution, then normalized to a sum of 1. 
#' The capture probabilities are stored in 'captureProb.'
#' 
#' The default values for meanlog & sdlog were derived from multiplexed whole-
#' exome capture of 16 human samples using the Agilent XT2 capture platform.
#' 
#' \code{cnvSimInterval} does not intend to make realistic exon sizes or 
#' spacing. 
#' All simulated intervals are 50 basepairs wide with start sites every 200 
#' basepairs.
#' The structure of the object is simple and users can easily create their own
#' interval objects containing capture probabilities if greater complexity is 
#' desired.
#' 
#' Setting the 'seed' parameter allows for reproducible simulation studies.
#' 
#' @seealso cnvSimCounts cnvSimPool
#' 
#' @importFrom stats rlnorm 
#' @import data.table
#' @export 

cnvSimInterval <- function(n = 100000, 
                           seqnames = "simSeq", 
                           meanlog = -12.36, 
                           sdlog = 0.7393,
                           seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  interval <- data.table(seqnames = "simSeq", start = (1:n) * 200L)
  interval[ , end := start + 50L]
  interval[ , captureProb := rlnorm(n, meanlog, sdlog)]
  interval[ , captureProb := captureProb/sum(captureProb)]
  interval[]
  
}
