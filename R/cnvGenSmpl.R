##----------------------------------------------------------------------------##
## cnvGenSmpl: generate a simulated sample
##----------------------------------------------------------------------------##

#' @name cnvGenSmpl
#' @title Generate a simulated sample
#' 
#' @param nr integer of length 1, the number of molecules for the sample
#' @param pe numeric, the multinomial probability at each exon
#' @param cs numeric, the possible copy states 
#' @param pc numeric, the prior probability for the copy state
#' @param cw integer of length 1, the width (number of exons) the cnv spans
#' @param seed numeric of length 1, the starting seed for the random number
#' generator, can be left NULL
#' 
#' @details 
#' `pe` essentiall provides the "genome" and needs to be consistent across
#' samples that will be compared. 
#' 
#' `cs` and `pc` must be the same length. The possible copy states assume the 
#' normal copy state (2 copies) is 1. Likewise, a single copy of an exon would
#' be denoted as 0.5. The values we suggest for the possible states and their
#' priors are given as the default values. 
#' 
#' @import utils
#' @import data.table
#' @export 


cnvGenSmpl <- function(nr, pe, cw, seed = NULL,
                       cs = c(0.001, 0.5, 1, 1.5, 2.0, 2.5, 3.0, 3.5, 4), 
                       pc = c(1e-6, 6.25e-4, 0.998747, 6.25e-4, 1e-6, 2.5e-7,
                              2.5e-7, 2.5e-7, 2.5e-7)) {
  
  if (length(pc) != length(cs)) stop("Invalid 'pc' and 'cs'.")
  if (!is.null(seed)) set.seed(seed)
  
  ne <- length(pe)
  states <- sample(cs, ne, TRUE, pc)
  if (cw > 1) {
    ind <- which(states != 1)
    states[sapply(ind, seq, length = cw)] <- rep(states[ind], each = cw)
  }
  adj_pe <- states*pe
  adj_pe <- adj_pe/sum(adj_pe)
  ref <- sample(ne, nr, TRUE, adj_pe)
  res <- data.table(ref)[ , .N, by = ref][order(ref)]
  setkey(res, ref)
  res <- res[J(1:ne), ]
  res[is.na(N), N := 0]
  res[ , act_cnvs := states]
  
}
