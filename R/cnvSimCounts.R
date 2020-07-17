##----------------------------------------------------------------------------##
## cnvSimCounts: generate simulated molecule counts
##----------------------------------------------------------------------------##

#' @name cnvSimCounts
#' @title Generate simulated molecule counts 
#' 
#' @param totalMolecules integer of length 1, the total number of molecules 
#' @param interval data.table [interval object][validObjects] with 
#' 'captureProb' field, see details
#' @param subject subject name/identifier 
#' @param variantWidth integer, gives the possible variant widths in contiguous
#' intervals, see details
#' @param CN numeric vector of possible copy numbers; 1.0 indicates diploid 
#' state, see details
#' @param cnProb numeric vector of probabilities corresponding to the possible
#' copy states in 'CN'
#' @param seed integer, passed to [set.seed()][base::set.seed()]
#' 
#' @details 
#' \code{cnvSimCounts} requires an [interval object][validObjects] with an 
#' added field, 'captureProb' defining the multinomial probability distribution
#' for interval coverage. 
#' In this multinomial distribution, a success at an interval indicates the 
#' interval was covered by a sequencing molecule.
#' 
#' \code{cnvSimCounts} will simulate variable-width copy number variants, with
#' the possible widths (number of contiguous intervals) given by the 
#' 'variantWidth' parameter.
#' All variant widths are simulated at an equal probability. 
#' 
#' The 'CN' parameter defines the possible copy states. 
#' To simply the computations, the \code{mcCNV} package defines 1.0 as the 
#' diploid state. 
#' The "actual" copies are given by multiplying 'CN' by 2. 
#' As such, all entries in the 'CN' parameter must be a multiple of 0.5.
#' 
#' We adjust the capture probabilities by multiplying the probability by the
#' simulated copy number.
#' For example, when the copy number is 1 (the diploid state), we do not wish
#' to adjust the probability. 
#' However, if say 3 copies of the interval are present, the probability of 
#' capturing that interval is increased by 1.5.
#' 
#' We have found, likely due to sequencing and mapping errors, even true 
#' homozygous deletions can have a few reads.
#' We account for this by using the multiplier 0.001 for intervals with complete
#' deletions (copy number is 0.0).
#' 
#' @seealso cnvSimCounts cnvSimPool
#' 
#' @import data.table
#' @export 


cnvSimCounts <- function(totalMolecules = 1e7L,
                         interval = cnvSimInterval(),
                         subject = "simulatedSubject",
                         variantWidth = 1L,
                         CN = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0), 
                         cnProb = c(1e-6, 6.25e-4, 0.998747, 6.25e-4, 1e-6, 
                                    2.5e-7, 2.5e-7, 2.5e-7, 2.5e-7),
                         seed = NULL) {
  
  ## Input checks; the sample function rescales the 'prob' parameter, and the
  ## capture probabilities are rescaled after multiplying the CN, so no checking
  ## for valid probabilities is required
  stopifnot(is.integer(totalMolecules))
  stopifnot(cnvValidInterval(interval))
  stopifnot(is.numeric(interval$captureProb))
  stopifnot(length(subject) == 1)
  stopifnot(!anyNA(cnProb))
  stopifnot(is.integer(variantWidth))
  stopifnot(length(CN) == length(cnProb))
  stopifnot(is.numeric(CN))
  stopifnot(is.numeric(cnProb))
  stopifnot(all(CN %% 0.5 == 0))
  if (!all(c(0.0, 0.5, 1.0, 1.5) %in% CN)) {
    warning("'CN' does not contain common copy-states, ",
            "e.g. 0.5 (1 copy), 1.0 (2 copies); see ?cnvSimCounts")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  cnts <- copy(interval)
  cnts[ , actCN := sample(x = CN, size = .N, replace = TRUE, prob = cnProb)]
  
  if (variantWidth > 1 || length(variantWidth > 1)) {
    setindex(cnts, actCN)
    if (cnts[ , sum(actCN != 1)] > 0) {
      cnts[ , vwid := 0L]
      cnts[actCN != 1, 
           vwid := sample(x = variantWidth, size = .N, replace = TRUE)]
      ## Expanding by seqnames here makes sure a variant doesn't spill off the
      ## end of one chromosome onto the beginning of another, or off the end of
      ## the genome
      ind <- cnts[ , 
                  .(stop = max(.I), 
                    ind = unlist(mapply(seq, .I, length = vwid))),
                  by = seqnames]
      ind[ , newCN := cnts[actCN != 1, unlist(mapply(rep, actCN, each = vwid))]]
      ind <- ind[ind <= stop]
      cnts[ind$ind, actCN := ind$newCN]
      cnts[ , vwid := NULL]
    }
  }
  
  cnts[ , 
       molCount := rmultinom(n = 1,
                             size = totalMolecules,
                             prob = captureProb*pmax(0.001, actCN))]
  cnts[ , subject := subject]
  cnts[ , nCoverMult := 0L]
  cnts[]
  
}
