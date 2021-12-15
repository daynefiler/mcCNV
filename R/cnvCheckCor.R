##----------------------------------------------------------------------------##
## cnvCheckCor: Calculate sample-sample & sample-pool correlations
##----------------------------------------------------------------------------##

#' @name cnvCheckCor
#' @title Calculate sample-sample & sample-pool correlations
#' 
#' @param counts data.table [counts object][validObjects]
#' 
#' @return 
#' List with two elements:
#' \describe{
#'   \item{sampleSample}{pair-wise correlations}
#'   \item{samplePool}{correlation of each sample to the sum counts
#'   of all other samples}
#' }
#' 
#' @import data.table
#' @importFrom stats cor
#' @export

cnvCheckCor <- function(counts) {
  stopifnot(cnvValidCounts(counts))
  mat <- cnvCountsToMatrix(counts)
  ## subject-subject correlations
  ss <- cor(mat)
  ## subject-pool correlations
  spCor <- function(s) cor(mat[ , s], rowSums(mat[ , colnames(mat) != s]))
  sp <- sapply(colnames(mat), spCor)
  list(sampleSample = ss, samplePool = sp)
}