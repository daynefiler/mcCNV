##----------------------------------------------------------------------------##
## cnvCountsToMatrix
##----------------------------------------------------------------------------##

#' @name cnvCountsToMatrix
#' @title Convert counts object to interval by sample matrix
#' @description \code{cnvCountsToMatrix} is a convenience function for casting
#' a long-format [counts object][validObjects] into an interval by sample 
#' counts matrix.
#' 
#' @param counts data.table [counts object][validObjects]
#' 
#' @return interval by sample count matrix
#' 
#' @import data.table
#' @export

cnvCountsToMatrix <- function(counts) {
  stopifnot(cnvValidCounts(counts))
  svec <- unique(counts$subject)
  counts[ , intName := sprintf("%s:%d-%d", seqnames, start, end)]
  mat <- dcast(counts, intName ~ subject, value.var = "molCount")
  mat <- as.matrix(mat, rownames = "intName")
}

