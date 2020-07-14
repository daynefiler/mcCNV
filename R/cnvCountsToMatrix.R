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
  counts <- copy(counts)
  svec <- unique(counts$subject)
  counts[ , intName := sprintf("%s:%d-%d", seqnames, start, end)]
  frm <- formula(intName + start + seqnames ~ subject)
  mat <- dcast(counts, frm, value.var = "molCount")
  setorder(mat, seqnames, start)
  mat[ , c("start", "seqnames") := NULL]
  mat <- as.matrix(mat, rownames = "intName")
  mat
}

