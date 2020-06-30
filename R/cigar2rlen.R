#' @title Convert cigar string to length
#' @description Converts cigar string to length
#' @useDynLib mcCNV, .registration = TRUE, .fixes = "C_"
#' @export

cigar2rlen <- function(x) {

  n <- length(x)
  .C(C_cigar2rlen, as.character(x), n, integer(n))[[3]]

}
