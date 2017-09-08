#' @title Convert cigar string to length
#' @description Converts cigar string to length
#' @useDynLib cnvR, .registration = TRUE
#' @export

cigar2rlen <- function(x) {

  n <- length(x)
  lens <- .C("cigar2rlen", as.character(x), as.integer(n), result = integer(n),
             PACKAGE = "cnvR")
  lens[["result"]]

}
