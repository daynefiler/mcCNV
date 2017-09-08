##----------------------------------------------------------------------------##
## cnvFilterReads
##----------------------------------------------------------------------------##

#' @name filterReads
#' @title Modifies given data.table
#' 
#' @description \code{cnvFilterReads} modifies the given data.table object, 
#' adding the 'ff' ("filter flag") field, specifying which reads to remove
#' from furhter analysis. 
#' 
#' @param x data.table object to be modified by reference

.filterReads <- function(x) {
  
  x[ , ff := 0L]
  
  ## Keep only reads paired to a read on the reverse strand of same reference
  x[!flag %in% c(99, 147, 163, 83, 97, 145, 161, 81), ff := ff + 1L]
  x[ , ff := ff + any(mapq < 20, na.rm = TRUE)*2L, by = qname]
  x[ , ff := ff + (!any(is.na(XA)))*4L, by = qname]
  x[rname != mrnm, ff := ff + 8L]
  
  ## Remove pairs where pos-strand position > neg-strand position
  setorder(x, qname, strand)
  x[ , ff := ff + (pos[1] > pos[2])*16L, by = qname]
  
  ## For the valid pairs, estimate insert size, and remove insert sizes +/-
  ## 5 deviations
  isize_break <- x[strand == "+", 
                   median(isize, na.rm = TRUE) + 5*mad(isize, na.rm = TRUE)]
  x[abs(isize) > isize_break, ff := ff + 32L]
  
  x[]
  
}
