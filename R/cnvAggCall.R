##----------------------------------------------------------------------------##
## cnvAggCall: Aggregate CN-state calls from different widths to single call
##----------------------------------------------------------------------------##

#' @name cnvAggCall
#' @title Aggregate CN-state calls from different widths to single call
#' 
#' @param dat the input data, see details
#' @inheritParams cnvCallCN
#' 
#' @details 
#' Need to add
#' 
#' C1 -- TRUE if any interval of the given width that overlaps the exon was
#' called a CNV
#' 
#' @import data.table

cnvAggCall <- function(calls) {
  
  stopifnot(cnvValidCalls(calls))
  
  cols <- paste0("int", seq(max(calls$intWidth)))
  calls[ , (cols) := tstrsplit(intName, split = ";", names = FALSE)]
  calls[ , intName := NULL]
  
  ind <- lapply(cols, function(x) calls[ , which(!is.na(get(x)))])
  calls <- calls[unlist(ind)]
  calls[ , loc := rep(1:length(ind), sapply(ind, length))]
  calls[ , tmp := paste0("int", loc)]
  rm(ind); gc()
  
  for (i in cols) calls[tmp == i, intName := get(i)]
  
  calls[ , (c(cols, "tmp")) := NULL]
  
  calls[ , CNV := CN != 1]
  setindexv(calls, c("subject", "intName"))
  calls[ , anyCNV := any(CNV, na.rm = TRUE), by = list(subject, intName)]
  
  if ("actCN" %in% names(calls)) {
    calls[ , actCN1 := actCN[intWidth == 1], by = .(subject, intName)]
  }
  
  setorder(calls, subject, intName, -CNV, -cnLogP)
  ind <- calls[ , .(ind = .I[1]), by = .(subject, intName)]$ind
  calls[ , CNV := NULL]
  calls <- calls[ind]
  calls[]
  
}