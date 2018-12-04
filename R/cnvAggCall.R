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
#' @importFrom tidyr separate 
#' @export

cnvAggCall <- function(dat, width = 5) {
  
  cols <- paste0("r", seq(width))
  dat <- separate(dat, 
                  col = "ref", 
                  into = cols, 
                  sep = ";", 
                  fill = "right", 
                  remove = FALSE)
  gc()
  
  ind <- lapply(cols, function(x) dat[ , which(!is.na(get(x)))])
  dat <- dat[unlist(ind)]
  dat[ , loc := rep(1:length(ind), sapply(ind, length))]
  dat[ , tmp := paste0("r", loc)]
  rm(ind); gc()
  
  for (i in cols) dat[tmp == i, sngl := get(i)]
  
  dat[ , (c(cols, "tmp")) := NULL]
  gc()
  dat[ , CNV := CN != 1]
  dat[ , C1 := any(CNV), by = list(sbj, sngl)]
  
  if (any(grepl("actCN", colnames(dat)))) {
    dat[ ,
         actCNSngl := actCN[width == 1], 
         by = list(sbj, sngl)]
  }
  
  
  setorder(dat, sbj, sngl, -CNV, -lp)
  ind <- dat[ , list(ind = .I[1]), by = list(sbj, sngl)]$ind
  dat[ , CNV := NULL]
  dat <- dat[ind]
  rm(ind); gc()
  
  dat[]
  
}