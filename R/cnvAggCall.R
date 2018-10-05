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
#' C2 -- TRUE only over the interval that had the greatest likelihood (this 
#' is a subset of C1, eg look where there were positive calls, then subset
#' down to the most likely interval within the C1 calls)
#' 
#' @import data.table
#' @importFrom tidyr separate 
#' @export

cnvAggCall <- function(dat, width = 5) {
  
  cols <- paste0("p", seq(width))
  dat <- separate(dat, 
                  col = "ref", 
                  into = cols, 
                  sep = ":", 
                  fill = "right", 
                  remove = FALSE)
  gc()
  
  ind <- lapply(cols, function(x) dat[ , which(!is.na(get(x)))])
  dat <- dat[unlist(ind)]
  dat[ , loc := rep(1:length(ind), sapply(ind, length))]
  dat[ , tmp := paste0("p", loc)]
  rm(ind); gc()
  for (i in cols) {
    dat[tmp == i, sngl := get(i)]
  }
  dat[ , (c(cols, "tmp")) := NULL]
  gc()
  dat[ , C1 := any(CN != 1), by = list(sbj, sngl)]
  if (any(grepl("actCN", colnames(dat)))) {
    dat[ ,
         actCNSngl := actCN[width == 1], 
         by = list(sbj, sngl)]
  }
  dat <- dat[!(C1 & CN == 1)]
  ## Initial versions had a weight parameter; when TRUE multiply the likelihood
  ## values by the width. 
  # if (weight)  dat[ , wtlk := lk*width]
  # lkcol <- ifelse(weight, "wtlk", "lk")
  lkcol <- "lk"
  setorderv(dat, c("sbj", "sngl", lkcol), c(1, 1, -1))
  ind <- dat[ , list(ind = .I[1]), by = list(sbj, sngl)]$ind
  dat <- dat[ind]
  rm(ind); gc()
  
  dat[ , refn := as.integer(sub("ref", "", sngl))]
  
  s1 <- dat[which(C1), list(sbj, refn, width, loc)]
  ind <- s1[ , rep(refn, width) + unlist(sapply(width, seq)) - rep(loc, width)]
  s1 <- s1[rep(seq(.N), width)]
  s1[ , ind := ind]
  rm(ind); gc()
  setkey(dat, sbj, refn)
  cols <- c("ref", lkcol, "C1", "loc", "width")
  cnms <- c("ref", lkcol, "C1", "altl", "altw")
  s1[ , (cnms) := dat[s1[ , list(sbj, ind)], .SD, .SDcols = cols]]
  
  s1 <- s1[!is.na(C1)][(C1)]
  setorderv(s1, c("sbj", "refn", lkcol), c(1, 1, -1))
  
  ind <- s1[ , list(ind = .I[1]), by = list(sbj, refn)]$ind
  s2 <- s1[ind]
  rm(ind); gc()
  s2[ , lbnd := ind + 1    - altl]
  s2[ , ubnd := ind + altw - altl]
  s2[ , C2 := refn >= lbnd & refn <= ubnd]
  
  setkey(dat, sbj, refn)
  setkey(s2, sbj, refn)
  
  dat <- s2[ , list(sbj, refn, C2)][dat]
  dat[is.na(C2), C2 := FALSE]
  
  dat[]
  
}