##----------------------------------------------------------------------------##
## cnvAggCall: Aggregate CN-state calls from different widths to single call
##----------------------------------------------------------------------------##

#' @name cnvAggCall
#' @title Aggregate CN-state calls from different widths to single call
#' 
#' @param dat the input data, see details
#' @param width the maximum window size used to analyze the data
#' 
#' @details 
#' Need to add
#' 
#' @import data.table
#' @export

cnvAggCall <- function(dat, width = 5) {
  
  refs <- strsplit(dat$ref, ":")
  refs <- lapply(seq(width), function(x) sapply(refs, "[", x))
  cols <- paste0("p", seq(width))
  dat[ , (cols) := refs]
  setkeyv(dat, cols = cols)
  rm(refs)
  
  dat[ , wd := width - rowSums(is.na(.SD)), .SDcols = cols]
  
  getSngl <- function(dat, col) {
    dat[ , sngl := get(col)]
    dat[!is.na(sngl)]
  }
  
  dat <- rbindlist(lapply(cols, getSngl, dat = dat))
  dat[ , (cols) := NULL]
  setorder(dat, sngl)
  dat[ , CNCall := any(CN != 1), by = sngl]
  if (any(grepl("actCN", colnames(dat)))) {
    dat[ ,
         actCNSngl := actCN[wd == 1], 
         by = sngl]
  }
  dat <- split(dat, by = "CNCall")
  setorder(dat[["FALSE"]], sngl,  lk)
  dat[["TRUE"]] <- dat[["TRUE"]][CN != 1]
  setorder(dat[["TRUE"]],  sngl, -lk)
  dat <- rbindlist(dat)
  ind <- dat[ , list(ind = .I[1]), by = sngl]
  dat <- dat[ind$ind]
  rm(ind)
  
  dat[]
  
}