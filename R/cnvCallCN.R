##----------------------------------------------------------------------------##
## cnvCallCN: Call CN-states
##----------------------------------------------------------------------------##

#' @name cnvCallCN
#' @title Call copy number-state for the given data
#' 
#' @param cnts the input data, see details
#' @param width integer of length 1, the maximum width (number of consecutive
#' exons) to evaluate simultaneously
#' @param outfile character of length 1, when not NULL, save .RDS file to the
#' file given
#' @inheritParams callCN
#' 
#' @details 
#' Need to add
#' 
#' @import data.table
#' @export

# cnts <- readRDS("~/Desktop/sub/sim_d100_w1_r0001.RDS"); prior <- 0.05; width <- 4; min.dlt <- 20; max.its <- 30
# cnts <- readRDS("d100/w1/sim_d100_w1_r0001.RDS"); prior <- 0.05; width <- 6; min.dlt <- 20; max.its <- 30
# cnts <- readRDS("simData/d050/w1/sim_d050_w1_r0001.RDS"); prior <- 0.0034; width <- 5; min.dlt <- 20; max.its <- 30; shrink = TRUE
cnvCallCN <- function(cnts, prior, width = 5, min.dlt = 20, max.its = 30, 
                      outfile = NULL, agg = FALSE, shrink = TRUE,
                      keep.cols = NULL, verbose = FALSE,
                      return.res = TRUE) {
  
  if (is.data.frame(cnts)) {
    if (!is.data.table(cnts)) cnts <- as.data.table(cnts)
  } else {
    if (verbose) cat("Reading file...")
    cnts <- readRDS(file = cnts)
    if (!is.data.table(cnts)) cnts <- as.data.table(cnts)
    if (verbose) cat("done.\n")
  }
  
  ## Need to add more data checks. 
  
  if (width > 1) {
    widths <- seq(width)
    if (verbose) cat("Collapsing exons...")
    cnts <- rbindlist(lapply(widths, .clpsExon, dat = cnts))
    if (verbose) cat("done.\n")
  } else {
    cnts[ , width := 1]
  }
  
  if (verbose) cat("Calling CNVs...")
  cnts <- .callCN(cnts = cnts, 
                  min.dlt = min.dlt, 
                  max.its = max.its, 
                  prior = prior,
                  shrink = shrink)
  if (verbose) cat("done.\n")
  
  its <- attr(cnts, "its")
  
  if (!is.null(keep.cols)) {
    keep.anyway <- c("sbj", "ref", "CN")
    if (agg) keep.anyway <- c(keep.anyway, "lk", "width", "loc")
    cols <- names(cnts)
    cols <- intersect(c(keep.cols, keep.anyway), cols)
    cnts <- cnts[ , .SD, .SDcols = cols]
  }
  
  if (agg & width > 1) cnts <- cnvAggCall(cnts, width = width)
  
  setattr(cnts, "its", its)
  
  if (!is.null(outfile)) saveRDS(cnts, file = outfile)
  
  if (return.res) return(cnts[])
  
  TRUE
  
}