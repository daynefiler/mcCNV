##----------------------------------------------------------------------------##
## cnvCallCN: Call CN-states
##----------------------------------------------------------------------------##

#' @name cnvCallCN
#' @title Call copy number-state for the given data
#' 
#' @param cnts the input data, see details
#' @param prior numeric of length 1, the prior probability of having a CNV
#' @param width integer of length 1, the maximum width (number of consecutive
#' exons) to evaluate simultaneously
#' @param min.dlt integer of length 1, the target number of changes in copy-
#' state to stop the alogorithm 
#' @param max.its integer of length 1, the maximum number of iterations
#' @param outfile character of length 1, when not NULL, save .RDS file to the
#' file given
#' 
#' @details 
#' Need to add
#' 
#' @import data.table
#' @importFrom parallel mclapply
#' @export

# cnts <- readRDS("~/Desktop/sub/sim_d100_w1_r0001.RDS"); prior <- 0.05; width <- 4; min.dlt <- 20; max.its <- 30
# cnts <- readRDS("d100/w1/sim_d100_w1_r0001.RDS"); prior <- 0.05; width <- 6; min.dlt <- 20; max.its <- 30
cnvCallCN <- function(cnts, prior, width = 5, min.dlt = 20, max.its = 30, 
                      outfile = NULL) {
  
  if (is.data.frame(cnts)) {
    if (!is.data.table(cnts)) cnts <- as.data.table(cnts)
  } else {
    cnts <- readRDS(file = cnts)
    if (!is.data.table(cnts)) cnts <- as.data.table(cnts)
  }
  
  ## Need to add more data checks. 
  
  widths <- seq(width)
  cnts <- mclapply(widths, .clpsExon, dat = cnts, mc.cores = width)
  cnts <- mclapply(cnts, .callCN, 
                   min.dlt = min.dlt, 
                   max.its = max.its, 
                   prior = prior, 
                   mc.cores = width)
  its <- sapply(cnts, attr, which = "its")
  cnts <- rbindlist(cnts)
  setattr(cnts, "its", its)
  
  if (is.null(outfile)) {
    return(cnts[])
  } else {
    cols <- c("ref", "sbj", "N", "CN", "lk")
    if (!is.null(cnts$actCN)) cols <- c(cols, "actCN")
    cnts <- cnts[ , .SD, .SDcols = cols]
    setattr(cnts, "its", its)
    saveRDS(cnts, file = outfile)
    return(TRUE)
  }
  
}