#' @name cnvGatherCounts
#' @title Gather saved counts into one data.table
#' 
#' @param input Character, directories and/or files to read in 
#' @inheritParams cnvGetCounts
#'  
#' @description This help file needs a lot of work
#' 
#' @import data.table
#' @importFrom utils file_test
#' @export

cnvGatherCounts <- function(input, outfile = NULL, results = TRUE) {
  
  ## Check input parameters
  if (is.null(outfile) & !results) {
    stop("'outfile' is NULL and 'results' is FALSE -- any results would be ",
         "lost at function termination.")
  }
  
  input <- path.expand(input)
  fls <- input[file_test("-f", input)]
  drs <- input[file_test("-d", input)]
  fls <- c(fls, sapply(drs, list.files, full.names = TRUE))
  
  cts <- lapply(fls, readRDS)
  cts <- rbindlist(cts)
  
  ## Write ouputs
  if (!is.null(outfile)) {
    saveRDS(cts, file = outfile)
  }
  
  if (results) return(cts[])
  NULL
  
}