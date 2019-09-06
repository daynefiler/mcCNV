#' @name cnvGatherCounts
#' @title Gather saved counts into one data.table
#' 
#' @param input Character, directories and/or files to read in 
#' @inheritParams cnvGetCounts
#' @param glob Logical, should the paths for input counts be expanded using 
#' wildcards?
#'  
#' @description This help file needs a lot of work
#' 
#' @import data.table
#' @importFrom utils file_test
#' @export

cnvGatherCounts <- function(input, outfile = NULL, results = TRUE, 
                            glob = FALSE) {
  
  ## Check input parameters
  if (is.null(outfile) & !results) {
    stop("'outfile' is NULL and 'results' is FALSE -- any results would be ",
         "lost at function termination.")
  }
  
  if (glob) input <- Sys.glob(input)
  input <- path.expand(input)
  fls <- input[file_test("-f", input)]
  drs <- input[file_test("-d", input)]
  fls <- c(fls, sapply(drs, list.files, full.names = TRUE))
  fls <- unlist(fls, use.names = FALSE)
  
  cts <- lapply(fls, readRDS)
  cts <- rbindlist(cts)
  
  ## Write ouputs
  if (!is.null(outfile)) {
    saveRDS(cts, file = outfile)
  }
  
  if (results) return(cts[])
  NULL
  
}