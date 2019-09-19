#' @name cnvGatherCounts
#' @title Gather saved counts into one data.table
#' 
#' @param input Character, directories and/or files to read in 
#' @inheritParams cnvGetCounts
#' @param glob Logical, should the paths for input counts be expanded using 
#' wildcards?
#' @param ... arguments passed to \code{\link[parallel]{mclapply}}
#' @param cast Logical of length 1, should the data be cast from long to wide 
#' (i.e. from a single count column to 1 sample per column?)
#'  
#' @description This help file needs a lot of work
#' 
#' @import data.table
#' @importFrom utils file_test
#' @importFrom parallel mclapply
#' @export

cnvGatherCounts <- function(input, outfile = NULL, results = TRUE, 
                            glob = FALSE, cast = FALSE,...) {
  
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
  
  cts <- mclapply(fls, readRDS, ...)
  cts <- rbindlist(cts)
  
  if (cast) {
    cts <- dcast(cts, ref + chr ~ sbj, value.var = "N")
  }
  
  ## Write ouputs
  if (!is.null(outfile)) {
    saveRDS(cts, file = outfile)
  }
  
  if (results) return(cts[])
  NULL
  
}