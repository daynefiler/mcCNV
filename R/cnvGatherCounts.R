##----------------------------------------------------------------------------##
## cnvGatherCounts
##----------------------------------------------------------------------------##

#' @name cnvGatherCounts
#' @title Gather saved count objects into one data.table
#'  
#' @description \code{cnvGatherCounts} reads [count objects][validObjects] 
#' stored in disk space, and combines them into a single multiple-sample count
#' object.
#' 
#' @param files Character, file paths pointing to count objects; see details
#' @param glob Logical, should the given 'files' parameter be passed to 
#' [Sys.glob()][base::Sys.glob()]? Defaults to TRUE
#' @param ... arguments passed to [mclapply][parallel::mclapply]
#' 
#' @details
#' The 'files' parameter must point to valid [count objects][validObjects] 
#' saved to disk using the [saveRDS()][base::saveRDS] function. When 'glob' is
#' \code{TRUE}, the default, any wildcards in 'files' will be expanded.
#' 
#' @import data.table
#' @importFrom utils file_test
#' @importFrom parallel mclapply
#' @export

cnvGatherCounts <- function(files, glob = TRUE, ...) {
  
  stopifnot(is.character(files))
  stopifnot(is.logical(glob) && length(glob) == 1)
  if (glob) files <- Sys.glob(files)
  filesExist <- sapply(files, file.exists)
  if (!all(filesExist)) stop("At least one file in 'files' does not exist.")
  cts <- mclapply(files, readRDS, ...)
  valid <- sapply(cts, cnvValidCounts)
  vMsg <- "At least one file in 'files' does not store a valid counts object."
  if (any(!valid)) stop(vMsg)
  cts <- rbindlist(cts)
  cts[]
  
}