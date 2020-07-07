#' @rdname validObjects
#' @name validObjects
#' @title Functions to check for appropriate data structures
#' @description For efficiency, the mcCNV package does not currently formalize
#' data structures using the S4 system. Rather, the mcCNV package utilizes
#' [data.table][data.table::data.table] objects, requiring specific fields. 

NULL

#' @rdname validObjects 
#' @section Interval objects:
#' Interval objects define the intervals for computing copy number; inspired
#' by the GenomicRanges package, they require 'seqnames', 'start', and 'end'.
#' The start and end fields must be integers and all end values
#' must be greater than all start values.
#' Finally, 'seqnames' cannot contain the semicolon (';'), colon (';'), or dash
#' ('-') characters. These characters are used in subsequent processing steps
#' to define and combine intervals, e.g. \code{"seq:1-10;seq:15-25"} would
#' represent a combined interval consisting of positions 1-10 and 15-25 on 
#' 'seq'.
#' 
#' Converting a [GRanges][GenomicRanges::GRanges] object with 
#' [as.data.table][data.table::as.data.table] will create a valid interval 
#' object. 
#' 
#' Of note, the \code{mcCNV} package uses 1-based positions (the standard for 
#' R programming) for both the start and end positions. 
#' [fread][data.table::fread] can typically load BED files directly, but the 
#' BED file specification uses 0-based start and 1-based end positions. Users
#' are responsible for ensuring the start and end positions are converted
#' appropriately, e.g. \code{interval[ , start := start + 1]}.
#' @import data.table
#' @export

cnvValidInterval <- function(x) {
  t1 <- is.data.table(x)
  t2 <- try(all(c("seqnames", "start", "end") %in% names(x)), silent = TRUE)
  t3 <- try(is.integer(x$start), silent = TRUE)
  t4 <- try(is.integer(x$end), silent = TRUE)
  t5 <- try(all(x$end > x$start), silent = TRUE)
  t6 <- try(!any(grepl("-|;|:", x$seqnames)))
  t1 && t2 && t3 && t4 && t5
}

#' @rdname validObjects 
#' @section Count objects:
#' Count objects list the molecule counts over the given intervals. Count 
#' objects require the same fields and specifications as interval objects 
#' (listed above). Additionally, count objects have integer fields 'molCount'
#' giving the number of overlapping molecules & 'nCoverMult' giving the number
#' of the overlapping molecules that overlapped more than one interval. 
#' Finally, count objects have a 'subject' field giving the subject represented.
#' Count objects are, by default, long and can be combined using the data.table
#' convention, [rbindlist][data.table::rbindlist()]. Count objects are created
#' by [cnvGetCounts].
#' 
#' We provide the [cnvGatherCounts] convenience function for reading 
#' [saved count objects][base::saveRDS], and combining them into a single 
#' multiple-subject count object required by [cnvCallCN].
#' 
#' @import data.table
#' @export

cnvValidCounts <- function(x) {
  t1 <- cnvValidInterval(x)
  flds <- c("molCount", "nCoverMult", "subject")
  t2 <- try(all(flds %in% names(x)), silent = TRUE)
  t3 <- try(is.integer(x$molCount), silent = TRUE)
  t4 <- try(is.integer(x$nCoverMult), silent = TRUE)
  t1 && t2 && t3 && t4
}
