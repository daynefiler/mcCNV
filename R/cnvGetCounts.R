##----------------------------------------------------------------------------##
## cnvGetCounts
##----------------------------------------------------------------------------##

#' @name cnvGetCounts
#' @title Count molecules from a bam file overlapping given intervals
#' 
#' @description \code{cnvGetCounts} takes a bam file and an 
#' [interval object][validObjects] and returns the number of 
#' overlapping molecules for each of the given intervals
#' 
#' @param bamfile Character of length 1, the file path to the .bam file
#' @param interval data.table object containing the intervals to count 
#' molecules over; see \code{?cnvValidInterval}
#' @param subject Character of length 1, the subject name for bam file; 
#' defaults to 'bamfile' parameter when NULL
#' @param verbose TRUE/FALSE
#' 
#' @details
#' 
#' @return data.table [counts object][validObjects]
#' 
#' @import data.table
#' @importFrom IRanges IRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom Rsamtools scanBamHeader scanBamFlag
#' @export 

cnvGetCounts <- function(bamfile, interval, subject = NULL, verbose = TRUE) {
  
  ## Check input parameters
  stopifnot(is.logical(verbose) && length(verbose) == 1)
  if (length(bamfile) > 1) stop("'bamfile' must be of length 1")
  if (!file.exists(bamfile)) stop("Given 'bamfile' does not exist")
  if (is.null(subject)) subject <- bamfile
  if (length(subject) > 1) stop("'subject' must be of length 1")
  if (!cnvValidInterval(interval)) {
    stop("Invalid 'interval'; see ?cnvValidInterval")
  }
  
  ## Get reference intervals from bamfile
  ref <- scanBamHeader(bamfile, what = "targets")[[1]][[1]]
  ref <- ref[as.character(unique(interval$seqnames))]
  ref <- GRanges(names(ref), IRanges(start = 1, width = ref))
  names(ref) <- seqnames(ref)
  ref <- ref[unique(interval$seqnames)]
  
  ## Load reads from bamfile that span the given reference 
  flds <- c("qname", "flag", "strand", "rname", "pos", "mapq", 
            "mrnm", "cigar", "isize")
  tags <- c("XA")
  flgs <- scanBamFlag(isProperPair = TRUE, 
                      isPaired = TRUE, 
                      isDuplicate = FALSE,
                      isSupplementaryAlignment = FALSE,
                      isSecondaryAlignment = FALSE)
  
  cts <- vector(mode = "list", length = length(ref))
  names(cts) <- names(ref)
  
  if (verbose) cat("Processing bamfile: ", bamfile, "\n")
  
  for (r in names(ref)) {
    
    if (verbose) cat("Reading ", r, "...")
    rds <- cnvReadBam(bamfile = bamfile, 
                      what = flds, 
                      tags = tags, 
                      flag = flgs,
                      which = ref[r],
                      keyby = "qname")
    setorder(rds, qname, strand)
    if (verbose) cat("done.\n")
    if (nrow(rds) == 0) {
      cts[[r]] <- data.table(interval[seqnames == r], 
                             molCount = NA_integer_,
                             ncCoverMult = NA_integer_)
      
    } else {
      cts[[r]] <- .procReads(rds, interval[seqnames == r], verbose = verbose)
    }
  }
  
  cts <- rbindlist(cts)
  cts[ , subject := subject]
  cts[]
  
}
