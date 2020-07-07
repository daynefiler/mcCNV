##----------------------------------------------------------------------------##
## cnvGetCounts & .procReads helper function
##----------------------------------------------------------------------------##

#' @name procReads
#' @title Process reads from a bam file 
#' 
#' @param rds data.table with reads 
#' @param int GRanges object with refs
#' @param verbose verbose
#' 
#' @description This help file needs a lot of work
#' 
#' @importFrom IRanges IRanges findOverlapPairs
#' @importClassesFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges ranges
#' @importFrom GenomeInfoDb seqnames
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom S4Vectors elementMetadata elementMetadata<-
#' @importFrom stats rbinom

.procReads <- function(rds, int, verbose) {
  
  ## Initialize column to store the filter flags
  if (verbose) cat("Filtering reads...")
  rds[ , ff := 0L]
  
  ## ff ('filter flag') bit-key:
  ## 1  - read not paired in the correct orientation
  ## 2  - read mapping quality <20 for either read
  ## 4  - both read mapping positions ambiguous
  ## 8  - pairs mapped to different references 
  ## 16 - pos-strand position > neg-strand position 
  ## 32 - insert size > median(isize) + 5*mad(isize)
  
  rds[!flag %in% c(99, 147, 83, 163), ff := ff + 1L]
  rds[ , ff := ff + any(mapq < 20, na.rm = TRUE)*2L, by = qname]
  rds[ , ff := ff + (!any(is.na(XA)))*4L, by = qname]
  rds[rname != mrnm, ff := ff + 8L]
  rds[!bitwAnd(ff, 1) & !bitwAnd(ff, 8), 
      ff := ff + (pos[1] > pos[2])*16L, by = qname]
  isize_break <- rds[strand == "+",
                     median(isize, na.rm = TRUE) + 5*mad(isize, na.rm = TRUE)]
  rds[abs(isize) > isize_break, ff := ff + 32L]
  if (verbose) cat("done.\n")
  
  ## Define molecules
  if (verbose) cat("Defining molecules...")
  mls <- rds[isize > 0 & ff == 0, 
             .(seqnames = rname, start = pos, end = pos + isize, qname = qname)]
  if (verbose) cat("done.\n")
  
  ## Find overlaps
  if (verbose) cat("Counting overlaps...")
  setkey(mls, seqnames, start, end)
  setkey(int, seqnames, start, end)
  ovlp <- foverlaps(mls, int, which = TRUE, nomatch = NULL)
  ovlp[ , w := 1/.N, by = xid]
  cts <- ovlp[ , .(molCount = sum(w == 1), nCoverMult = sum(w != 1)), by = yid]
  cts[ , add := 0L]
  cts[nCoverMult > 0, add := sapply(nCoverMult, rbinom, n = 1, prob = 0.5)]
  cts[!is.na(add), molCount := molCount + add]
  with(cts, int[yid, c("molCount", "nCoverMult") := .(molCount, nCoverMult)])
  int[is.na(molCount), c("molCount", "nCoverMult") := 0L]
  cat("done.\n")
  int[]
  
}

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
  if (missing(bamfile)) stop("Must provide 'bamfile'; see ?cnvGetCounts")
  if (missing(interval)) stop("Must provide 'interval'; see ?cnvGetCounts")
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
