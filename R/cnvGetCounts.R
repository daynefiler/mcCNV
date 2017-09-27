##----------------------------------------------------------------------------##
## cnvGetCounts: given bam and int files, return molecule counts per interval
##----------------------------------------------------------------------------##

#' @name cnvGetCounts
#' @title Given bam and int files, return molecule counts per interval
#' 
#' @description \code{cnvGetCounts} takes a bam file and an interval file and
#' returns the overlapping molecules for each interval in the interval file 
#' that span a given reference sequence. 
#' 
#' @param bamfile Character of length 1, the file path to the .bam file
#' @param intfile Character of length 1, the file path ot the interval file
#' @param refname Character of length 1, the name of the reference sequence
#' to evaluate
#' @param outfile Character of length 1, the file path to write results (.csv),
#' no results written if NULL
#' @param return TRUE/FALSE, should the results be returned to the console?
#' Should not be FALSE if outfile is NULL
#' @param verbose TRUE/FALSE, if TRUE progress messages will be sent to the 
#' console
#' 
#' @return IRanges object with reference intervals
#' 
#' @import utils
#' @import Rsamtools 
#' @import data.table
#' @import IRanges 
#' @import GenomicRanges
#' @export 


cnvGetCounts <- function(bamfile, intfile, refname, 
                         outfile = NULL, results = TRUE, verbose = FALSE) {
  
  ## Check input parameters
  if (is.null(outfile) & !results) {
    stop("'outfile' is NULL and 'results' is FALSE -- any results would be ",
         "lost at function termination.")
  }
  
  if (missing(refname)) stop("Must provide 'refname'; see ?cnvGetCounts.")
  if (missing(bamfile)) stop("Must provide 'bamfile'; see ?cnvGetCounts.")
  if (missing(intfile)) stop("Must provide 'intfile'; see ?cnvGetCounts.")
  
  if (is.character(bamfile) & length(bamfile) == 1) {
    if (!file.exists(bamfile)) {
      stop("Given 'bamfile' does not exist.")
    }
  } else {
    stop("'bamfile' must be a character of length 1 and point to a valid file.")
  }
  
  if (is.character(refname) & length(refname) == 1) {
    read_header <- TRUE
  } else {
    if (is(refname, "IRanges") & length(refname) == 1) {
      read_header <- FALSE
    } else {
      stop("'refname' must be valid character or IRanges object of length 1.")
    }
  }
  
  if (is.character(intfile)) {
    read_intfile <- TRUE
    if (length(intfile) > 1) stop("'intfile' file path must be length 1.")
    if (!file.exists(intfile)) stop("Given 'intfile' does not exist.")
  } else {
    if (is(intfile, "IRanges")) {
      read_intfile <- FALSE
    } else {
      stop("'intfile' must be valid file or IRanges object.")
    }
  }
  
  ## Load and format the interval file 
  if (read_intfile) {
    if (verbose) cat("Reading interval file..")
    intStr <- fread(intfile, header = FALSE, sep = "\t", col.names = "str")$str
    intTbl <- strcapture(pattern = "(.*?):([[:digit:]]+)-([[:digit:]]+)",
                         x = intStr,
                         proto = data.frame(chr = character(), 
                                            start = integer(),
                                            end = integer()))
    intTbl$str <- intStr
    intRng <- with(intTbl[intTbl$chr == refname, ], 
                   IRanges(start = start, end = end, names = str))
    if (verbose) cat("done.\n")
  }
  
  ## Get reference intervals from bamfile 
  if (read_header) {
    if (verbose) cat("Reading bam header...")
    refList <- cnvReadRefs(bamfile = bamfile)
    if (!refname %in% names(refList)) {
      stop("Given 'refname' not in bam header.")
    }
    ref <- refList[refname]
    if (verbose) cat("done.\n")
  } else {
    ref <- refname
  }
  
  ## Load reads from bamfile that span the given reference 
  flds <- c("qname", "flag", "strand", "rname", "pos", "mapq", 
            "mrnm", "cigar", "isize")
  tags <- c("XA")
  if (verbose) cat("Reading bam file...")
  rds <- cnvReadBam(bamfile = bamfile, 
                    what = flds, 
                    tags = tags, 
                    which = ref,
                    keyby = "qname")
  setorder(rds, qname, strand)
  if (verbose) cat("done.\n")
  
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
  
  rds[!flag %in% c(99, 147, 163, 83, 97, 145, 161, 81), ff := ff + 1L]
  rds[ , ff := ff + any(mapq < 20, na.rm = TRUE)*2L, by = qname]
  rds[ , ff := ff + (!any(is.na(XA)))*4L, by = qname]
  rds[rname != mrnm, ff := ff + 8L]
  rds[ , ff := ff + (pos[1] > pos[2])*16L, by = qname]
  isize_break <- rds[strand == "+", 
                   median(isize, na.rm = TRUE) + 5*mad(isize, na.rm = TRUE)]
  rds[abs(isize) > isize_break, ff := ff + 32L]
  if (verbose) cat("done.\n")
  
  ## Define molecules
  if (verbose) cat("Defining molecules...")
  rds[ , pos2 := pos + cigar2rlen(cigar) - 1]
  mls <- rds[!as.logical(ff), 
             list(q1 = pos[1], q2 = pos2[1], q3 = pos[2], q4 = pos2[2]), 
             by = qname]
  
  molRng <- with(mls, IRanges(start = q1, end = q4, names = qname))
  elementMetadata(molRng) <- mls[ , list(q2, q3)]
  if (verbose) cat("done.\n")
  
  ## Find overlaps
  if (verbose) cat("Getting overlaps and collapsing by interval...")
  prs <- findOverlapPairs(intRng, molRng)
  keep <- c("first.names", "first.start", "first.end", 
            "second.X.names", "second.q2", "second.q3")
  prs <- as.data.table(prs)[ , .SD,.SDcols = keep]
  setnames(prs, c("ref", "r1", "r2", "mol", "q2", "q3"))
  
  ## lf ('location flag') bit-key
  ## 1 - pos-strand mol overlaps
  ## 2 - neg-strand mol overlaps
  prs[ , lf := 0L]
  prs[q2 > r1, lf := lf + 1L]
  prs[q3 < r2, lf := lf + 2L]
  
  prs[ , wt := 1/.N, by = mol]
  cts <- prs[ , 
              list(N = sum(wt), 
                   hang = length(which(lf != 3)), 
                   mult = length(which(wt < 1))),
              by = ref]
  setkey(cts, ref)
  cts <- cts[names(intRng)]
  cts[is.na(N), N:= 0]
  if (verbose) cat("done.\n")
  
  ## Write ouputs
  if (!is.null(outfile)) {
    if (verbose) cat("Writing output file...")
    fwrite(cts, file = outfile)
    if (verbose) cat("done.\n")
  }
  
  if (results) return(cts[])
  
  NULL
  
}
