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
#' @param int GRanges object containing the intervals (exons) to count reads
#' into
#' @param sbj Character or numeric of length 1, the subject name for bam file; 
#' defaults to the input filename when NULL
#' @param outfile Character of length 1, the file path to write results (.RDS),
#' no results written if NULL
#' @param return TRUE/FALSE, should the results be returned to the console?
#' Should not be FALSE if outfile is NULL
#' @param verbose TRUE/FALSE, if TRUE progress messages will be sent to the 
#' console
#' 
#' @return data.table object with counts for the intervals given in int
#' 
#' @import data.table
#' @importFrom IRanges IRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom Rsamtools scanBamHeader
#' @export 

cnvGetCounts <- function(bamfile, int, sbj = NULL, outfile = NULL, 
                         results = TRUE, verbose = FALSE) {
  
  ## Check input parameters
  if (is.null(outfile) & !results) {
    stop("'outfile' is NULL and 'results' is FALSE -- any results would be ",
         "lost at function termination.")
  }
  
  if (missing(bamfile)) stop("Must provide 'bamfile'; see ?cnvGetCounts.")
  if (missing(int)) stop("Must provide 'int'; see ?cnvGetCounts.")
  
  if (is.character(bamfile) & length(bamfile) == 1) {
    if (!file.exists(bamfile)) {
      stop("Given 'bamfile' does not exist.")
    }
  } else {
    stop("'bamfile' must be a character of length 1 and point to a valid file.")
  }
  
  if (is.null(sbj)) {
    sbj <- sub(".bam$", "", bamfile)
  } else {
    if (length(sbj) > 1) {
      warning("length(sbj) > 1, only using the first element")
      sbj <- sbj[1]
    }
    if (!is.character(sbj)) sbj <- as.character(sbj)
  }
  
  ## Get reference intervals from bamfile 
  ref <- scanBamHeader(bamfile, what = "targets")[[1]][[1]]
  ref <- ref[as.character(unique(seqnames(int)))]
  ref <- GRanges(names(ref), IRanges(start = 1, width = ref))
  
  ## Load reads from bamfile that span the given reference 
  flds <- c("qname", "flag", "strand", "rname", "pos", "mapq", 
            "mrnm", "cigar", "isize")
  tags <- c("XA")
  
  cts <- vector(mode = "list", length = length(ref))
  names(cts) <- seqnames(ref)
  
  for (r in as.character(seqnames(ref))) {
    
    if (verbose) cat("Reading ", r, "...")
    rds <- cnvReadBam(bamfile = bamfile, 
                      what = flds, 
                      tags = tags, 
                      which = GRanges(ref[seqnames(ref) == r]),
                      keyby = "qname")
    setorder(rds, qname, strand)
    if (verbose) cat("done.\n")
    if (nrow(rds) == 0) {
      cts[[r]] <- data.table(ref = names(int[seqnames(int) == r]), 
                             N = 0L,
                             hang = NA_integer_, 
                             mult = NA_integer_, 
                             chr = r)
    } else {
      cts[[r]] <- .procReads(rds, int[seqnames(int) == r], verbose = verbose)
      set(cts[[r]], j = "chr", value = r)
    }
  }
  
  cts <- rbindlist(cts)
  cts[ , sbj := sbj]
  
  ## Write ouputs
  if (!is.null(outfile)) {
    if (verbose) cat("Writing output file...")
    saveRDS(cts, file = outfile)
    if (verbose) cat("done.\n")
  }
  
  if (results) return(cts[])
  
  NULL
  
}
