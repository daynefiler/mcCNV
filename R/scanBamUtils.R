##----------------------------------------------------------------------------##
## Utility functions for unpacking Rsamtools::scanBam objects
##----------------------------------------------------------------------------##

#' @name unpackTag 
#' @title Unpack the tag values from object returned by Rsamtools::scanBam
#' 
#' @description \code{.unpackTag} simply extracts the tag values from the list
#' returned by Rsamtools:scanBam. This allows for collapsing the list into
#' a data.table object. 
#' 
#' @param l List object returned by Rsamtools::scanBam
#' @param tags Character, the tag names to unpack
#' 
#' @return List object that resembles the output from Rsamtools:scanBam

.unpackTag <- function(l, tags) {
  
  tmp <- function(x) {
    x[tags] <- x$tag[tags]
    n <- length(x[[1]])
    ind <- sapply(x[tags], is.null)
    if (n > 0 & any(ind)) x[tags[ind]] <- list(rep(NA, n))
    x$tag <- NULL
    x
  }
  l <- lapply(l, tmp)
  l
}

#' @name unpackSeq
#' @title Unpack the seq values from object returned by Rsamtools::scanBam
#' 
#' @description \code{.unpackSeq} converts the seq values to simple character 
#' strings from the list returned by Rsamtools:scanBam. This allows for 
#' collapsing the list into a data.table object. 
#' 
#' @param l List object returned by Rsamtools::scanBam
#' 
#' @return List object that resembles the output from Rsamtools:scanBam

.unpackSeq <- function(l) {
  
  tmp <- function(x) {
    x$seq <- as.character(l$seq)
    x
  }
  
  l <- lapply(l, tmp)
  l
  
}

#' @name padGRanges
#' @title Add padding to GRanges object
#' 
#' @description \code{cnvReadBam} serves as a wrapper to call Rsamtools:scanBam
#' and processes the list result into a data.table. 
#' 
#' @param g GRanges object
#' @param pad Integer of length 1, the amount of padding to add
#' @param keyby Character, the field to key the data.table by
#' 
#' @return GRanges object with padding added
#' 
#' @importFrom GenomicRanges start end

.padGRanges <- function(g, pad) {
  if (length(pad) > 1 )
  pad <- abs(pad)
  start(g) <- start(g) - pad
  end(g) <- end(g) + pad
  g
}

#' @name cnvReadBam
#' @title Wrapper to Rsamtools:scanBam to load .bam file as a data.table object
#' 
#' @description \code{cnvReadBam} serves as a wrapper to call Rsamtools:scanBam
#' and processes the list result into a data.table. 
#' 
#' @param bamfile Character of length 1, the file path to the .bam file
#' @param which IRanges object, the range of sequences to load, passed to 
#' Rsamtools:ScanBamParam 'which' parameter
#' @param what Character, the fields to extract from the bam file, passed to 
#' Rsamtools:ScanBamParam 'what' parameter
#' @param tags Character, the tag fields to extract from the bam file, passed to 
#' Rsamtools:ScanBamParam 'tags' parameter
#' @param keyby Character, the field to key the data.table by
#' 
#' @return data.table object with specified fields & tags extracted from 
#' .bam file
#' 
#' @importFrom Rsamtools ScanBamParam scanBam
#' @export 

cnvReadBam <- function(bamfile, which, what, tags = NULL, keyby = NULL) {
  
  p <- ScanBamParam(which = which, what = what, tag = tags)
  s <- scanBam(file = bamfile, index = bamfile, param = p)
  if (!is.null(tags)) s <- .unpackTag(l = s, tags = tags)
  if ("seq" %in% what) s <- .unpackSeq(l = s)
  if (length(s[[1]][[1]]) == 0) {
    return(structure(s[[1]], class = c("data.table", "data.frame")))
  }
  s <- rbindlist(s)
  if (!is.null(keyby)) setkeyv(s, keyby)
  s[]
  
}


#' @name cnvReadRefs
#' @title Wrapper to Rsamtools:scanBamHeader to extract reference names & 
#' intervals from bam files. 
#' 
#' @description \code{cnvReadRefs} serves as a wrapper to call 
#' Rsamtools:scanBamHeader and return IRanges object with reference intervals. 
#' 
#' @param bamfile Character of length 1, the file path to the .bam file
#' 
#' @return IRanges object with reference intervals
#' 
#' @importFrom Rsamtools scanBamHeader
#' @importFrom IRanges IRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges GRanges
#' @export 

cnvReadRefs <- function(bamfile) {
  
  h <- scanBamHeader(bamfile, what = "targets")[[1]]
  r <- IRanges(start = 1, width = h$targets)
  GRanges(seqnames = names(h$targets), ranges = r)
  
}


#' @name cnvInt2GRanges
#' @title Read in GATK-style .interval file to GRanges. 
#' 
#' @description \code{cnvInt2GRanges} reads a GATK-style .interval file and 
#' returns a GRanges object for the intervals in the file. 
#' 
#' @param intfile Character of length 1, the file path to the .interval file
#' @param pad Integer of length 1, the number of basepairs to pad around the 
#' bed file intervals
#' 
#' @note 
#' 
#' @return GRanges object with reference intervals
#' 
#' @import data.table
#' @importFrom GenomicRanges GRanges 
#' @importFrom GenomeInfoDb seqnames
#' @importClassesFrom GenomicRanges GRanges
#' @export 

cnvInt2GRange <- function(intfile, pad = NULL, skip = 0) {
  int <- fread(intfile, header = FALSE, sep = "\t", skip = skip)[[1]]
  int <- GRanges(int)
  if (!is.null(pad)) int <- .padGRanges(int)
  names(int) <- paste(seqnames(int), ranges(int), sep = ":")
  int
}

#' @name cnvBed2GRanges
#' @title Read in .bed file to GRanges. 
#' 
#' @description \code{cnvInt2GRanges} reads a GATK-style .interval file and 
#' returns a GRanges object for the intervals in the file. 
#' 
#' @param intfile Character of length 1, the file path to the .interval file
#' @param pad Integer of length 1, the number of basepairs to pad around the 
#' bed file intervals
#' @param skip Integer of length 1, the number of lines to skip when reading
#' the BED file, passed to [data.table::fread]. 
#' @param adj Logical of length 1, when TRUE adjust the start positions by 
#' adding 1 to convert 0-based BED coordinates to 1-based
#' 
#' @note BED files have 0-based coordinates, and this function defaults to 
#' adding 1 to all start positions. See this 
#' [Biostars post](https://www.biostars.org/p/84686/) for more information.
#' 
#' @return GRanges object with reference intervals
#' 
#' @importFrom IRanges IRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges ranges
#' @importFrom GenomeInfoDb seqnames
#' @importClassesFrom GenomicRanges GRanges
#' @export 

cnvBed2GRange <- function(bedfile, pad = NULL, skip = 0, adj = TRUE) {
  int <- fread(bedfile, skip = skip)
  if (ncol(int) > 3) {
    nms <- int[[4]]
    hasNames <- TRUE
  } else {
    hasNames <- FALSE
  }
  int <- GRanges(seqnames = int[[1]], 
                 IRanges(start = int[[2]] + 1, end = int[[3]]))
  if (hasNames) {
    names(int) <- nms
  } else {
    names(int) <- paste(seqnames(int), ranges(int), sep = ":")
  }
  if (!is.null(pad)) int <- .padGRanges(int)
  int
}


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
  
  rds[!flag %in% c(99, 147, 163, 83, 97, 145, 161, 81), ff := ff + 1L]
  rds[ , ff := ff + any(mapq < 20, na.rm = TRUE)*2L, by = qname]
  rds[ , ff := ff + (!any(is.na(XA)))*4L, by = qname]
  rds[rname != mrnm, ff := ff + 8L]
  rds[!bitwAnd(ff, 1) & ! bitwAnd(ff, 8), 
      ff := ff + (pos[1] > pos[2])*16L, by = qname]
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
  
  rds <- with(mls, IRanges(start = q1, end = q4, names = qname))
  elementMetadata(rds) <- mls[ , list(q2, q3)]
  rm(mls); gc()
  if (verbose) cat("done.\n")
  
  ## Find overlaps
  if (verbose) cat("Getting overlaps and collapsing by interval...")
  prs <- findOverlapPairs(ranges(int), rds)
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
              list(N = sum(wt == 1), 
                   hang = length(which(lf != 3)), 
                   mult = length(which(wt < 1))),
              by = ref]
  cts[mult > 0, N := N + sapply(mult, rbinom, n = 1, prob = 0.5)]
  setkey(cts, ref)
  cts <- cts[names(int)]
  cts[is.na(N), N:= 0]
  if (verbose) cat("done.\n")
  
  cts[]
  
}

# .parseXa <- function(x) {
#   
#   proto <- data.frame(rname = character(), 
#                       strand = character(), 
#                       pos = numeric(), 
#                       cigar = character(), 
#                       nm = numeric())
#   pattern <- "(.*?),([[:punct:]])([[:digit:]]+),([[:alnum:]]+),([[:digit:]]+)"
#   
#   x <- strsplit(x, ";")
#   lapply(x, strcapture, pattern = pattern, proto = proto)
#   
# }