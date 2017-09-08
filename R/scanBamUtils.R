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
#' @import Rsamtools 
#' @export 

cnvReadBam <- function(bamfile, which, what, tags = NULL, keyby = NULL) {
  
  p <- ScanBamParam(which = which, what = what, tag = tags)
  s <- scanBam(file = bamfile, index = bamfile, param = p)
  if (!is.null(tags)) s <- .unpackTag(l = s, tags = tags)
  if ("seq" %in% what) s <- .unpackSeq(l = s)
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
#' @import Rsamtools 
#' @export 

cnvReadRefs <- function(bamfile) {
  
  h <- scanBamHeader(bamfile)[[1]]
  refRng <- IRanges(start = 0, width = h$targets, names = names(h$targets))
  refRng
  
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