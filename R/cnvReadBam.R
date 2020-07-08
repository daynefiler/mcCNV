##----------------------------------------------------------------------------##
## Utility functions for unpacking Rsamtools::scanBam objects
##----------------------------------------------------------------------------##

#' @name cnvReadBam
#' @title Wrapper to Rsamtools:scanBam to load .bam file as a data.table object
#' 
#' @description \code{cnvReadBam} serves as a wrapper to call Rsamtools:scanBam
#' and processes the list result into a data.table. 
#' 
#' @param bamfile Character of length 1, the file path to the .bam file
#' @param which IRangesList or GRanges object, the range of sequences to load, 
#' passed to [Rsamtools::ScanBamParam] 'which' parameter; loads all sequences
#' when omitted
#' @param what Character, the fields to extract from the bam file, passed to 
#' [Rsamtools::ScanBamParam] 'what' parameter
#' @param tags Character, the tag fields to extract from the bam file, passed to 
#' [Rsamtools::ScanBamParam] 'tags' parameter
#' @param keyby Character, the field to key the data.table by
#' 
#' @return data.table object with specified fields & tags extracted from 
#' .bam file
#' 
#' @importFrom Rsamtools ScanBamParam scanBam bamTag<- scanBamFlag
#' @export 

cnvReadBam <- function(bamfile, 
                       which, 
                       what = scanBamWhat(),
                       flag = scanBamFlag(),
                       tags = NULL, 
                       keyby = NULL) {
  
  if (missing(which)) {
    p <- ScanBamParam(what = what, flag = flag)
  } else {
    p <- ScanBamParam(which = which, what = what, flag = flag)
  }
  if (!is.null(tags)) bamTag(p) <- tags
  s <- scanBam(file = bamfile, index = bamfile, param = p)
  if (!is.null(tags)) s <- .unpackTag(l = s, tags = tags)
  if ("seq" %in% what) s <- .unpackSeq(l = s)
  if ("qual" %in% what) s <- .unpackQual(l = s)
  if (length(s[[1]][[1]]) == 0) {
    return(structure(s[[1]], class = c("data.table", "data.frame")))
  }
  s <- rbindlist(s)
  if (!is.null(keyby)) setkeyv(s, keyby)
  s[]
  
}