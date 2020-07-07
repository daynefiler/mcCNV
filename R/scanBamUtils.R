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
  lapply(l, tmp)
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
  lapply(l, function(x) {x$seq <- as.character(x$seq); x})
}

#' @name unpackQual
#' @title Unpack the qual values from object returned by Rsamtools::scanBam
#' 
#' @description \code{.unpackQual} converts the qual values to simple character 
#' strings from the list returned by Rsamtools:scanBam. This allows for 
#' collapsing the list into a data.table object. 
#' 
#' @param l List object returned by Rsamtools::scanBam
#' 
#' @return List object that resembles the output from Rsamtools:scanBam

.unpackQual <- function(l) {
  lapply(l, function(x) {x$qual <- as.character(x$qual); x})
}

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

