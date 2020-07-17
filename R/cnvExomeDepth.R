#' @title Run ExomeDepth on mcCNV counts object
#' @description Wrapper to run the ExomeDepth algorithm on a 
#' [counts object][validObjects]
#' @param counts data.table [counts object][validObjects]
#' @details 
#' Runs ExomeDepth using the default parameters, then maps the call information 
#' back to the counts object.
#' Will use 'mclapply' to parallelize the computation when available.
#' 
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom ExomeDepth select.reference.set CallCNVs
#' @importClassesFrom ExomeDepth ExomeDepth
#' @export 

cnvExomeDepth <- function(counts) {
  
  cmat <- cnvCountsToMatrix(counts)
  int <- unique(counts[ , .(seqnames, start, end)])
  int[ , intName := sprintf("%s:%d-%d", seqnames, start, end)]
  setkey(int, seqnames, start, end)
  sbjVec <- colnames(cmat)
  
  getRef <- function(sbj) {
    select.reference.set(test.counts = cmat[ , sbj],
                         reference.counts = cmat[ , setdiff(sbjVec, sbj)],
                         bin.length = int$end - int$start + 1,
                         n.bins.reduced = min(1e4, nrow(cmat)))
  }
  
  refList <- mclapply(sbjVec, getRef)
  names(refList) <- sbjVec
  
  calcCN <- function(sbj) {
    cn <- new('ExomeDepth', 
              test = cmat[ , sbj],
              reference = rowSums(cmat[ , refList[[sbj]]$reference.choice]), 
              formula = 'cbind(test, reference) ~ 1')
    cn <- CallCNVs(x = cn, 
                   chromosome = int$seqnames, 
                   start = int$start,
                   end = int$end,
                   name = int$intName)
    cn
  }
  
  cnList <- mclapply(sbjVec, calcCN)
  
  xpndCNV <- function(x, s) {
    d <- as.data.table(x@CNV.calls)
    i1 <- d[ , unlist(mapply(seq, start.p, end.p, SIMPLIFY = FALSE))]
    d <- d[d[ , rep(.I, nexons)]]
    d <- d[ , .(type, nexons, BF, reads.expected, reads.observed)]
    d <- cbind(int[i1], d)
    d[ , subject := s]
    d[]
  }
  
  calls <- mcmapply(xpndCNV, x = cnList, s = sbjVec, SIMPLIFY = FALSE)
  calls <- rbindlist(calls)
  setkey(calls, subject, seqnames, start, end)
  setkey(counts, subject, seqnames, start, end)
  calls <- calls[counts]
  calls[]
  
}

