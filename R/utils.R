##----------------------------------------------------------------------------##
## Utility functions for unpacking Rsamtools::scanBam objects
##----------------------------------------------------------------------------##

.unpackTag <- function(l, tags) {
  
  # @description \code{.unpackTag} simply extracts the tag values from the list
  # returned by Rsamtools:scanBam. This allows for collapsing the list into
  # a data.table object. 
  # 
  # @param l List object returned by Rsamtools::scanBam
  # @param tags Character, the tag names to unpack
  
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

.unpackSeq <- function(l) {
  # @description \code{.unpackSeq} converts the seq values to simple character 
  # strings from the list returned by Rsamtools:scanBam. This allows for 
  # collapsing the list into a data.table object. 
  # 
  # @param l List object returned by Rsamtools::scanBam
  lapply(l, function(x) {x$seq <- as.character(x$seq); x})
}

.unpackQual <- function(l) {
  # @description \code{.unpackQual} converts the qual values to simple character 
  # strings from the list returned by Rsamtools:scanBam. This allows for 
  # collapsing the list into a data.table object. 
  # 
  # @param l List object returned by Rsamtools::scanBam
  lapply(l, function(x) {x$qual <- as.character(x$qual); x})
}

##----------------------------------------------------------------------------##
## .procReads: process reads from bam file, used in cnvGetCounts
##----------------------------------------------------------------------------##

#' @import data.table
#' @importFrom stats rbinom mad median

.procReads <- function(rds, int, verbose) {
  
  # @param rds data.table with reads 
  # @param int GRanges object with refs
  # @param verbose verbose
  
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

##----------------------------------------------------------------------------##
## .calcShrPhi: shrink the phi values
##----------------------------------------------------------------------------##

.calcShrPhi <- function(phi) {
  
  if (any(is.na(phi))) stop("Error in .calcShrPhi due to missing values.")
  nPhi <- length(phi)
  xi <- seq(0, 1, length.out = 1000)
  asd <- sapply(xi, function(x) 1/sum((phi - x)^2))
  xi <- xi[which.min(diff(asd)/diff(xi))]
  dlt <- (sum((phi - mean(phi))^2)/(nPhi - 1))/(sum((phi - xi)^2)/(nPhi - 2))
  sPhi <- (1 - dlt)*phi + dlt*xi
  sPhi
  
}

##----------------------------------------------------------------------------##
## .clpsExon: Collapse exons by specified width 
##----------------------------------------------------------------------------##

#' @import data.table

.clpsExon <- function(dat, cw) {
  
  # @param dat data.table count object, plus the intName column
  # @param cw integer of length 1, the width to collapse by ('collapse width')
  
  ## Set the order based on the chromosome & start position, or the ref 
  ## column if chr and spos not given (helpful for simulated data)
  setorder(dat, seqnames, start, end)
  
  aw <- -1*(cw - 1)
  shiftby <- c("subject", "seqnames")
  scols <- intersect(c("intName", "molCount", "actCN"), names(dat))
  dat <- dat[ , 
             shift(x = .SD, n = 0:aw, type = "shift", give.names = TRUE), 
             by = shiftby,
             .SDcols = scols]
  countCols <- grep("^molCount_shift", names(dat), value = TRUE)
  dat[ , molCount := rowSums(.SD), .SDcols = countCols]
  dat[ , (countCols) := NULL]
  intNmCols <- grep("^intName_shift", names(dat), value = TRUE)
  dat[ , intName := do.call(paste, c(.SD, sep = ";")), .SDcols = intNmCols]
  dat[ , (intNmCols) := NULL]
  if ("actCN" %in% scols) {
    cnCols <- grep("^actCN_shift", names(dat), value = TRUE)
    dat[ , actCN := do.call(paste, c(.SD, sep = ";")), .SDcols = cnCols]
    dat[ , (cnCols) := NULL]
  }
  dat <- dat[!is.na(molCount)]
  dat[ , molCount := as.integer(molCount)]
  dat[ , width := as.integer(cw)]
  dat[]
  
}

##----------------------------------------------------------------------------##
## .callCN: Call CNS
##----------------------------------------------------------------------------##

#' @importFrom Rfast rowMaxs
#' @importFrom matrixStats rowLogSumExps
#' @importFrom stats dnbinom
#' @import data.table

.callCN <- function(cnts, delta, iterations, prior, shrink = TRUE) {
  
  # @param cnts data.table object containing the counts
  # @param prior numeric of length 1, the prior probability of having a CNV
  # @param delta integer of length 1, the target number of changes in copy-
  # state to stop the alogorithm 
  # @param iterations integer of length 1, the maximum number of iterations
  # @param shrink logical of length 1, shrinkage applied to phi when TRUE
  
  cnts[ , CN := 1]
  
  cs <- c(0.001, 0.5, 1, 1.5, 2.0, 2.5)
  nstates <- length(cs)
  if (prior > 1/nstates || prior < 0) {
    if (prior < 0) {
      stop("Invalid prior; must be greater than 0 and less than 1/", nstates)
    }
    prior <- 1/nstates
    warning("Prior > 1/", nstates, "; using 1/", nstates, ".")
  }
  cp <- rep(prior, nstates); cp[which(cs == 1)] <- 1 - prior*(nstates - 1)
  chngVec <- rep(NA_integer_, iterations)
  it <- 1
  repeat {
    
    cnts[ , adjN := molCount/CN]
    cnts[ , use := CN > 0.001 & adjN > 10]
    cnts[ , geomn := exp(mean(log(adjN[use]))), by = intName]
    cnts[ , sf := median(adjN/geomn, na.rm = TRUE), by = list(subject, width)]
    cnts[ , mn := mean(adjN[use]/sf[use]), by = intName]
    cnts[ , vr :=  var(adjN[use]/sf[use]), by = intName]
    
    shrPhi <- cnts[ , 
                    list(n = .N, mn = mn[1], vr = vr[1], isf = sum(1/sf)), 
                    by = list(intName, width)]
    # shrPhi[is.na(vr), vr := mn]
    shrPhi[ , phi := (n*vr + mn*isf)/(mn^2*isf), by = intName]
    shrPhi[(phi < 0), phi := 0]
    if (shrink) shrPhi[!is.na(phi), phi := .calcShrPhi(phi), by = width]
    setkey(cnts, intName)
    setkey(shrPhi, intName)
    cnts <- shrPhi[ , list(intName, phi)][cnts]
    setkey(cnts, subject, intName)
    cnts[ , oldCN := CN]
    calcProb <- function(x) {
      cnts[ , dnbinom(molCount, sf/phi, 1/(mn*phi*x + 1), log = TRUE)]
    }
    probMat <- do.call(cbind, lapply(cs, calcProb))
    probMat <- sweep(probMat, 2, log(cp), "+")
    psum <- rowLogSumExps(probMat, na.rm = TRUE)
    cnts[ , CN := cs[rowMaxs(probMat)]]
    cnts[ , llk := rowMaxs(probMat, value = TRUE)]
    cnts[ , lp  := llk - psum]
    cnts[ , lp1 := probMat[ , cs == 1] - psum]
    rm(probMat); gc()
    nchng <- cnts[oldCN != CN, .N]
    chngVec[it] <- nchng
    
    if (nchng < delta | it == iterations) break
    
    it <- it + 1
    cnts[ , phi := NULL]
    
  }
  
  cnts[ , oldCN := NULL]
  setattr(cnts, "delta", chngVec)
  
}
