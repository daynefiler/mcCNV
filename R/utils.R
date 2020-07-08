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
  
  # @param dat data.table object containing the counts
  # @param cw integer of length 1, the width to collapse by
  
  ## Set the order based on the chromosome & start position, or the ref 
  ## column if chr and spos not given (helpful for simulated data)
  setorderv(dat, intersect(c("chr", "spos", "ref"), names(dat)))
  
  aw <- -1*(cw - 1)
  shiftby <- if ("chr" %in% names(dat)) c("sbj", "chr") else "sbj"
  scols <- c("ref", "N")
  if ("actCN" %in% names(dat)) scols <- c(scols, "actCN")
  dat <- dat[ , 
              shift(x = .SD, n = 0:aw, type = "shift", give.names = TRUE), 
              by = shiftby,
              .SDcols = scols]
  dat[ , N := rowSums(.SD), .SDcols = grep("^N_shift", colnames(dat))]
  dat[ , 
       ref := do.call(paste, c(.SD, sep = ";")), 
       .SDcols = grep("ref", colnames(dat))]
  if (any(grepl("actCN", colnames(dat)))) {
    dat[ ,
         actCN := do.call(paste, c(.SD, sep = ":")), 
         .SDcols = grep("actCN", colnames(dat))]
  }
  dat[ , grep("_shift_", colnames(dat), value = TRUE) := NULL]
  dat <- dat[!is.na(N)]
  dat[ , width := cw]
  dat[]
  
}

##----------------------------------------------------------------------------##
## .callCN: Call CNS
##----------------------------------------------------------------------------##

#' @importFrom Rfast rowMaxs
#' @importFrom matrixStats rowLogSumExps
#' @importFrom stats dnbinom
#' @import data.table

.callCN <- function(cnts, min.dlt, max.its, prior, shrink = TRUE) {
  
  # @param cnts data.table object containing the counts
  # @param prior numeric of length 1, the prior probability of having a CNV
  # @param min.dlt integer of length 1, the target number of changes in copy-
  # state to stop the alogorithm 
  # @param max.its integer of length 1, the maximum number of iterations
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
  
  it <- 1
  repeat {
    
    cnts[ , adjN := N/CN]
    cnts[ , use := CN > 0.001 & adjN > 10]
    cnts[ , geomn := exp(mean(log(adjN[use]))), by = ref]
    cnts[ , sf := median(adjN/geomn, na.rm = TRUE), by = list(sbj, width)]
    cnts[ , mn := mean(adjN[use]/sf[use]), by = ref]
    cnts[ , vr :=  var(adjN[use]/sf[use]), by = ref]
    
    shrPhi <- cnts[ , 
                    list(n = .N, mn = mn[1], vr = vr[1], isf = sum(1/sf)), 
                    by = list(ref, width)]
    # shrPhi[is.na(vr), vr := mn]
    shrPhi[ , phi := (n*vr + mn*isf)/(mn^2*isf), by = ref]
    shrPhi[(phi < 0), phi := 0]
    if (shrink) shrPhi[!is.na(phi), phi := .calcShrPhi(phi), by = width]
    setkey(cnts, ref)
    setkey(shrPhi, ref)
    cnts <- shrPhi[ , list(ref, phi)][cnts]
    setkey(cnts, sbj, ref)
    cnts[ , oldCN := CN]
    calcProb <- function(x) {
      cnts[ , dnbinom(N, sf/phi, 1/(mn*phi*x + 1), log = TRUE)]
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
    
    if (nchng < min.dlt | it == max.its) break
    
    it <- it + 1
    cnts[ , phi := NULL]
    
  }
  
  cnts[ , oldCN := NULL]
  setattr(cnts, "its", it)
  cnts[]
  
}
