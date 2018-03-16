##----------------------------------------------------------------------------##
## .calcShrPhi: shrink the phi values
##----------------------------------------------------------------------------##

#' @name calcShrPhi
#' @title Shrink phi values
#' 
#' @param phi numeric, the phi values to shrink
#' 
#' @details 
#' Need to add.
#' 

.calcShrPhi <- function(phi) {
  
  nPhi <- length(phi)
  xi <- seq(0, 1, length.out = 1000)
  asd <- sapply(xi, function(x) 0.05/sum((phi - x)^2))
  xi <- xi[which.min(diff(asd)/diff(xi))]
  dlt <- (sum((phi - mean(phi))^2)/(nPhi - 1))/(sum((phi - xi)^2)/(nPhi - 2))
  sPhi <- (1 - dlt)*phi + dlt*xi
  sPhi
  
}

##----------------------------------------------------------------------------##
## .calcShrPhi: Calculate the MLCN-state for a subject
##----------------------------------------------------------------------------##

#' @name calcMLCN
#' @title Calculate the maximum likelihood copy number-state for a subject
#' 
#' @param N integer of length 1, the molecule counts for the subject
#' @param mu numeric of length 1, the locus mean counts for the group 
#' @param sf numeric of length 1, the size-factor adjustment for the subject
#' @param phi numeric of length 1, the locus phi/dispersion value for the group 
#' @param prior numeric of length 1, the prior probabilty of having a CNV
#' 
#' @details 
#' Need to add
#' 
#' @return list, first element gives the most likely copy number-state and 
#' the second element gives the likelihood value
#' 

.calcMLCN <- function(N, mu, sf, phi, prior) {
  
  cs <- c(0.001, 0.5, 1, 1.5, 2.0, 2.5, 3.0, 3.5, 4)
  cp <- rep(prior, length(cs)); cp[which(cs == 1)] <- 1 - prior
  pr <- lapply(cs, function(x) list(p = 1/(mu*phi*x + 1), r = sf/phi))
  cdf <- sapply(pr, function(x) pnbinom(N, size = x$r, prob = x$p))
  lk <- ifelse(cdf >= 0.5, 2*(1 - cdf), 2*cdf)
  list(cs[which.max(lk*cp)], max(lk*cp))
  
}

##----------------------------------------------------------------------------##
## .clpsExon: Collapse exons by specified width 
##----------------------------------------------------------------------------##

#' @name clpsExon
#' @title Collapse exons by specified width 
#' 
#' @param dat data.table object containing the counts
#' @param cw integer of length 1, the width to collapse by
#' 
#' @details 
#' Need to add
#' 
#' @import data.table

.clpsExon <- function(dat, cw) {
  
  dat <- dat[ , 
              shift(.SD, 1:cw - 1, type = "lead", give.names = TRUE), 
              by = sbj]
  dat[ , N := rowSums(.SD), .SDcols = grep("N", colnames(dat))]
  dat[ , 
       ref := do.call(paste, c(.SD, sep = ":")), 
       .SDcols = grep("ref", colnames(dat))]
  if (any(grepl("actCN", colnames(dat)))) {
    dat[ ,
         actCN := do.call(paste, c(.SD, sep = ":")), 
         .SDcols = grep("actCN", colnames(dat))]
  }
  dat[ , grep("_lead_", colnames(dat), value = TRUE) := NULL]
  dat <- dat[!is.na(N)]
  dat[ , width := cw]
  dat[]
  
}

##----------------------------------------------------------------------------##
## .callCN: Call CNS
##----------------------------------------------------------------------------##

#' @name callCN
#' @title Call copy number-state
#' 
#' @param cnts data.table object containing the counts
#' @param prior numeric of length 1, the prior probability of having a CNV
#' @param min.dlt integer of length 1, the target number of changes in copy-
#' state to stop the alogorithm 
#' @param max.its integer of length 1, the maximum number of iterations
#' @param shrink logical of length 1, shrinkage applied to phi when TRUE
#' 
#' @details 
#' Need to add
#' 
#' @importFrom Rfast rowMaxs
#' @import data.table

.callCN <- function(cnts, min.dlt, max.its, prior, shrink = TRUE) {
  
  cnts[ , CN := 1]
  
  cs <- c(0.001, 0.5, 1, 1.5, 2.0, 2.5, 3.0, 3.5, 4)
  cp <- rep(prior, length(cs)); cp[which(cs == 1)] <- 1 - prior
  
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
    shrPhi[is.na(vr), vr := mn]
    shrPhi[ , phi := (n*vr + mn*isf)/(mn^2*isf), by = ref]
    if (shrink) shrPhi[ , phi := .calcShrPhi(phi), by = width]
    setkey(cnts, ref)
    setkey(shrPhi, ref)
    cnts <- shrPhi[ , list(ref, phi)][cnts]
    setkey(cnts, sbj, ref)
    cnts[ , oldCN := CN]
    calcProb <- function(x) cnts[ , pnbinom(N, sf/phi, 1/(mn*phi*x + 1))]
    calcExpt <- function(x) cnts[ , N > sf*mn*x]
    probMat <- do.call(cbind, lapply(cs, calcProb))
    exptMat <- do.call(cbind, lapply(cs, calcExpt))
    exptMat[is.na(exptMat)] <- FALSE
    probMat[exptMat] <- 1 - probMat[exptMat]
    rm(exptMat); gc()
    probMat <- sweep(probMat, 2, cp, "*")
    cnts[ , CN := cs[rowMaxs(probMat)]]
    cnts[ , lk := rowMaxs(probMat, value = TRUE)]
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

