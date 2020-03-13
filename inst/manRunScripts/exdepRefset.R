library(data.table)
library(ExomeDepth)
library(rslurm)
library(dlfUtils)

dfil <- "/projects/sequence_analysis/vol5/dfiler"
cnvd <- file.path(dfil, "CNV", "realData37")

cts <- readRDS(file.path(cnvd, "allCounts.counts"))
smpls <- data.table(smpl = grep("_", names(cts), value = TRUE))

getRefSet <- function(smpl) {
  
  dfil <- "/projects/sequence_analysis/vol5/dfiler"
  cnvd <- file.path(dfil, "CNV", "realData37")
  cts <- readRDS(file.path(cnvd, "allCounts.counts"))
  mat <- as.matrix(cts[ , .SD, .SDcols = grep("_", names(cts), value = TRUE)])
  mat[is.na(mat)] <- 0
  wids <- cts$wid
  oth <- setdiff(colnames(mat), smpl)
  sel <- select.reference.set(mat[ , smpl], mat[ , oth], wids, smpl)
  res <- as.data.table(sel[[2]])
  setkey(res, ref.samples)
  res[ , slct := FALSE]
  res[J(sel[[1]]), slct := TRUE]
  res[ , smpl := smpl]
  res[]
  
}

slurmArray(getRefSet, 
           params = smpls, 
           jobname = "exdepRefset")



