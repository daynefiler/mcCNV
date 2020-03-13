##----------------------------------------------------------------------------##
## Script to analyze the real data
##----------------------------------------------------------------------------##

library(data.table)
library(rslurm)
library(cnvR)
library(tidyr)
library(dlfUtils)

cnvDir <- "/projects/sequence_analysis/vol5/dfiler/CNV"
jbnm <- "idtSubset"
rdDir <- file.path(cnvDir, "realData37")
resDir <- file.path(rdDir, "idtSubset")
slrmDir <- file.path(cnvDir, paste0("_rslurm_", jbnm))
if (dir.exists(slrmDir)) unlink(slrmDir, recursive = TRUE)
if (dir.exists(resDir)) unlink(resDir, recursive = TRUE)
dir.create(resDir)
countDirs <- Sys.glob(file.path(rdDir, "*", "counts"))
allFls <- data.table(path = list.files(countDirs, full.names = TRUE))
allFls[ , proj := gsub(paste0(rdDir, "|/|counts.+$"), "", path)]
allFls <- allFls[!grepl("oto", proj)]
allFls[ , smpl := sub(".counts", "", basename(path))]

cts <- readRDS(file.path(rdDir, "allCounts.counts"))
metaCols <- c("ref", "GC", "probeCount", "chr", "start", "end", "wid")
smplCts <- cts[ , 
                lapply(.SD, sum, na.rm = TRUE), 
                .SDcols = setdiff(names(cts), metaCols)]
smplCts <- melt(smplCts, 
                variable.name = "projSmpl", 
                value.name = "N", 
                measure.vars = names(smplCts))
smplCts <- separate(smplCts, "projSmpl", 
                    c("proj", "smpl"), 
                    sep = "_", 
                    extra = "merge", 
                    remove = FALSE)
setkey(smplCts, proj, smpl)
setkey(allFls, proj, smpl)

allFls <- smplCts[allFls]
bnd <- allFls[ , {x = reduce3mad(N, TRUE); .(l = x[1], h = x[2])}, by = proj]
allFls <- merge(allFls,  bnd)
allFls <- allFls[N >= l & N <= h]

setkey(1234)
smplPath <- function(x, path) sample(path, size = 9)
pars <- allFls[ , .(sbj = lapply(1:20, smplPath, path = path)), by = proj]
pars[ , ind := 1:.N, by = proj]
pars <- pars[ind == 1 | proj != "IDT-MC"]
pars[ , proj := sprintf("red-%s-%02d", proj, ind)]
pars[ , wdir := resDir]
pars[ , prior := 0.060]
pars[grepl("WGS|IDT", proj), prior := 0.006]
pars[ , ind := NULL]

runCalc <- function(sbj, proj, wdir, prior) {
  sbj <- unlist(sbj)
  outFile <- file.path(wdir, sprintf("%s.results", proj))
  kp <- c("ref", "sbj", "N", "mn", "phi", "width", "CN", "lp", "lp1", "vr")
  cts <- cnvGatherCounts(sbj)
  res <- try(cnvCallCN(cnts = cts, 
                       prior = prior, 
                       keep.cols = kp, 
                       verbose = TRUE, 
                       outfile = outFile, 
                       return.res = FALSE))
  if (is(res, "try-error")) return(res)
  TRUE
}

slurm_apply(f = runCalc, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = jbnm,
            slurm_options = list(mem = 20000,
                                 array = sprintf("0-%d", nrow(pars) - 1),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "8-00:00:00"))

