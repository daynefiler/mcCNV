##----------------------------------------------------------------------------##
## Script to analyze the real data
##----------------------------------------------------------------------------##

library(data.table)
library(rslurm)
library(mcCNV)

cnvDir <- "/projects/sequence_analysis/vol5/dfiler/CNV"
jbnm <- "realDataAnalysis"
rdDir <- file.path(cnvDir, "realData37")
resDir <- file.path(rdDir, "results")
slrmDir <- file.path(cnvDir, paste0("_rslurm_", jbnm))
if (dir.exists(slrmDir)) unlink(slrmDir, recursive = TRUE)
if (dir.exists(resDir)) unlink(resDir, recursive = TRUE)
dir.create(resDir)
countDirs <- Sys.glob(file.path(rdDir, "*", "counts"))
allFls <- data.table(path = list.files(countDirs, full.names = TRUE))
allFls[ , cap := gsub(paste0(rdDir, "|/|counts.+$"), "", path)]
## Do not include IDT-IC in 'ind' even though independently captured. Here, the 
## 'ind' flag is really used for sampling pools and we only have 16 samples.
allFls[ , ind := cap %in% c("otoIC", "NCGENES")]
## Exclude C047_B096 due to very low reads
allFls <- allFls[!grepl("C047_B096", path)]

getSbj <- function(y) allFls[cap == y, unique(path)]
otoIC <- getSbj("otoIC")
ncgIC <- getSbj("NCGENES")
## Set seed twice to match the dirichlet pools
set.seed(5678)
pools <- lapply(1:30, function(x) sample(otoIC, size = 16))
pools <- c(pools, lapply(1:30, function(x) sample(otoIC, size = 12)))
set.seed(5678)
pools <- c(pools, lapply(1:30, function(x) sample(ncgIC, size = 16)))
pools <- c(pools, lapply(allFls[!(ind), unique(cap)], getSbj))
proj <- c(sprintf("otoIC16-%02d", 1:30),
          sprintf("otoIC12-%02d", 1:30),
          sprintf("NCGENES-%02d", 1:30),
          allFls[!(ind), unique(cap)])
pars <- data.table(sbj = pools, proj = proj)
pars[ , wdir := file.path(rdDir, "results")]
pars[             , prior := 0.060]
pars[grepl("WGS|IDT", proj), prior := 0.006]

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


##----------------------------------------------------------------------------##
## Aggregate calls for each result
##----------------------------------------------------------------------------##

library(data.table)
library(rslurm)
library(mcCNV)

cnvDir <- "/projects/sequence_analysis/vol5/dfiler/CNV"
resDir <- file.path(cnvDir, "realData37", "results")
resFls <- list.files(resDir, full.names = TRUE, pattern = ".results$")
resFls <- data.table(fl = resFls)

aggRes <- function(fl) {
  dat <- readRDS(fl)
  agg <- try(cnvAggCall(dat, width = 5))
  if (is(agg, "try-error")) return(agg)
  saveRDS(agg, paste0(fl, ".agg"))
  TRUE
}

slurm_apply(f = aggRes, 
            params = resFls, 
            nodes = nrow(resFls),
            cpus_per_node = 1,
            jobname = "aggRealRes",
            slurm_options = list(mem = 20000,
                                 array = sprintf("0-%d", nrow(resFls) - 1),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "8-00:00:00"))

