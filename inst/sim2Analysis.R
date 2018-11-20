##----------------------------------------------------------------------------##
## Script to do the analysis on the whole simulation 2 at prior of 0.03
##----------------------------------------------------------------------------##

library(cnvR)
library(rslurm)
library(data.table)

## Directory with simulated data
wdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
                  "dfiler", "CNV")
idir <- file.path(wdir, "sim2Data")
deps <- seq(5, 100, 5) ## Sequencing depths 


## Set up the file system
odir <- file.path(wdir, "sim2Analysis")
depDir <- sprintf("d%0.3d", deps)
mkdirs <- sapply(file.path(odir, depDir), dir.create, recursive = TRUE)
all(mkdirs)

pars <- expand.grid(dep = deps, rep = seq(1000))
pars <- as.data.table(pars)
pars[ , wdir := wdir]
pars[ , prior := 0.06]

doCalc <- function(prior, dep, rep, wdir) {
  ifile <- sprintf("sim_d%0.3d_r%0.4d.RDS", dep, rep)
  idir <- file.path(sprintf("d%0.3d", dep))
  dat <- readRDS(file.path(wdir, "sim2Data", idir, ifile))
  priorFmt <- sub("0.", "", sprintf("p%0.4f", prior))
  ofile <- sprintf("sRes_%s_d%0.3d_r%0.4d.RDS", priorFmt, dep, rep)
  out <- file.path(wdir, "sim2Analysis", idir, ofile)
  kp <- c("ref", "sbj", "N", "actCN", "mn", "phi", "width", "CN", "lk", "lk1")
  smpls <- try(cnvCallCN(cnts = dat, 
                         prior = prior, 
                         outfile = out,
                         agg = TRUE, 
                         shrink = TRUE,
                         keep.cols = kp, 
                         width = 5,
                         verbose = TRUE))
  if (!is(smpls, 'try-error')) {
    smpls[ , ACT := actCNSngl != 1]
    makeNull <- function(x) paste(rep(1, x), collapse = ":")
    smpls[ , null := sapply(width, makeNull)]
    smpls[ , PRO := actCN != null]
    res <- smpls[ , .N, by = list(ACT, C1, C2, PRO)]
    setorder(res, -ACT, -C1, -C2, -PRO)
    res[ , prior := prior]
    res[ , dep := dep]
    res[ , rep := rep]
    return(res[])
  } else {
    return(FALSE)
  }
}

slurm_apply(f = doCalc, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "sim2Analysis", 
            slurm_options = list(mem = 16000,
                                 array = sprintf("0-%d%%%d", 
                                                 nrow(pars) - 1, 
                                                 300),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "10-00:00:00"))


##----------------------------------------------------------------------------##
## After the job finishes...
##----------------------------------------------------------------------------##

library(cnvR)
library(rslurm)
library(data.table)
library(lattice)
library(latticeExtra)
library(RColorBrewer)

col2alpha <- function(col, alpha = 0.5) {
  tmp <- col2rgb(col)
  rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha)  
}

calcMCC <- function(tp, tn, fp, fn) {
  tp <- as.numeric(tp)
  tn <- as.numeric(tn)
  fp <- as.numeric(fp)
  fn <- as.numeric(fn)
  (tp*tn - fp*fn)/(sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)))
}

# sjob <- slurm_job("sim2Analysis", 20000)
# res <- get_slurm_out(sjob, ncores = 60)
# res <- rbindlist(res)
# saveRDS(res, "sim2Analysis.RDS")
res <- procRes(readRDS("sim2Analysis.RDS"))

pltStat(res[[1]], res[[2]], "tpr", "TPR", 0.95)
pltStat(res[[1]], res[[2]], "pf", "FP/(TP + FP)", 0.05)
pltStat(res[[1]], res[[2]], "mcc", "MCC", 0.95)

