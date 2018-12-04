##----------------------------------------------------------------------------##
## Script to do the analysis on the whole simulation at prior of 0.06
##----------------------------------------------------------------------------##

library(cnvR)
library(rslurm)
library(data.table)

## Directory with simulated data
wdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
                  "dfiler", "CNV")
idir <- file.path(wdir, "sim1Data")
deps <- seq(5, 100, 5) ## Sequencing depths 
cws <- 1:5 ## Sizes of cnvs (number of exons spanned) 

## Set up the file system
odir <- file.path(wdir, "sim1Analysis")
depDir <- with(expand.grid(deps, cws),
               sprintf("d%0.3d/w%0.1d", Var1, Var2))
mkdirs <- sapply(file.path(odir, depDir), dir.create, recursive = TRUE)
all(mkdirs)

pars <- expand.grid(dep = deps, cw = cws, rep = seq(200))
pars <- as.data.table(pars)
pars[ , wdir := wdir]
pars[ , prior := 0.06]

doCalc <- function(prior, dep, cw, rep, wdir) {
  ifile <- sprintf("sim_d%0.3d_w%0.1d_r%0.4d.RDS", dep, cw, rep)
  idir <- file.path(sprintf("d%0.3d", dep), sprintf("w%0.1d", cw))
  dat <- readRDS(file.path(wdir, "sim1Data", idir, ifile))
  priorFmt <- sub("0.", "", sprintf("p%0.4f", prior))
  ofile <- sprintf("sRes_%s_d%0.3d_w%0.1d_r%0.4d.RDS", priorFmt, dep, cw, rep)
  out <- file.path(wdir, "sim1Analysis", idir, ofile)
  kp <- c("ref", "sbj", "N", "actCN", "mn", "phi", "width", "CN", "lp", "lp1")
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
    res <- smpls[ , .N, by = list(ACT, C1, PRO)]
    setorder(res, -ACT, -C1, -PRO)
    res[ , prior := prior]
    res[ , dep := dep]
    res[ , width := cw]
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
            jobname = "sim1Analysis", 
            slurm_options = list(mem = 16000,
                                 array = sprintf("0-%d%%%d", 
                                                 nrow(pars) - 1, 
                                                 500),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "10-00:00:00"))







