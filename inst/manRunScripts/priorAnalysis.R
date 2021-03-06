##----------------------------------------------------------------------------##
## Script to do the prior selection analysis
##----------------------------------------------------------------------------##

library(mcCNV)
library(rslurm)
library(data.table)

## Create vector or priors
priors <- numeric(25)
priors[1] <- 1/6
for (i in seq_along(priors)) {
  priors[i + 1] <- priors[i]/1.3
}
priors <- round(priors, 4)

## Directory with simulated data
wdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
                  "dfiler", "CNV")
idir <- file.path(wdir, "sim1Data")
deps <- seq(5, 100, 5) ## Sequencing depths 
cws <- 1:5 ## Sizes of cnvs (number of exons spanned) 

## Set up the file system
odir <- file.path(wdir, "priorAnalysis")
depDir <- with(expand.grid(priors, deps, cws),
               sprintf("p%0.04f/d%0.3d/w%0.1d", Var1, Var2, Var3))
depDir <- sub("0.", "", depDir)
mkdirs <- sapply(file.path(odir, depDir), dir.create, recursive = TRUE)
all(mkdirs)

pars <- expand.grid(prior = priors, dep = deps, cw = cws, rep = seq(200))
pars <- as.data.table(pars)
set.seed(1234)
pars <- pars[pars[ , .I[sample(.N, 20)], by = list(prior, dep, cw)]$V1]
pars[ , wdir := wdir]

doCalc <- function(prior, dep, cw, rep, wdir) {
  ifile <- sprintf("sim_d%0.3d_w%0.1d_r%0.4d.RDS", dep, cw, rep)
  idir <- file.path(sprintf("d%0.3d", dep), sprintf("w%0.1d", cw))
  dat <- readRDS(file.path(wdir, "sim1Data", idir, ifile))
  priorFmt <- sub("0.", "", sprintf("p%0.4f", prior))
  ofile <- sprintf("sRes_%s_d%0.3d_w%0.1d_r%0.4d.RDS", priorFmt, dep, cw, rep)
  odir <- file.path(priorFmt, idir)
  out <- file.path(wdir, "priorAnalysis", odir, ofile)
  kp <- c("ref", "sbj", "N", "actCN", "mn", "phi", "width", "CN", "lp", "lp1")
  smpls <- try(cnvCallCN(cnts = dat, 
                         prior = prior, 
                         outfile = out,
                         agg = TRUE, 
                         shrink = TRUE,
                         keep.cols = kp, 
                         width = 1:5,
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
            jobname = "priorAnalysis", 
            slurm_options = list(mem = 16000,
                                 array = sprintf("0-%d%%%d", 
                                                 nrow(pars) - 1, 
                                                 500),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "10-00:00:00"))




