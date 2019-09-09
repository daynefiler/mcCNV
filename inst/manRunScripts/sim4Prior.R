##----------------------------------------------------------------------------##
## Script to do the prior selection analysis -- take 2
##----------------------------------------------------------------------------##

library(cnvR)
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
idir <- file.path(wdir, "sim4Data")

## Set up the file system
odir <- file.path(wdir, "sim4Prior")
pdirs <- sub("0.", "", sprintf("p%0.04f", priors))
mkdirs <- sapply(file.path(odir, pdirs), dir.create, recursive = TRUE)
all(mkdirs)

pars <- expand.grid(prior = priors, rep = seq(500))
pars <- as.data.table(pars)
pars[ , wdir := wdir]

doCalc <- function(prior, rep, wdir) {
  dat <- readRDS(file.path(wdir, "sim4Data", sprintf("sim_r%0.4d.RDS", rep)))
  priorFmt <- sub("0.", "", sprintf("p%0.4f", prior))
  ofile <- sprintf("sRes_%s_r%0.4d.RDS", priorFmt, rep)
  out <- file.path(wdir, "sim4Prior", priorFmt, ofile)
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
            jobname = "sim4Prior", 
            slurm_options = list(mem = 16000,
                                 array = sprintf("0-%d%%%d", 
                                                 nrow(pars) - 1, 
                                                 500),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "10-00:00:00"))

