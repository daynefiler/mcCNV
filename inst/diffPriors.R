##----------------------------------------------------------------------------##
## Script to do the prior selection analysis
##----------------------------------------------------------------------------##

library(cnvR)
library(rslurm)
library(data.table)

## Create vector or priors
priors <- numeric(25)
priors[1] <- 0.5
for (i in seq_along(priors)) {
  priors[i + 1] <- priors[i]/1.3
}
priors <- round(priors, 4)

## Directory with simulated data
wdir <- file.path("/", "projects", "sequence_analysis", "vol4", 
                  "dfiler", "cnvR")
idir <- file.path(wdir, "simData")
deps <- seq(5, 100, 5) ## Sequencing depths 
cws <- 1:5 ## Sizes of cnvs (number of exons spanned) 

## Set up the file system
# odir <- file.path(wdir, "priorAnalysis")
# dir.create(odir)
# depDir <- with(expand.grid(priors, deps, cws),
#                sprintf("p%0.04f/d%0.3d/w%0.1d", Var1, Var2, Var3))
# depDir <- sub("0.", "", depDir)
# sapply(file.path(odir, depDir), dir.create, recursive = TRUE)

pars <- expand.grid(prior = priors, dep = deps, cw = cws, rep = seq(1000))
pars <- as.data.table(pars)
set.seed(1234)
pars <- pars[pars[ , .I[sample(.N, 50)], by = list(prior, dep, cw)]$V1]
pars[ , wdir := wdir]

doCalc <- function(prior, dep, cw, rep, wdir) {
  ifile <- sprintf("sim_d%0.3d_w%0.1d_r%0.4d.RDS", dep, cw, rep)
  idir <- file.path(sprintf("d%0.3d", dep), sprintf("w%0.1d", cw))
  dat <- readRDS(file.path(wdir, "simData", idir, ifile))
  priorFmt <- sub("0.", "", sprintf("p%0.4f", prior))
  ofile <- sprintf("sRes_%s_d%0.3d_w%0.1d_r%0.4d.RDS", priorFmt, dep, cw, rep)
  odir <- file.path(priorFmt, idir)
  out <- file.path(wdir, "priorAnalysis", odir, ofile)
  smpls <- try(cnvCallCN(cnts = dat, prior = prior, outfile = out))
  !is(smpls, 'try-error')
}

res <- slurm_apply(f = doCalc, 
                   params = pars, 
                   nodes = 400,
                   cpus_per_node = 5,
                   jobname = "cnvDiffPrior", 
                   slurm_options = list(mem = 16000,
                                        'cpus-per-task' = 1,
                                        error = "cnvDiffPrior-err.txt",
                                        time = "96:00:00"))


