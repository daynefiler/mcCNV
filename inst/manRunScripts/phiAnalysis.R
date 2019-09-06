##----------------------------------------------------------------------------##
## Script to do the analysis on the whole simulation at prior of 0.03
##----------------------------------------------------------------------------##

library(cnvR)
library(rslurm)
library(data.table)

## Directory with simulated data
wdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
                  "dfiler", "CNV")
idir <- file.path(wdir, "sim1Analysis")
deps <- seq(5, 100, 5) ## Sequencing depths 
cws <- 1:5 ## Sizes of cnvs (number of exons spanned) 

pars <- expand.grid(dep = deps, cw = cws, rep = seq(200))
pars <- as.data.table(pars)
pars[ , wdir := wdir]

getPhi <- function(dep, cw, rep, wdir) {
  ifile <- sprintf("sRes_p0600_d%0.3d_w%0.1d_r%0.4d.RDS", dep, cw, rep)
  idir <- file.path(sprintf("d%0.3d", dep), sprintf("w%0.1d", cw))
  smpls <- readRDS(file.path(wdir, "sim1Analysis", idir, ifile))
  if (!is(smpls, 'try-error')) {
    smpls <- smpls[ , list(ref, width, mn, phi)]
    smpls <- smpls[ ,
                   list(mn = list(summary(mn)), phi = list(summary(phi))),
                   by = width]
    smpls[ , dep := dep]
    smpls[ , width := cw]
    smpls[ , rep := rep]
    return(smpls[])
  } else {
    return(FALSE)
  }
}

slurm_apply(f = getPhi, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "phiAnalysis", 
            slurm_options = list(mem = 4000,
                                 array = sprintf("0-%d%%%d", 
                                                 nrow(pars) - 1, 
                                                 1000),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "96:00:00"))


