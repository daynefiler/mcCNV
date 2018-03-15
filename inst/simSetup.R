##----------------------------------------------------------------------------##
## Script to create the simulated dataset
##----------------------------------------------------------------------------##

library(cnvR)
library(rslurm)
library(data.table)

## Setup simulation
## Directory for storing the data
fdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
                  "dfiler", "cnvR", "simData")
deps <- seq(5, 100, 5) ## Sequencing depths 
cw <- 1:5 ## Sizes of cnvs (number of exons spanned) 
reps <- 200 ## Number of repetitions 
Ne <- 172000 ## Number of exons
Ns <- 16 ## Number of samples

## Set up the file system
if (dir.exists(fdir)) unlink(fdir, recursive = TRUE, force = TRUE)
dir.create(fdir)
depDir <- sprintf("d%0.3d/w%0.1d", rep(deps, each = length(cw)), cw)
sapply(file.path(fdir, depDir), dir.create, recursive = TRUE)

pars <- expand.grid(Ns = Ns, Ne = Ne, dep = deps, 
                    cw = cw, rep = seq(reps), fdir = fdir)
pars <- as.data.table(pars)
set.seed(1234)
pars[ , seed := sample(1e6, .N)]

savePool <- function(Ns, Ne, dep, cw, rep, seed, fdir) {
  wndw <- c(dep - 3, dep + 3)*1e6
  smpls <- try(cnvGenPool(ns = Ns, ne = Ne, wndw = wndw, cw = cw, seed = seed))
  odir <- file.path(fdir, sprintf("d%0.3d", dep), sprintf("w%0.1d", cw))
  fname <- sprintf("sim_d%0.3d_w%0.1d_r%0.4d.RDS", dep, cw, rep)
  saveRDS(smpls, file = file.path(odir, fname))
  !is(smpls, 'try-error')
}

slurm_apply(f = savePool, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "genSmpls",
            slurm_options = list(mem = 12000,
                                 array = sprintf("0-%d%%%d",
                                                 nrow(pars) - 1,
                                                 1000),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "96:00:00"))



