##----------------------------------------------------------------------------##
## Script to create the simulated dataset with no CNVs
##----------------------------------------------------------------------------##

library(cnvR)
library(rslurm)
library(data.table)

## Setup simulation
## Directory for storing the data
fdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
                  "dfiler", "CNV", "sim3Data")
deps <- seq(5, 100, 5) ## Sequencing depths 
reps <- 200 ## Number of repetitions 
Ne <- 172000 ## Number of exons
Ns <- 16 ## Number of samples

## Set up the file system
if (dir.exists(fdir)) unlink(fdir, recursive = TRUE, force = TRUE)
dir.create(fdir)
depDir <- sprintf("d%0.3d", deps)
sapply(file.path(fdir, depDir), dir.create, recursive = TRUE)
mkdirs <- sapply(file.path(odir, depDir), dir.create, recursive = TRUE)
all(mkdirs)

pars <- expand.grid(Ns = Ns, Ne = Ne, dep = deps, rep = seq(reps), fdir = fdir)
pars <- as.data.table(pars)
set.seed(9101)
pars[ , seed := sample(1e6, .N)]

savePool <- function(Ns, Ne, dep, rep, seed, fdir) {
  wndw <- c(dep - 3, dep + 3)*1e6
  smpls <- try(cnvGenPool(ns = Ns, ne = Ne, wndw = wndw, seed = seed,
                          cw = 1, pc = 1, cs = 1))
  odir <- file.path(fdir, sprintf("d%0.3d", dep))
  fname <- sprintf("sim_d%0.3d_r%0.4d.RDS", dep, rep)
  saveRDS(smpls, file = file.path(odir, fname))
  !is(smpls, 'try-error')
}

slurm_apply(f = savePool, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "genSmpls3",
            slurm_options = list(mem = 12000,
                                 array = sprintf("0-%d%%%d",
                                                 nrow(pars) - 1,
                                                 1000),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "96:00:00"))

