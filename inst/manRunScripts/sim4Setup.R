library(cnvR)
library(rslurm)
library(data.table)

## Setup simulation
## Directory for storing the data
fdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
                  "dfiler", "CNV", "sim4Data")

## Set up the file system
if (dir.exists(fdir)) unlink(fdir, recursive = TRUE, force = TRUE)
dir.create(fdir)

pars <- data.table(rep = 1:500, fdir = fdir)
set.seed(2222)
pars[ , seed := sample(1e6, .N)]

savePool <- function(rep, seed, fdir) {
  wndw <- c(35031844, 62614330)
  smpls <- try(cnvGenPool(ns = 16, ne = 172736, wndw = wndw, seed = seed,
                          cw = 5, meanlog = -12.27, sdlog = 0.6588))
  fname <- sprintf("sim_r%0.4d.RDS", rep)
  saveRDS(smpls, file = file.path(fdir, fname))
  !is(smpls, 'try-error')
}

slurm_apply(f = savePool, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "genSmpls4",
            slurm_options = list(mem = 12000,
                                 array = sprintf("0-%d%%%d",
                                                 nrow(pars) - 1,
                                                 1000),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "96:00:00"))

