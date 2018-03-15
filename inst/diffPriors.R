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
wdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
                  "dfiler", "cnvR")
idir <- file.path(wdir, "simData")
deps <- seq(5, 100, 5) ## Sequencing depths 
cws <- 1:5 ## Sizes of cnvs (number of exons spanned) 

## Set up the file system
odir <- file.path(wdir, "priorAnalysis")
dir.create(odir)
depDir <- with(expand.grid(priors, deps, cws),
               sprintf("p%0.04f/d%0.3d/w%0.1d", Var1, Var2, Var3))
depDir <- sub("0.", "", depDir)
sapply(file.path(odir, depDir), dir.create, recursive = TRUE)

pars <- expand.grid(prior = priors, dep = deps, cw = cws, rep = seq(1000))
pars <- as.data.table(pars)
set.seed(1234)
pars <- pars[pars[ , .I[sample(.N, 20)], by = list(prior, dep, cw)]$V1]
pars[ , wdir := wdir]

doCalc <- function(prior, dep, cw, rep, wdir) {
  ifile <- sprintf("sim_d%0.3d_w%0.1d_r%0.4d.RDS", dep, cw, rep)
  idir <- file.path(sprintf("d%0.3d", dep), sprintf("w%0.1d", cw))
  dat <- readRDS(file.path(wdir, "simData", idir, ifile))
  priorFmt <- sub("0.", "", sprintf("p%0.4f", prior))
  ofile <- sprintf("sRes_%s_d%0.3d_w%0.1d_r%0.4d.RDS", priorFmt, dep, cw, rep)
  odir <- file.path(priorFmt, idir)
  out <- file.path(wdir, "priorAnalysis", odir, ofile)
  smpls <- try(cnvCallCN(cnts = dat, 
                         prior = prior, 
                         outfile = out,
                         agg = TRUE, 
                         shrink = TRUE,
                         keep.cols = c("ref", "sbj", "N", "actCN", 
                                       "width", "CN", "lk"), 
                         width = 5, 
                         verbose = TRUE))
  !is(smpls, 'try-error')
}

res <- slurm_apply(f = doCalc, 
                   params = pars, 
                   nodes = nrow(pars),
                   cpus_per_node = 1,
                   jobname = "priorAnalysis", 
                   slurm_options = list(mem = 20000,
                                        array = sprintf("0-%d%%%d", 
                                                        nrow(pars) - 1, 
                                                        1000),
                                        'cpus-per-task' = 1,
                                        error =  "%A_%a.err",
                                        output = "%A_%a.out",
                                        time = "96:00:00"))




##----------------------------------------------------------------------------##
## Script to rerun failed jobs
##----------------------------------------------------------------------------##

# library(cnvR)
# library(rslurm)
# library(data.table)
# 
# ## Create vector or priors
# priors <- numeric(25)
# priors[1] <- 0.5
# for (i in seq_along(priors)) {
#   priors[i + 1] <- priors[i]/1.3
# }
# priors <- round(priors, 4)
# 
# ## Directory with simulated data
# wdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
#                   "dfiler", "cnvR")
# idir <- file.path(wdir, "simData")
# deps <- seq(5, 100, 5) ## Sequencing depths 
# cws <- 1:5 ## Sizes of cnvs (number of exons spanned) 
# odir <- file.path(wdir, "priorAnalysis")
# 
# pars <- expand.grid(prior = priors, dep = deps, cw = cws, rep = seq(1000))
# pars <- as.data.table(pars)
# set.seed(1234)
# pars <- pars[pars[ , .I[sample(.N, 50)], by = list(prior, dep, cw)]$V1]
# pars[ , wdir := wdir]
# 
# fls <- list.files(odir, recursive = TRUE)
# fls <- basename(fls)
# pars[ , pf := sub("0.", "", sprintf("p%0.4f", prior))]
# pars[ , of := sprintf("sRes_%s_d%0.3d_w%0.1d_r%0.4d.RDS", pf, dep, cw, rep)]
# pars[ , rerun := !of %in% fls]
# pars <- pars[(rerun)]
# pars <- pars[ , list(prior, dep, cw, rep, wdir)]
# 
# doCalc <- function(prior, dep, cw, rep, wdir) {
#   ifile <- sprintf("sim_d%0.3d_w%0.1d_r%0.4d.RDS", dep, cw, rep)
#   idir <- file.path(sprintf("d%0.3d", dep), sprintf("w%0.1d", cw))
#   dat <- readRDS(file.path(wdir, "simData", idir, ifile))
#   priorFmt <- sub("0.", "", sprintf("p%0.4f", prior))
#   ofile <- sprintf("sRes_%s_d%0.3d_w%0.1d_r%0.4d.RDS", priorFmt, dep, cw, rep)
#   odir <- file.path(priorFmt, idir)
#   out <- file.path(wdir, "priorAnalysis", odir, ofile)
#   smpls <- try(cnvCallCN(cnts = dat, prior = prior, outfile = out))
#   !is(smpls, 'try-error')
# }
# 
# res <- slurm_apply(f = doCalc, 
#                    params = pars, 
#                    nodes = nrow(pars),
#                    cpus_per_node = 1,
#                    jobname = "priorRerun2", 
#                    slurm_options = list(mem = 32000,
#                                         array = sprintf("0-%d%%%d", 
#                                                         nrow(pars) - 1, 
#                                                         400),
#                                         'cpus-per-task' = 1,
#                                         error =  "%A_%a.err",
#                                         output = "%A_%a.out",
#                                         time = "96:00:00"))
# 
# 
# ##----------------------------------------------------------------------------##
# ## Script to analyze the results
# ##----------------------------------------------------------------------------##
# 
# library(cnvR)
# library(rslurm)
# library(data.table)
# library(stringr)
# 
# wdir <- file.path("/", "projects", "sequence_analysis", "vol5", 
#                   "dfiler", "cnvR")
# odir <- file.path(wdir, "priorAnalysis")
# pars <- list.files(odir, recursive = TRUE, full.names = TRUE)
# pars <- data.table(objFile = pars)
# 
# priorAnalysis <- function(objFile) {
#   
#   prior <- as.numeric(sub("p", "", str_extract(objFile, "p[0-9]+")))/10000
#   width <- as.numeric(sub("w", "", str_extract(objFile, "w[0-9]+")))
#   depth <- as.numeric(sub("d", "", str_extract(objFile, "d[0-9]+")))
#   
#   dat <- readRDS(objFile)
#   dat <- cnvAggCall(dat)
#   out <- data.table(objFile = objFile,
#                     prior = prior,
#                     width = width,
#                     depth = depth,
#                     tn = dat[actCNSngl == 1 & CN == 1, .N],
#                     tp = dat[actCNSngl != 1 & CN != 1, .N],
#                     fp = dat[actCNSngl == 1 & CN != 1, .N],
#                     fn = dat[actCNSngl != 1 & CN == 1, .N],
#                     fc = dat[actCNSngl != CN, .N],
#                     calls = list(dat[ , .N, by = list(actCNSngl, CN)]))
#   
#   out[]
#   
# }
# 
# jobSetup <- function(pars, name) {
#   slurm_apply(f = priorAnalysis, 
#               params = pars, 
#               nodes = nrow(pars),
#               cpus_per_node = 1,
#               jobname = name, 
#               slurm_options = list(mem = 12228,
#                                    array = sprintf("0-%d%%%d", 
#                                                    nrow(pars) - 1, 
#                                                    800),
#                                    'cpus-per-task' = 1,
#                                    error =  "%A_%a.err",
#                                    output = "%A_%a.out",
#                                    time = "96:00:00"))
# }
# 
# p1 <- pars[1:60000]
# p2 <- pars[60001:120000]
# p3 <- pars[120001:nrow(pars)]
# 
# jobSetup(p1, "priorAnalysis1")
# jobSetup(p2, "priorAnalysis2")
# jobSetup(p3, "priorAnalysis3")

