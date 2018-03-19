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

pars <- expand.grid(prior = priors, dep = deps, cw = cws, rep = seq(200))
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
                         weight = TRUE,
                         verbose = TRUE))
  if (!is(smpls, 'try-error')) {
    smpls[ , ACT := actCNSngl != 1]
    smpls[is.na(C2), C2 := FALSE]
    res <- smpls[ , .N, by = list(ACT, C1, C2)][order(-ACT, -C1, -C2)]
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
            slurm_options = list(mem = 20000,
                                 array = sprintf("0-%d%%%d", 
                                                 nrow(pars) - 1, 
                                                 1000),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "96:00:00"))


##----------------------------------------------------------------------------##
## Script to analyze the results
##----------------------------------------------------------------------------##

library(cnvR)
library(rslurm)
library(data.table)
library(lattice)
library(latticeExtra)

sjob <- slurm_job("priorAnalysis", 52000)

res <- get_slurm_out(sjob)
saveRDS(res, "~/Desktop/priorRes.RDS")
res <- rbindlist(res)

res[ , fdep := as.factor(dep)]
res[ , fwid := as.factor(width)]

resSmry <- res[ , 
               list(tp1 = sum(N[ ACT &  C1]),
                    fp1 = sum(N[!ACT &  C1]),
                    tn1 = sum(N[!ACT & !C1]),
                    fn1 = sum(N[ ACT & !C1]),
                    tp2 = sum(N[ ACT &  C2]),
                    fp2 = sum(N[!ACT &  C2]),
                    tn2 = sum(N[!ACT & !C2]),
                    fn2 = sum(N[ ACT & !C2])),
               by = list(prior, fdep, fwid, rep)]

resSmry[ , fpr1 := fp1/(fp1 + tn1)]
resSmry[ , tpr1 := tp1/(tp1 + fn1)]
resSmry[ , spc1 := tn1/(fp1 + tn1)]
resSmry[ , pf1  := fp1/(fp1 + tp1)]
resSmry[ , fpr2 := fp2/(fp2 + tn2)]
resSmry[ , tpr2 := tp2/(tp2 + fn2)]
resSmry[ , spc2 := tn2/(fp2 + tn2)]
resSmry[ , pf2  := fp2/(fp2 + tp2)]

resMn <- resSmry[ , lapply(.SD, mean), by = list(prior, fdep, fwid)]
resSD <- resSmry[ , lapply(.SD, sd),   by = list(prior, fdep, fwid)]

plt <- xyplot(tpr1 ~ pf1 | fdep + fwid, data = resMn, xlim = c(0, 1),
              xlab = "FP/(TP + FP)",
              ylab = "TPR",
              panel = function(x, y, ...) {
                trellis.par.set(pty = "s")
                panel.xyplot(x, y, type = "l", lwd = 2)
                # panel.abline(h = 1, lty = "dashed")
                panel.abline(v = 0.05, lty = "dashed")
                panel.abline(a = 0, b = 1, lty = "dashed")
                panel.abline(h = 0.95, lty = "dashed")
              })

plt <- plt + as.layer(xyplot(tpr2 ~ pf2 | fdep + fwid, data = resMn, type = "l", col = "red", lwd = 2))

pdf("~/Github/cnvR/inst/junk/priorRes.pdf", width = 11, height = 8.5)
plt
graphics.off()







