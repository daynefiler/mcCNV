##----------------------------------------------------------------------------##
## Script to fit the dirichlet distribution
##----------------------------------------------------------------------------##

library(data.table)
library(rslurm)
library(Rfast)

# cts <- readRDS("realData37/AllSamples.counts")
# cts[ , proj := dirname(dirname(sbj))]
# saveRDS(cts[N < 5, .N, by = list(sbj, N)], file = "realData37/lowCounts.RDS")
# saveRDS(cts[N > 2000, .N, by = list(sbj)], file = "realData37/highCounts.RDS")
# nExonsLow <- cts[N < 5, length(unique(ref))]
# getUniqueLow <- function(i) {
#   nExonsLow - cts[N < 5 & sbj != i, length(unique(ref))]
# }
# nExonsHigh <- cts[N > 2000, length(unique(ref))]
# getUniqueHigh <- function(i) {
#   nExonsHigh - cts[N > 2000 & sbj != i, length(unique(ref))]
# }
# unqExonRmLow <- sapply(cts[ , unique(sbj)], getUniqueLow)
# saveRDS(unqExonRmLow, file = "realData37/uniqueExonLessThan5.RDS")
# unqExonRmHigh <- sapply(cts[ , unique(sbj)], getUniqueHigh)
# saveRDS(unqExonRmHigh, file = "realData37/uniqueExonGreaterThan2000.RDS")
# cts <- cts[!ref %in% cts[N < 5, unique(ref)]]
# cts[ , p := N/sum(N), by = sbj]
# cts <- cts[ , list(ref, sbj, proj, p, N)]
# saveRDS(cts, "realData37/AllAtLeast5.counts")
## Now remove refs with greater than 2000, but first remove NCG_00211, NCG_00106
## and NCG_00068 as these samples have much greater coverage than others and 
## including them would result in the removal of many more exons ~ 500 vs ~2400
# rmsmpls <- paste0("NCGENES/markdup/NCG_00", c("211", "106", "068"), ".markdup")
# cts <- cts[!sbj %in% rmsmpls]
# cts <- cts[!ref %in% cts[N > 2000, unique(ref)]]
# cts[ , p := N/sum(N), by = sbj]
# saveRDS(cts, "realData37/AllAtLeast5AndLessThan2000.counts")

getAlpha0 <- function(sbj, fname) {
  
  its <- 500
  hm <- file.path("/", "projects", "sequence_analysis", "vol5", "dfiler")
  x <- readRDS(file.path(hm, "CNV", "realData37", fname))
  sbs <- unlist(sbj)
  x <- x[sbj %in% sbs]
  setorder(x, sbj, ref)
  N <- x[ , list(N = sum(N)), by = sbj]
  x <- dcast(x, sbj ~ ref, value.var = "p")
  sbj <- x[ , sbj]
  x[ , sbj := NULL]
  x <- as.matrix(x)
  rownames(x) <- sbj
  dm <- dim(x)
  n <- dm[1]
  p <- dm[2]
  m <- colmeans(x)
  zx <- t(Log(x))
  down <- -sum(m * (rowmeans(zx) - log(m)))
  sa <- 0.5 * (p - 1)/down
  a1 <- sa * m
  gm <- rowsums(zx)
  z <- n * Digamma(sa)
  g <- z - n * Digamma(a1) + gm
  qk <- -n * Trigamma(a1)
  b <- sum(g/qk)/(1/z - sum(1/qk))
  a2 <- a1 - (g - b)/qk
  i <- 1L
  tl <- numeric(length = its)
  a0 <- numeric(length = its)
  ll <- numeric(length = its)
  while (i <= its) {
    a1 <- a2
    z <- n * digamma(sum(a1))
    g <- z - n * Digamma(a1) + gm
    qk <- -n * Trigamma(a1)
    b <- sum(g/qk)/(1/z - sum(1/qk))
    a2 <- a1 - (g - b)/qk
    cat("iteration:", i, "\n")
    tl[i] <- sum(abs(a2 - a1))
    a0[i] <- sum(a2)
    ll[i] <- n*Lgamma(a0[i]) - n*sum(Lgamma(a2)) + sum(zx*(a2 - 1))
    i <- i + 1L
  }
  if (is.null(colnames(x))) {
    names(a2) <- paste("X", 1:p, sep = "")
  }
  else names(a2) <- colnames(x)
  res <- list(a = a2, a0 = a0[its], a0vec = a0, ll = ll, tl = tl, tc = N)
  res
  
}

cts <- readRDS("realData37/AllAtLeast5.counts")

set.seed(5678)
getSbj <- function(y) cts[proj == y, unique(sbj)]
ncg <- getSbj("NCGENES")
pools <- lapply(1:30, function(x) sample(ncg, size = 16))
pools <- c(pools, lapply(c("Pool1", "Pool2", "SMA"), getSbj))

pars <- data.table(sbj = pools, fname = "AllAtLeast5.counts")

slurm_apply(f = getAlpha0, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "getAlpha0",
            slurm_options = list(mem = 20000,
                                 array = sprintf("0-%d", nrow(pars) - 1),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "8-00:00:00"))


cts <- readRDS("realData37/AllAtLeast5AndLessThan2000.counts")

set.seed(5678)
getSbj <- function(y) cts[proj == y, unique(sbj)]
ncg <- getSbj("NCGENES")
pools <- lapply(1:30, function(x) sample(ncg, size = 16))
pools <- c(pools, lapply(c("Pool1", "Pool2", "SMA"), getSbj))

pars <- data.table(sbj = pools, fname = "AllAtLeast5AndLessThan2000.counts")

slurm_apply(f = getAlpha0, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "getAlpha0_new",
            slurm_options = list(mem = 20000,
                                 array = sprintf("0-%d", nrow(pars) - 1),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "8-00:00:00"))



