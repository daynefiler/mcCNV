##----------------------------------------------------------------------------##
## Script to fit the dirichlet distribution
##----------------------------------------------------------------------------##

library(data.table)
library(rslurm)
library(Rfast)

getAlpha0 <- function(sbj, fname) {
  
  its <- 2500
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

# cts <- readRDS("realData37/otoAll.counts")
# cts[ , proj := dirname(dirname(sbj))]
# saveRDS(cts[N < 5, list(exonN = .N), by = list(proj, sbj, readN = N)],
#         file = "realData37/otolowCounts.RDS")
# saveRDS(cts[N > 2000, .N, by = list(proj, sbj)],
#         file = "realData37/otohighCounts.RDS")
# ## Using 2 as the low cutoff, due to the decreased number of exons
# nExonsLow <- cts[N < 2, length(unique(ref))]
# getUniqueLow <- function(i) {
#   nExonsLow - cts[N < 2 & sbj != i, length(unique(ref))]
# }
# nExonsHigh <- cts[N > 2000, length(unique(ref))]
# getUniqueHigh <- function(i) {
#   nExonsHigh - cts[N > 2000 & sbj != i, length(unique(ref))]
# }
# unqExonRmLow <- sapply(cts[ , unique(sbj)], getUniqueLow)
# saveRDS(unqExonRmLow, file = "realData37/otoUniqueExonLessThan5.RDS")
# unqExonRmHigh <- sapply(cts[ , unique(sbj)], getUniqueHigh)
# saveRDS(unqExonRmHigh, file = "realData37/uniqueExonGreaterThan2000.RDS")
# 
# ## Remove sample otoIC/markdup/C047_B096.markdup due to very low reads and 
# ## 1126 unique exons with less than 2 reads
# cts <- cts[sbj != "otoIC/markdup/C047_B096.markdup"]
# 
# cts <- cts[!ref %in% cts[N < 2, unique(ref)]]
# cts[ , p := N/sum(N), by = sbj]
# cts <- cts[ , list(ref, sbj, proj, p, N)]
# saveRDS(cts, "realData37/otoAllAtLeast2.counts")
# cts <- cts[!ref %in% cts[N > 2000, unique(ref)]]
# cts[ , p := N/sum(N), by = sbj]
# saveRDS(cts, "realData37/otoAllAtLeast2AndLessThan2000.counts")

cts <- readRDS("realData37/otoAllAtLeast2.counts")

set.seed(5678)
getSbj <- function(y) cts[proj == y, unique(sbj)]
ic <- getSbj("otoIC")
pools <- lapply(1:30, function(x) sample(ic, size = 16))
pools <- c(pools, lapply(1:30, function(x) sample(ic, size = 12)))
pools <- c(pools, lapply(cts[proj != "otoIC", unique(proj)], getSbj))

pars <- data.table(sbj = pools, fname = "otoAllAtLeast2.counts")

slurm_apply(f = getAlpha0, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "getAlpha0-oto",
            slurm_options = list(mem = 20000,
                                 array = sprintf("0-%d", nrow(pars) - 1),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "8-00:00:00"))

cts <- readRDS("realData37/otoAllAtLeast2AndLessThan2000.counts")

set.seed(5678)
getSbj <- function(y) cts[proj == y, unique(sbj)]
ic <- getSbj("otoIC")
pools <- lapply(1:30, function(x) sample(ic, size = 16))
pools <- c(pools, lapply(1:30, function(x) sample(ic, size = 12)))
pools <- c(pools, lapply(cts[proj != "otoIC", unique(proj)], getSbj))

pars <- data.table(sbj = pools, fname = "otoAllAtLeast2AndLessThan2000.counts")

slurm_apply(f = getAlpha0, 
            params = pars, 
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "getAlpha0-oto_new",
            slurm_options = list(mem = 20000,
                                 array = sprintf("0-%d", nrow(pars) - 1),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "8-00:00:00"))


# library(rslurm)
# res1 <- get_slurm_out(slurm_job("getAlpha0oto", 74))
# saveRDS(res1, file = "results_alpha0_oto.RDS")
# res2 <- get_slurm_out(slurm_job("getAlpha0oto_new", 74))
# saveRDS(res2, file = "results_alpha0_oto_new.RDS")
