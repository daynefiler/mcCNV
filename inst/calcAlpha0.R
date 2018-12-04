##----------------------------------------------------------------------------##
## Script to fit the dirichlet distribution
##----------------------------------------------------------------------------##

library(data.table)
library(rslurm)
library(Rfast)

# cts <- readRDS("realData37/AllSamples.counts")
# cts[ , proj := dirname(dirname(sbj))]
# saveRDS(cts[N < 5, .N, by = list(sbj, N)], file = "realData37/lowCounts.RDS")
# nExons <- cts[N < 5, length(unique(ref))]
# getUniqueLow <- function(i) {
#   nExons - cts[N < 5 & sbj != i, length(unique(ref))]
# }
# unqExonRm <- sapply(cts[ , unique(sbj)], getUniqueLow)
# saveRDS(unqExonRm, file = "realData37/uniqueExonLessThan5.RDS")
# cts <- cts[!ref %in% cts[N < 5, unique(ref)]]
# cts[ , p := N/sum(N), by = sbj]
# cts <- cts[ , list(ref, sbj, proj, p, N)]
# saveRDS(cts, "realData37/AllAtLeast5.counts")
cts <- readRDS("realData37/AllAtLeast5.counts")

set.seed(5678)
getSbj <- function(y) cts[proj == y, unique(sbj)]
ncg <- getSbj("NCGENES")
pools <- lapply(1:30, function(x) sample(ncg, size = 16))
pools <- c(pools, lapply(c("Pool1", "Pool2", "SMA"), getSbj))

pars <- data.table(sbj = pools)

getAlpha0 <- function(sbj) {
  
  its <- 500
  hm <- file.path("/", "projects", "sequence_analysis", "vol5", "dfiler")
  x <- readRDS(file.path(hm, "CNV", "realData37", "AllAtLeast5.counts"))
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





