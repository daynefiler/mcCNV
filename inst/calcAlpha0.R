library(data.table)
library(rslurm)
library(Rfast)

# cts <- readRDS("realData37/AllSubjects.counts")
# cts[ , proj := dirname(dirname(sbj))]
# cts <- cts[!ref %in% cts[N == 0, unique(ref)]]
# cts[ , p := N/sum(N), by = sbj]
# cts <- cts[ , list(ref, sbj, proj, p)]
# saveRDS(cts, "realData37/AllNonZero.counts")
cts <- readRDS("realData37/AllNonZero.counts")

set.seed(1234)
getSbj <- function(y) cts[proj == y, unique(sbj)]
ncg <- getSbj("NCGENES")
pools <- lapply(1:30, function(x) sample(ncg, size = 16))
pools <- c(pools, lapply(c("Pool1", "Pool2", "SMA"), getSbj))

pars <- data.table(sbj = pools)

getAlpha0 <- function(sbj) {
  
  hm <- file.path("/", "projects", "sequence_analysis", "vol5", "dfiler")
  cts <- readRDS(file.path(hm, "CNV", "realData37", "AllNonZero.counts"))
  sbs <- unlist(sbj)
  cts <- cts[sbj %in% unlist(sbs)]
  setorder(cts, sbj, ref)
  cts <- dcast(cts, sbj ~ ref, value.var = "p")
  sbj <- cts[ , sbj]
  cts[ , sbj := NULL]
  cts <- as.matrix(cts)
  rownames(cts) <- sbj
  fit <- diri.nr2(cts)
  fit
  
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









