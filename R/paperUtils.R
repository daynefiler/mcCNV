#' @import dlfUtils

.calcMCC <- function(tp, tn, fp, fn) {
  tp <- as.numeric(tp)
  tn <- as.numeric(tn)
  fp <- as.numeric(fp)
  fn <- as.numeric(fn)
  (tp*tn - fp*fn)/(sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)))
}

.procRes <- function(x) {
  if (!"width" %in% names(x)) x[ , width := NA]
  if (!"dep"   %in% names(x)) x[ ,   dep := NA]
  x[ , fdep := as.factor(dep)]
  x[ , fwid := as.factor(width)]
  x[ , N := as.numeric(N)]
  xSmry <- x[ , 
              list(tp = sum(N[ PRO &  C1]),
                   fp = sum(N[!ACT &  C1 & !PRO]),
                   tn = sum(N[!ACT & !C1]),
                   fn = sum(N[ ACT & !C1])),
              by = list(prior, fdep, fwid, rep)]
  xSmry[ , fpr := fp/(fp + tn)]
  xSmry[ , tpr := tp/(tp + fn)]
  xSmry[ , spc := tn/(fp + tn)]
  xSmry[ , pf  := fp/(fp + tp)]
  xSmry[ , mcc := .calcMCC(tp, tn, fp, fn)]
  xMn <- xSmry[ , lapply(.SD, mean), by = list(prior, fdep, fwid)]
  xSD <- xSmry[ , lapply(.SD, sd),   by = list(prior, fdep, fwid)]
  list(xMn, xSD)
}

.pltStat <- function(mndat, sddat, stat, ylab, h) {
  plot.new()
  plot.window(ylim = c(0, 1), xlim = c(5, 100))
  abline(h = h, lty = "dashed", col = "darkgrey")
  abline(v = seq(5, 100, 5), lty = "dotted", col = "lightgrey")
  for (i in 1:length(unique(mndat$fwid))) {
    fdep <- mndat[order(fdep)][fwid == i, as.numeric(as.character(fdep))]
    sval <- mndat[order(fdep)][fwid == i, get(stat)]
    SD   <- sddat[order(fdep)][fwid == i, get(stat)]
    blues <- tableau_seq_gradient_pal(palette = "Blue")(1:5/5)
    cols <- c(blues, "#fc7d0b")[i]
    lines(fdep, sval, col = cols, lwd = 2)
    polygon(x = c(fdep, rev(fdep)), 
            y = c(sval + 3*SD, rev(sval - 3*SD)), 
            col = col2alpha(cols, alpha = 0.25), 
            border = NA)
  }
  C <- seq(40, 250, 30)
  axis(side = 1, at = C*76002546/200E6, labels = paste0(C, "x"), 
       tcl = -0.8, mgp = c(3, 2, 0), col = "gray50", col.axis = "gray50")
  axis(side = 1, at = seq(10, 100, by = 10), tcl = -0.4)
  axis(side = 2)
  mtext("Depth (millions of mapped read pairs; approx. coverage)", 
        side = 1, line = 4)
  mtext(ylab, side = 2, line = 3)
}


.pltMnVrCont <- function(sub, nbin = 100, col1 = "darkblue", 
                         col2 = "darkorange", ncontour = 20, 
                         grpVec = sub$cap == "ic",
                         grp1 = "IC", grp2 = "MC",
                         useHeavy = FALSE) {
  mvr <- with(sub[!is.na(vr)], 
              c(quantile(mn, c(0.001, 0.999)), 
                quantile(vr, c(0.001, 0.999))))
  d1 <- with(sub[(grpVec) & !is.na(vr)], 
             log2dDen(mn, vr, lims = mvr, nbin = nbin))
  d2 <- with(sub[!(grpVec) & !is.na(vr)], 
             log2dDen(mn, vr, lims = mvr, nbin = nbin))
  aseq <- seq(0.6, 1, length.out = ncontour)
  c1pal <- sapply(col2alpha(col1, aseq), alpha2opaque)
  c2pal <- sapply(col2alpha(col2, aseq), alpha2opaque)
  layout(matrix(c(2, 1, 4, 3), ncol = 2), heights = c(1, 5), widths = c(5, 1))
  par(mar = c(4, 4, 0.1, 0.1))
  contour(d1, col = c1pal, nlevels = ncontour, axes = FALSE, drawlabels = FALSE)
  contour(d2, add = TRUE, col = c2pal, nlevels = ncontour, drawlabels = FALSE)
  par(usr = log10(mvr))
  d1lm <- with(sub[(grpVec) & !is.na(vr)], lm(log10(vr) ~ log10(mn)))
  abline(d1lm, col = col1, lwd = 2, lty = "dotted")
  d2lm <- with(sub[!(grpVec) & !is.na(vr)], lm(log10(vr) ~ log10(mn)))
  abline(d2lm, col = col2, lwd = 2, lty = "dotted")
  xticks <- axisTicks(usr = log10(mvr[1:2]), log = TRUE)
  yticks <- axisTicks(usr = log10(mvr[3:4]), log = TRUE)
  axis(side = 1, at = log10(xticks), labels = xticks)
  axis(side = 2, at = log10(yticks), labels = yticks)
  title(xlab = "Mean", ylab = "Variance")
  par(mar = c(0.1, 4, 0.1, 0.1))
  abline(0, 1, lty = "dashed", col = "darkgray")
  plot.new()
  mnDen1 <- density(sub[(grpVec) & !is.na(vr), log10(mn)], 
                    from = log10(mvr[1]), to = log10(mvr[2]))
  mnDen2 <- density(sub[!(grpVec) & !is.na(vr), log10(mn)], 
                    from = log10(mvr[1]), to = log10(mvr[2]))
  par(usr = c(log10(mvr)[1:2], -0.1, max(c(mnDen1$y, mnDen2$y))*1.04))
  lines(mnDen1, col = col1, lwd = 2)
  lines(mnDen2, col = col2, lwd = 2)
  par(mar = c(4, 0.1, 0.1, 0.1))
  plot.new()
  vrDen1 <- density(sub[(grpVec) & !is.na(vr), log10(vr)], 
                    from = log10(mvr[3]), to = log10(mvr[4]))
  vrDen2 <- density(sub[!(grpVec) & !is.na(vr), log10(vr)], 
                    from = log10(mvr[3]), to = log10(mvr[4]))
  par(usr = c(-0.1, max(c(vrDen1$y, vrDen2$y))*1.04, log10(mvr)[3:4]))
  lines(x = vrDen1$y, y = vrDen1$x, col = col1, lwd = 2)
  lines(x = vrDen2$y, y = vrDen2$x, col = col2, lwd = 2)
  par(mar = rep(0, 4))
  plot.new()
  legend(x = "center", 
         lwd = 4, 
         col = c(col1, col2), 
         legend = c(grp1, grp2), 
         bty = "n")
}

getQuantSub <- function(sub) {
  t1 <- !is.na(sub$vr)
  t2 <- sub[ , vr < quantile(vr, 0.999, na.rm = TRUE)]
  t3 <- sub[ , vr > quantile(vr, 0.001, na.rm = TRUE)]
  t4 <- sub[ , mn < quantile(mn, 0.999, na.rm = TRUE)]
  t5 <- sub[ , mn > quantile(mn, 0.001, na.rm = TRUE)]
  t1 & t2 & t3 & t4 & t5
}

# ncgMnVr <- readRDS(file.path("..", "Data", "ncgMeanVar.results"))
# ncgSub <- c(sprintf("NCGENES-%d", 10:12), "Pool1", "Pool2", "SMA")
# ncgMnVr <- ncgMnVr[proj %in% ncgICSub | ic == 2]
# pdf(paste0(prefix, "ncgMnVr.pdf"), width = 4.5, height = 3, pointsize = 10)
# pltMnVrCont(ncgMnVr, nbin = 100)
# graphics.off()
# 
# .pltMnVrCont(mv[grepl("IDT", proj)], 50, "blue", "orange", 20)
# 
# 
# subNcgIdt <- c("WGS", "IDT-MC")
# subNcgIdt <- c("NCGENES-10", "IDT-IC")
# mvSub <- mv[proj %in% subNcgIdt]
# mvSub[grepl("IDT", proj), cap := "mc"]
# .pltMnVrCont(mvSub, 50, "blue", "orange", 20)
# 
# cts <- melt(cts, 
#             id = c("ref", "GC", "probeCount", "chr", "start", "end", "wid"), 
#             value.name = "N", variable.name = "sbj")
# 
# metaCols <- c("ref", "GC", "probeCount", "chr", "start", "end", "wid")
# smplCts <- cts[ , 
#                lapply(.SD, sum, na.rm = TRUE), 
#                .SDcols = setdiff(names(cts), metaCols)]
# smplCts <- melt(smplCts, 
#                 variable.name = "projSmpl", 
#                 value.name = "N", 
#                 measure.vars = names(smplCts))
# smplCts <- separate(smplCts, "projSmpl", 
#                     c("proj", "smpl"), 
#                     sep = "_", 
#                     extra = "merge", 
#                     remove = FALSE)
# smplCts <- smplCts[!grepl("oto", proj)]
# bs <- beeswarm(N ~ as.factor(proj), data = smplCts, cex = 0.8, do.plot = FALSE)
# smplCts[ , c("x", "y") := bs[ , c("x", "y")]]
# smplCts[ , col := ifelse(grepl("IDT", proj), "darkblue", "darkorange")]
# smplCts[ , pch := ifelse(proj %in% c("IDT-IC", "NCGENES"), 16, 17)]
# pdf("readCounts.pdf", width = 10, height = 6)
# par(mar = c(2, 4, 1, 1) + 0.1)
# plot(y ~ x, data = smplCts, 
#      ann = FALSE, 
#      xaxt = "n", 
#      bty = "n", 
#      pch = pch, 
#      col = col, 
#      cex = 0.8)
# boxplot(N ~ as.factor(proj), 
#         data = smplCts, 
#         frame = FALSE, 
#         add = TRUE,
#         ann = FALSE,
#         axes = FALSE,
#         col = "transparent",
#         border = "gray30",
#         outline = FALSE,
#         boxwex = 0.4)
# axis(side = 1, at = 1:7, tick = FALSE, labels = unique(smplCts$proj))
# title(ylab = "Read pairs per sample")
# legend("topright", 
#        legend = c("Independent captures", "Mulitplexed capture", 
#                   "IDT Human Exome", "Agilent SureSelect All Exon"), 
#        pch = c(16, 18, NA, NA), 
#        fill = c(NA, NA, "darkblue", "darkorange"), 
#        bty = "n", 
#        border = NA)
# graphics.off()



