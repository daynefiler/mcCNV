.line2user <- function(line, side, outer = FALSE) {
  unit <- if (outer) "nic" else "npc"
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(c(0, lh), 'inches', unit))
  y_off <- diff(grconvertY(c(0, lh), 'inches', unit))
  switch(side,
         `1` = grconvertY(-line * y_off, unit, 'user'),
         `2` = grconvertX(-line * x_off, unit, 'user'),
         `3` = grconvertY(1 + line * y_off, unit, 'user'),
         `4` = grconvertX(1 + line * x_off, unit, 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

.col2alpha <- function(col, alpha = 0.5) {
  tmp <- col2rgb(col)
  rgb(tmp[1]/255, tmp[2]/255, tmp[3]/255, alpha)
}

.calcMCC <- function(tp, tn, fp, fn) {
  tp <- as.numeric(tp)
  tn <- as.numeric(tn)
  fp <- as.numeric(fp)
  fn <- as.numeric(fn)
  (tp*tn - fp*fn)/(sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)))
}

.procRes <- function(x) {
  if (!"width" %in% names(x)) x[ , width := 1L]
  x[ , fdep := as.factor(dep)]
  x[ , fwid := as.factor(width)]
  x[ , N := as.numeric(N)]
  xSmry <- x[ , 
              list(tp = sum(N[ ACT &  C1]),
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
    cols <- c(brewer.pal(9, "Blues")[4:8], "darkorange")[i]
    lines(fdep, sval, col = cols, lwd = 2)
    polygon(x = c(fdep, rev(fdep)), 
            y = c(sval + 3*SD, rev(sval - 3*SD)), 
            col = .col2alpha(cols, alpha = 0.25), 
            border = NA)
  }
  axis(side = 1, at = seq(10, 100, by = 10))
  axis(side = 2)
  mtext("Depth", side = 1, line = 3)
  mtext(ylab, side = 2, line = 3)
}

#-------------------------------------------------------------------------------
# addfiglab: Add a figure label
#-------------------------------------------------------------------------------

#' @name addfiglab
#' @title Add a figure label
#' @description Add a figure label based on x and y line values
#'
#' @export

.addfiglab <- function(lab, xl = par()$mar[2], yl = par()$mar[3], cex = 1) {
  
  text(x = .line2user(xl, 2), y = .line2user(yl, 3), 
       lab, xpd = NA, font = 2, cex = cex, adj = c(0, 1))
  
}

#-------------------------------------------------------------------------------


