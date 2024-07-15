###==========================================================================###
# This function is a modified version of mdsPlot() from minfi                  #
# Modified by Roy Lardenoije                                                   #
###==========================================================================###

# Modified mdsPlot function from minfi to have the legend outside the plot 
# area
mdsPlot2 <- function (dat, numPositions = 1000, PCs = c(1,2), sampNames = NULL, sampGroups = NULL, xlim, ylim, 
                      pch = 19, pal = brewer.pal(8, "Dark2"), legendPos = legendPos, 
                      legendNCol = 2, legendTitle = "", inset = c(-0.3, 0), main = NULL){
  if (is(dat, "MethylSet") || is(dat, "RGChannelSet")) {
    b <- getBeta(dat)
  }
  else if (is(dat, "matrix")) {
    b <- dat
  }
  else {
    stop("dat must be an 'MethylSet', 'RGChannelSet', or 'matrix'.")
  }
  if (is.null(main)) {
    main <- sprintf("Beta MDS\n%d most variable positions", numPositions)
  }
  o <- order(rowVars(b), decreasing = TRUE)[seq_len(numPositions)]
  d <- dist(t(b[o, ]))
  fit <- cmdscale(d, k = max(PCs))
  if (missing(xlim)) 
    xlim <- range(fit[, PCs[1]]) * 1.2
  if (missing(ylim)) 
    ylim <- range(fit[, PCs[2]]) * 1.2
  if (is.null(sampGroups)) 
    sampGroups <- rep(1, numPositions)
  if(is.logical(sampGroups)){
    sampGroups <- factor(sampGroups, levels = c("TRUE", "FALSE"))
  } else {
    sampGroups <- as.factor(sampGroups)
  }
  col <- pal[sampGroups]
  par(mar = c(5.1, 4.1, 4.1, 8.1))
  if (is.null(sampNames)) {
    plot(x = fit[, PCs[1]], y = fit[, PCs[2]], col = col, pch = pch, xlim = xlim, ylim = ylim, 
         xlab = paste0("PC", PCs[1]), ylab = paste0("PC", PCs[2]), main = main)
  }
  else {
    plot(x = 0, y = 0, type = "n", xlim = xlim, ylim = ylim, xlab = paste0("PC", PCs[1]), ylab = paste0("PC", PCs[2]), 
         main = main)
    text(x = fit[, PCs[1]], y = fit[, PCs[2]], sampNames, col = col)
  }
  numGroups <- length(unique(sampGroups))
  # if (missing(legendNCol)) 
  #     legendNCol <- numGroups
  if (numGroups > 1) {
    legend(x = legendPos, legend = levels(sampGroups), ncol = legendNCol, 
           text.col = pal[seq_len(numGroups)], inset = inset, xpd = TRUE, cex = 0.75, 
           title = legendTitle)
  }
}