###==============================================================================================###
# This function is a modified version of champ.QC()                                                #
# Modified by Roy Lardenoije                                                                       #
# Version 1, last updated 2020-06-03                                                               #
###==============================================================================================###

# This function is a modified version of champ.QC(), using mdsPlot2, and a modified dendrogram (e.g. 
# legend outside plot)
custom.QC <- function (beta, pheno, mdsPlot = TRUE, densityPlot = TRUE, dendrogram = TRUE, 
                       PDFplot = TRUE, Rplot = TRUE, Feature.sel = "None", main = "", 
                       leaves_cex = 0.4, labels_cex = 0.6, sampNames = TRUE, legendPos = "topright", 
                       legendTitle = "", legendInset =  c(-0.3, 0), legendNCol = 2, resultsDir = "./QCimages/", density_xlab = ""){
    message("[===========================]")
    message("[<<<<< custom.QC START >>>>>>]")
    message("-----------------------------")
    require(minfi)
    require(dendextend)
    
    # Modified mdsPlot function from minfi to have the legend outside the plot area
    mdsPlot2 <- function (dat, numPositions = 1000, sampNames = NULL, sampGroups = NULL, xlim, ylim, 
                          pch = 19, pal = rainbow(length(unique(sampGroups)), alpha = 0.4), legendPos = legendPos, 
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
        fit <- cmdscale(d)
        if (missing(xlim)) 
            xlim <- range(fit[, 1]) * 1.2
        if (missing(ylim)) 
            ylim <- range(fit[, 2]) * 1.2
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
            plot(x = fit[, 1], y = fit[, 2], col = col, pch = pch, xlim = xlim, ylim = ylim, 
                 xlab = "PC1", ylab = "PC2", main = main)
        }
        else {
            plot(x = 0, y = 0, type = "n", xlim = xlim, ylim = ylim, xlab = "PC1", ylab = "PC2", 
                 main = main)
            text(x = fit[, 1], y = fit[, 2], sampNames, col = col)
        }
        numGroups <- length(unique(sampGroups))
        # if (missing(legendNCol)) 
        #     legendNCol <- numGroups
        if (numGroups > 1) {
            legend(x = legendPos, legend = levels(sampGroups), ncol = legendNCol, pch = 19,
                   col = pal[seq_len(numGroups)], inset = inset, xpd = TRUE, cex = 0.75, 
                   title = legendTitle)
        }
    }
    
    # Modified densityPlot function from minfi to have more control of the legend
    densityPlot2 <- function (dat, sampGroups = NULL, main = "", xlab = "BETA", 
                              pal = rainbow(length(unique(sampGroups)), alpha = 0.4), xlim, ylim, add = TRUE, legend = TRUE, 
                              legendTitle = "", ...) 
    {
      if (is(dat, "RGChannelSet") || is(dat, "MethylSet")) {
        b <- getBeta(dat)
      }
      else if (is(dat, "matrix")) {
        b <- dat
      }
      else {
        stop("argument 'dat' must be an 'RGChannelSet', a 'MethylSet' or ", 
             "matrix.")
      }
      d <- apply(b, 2, function(x) density(as.vector(x), na.rm = TRUE))
      if (missing(ylim)) 
        ylim <- range(sapply(d, function(i) range(i$y)))
      if (missing(xlim)) 
        xlim <- range(sapply(d, function(i) range(i$x)))
      if (is.null(sampGroups)) {
        sampGroups <- rep(1, ncol(b))
      }
      else if (length(sampGroups) == 1) {
        sampGroups <- rep(sampGroups, ncol(b))
      }
      if(is.logical(sampGroups)){
        sampGroups <- factor(sampGroups, levels = c("TRUE", "FALSE"))
      } else {
        sampGroups <- as.factor(sampGroups)
      }
      if (add) {
        plot(x = 0, type = "n", ylim = ylim, xlim = xlim, ylab = "Density", 
             xlab = xlab, main = main, ...)
        abline(h = 0, col = "grey80")
      }
      for (i in seq_along(d)) {
        lines(d[[i]], col = pal[sampGroups[i]])
      }
      if (legend & length(levels(sampGroups)) > 1) {
        legend("topright", legend = levels(sampGroups), col = pal[seq_len(numGroups)], fill = pal, 
               title = legendTitle)
      }
    }
    if (!file.exists(resultsDir)) 
        dir.create(resultsDir)
    if(PDFplot) message("custom.QC Results will be saved in ", resultsDir)
    message("[QC plots will proceed with ", dim(beta)[1], " probes and ", dim(beta)[2], " samples.]\n")
    if (min(beta, na.rm = TRUE) == 0) {
        beta[beta == 0] <- 1e-06
        message("[", length(which(beta == 0)), " zeros detected in your dataset, will be replaced with 0.000001]\n")
    }
    if (ncol(beta) != length(pheno)) 
        stop("Dimensions of dataset samples, and phenotype samples must match. Please check your input.")
    message("<< Prepare Data. >>")
    if (mdsPlot) {
        if (sampNames) {
            sampNames <- colnames(beta)
        } else {
            sampNames <- NULL
        }
        if (Rplot) 
            mdsPlot2(beta, numPositions = 1000, sampGroups = pheno, sampNames = sampNames, 
                     legendPos = legendPos, legendTitle = legendTitle, inset = legendInset,
                     legendNCol = legendNCol)
        if (PDFplot) {
            pdf(paste(resultsDir, "mdsPlot.pdf", sep = "/"), width = 6, height = 4)
            mdsPlot2(beta, numPositions = 1000, sampGroups = pheno, sampNames = sampNames, 
                     legendPos = legendPos, legendTitle = legendTitle, inset = legendInset,
                     legendNCol = legendNCol)
            dev.off()
        }
        message("<< MDS plot done. >>\n")
    }
    if (densityPlot) {
        if (Rplot) {
          densityPlot2(beta, sampGroups = pheno, 
                      main = paste("Density plot (", nrow(beta), " probes)", sep = ""), 
                      xlab = density_xlab, legendTitle = legendTitle)
        }
        if (PDFplot) {
            pdf(paste(resultsDir, "densityPlot.pdf", sep = "/"), width = 6, height = 4)
            densityPlot2(beta, sampGroups = pheno, 
                        main = paste("Density plot (", nrow(beta), " probes)", sep = ""), 
                        xlab = density_xlab, legendTitle = legendTitle)
            dev.off()
        }
        message("<< Density plot done. >>\n")
    }
    if (dendrogram) {
        if (Feature.sel == "None") {
            message("< Dendrogram feature selection method >: No selection, directly use all CpGs to calculate distance matrix.")
            hc <- hclust(dist(t(beta)))
        }
        else if (Feature.sel == "SVD") {
            message("< Dendrogram feature selection method >: Use top SVD CpGs to calculate distance matrix.")
            SVD <- svd(beta)
            rmt.o <- EstDimRMT(beta - rowMeans(beta))
            M <- SVD$v[, 1:rmt.o$dim]
            rownames(M) <- colnames(beta)
            colnames(M) <- paste("Component", c(1:rmt.o$dim))
            hc <- hclust(dist(M))
        }
        dend <- as.dendrogram(hc)
        if(is.logical(pheno)){
          pheno <- factor(pheno, levels = c("TRUE", "FALSE"))
        }
        MyColor <- rainbow(length(table(pheno)), alpha = 0.4)
        names(MyColor) <- names(table(pheno))
        labels_colors(dend) <- MyColor[as.character(pheno[order.dendrogram(dend)])]
        dend <- dend %>% set("labels_cex", labels_cex)
        dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex", leaves_cex) %>% set("leaves_col", MyColor[as.character(pheno[order.dendrogram(dend)])])
        if (Rplot) {
            plot(dend, center = TRUE, main = paste(main, " (", nrow(beta), " probes)", sep = ""))
            legend("topright", fill = MyColor, legend = names(MyColor), inset = c(-0.05, -0.1), xpd = TRUE, cex = 0.75, title = legendTitle)
        }
        if (PDFplot) {
            pdf(paste(resultsDir, "sampleCluster.pdf", sep = "/"), width = floor(log(ncol(beta)) * 3), height = floor(log(ncol(beta)) * 2))
            plot(dend, center = TRUE, main = paste(main, " (", nrow(beta), " probes)", sep = ""))
            legend("topright", fill = MyColor, legend = names(MyColor), inset = c(-0.05, -0.1), xpd = TRUE, cex = 0.75, title = legendTitle)
            dev.off()
        }
        message("<< Dendrogram plot done. >>\n")
    }
    message("[<<<<<< custom.QC END >>>>>>>]")
    message("[===========================]")
}