###==========================================================================###
# This function is a modified version of predictSex() from wateRmelon          #
# Modified by Roy Lardenoije                                                   #
###==========================================================================###

# This modified version of predictSex() also makes a plot colored by reported sex
predictSex2 <- function (x, x.probes = NULL, pc = 2, plot = TRUE, irlba = TRUE, 
                         center = FALSE, scale. = FALSE, reported.sex, sampNames = NULL, 
                         axis.scale = 1.2){
    if (!is.null(x.probes)) {
        x <- x[x.probes, ]  
    }
    x <- na.omit(x)
    message(paste0("There are ", nrow(x), " X chromosome probes without NA values that will be used."))
    cm <- colMedians(x)
    message("Performing PCA")
    if (irlba & requireNamespace("irlba")) {
        if(nrow(x) < 3) {
            stop("There are less than 3 X chromosome probes without NA values! More are required by irlba.")
        }
        pca <- irlba::prcomp_irlba(x, n = pc, center = center, retx = F, scale. = scale.)$rot
    }
    else {
        pca <- prcomp(x, center = center, retx = F, scale. = scale.)$rot[, seq_len(pc)]
    }
    if (any(table(cm >= 0.48) >= (ncol(x) * 0.95))) {
        message("Data likely consists of single gender")
    }
    message("Generating Clusters")
    c1 <- kmeans(cm, centers = fivenum(cm)[c(2, 4)])
    c2 <- kmeans(pca[, pc], centers = 2)
    n <- 0
    while (sum((s <- c1$cluster + c2$cluster) == 3) > (length(c2$cluster)/2)) {
        c2 <- kmeans(pca[, pc], centers = 2)
        n <- n + 1
        if (n > 10) 
            stop("Cannot find meaningful clusters for data")
    }
    is.na(s[s == 3]) <- TRUE
    o <- order(t <- tapply(cm, s, mean))
    s[s == names(t)[o[1]]] <- "Male"
    s[s == names(t)[o[2]]] <- "Female"
    s[is.na(s)] <- "Undefined"
    if (plot) {
        col <- rainbow(nlevels(factor(s)))
        xlim <- range(pca[, 1]) * axis.scale    
        ylim <- range(pca[, pc]) * axis.scale
        par(mar = c(4.5, 4.5, 1, 8), mfrow = c(2, 1))
        if (is.null(sampNames)) {
            plot(x = pca[, 1], y = pca[, pc], col = col[factor(s)], xlab = "PC 1", 
                 ylab = paste("PC", pc))
            legend("right", legend = levels(factor(s)), col = col, 
                   pch = 1, inset = c(-0.3, 0), xpd = T, title = "Predicted sex")
            plot(x = pca[, 1], y = pca[, pc], col = col[factor(reported.sex)], xlab = "PC 1", 
                 ylab = paste("PC", pc))
            legend("right", legend = levels(factor(reported.sex)), col = col, 
                   pch = 1, inset = c(-0.3, 0), xpd = T, title = "Reported sex")
        }
        else {
            plot(x = 0, y = 0, type = "n", xlab = "PC 1", ylab = paste("PC", pc), xlim = xlim, ylim = ylim)
            text(x = pca[, 1], y = pca[, pc], sampNames, col = col[factor(s)], cex = 0.5)
            legend("right", legend = levels(factor(s)), col = col, 
                   pch = 1, inset = c(-0.3, 0), xpd = T, title = "Predicted sex")
            plot(x = 0, y = 0, type = "n", xlab = "PC 1", ylab = paste("PC", pc), xlim = xlim, ylim = ylim)
            text(x = pca[, 1], y = pca[, pc], sampNames, col = col[factor(reported.sex)], cex = 0.5)
            legend("right", legend = levels(factor(reported.sex)), col = col, 
                   pch = 1, inset = c(-0.3, 0), xpd = T, title = "Reported sex")
        }
    }
    return(s)
}