expectedHists <- function(expDists, binSize, obsDist) {

    # find the breaks according to bin size
    expMins <- sapply(append(expDists, list(obsDist)), function(.dt) min(.dt$PCC))
    expMaxs <- sapply(append(expDists, list(obsDist)), function(.dt) max(.dt$PCC))
    minpcc <- floor(min(expMins) / binSize) * binSize
    maxpcc <- ceiling(max(expMaxs) / binSize) * binSize
    brks <- seq(minpcc, maxpcc, binSize)

    H <- lapply(expDists, function(d) plotrix::weighted.hist(d$PCC, d$transitProb, brks, plot = FALSE))
    H <- lapply(H, function(h) {
        h$density <- h$counts / sum(h$counts)
        class(h) <- "histogram"
        h$breaks <- brks
        h
    })

    return(H)

}

#' @export
plotModelsHist <- function(rootDirs, binSize, xlim, plot = TRUE) {

    # load expected distributions from each model
    expDistPaths <- stringi::stri_join(rootDirs, "expDistribution.tsv")
    expDists <- lapply(expDistPaths, function(p) setNames(data.table::fread(p, header = FALSE), c("PCC", "transitProb")))

    # load the observed distribution (assumes it is the same for all models, so we use the first rootDir)
    obsDistPath <- stringi::stri_join(rootDirs[1], "obsDistribution.txt")
    obsDist <- setNames(data.table::fread(obsDistPath, header = FALSE), "PCC")

    # compute expected histograms for each model
    H <- setNames(expectedHists(expDists, binSize, obsDist), basename(rootDirs))
    brks <- H[[1]]$breaks

    # compute observed histogram
    H$Observed <- hist(obsDist$PCC, breaks = brks, plot = FALSE)
    H$Observed$density <- H$Observed$counts / sum(H$Observed$counts)

    # compute a naive histogram of expected with no weights (uniform), it is assumed the expected PCCs for all models are the same, so we just take the first dist
    H$Uniform <- hist(expDists[[1]]$PCC, breaks = brks, plot = FALSE)
    H$Uniform$density <- H$Uniform$counts / sum(H$Uniform$counts)

    # plot 
    if (missing(xlim)) xlim <- c(min(brks), max(brks))
    if (plot) plotHists(H, xlim)

    return(H)

}

#' @export
plotHists <- function(hists, xlim) {

    brks <- hists[[1]]$breaks
    if (missing(xlim)) xlim <- c(min(brks), max(brks))

    binIsXlim <- sapply(1:(length(brks) - 1), function(i) brks[i] >= xlim[1] & brks[i + 1] <= xlim[2])
    ylim.max <- max(sapply(hists, function(h) max(h$density[binIsXlim])))
    par(las = 2)
    plot(NULL, xlim = xlim, ylim = c(0, ylim.max), xlab = "PCC", ylab = "Density", xaxt = "n")
    axis(side = 1, at = brks, labels = brks)

    for (i in 1:length(hists)) {
        plot(hists[[i]], freq = FALSE, add = TRUE, col = NA, border = i, lty = i + 1)
    }
    legend("topleft", legend = names(hists), col = 1:length(hists), lty = (1:length(hists)) + 1)

}

#' @export
plotEnsemblsHist <- function(rootDirsList, binSize, xlim, plot = TRUE) {

    # load expected distributions from each ensembl of models
    expDistPathsList <- lapply(rootDirsList, function(rootDirs) stringi::stri_join(rootDirs, "expDistribution.tsv"))
    expEnsemblsDist <- lapply(expDistPathsList, function(expDistPaths) {
        expDists <- lapply(expDistPaths, function(p) setNames(data.table::fread(p, header = FALSE), c("PCC", "transitProb")))
        data.table::rbindlist(expDists)
    })

    # load the observed distribution of the ensembls (assumes it is the same for all ensembls, so we use the first ensembl)
    obsDistPaths <- stringi::stri_join(rootDirsList[[1]], "obsDistribution.txt")
    obsEnsemblDist <- lapply(obsDistPaths, function(p) setNames(data.table::fread(p, header = FALSE), "PCC"))
    obsEnsemblDist <- data.table::rbindlist(obsEnsemblDist)

    # compute expected histograms for each ensembl of models
    H <- setNames(expectedHists(expEnsemblsDist, binSize, obsEnsemblDist), names(rootDirsList))
    brks <- H[[1]]$breaks

    # compute observed histogram
    H$Observed <- hist(obsEnsemblDist$PCC, breaks = brks, plot = FALSE)
    H$Observed$density <- H$Observed$counts / sum(H$Observed$counts)

    # compute a naive histogram of expected with no weights (uniform), it is assumed the expected PCCs for all ensembls are the same, so we just take the first dist
    H$Uniform <- hist(expEnsemblsDist[[1]]$PCC, breaks = brks, plot = FALSE)
    H$Uniform$density <- H$Uniform$counts / sum(H$Uniform$counts)

    # plot
    if (missing(xlim)) xlim <- c(min(brks), max(brks))
    if (plot) plotHists(H, xlim)

    return(H)

}

#' @export
evaluateModels <- function(rootDirs, binSize) {

    H <- plotModelsHist(rootDirs, binSize, plot = FALSE)
    
    .o <- H$Observed$density
    H$Observed <- NULL
    .E <- lapply(H, function(h) h$density)

    return(sapply(.E, cosineSimilarity, .o))

}

cosineSimilarity <- function(x, y) {

    return(sum(x * y) / (sqrt(sum(x ^ 2)) * sqrt(sum(y ^ 2))))

}

#' @export
predictionCurves <- function(rootdir, gtSet, plot = TRUE, ids = NULL) {

    # load results for the rootdir
    resultsPath <- stringi::stri_join(rootdir, "results.tsv")
    results <- data.table::fread(resultsPath)

    # if ids provided for each genomic element, use that as id column
    if (!is.null(ids)) results$Id <- ids

    # remove missing values
    results <- results[!is.na(results$Pvalue)]

    # determine data labels according to ground truth
    results$Label <- ifelse(results$Id %in% gtSet, 1, 0)

    # compute prediction curves
    pcurveSet <- list()
    pcurveSet$obsROC <- PRROC::roc.curve(scores.class0 = -results$Pvalue, weights.class0 = results$Label, 
        curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE
    )
    pcurveSet$obsPRC <- PRROC::pr.curve(scores.class0 = -results$Pvalue, weights.class0 = results$Label,
        curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE
    )

    # plot if specified
    if (plot) {

        title <- basename(rootdir)
        par(mfrow = c(1, 2), pty = "s")
        plotPRROCcurves(pcurveSet["obsROC"], title)
        plotPRROCcurves(pcurveSet["obsPRC"], title)

    }

    return(pcurveSet)

}

#' @export
plotPRROCcurves = function(curves, title) {

    n <- length(curves)
    cols <- 1:n
    auc <- sapply(curves, function(.c) ifelse(.c$type == "ROC", .c$auc, .c$auc.integral))
    lgnd <- paste(names(curves), signif(auc, 3))

    plot(curves[[1]], rand.plot = TRUE, col = cols[1], auc.main = FALSE, main = title, lwd = 1)
    curves[[1]] <- NULL
    cols <- cols[-1]

    while (length(curves) > 0) {

        plot(curves[[1]], add = TRUE, rand.plot = TRUE, col = cols[1], lwd = 1)
        curves[[1]] <- NULL
        cols <- cols[-1]

    }
    legend("bottomright", legend = lgnd, lty = 1, lwd = 2, col = 1:n, cex = 0.7, bty = "n")

}
