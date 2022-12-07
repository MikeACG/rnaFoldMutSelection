#' @import data.table
#' @import GenomicRanges
#' @import arrow

#' @export
modelName <- function(cohort, k, x) {

    rootdir <- paste0(cohort, "_", k, "mer_", paste(x, collapse = "-"))

    return(rootdir)

}

#' @export
makeRoot <- function(rootdir, overwrite = FALSE) {

    rootis <- file.exists(rootdir)
    if (!rootis) {
        dir.create(rootdir)
    } else if (overwrite) {
        unlink(rootdir, recursive = TRUE)
        dir.create(rootdir)
    } else {
        stop(paste0("Directory ", rootdir, " already exists"))
    }

}

#' @export
loadKmerData <- function(sites, kmer, x, kmercol) {

    # load sites with kmer into memory
    cols2get <- c("seqnames", "start", "end", "isminus", x)
    scan_builder <- sites$NewScan()
    scan_builder$Filter(arrow::Expression$field_ref(kmercol) == kmer)
    scan_builder$Project(cols2get)
    scanner <- scan_builder$Finish()
    D <- data.table::data.table(as.data.frame(scanner$ToTable()))

    return(D)

}

#' @export
kmerifyMaf <- function(mafRanges, siteRanges, D, kmer) {

    # get subset of maf that hits sites
    mafSiteOv <- GenomicRanges::findOverlaps(mafRanges, siteRanges, ignore.strand = TRUE)
    m <- mafRanges[S4Vectors::queryHits(mafSiteOv)]
    m$ov <- mafSiteOv

    # assign strand of site
    m$isminus <- D$isminus[S4Vectors::subjectHits(mafSiteOv)]

    # change strand of maf if neccessary
    m$initialState[m$isminus] <- chartr(
        "ACGT",
        "TGCA",
        m$initialState[m$isminus]
    )
    m$finalState[m$isminus] <- chartr(
        "ACGT",
        "TGCA",
        m$finalState[m$isminus]
    )

    # assign mutation type
    if (length(m) > 0) {
        m$mutation <- paste0(kmer, ">", m$finalState)
    } else {
        m$mutation <- character(0)
    }

    return(m)

}

#' @export
fitKmer <- function(D, m, x, kmer, nuc) {

    # fit model for each type of mutation
    R <- list()
    mutTypes <- mtypes(kmer, nuc)
    for (mt in mutTypes) {

        msites <- table(S4Vectors::subjectHits(m$ov[m$mutation == mt]))
        D$mutated <- rep(0L, nrow(D))
        D$mutated[as.integer(names(msites))] <- as.integer(msites)
        n <- sum(D$mutated)
        if (length(x) > 0 & n > 0) {

            R[[mt]] <- simplifyLm(lm(reformulate(x, "mutated"), data = D))

        } else {

            R[[mt]] <- structure(
                list(cellval = n / nrow(D)),
                class = "mutMatrixCell"
            )

        }
        

    }

    return(R)

}

#' @export
sitifyGtf <- function(gtfSites, siteRanges, D, x, kmer) {

    # get target sites in siteRanges
    gtfSiteOv <- GenomicRanges::findOverlaps(gtfSites, siteRanges, ignore.strand = TRUE)
    g <- gtfSites[S4Vectors::queryHits(gtfSiteOv)]

    # add covariates and kmer
    GenomicRanges::mcols(g) <- cbind(
        D[S4Vectors::subjectHits(gtfSiteOv), ..x],
        data.table::as.data.table(GenomicRanges::mcols(g))
    )
    g$kmer <- rep(kmer, length(g))

    return(g)

}

#' @export
gtfyMaf <- function(m, g) {

    # compute overlap between maf and gtf sites
    mgOv <- GenomicRanges::findOverlaps(m, g, ignore.strand = TRUE)

    # if one mutations hits multiple sites, duplicate it
    mIdxs <- S4Vectors::queryHits(mgOv)
    m2 <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(m)[mIdxs],
        ranges = IRanges::IRanges(
            GenomicRanges::start(m)[mIdxs],
            GenomicRanges::end(m)[mIdxs]
        )
    )
    m$ov <- NULL
    GenomicRanges::mcols(m2) <- GenomicRanges::mcols(m)[mIdxs, ]

    # add gtf site id to mutations
    m2$id <- g$id[S4Vectors::subjectHits(mgOv)]

    return(m2)

}

#' @export
possibleKmers <- function(k) {

    flanks <- gtools::permutations(4, k - 1, c("A", "C", "G", "T"), repeats.allowed = TRUE)

    w <- ncol(flanks) / 2 
    left <- apply(flanks[, 1:w, drop = FALSE], 1, paste, collapse = "")
    right <- apply(flanks[, (w + 1):ncol(flanks), drop = FALSE], 1, paste, collapse = "")

    cKmers <- paste0(left, "C", right)
    tKmers <- paste0(left, "T", right)
    pKmers <- c(cKmers, tKmers)

    return(pKmers)

}

#' @export
getGtfSites <- function(gtf) {

    gtfRanges <- GenomicRanges::GRanges(
        seqnames = gtf$Chromosome,
        ranges = IRanges::IRanges(gtf$Start_Position, gtf$End_Position),
        strand = gtf$Strand
    )

    gtfSites <- GenomicRanges::slidingWindows(gtfRanges, 1)
    # ensure sites are ordered by transcription direction
    isMinus <- as.character(GenomicRanges::strand(gtfRanges)) == "-"
    gtfSites[isMinus] <- S4Vectors::endoapply(gtfSites[isMinus], rev)
    gtfSites <- unlist(gtfSites)
    gtfSites$id <- unlist(mapply(function(x, .t) rep(x, .t), gtf$id, GenomicRanges::width(gtfRanges)))

    # add a key
    gtfSites$key <- 1:length(gtfSites)

    return(gtfSites)

}

