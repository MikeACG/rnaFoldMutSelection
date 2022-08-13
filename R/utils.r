#' @export
ranges2track <- function(ranges, genome) {

    track <- Biostrings::BStringSet(lapply(genome, function(chr) Biostrings::BString(strrep("N", length(chr)))))
    
    for (i in 1:length(genome)) {
        
        chr <- names(genome)[i]
        chrRanges <- ranges[ranges$Chromosome == chr]

        iranges <- IRanges::disjoin(IRanges::IRanges(chrRanges$Start_Position, chrRanges$End_Position))
        value <- sapply(IRanges::width(iranges), function(w) strrep("Y", w))
        track[[chr]] <- Biostrings::replaceAt(track[[chr]], iranges, value)

    }

    return(track)

}

resolveDisjoinCvar <- function(o, d, s) {

    overlap <- IRanges::findOverlaps(o, d)
    if (length(d) == 0) return(s)
    drangesIdxs <- 1:length(d)

    set.seed(123)
    resolved <- sapply(drangesIdxs, function(j) {
        isIdx <- S4Vectors::subjectHits(overlap) == j
        orangesIdxs <- S4Vectors::queryHits(overlap)[isIdx]
        cvars <- unique(s[orangesIdxs])
        if (length(cvars) == 1) return(cvars)
        sample(cvars, 1)
    })

    return(resolved)
}

#' @export
cvar2track <- function(cvar, ranges, genome) {

    track <- Biostrings::BStringSet(lapply(genome, function(chr) Biostrings::BString(strrep("N", length(chr)))))
    
    for (i in 1:length(genome)) {
        
        chr <- names(genome)[i]
        isChr <- ranges$Chromosome == chr
        chrRanges <- ranges[isChr, ]
        chrCvar <- cvar[isChr]

        iranges <- IRanges::IRanges(chrRanges$Start_Position, chrRanges$End_Position)
        dranges <- IRanges::disjoin(iranges)
        chrCvar <- resolveDisjoinCvar(iranges, dranges, chrCvar)
        value <- mapply(strrep, chrCvar, IRanges::width(dranges))
        track[[chr]] <- Biostrings::replaceAt(track[[chr]], dranges, value)

    }

    ulvls <- unique(cvar)
    attr(track, "lvls") <- ulvls
    if (length(ulvls) == 1) attr(track, "lvls") <- c(ulvls, "N")
    return(track)

}

possibleKmers <- function(k, cvars, cgenomes) {

    flanks <- gtools::permutations(4, k - 1, c("A", "C", "G", "T"), repeats.allowed = TRUE)

    w <- ncol(flanks) / 2 
    left <- apply(flanks[, 1:w, drop = FALSE], 1, stringi::stri_join, collapse = "")
    right <- apply(flanks[, (w + 1):ncol(flanks), drop = FALSE], 1, stringi::stri_join, collapse = "")

    cKmers <- stringi::stri_join(left, "C", right)
    tKmers <- stringi::stri_join(left, "T", right)

    pKmers <- c(cKmers, tKmers)
    if (length(cvars) > 0) {

        lvls <- expand.grid(lapply(cgenomes[cvars], attr, "lvls"))
        lvls <- apply(as.matrix(lvls), 1, stringi::stri_join, collapse = "_")
        pKmers <- unlist(lapply(lvls, function(l) stringi::stri_join(pKmers, l, sep = "_")))

    }

    return(pKmers)

}

isPurine <- function(oddseqs, w) {

    centers <- unlist(Biostrings::extractAt(oddseqs, IRanges::IRanges(w + 1, w + 1)))
    isPurine <- as.character(centers) %in% c("A", "G")

    return(isPurine)

}

cvarKmer <- function(cgenome, chr, chrSiteRanges, kmerIsPurine) {

    # get the covariate at sites of kmers
    kmerCvars <- as.character(Biostrings::extractAt(cgenome[[chr]], chrSiteRanges))

    # determine the type of covariate
    psiteCvars <- attr(cgenome, "lvls")
    isAlwaysPresent <- !("N" %in% psiteCvars)

    if (isAlwaysPresent) {

        # compare the strand in which we are counting to the strand at the site and assign same or different
        nkmers <- length(chrSiteRanges)
        countStrand <- ifelse(kmerIsPurine, rep("-", nkmers), rep("+", nkmers))
        kmerCvars <- ifelse(countStrand == kmerCvars, "+", "-")

    }

    return(kmerCvars)

}

mat2hash <- function(M, w) {

    r <- rownames(M)
    .c <- colnames(M)
    v <- as.vector(M)

    combs <- as.matrix(expand.grid(Context = r, Mutation = .c))
    contextCenter <- substr(combs[, "Context"], w + 1, w + 1)
    vIsMut <- contextCenter != combs[, "Mutation"]
    keys <- apply(combs[vIsMut, ], 1, stringi::stri_join, collapse = ">")

    return(hash::hash(keys, v[vIsMut]))

}

cvarKmer2 <- function(cgenome, chr, chrSiteRanges, countStrand) {

    # get the covariate at sites of kmers
    kmerCvars <- as.character(Biostrings::extractAt(cgenome[[chr]], chrSiteRanges))

    # determine the type of covariate
    psiteCvars <- attr(cgenome, "lvls")
    isAlwaysPresent <- !("N" %in% psiteCvars)

    if (isAlwaysPresent) {

        # compare the strand in which we are counting to the strand at the site and assign same or different
        nkmers <- length(chrSiteRanges)
        kmerCvars <- ifelse(countStrand == kmerCvars, rep("+", nkmers), rep("-", nkmers))

    }

    return(kmerCvars)

}

loadRData <- function(fileName){

    load(fileName)
    get(ls()[ls() != "fileName"])

}

#' @export
analyzeFromSchedule <- function(idx, scheduleFile, gElmsPath, genomePath, mafsPath, snpfoldDir, rng, outdir, ...) {

    idx <- as.integer(idx)
    schedule <- data.table::fread(scheduleFile, sep = "\t")

    cohort <- as.character(schedule$cohort[idx])
    k <- as.integer(schedule$k[idx])
    cgenomes <- list()
    if (schedule$cgenomes[idx] != "simple") {

        cgenomesPaths <- unlist(strsplit(schedule$cgenomes[idx], ","))
        cgenomes <- setNames(lapply(cgenomesPaths, loadRData), gsub("[.].*", "", basename(cgenomesPaths)))
    
    }

    gElms <- loadRData(as.character(gElmsPath))
    genome <- loadRData(as.character(genomePath))
    mafs <- loadRData(as.character(mafsPath))
    snpfoldDir <- as.character(snpfoldDir)
    rng <- as.integer(rng)
    outdir <- as.character(outdir)

    rnaFoldMutSelection(cohort, k, cgenomes, gElms, genome, mafs, snpfoldDir, rng, outdir, ...)

}

#' @export
makeSchedule <- function(cohorts, ks, cgenomesPaths) {

    n <- length(cgenomesPaths)
    combs <- lapply(1:n, function(i) gtools::combinations(n, i, cgenomesPaths))
    combs <- lapply(combs, function(.c) apply(.c, 1, stringi::stri_join, collapse = ","))
    
    schedule <- expand.grid(cohort = cohorts, k = ks, cgenomes = c("simple", unlist(combs)))

    return(data.table::data.table(schedule))

}
