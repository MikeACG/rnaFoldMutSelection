#' @export
rnaFoldMutSelection <- function(cohort, k, cgenomes, gElms, genome, mafs, snpfoldDir, rng, outdir, runId) {

    cat("BEGIN\n")
    set.seed(rng)
    if (missing(runId)) {

        runId <- stringi::stri_join(
            cohort,
            "_",
            k,
            "mers_",
            ifelse(length(names(cgenomes)) > 0, stringi::stri_join(names(cgenomes), collapse = "-"), "simple")
        )

    }

    # create structure for output files
    outRoot <- stringi::stri_join(outdir, runId, "/")
    if (file.exists(outRoot)) stop(stringi::stri_join("Directory for results already exists at", outRoot))
    if (!file.exists(outdir)) dir.create(outdir)
    dir.create(outRoot)
    Mpath <- stringi::stri_join(outRoot, "mutationMatrix.RData")
    pfilesPath <- stringi::stri_join(outRoot, "pfiles/")
    dir.create(pfilesPath)
    expDistPath <- stringi::stri_join(outRoot, "expDistribution.tsv")
    file.create(expDistPath)
    obsDistPath <- stringi::stri_join(outRoot, "obsDistribution.txt")
    file.create(obsDistPath)

    # make binary genome based on regions of interest
    gtf <- data.table::rbindlist(gElms$gtf)
    bgenome <- ranges2track(gtf, genome)

    # calculate abundances
    cat("computing abundances...\n")
    abundance <- countAbundance(k, bgenome, cgenomes, genome)

    # calculate mutational matrix
    cat("computing mutational matrix...\n")
    M <- mutmat(abundance, cohort, mafs, bgenome, cgenomes, genome)
    save(M, file = Mpath)

    # clean some memory space
    bgenome <- gtf <- NULL; gc()

    # for each genomic element
    trSNPfoldPaths <- stringi::stri_join(snpfoldDir, gElms$Id, sep = "/")
    pPaths <- stringi::stri_join(pfilesPath, gElms$Id, ".tsv", sep = "")
    n <- nrow(gElms)
    pvalues <- numeric(n)
    nmuts <- integer(n)
    cat("---inference routine---\n")
    for (i in 1:n) {

        if (i %% 100 == 0) cat(stringi::stri_join("\t", i, "/", n, "...\n", sep = ""))

        # load SNPfold results
        trSNPfoldDt <- data.table::fread(trSNPfoldPaths[i])

        # assign transition probabilities and observed mutations at each site
        trTransit <- mutmatApplyTr(M, gElms$gtf[[i]], trSNPfoldDt, genome, cgenomes)
        trObs <- trCountMuts(gElms$gtf[[i]], trSNPfoldDt, mafs, cohort)
        P <- data.table::data.table(transitProb = trTransit, mutObs = trObs)

        # record PCC values for histogram construction purposes
        expDist <- data.table::data.table(PCC = trSNPfoldDt$PCC, transitProb = P$transitProb)
        obsDist <- data.table::data.table(PCC = rep(trSNPfoldDt$PCC, P$mutObs))

        # conduct test for selection
        pvalues[i] <- selectionTest(P, trSNPfoldDt)
        nmuts[i] <- sum(P$mutObs)

        # save intermediate data
        data.table::fwrite(P, pPaths[i], sep = "\t")
        if (nmuts[i] > 0) {
            data.table::fwrite(expDist, expDistPath, sep = "\t", append = TRUE)
            data.table::fwrite(obsDist, obsDistPath, append = TRUE)
        }

    }

    # save results
    R <- data.table::data.table(Id = gElms$Id, Name = gElms$Name, Mutations = nmuts, Pvalue = pvalues)
    resultsPath <- stringi::stri_join(outRoot, "results.tsv")
    data.table::fwrite(R, resultsPath, sep = "\t", na = "NA")
    cat("DONE\n")

}

countAbundance <- function(k, bgenome, cgenomes, genome) {

    cvars <- names(cgenomes)
    pKmers <- possibleKmers(k, cvars, cgenomes)
    w <- (k - 1) / 2

    abundance <- list()
    for (i in 1:length(genome)) {

        # get kmers of positions of interest
        chr <- names(genome)[i]
        chrPositions <- Biostrings::matchPattern("Y", bgenome[[chr]])
        chrKmerRanges <- IRanges::IRanges(IRanges::start(chrPositions) - w, IRanges::end(chrPositions) + w)
        kmers <- Biostrings::extractAt(genome[[chr]], chrKmerRanges)

        # orient by pyrimidine
        kmerIsPurine <- isPurine(kmers, w)
        kmers[kmerIsPurine] <- Biostrings::reverseComplement(kmers[kmerIsPurine])

        # if a genome covariate is used, modify kmers accordignly
        if (length(cvars) > 0) {

            chrSiteRanges <- IRanges::IRanges(IRanges::start(chrPositions), IRanges::end(chrPositions))
            cvarsKmers <- lapply(cgenomes[cvars], cvarKmer, chr, chrSiteRanges, kmerIsPurine)
            cvarsKmers <- apply(do.call(cbind, cvarsKmers), 1, stringi::stri_join, collapse = "_")
            kmers <- stringi::stri_join(kmers, cvarsKmers, sep = "_")

        }

        # count
        abundance[[i]] <- table(factor(kmers, levels = pKmers))

    }
    abundance <- do.call(rbind, abundance)

    abundance <- apply(abundance, 2, sum)
    abundanceProb <- abundance / sum(abundance)
    attr(abundanceProb, "abundanceCount") <- abundance
    attr(abundanceProb, "k") <- k
    attr(abundanceProb, "cvars") <- cvars
    if (length(cvars) == 0) attr(abundanceProb, "cvars") <- "simple"
    return(abundanceProb)

}

mutmat <- function(abundance, cohort, mafs, bgenome, cgenomes, genome) {

    pContexts <- names(abundance)
    cvars <- attr(abundance, "cvars")
    k <- attr(abundance, "k")
    w <- (k - 1) / 2

    chrMats <- list()
    for (i in 1:length(genome)) { # this correctly gives a matrix of zeros when there are no mutations

        # get data for this chromosome
        chr <- names(genome)[i]
        chrPositions <- Biostrings::matchPattern("Y", bgenome[[chr]])
        chrMuts <- mafs[[chr]]
        chrMuts <- chrMuts[chrMuts$Cohort == cohort]

        # get mutations that are in regions of interest
        chrMutsRanges <- IRanges::IRanges(chrMuts$Start_Position, chrMuts$End_Position)
        chrMutsHits <- IRanges::countOverlaps(chrMutsRanges, chrPositions)
        chrMutsThatHit <- chrMuts[chrMutsHits > 0, ]

        # get contexts of mutations that hit and orient by pyrimidine
        chrHitsContextRanges <- IRanges::IRanges(chrMutsThatHit$Start_Position - w, chrMutsThatHit$End_Position + w)
        chrMutsContext <- Biostrings::extractAt(genome[[chr]], chrHitsContextRanges)
        chrMutsContextIsPurine <- isPurine(chrMutsContext, w)
        chrMutsContext[chrMutsContextIsPurine] <- Biostrings::reverseComplement(chrMutsContext[chrMutsContextIsPurine])
        chrMutsThatHit$Tumor_Seq_Allele2[chrMutsContextIsPurine] <- chartr("ACGT", "TGCA", chrMutsThatHit$Tumor_Seq_Allele2[chrMutsContextIsPurine])

        # if a genome covariate is used, modify contexts accordignly
        if (cvars[1] != "simple") {

            chrSiteRanges <- IRanges::IRanges(chrMutsThatHit$Start_Position, chrMutsThatHit$End_Position)
            cvarsContexts <- lapply(cgenomes[cvars], cvarKmer, chr, chrSiteRanges, chrMutsContextIsPurine)
            cvarsContexts <- apply(do.call(cbind, cvarsContexts), 1, stringi::stri_join, collapse = "_")
            chrMutsContext <- stringi::stri_join(chrMutsContext, cvarsContexts, sep = "_")

        }

        # count the substitution type for the chromosome
        chrMats[[i]] <- table(
            factor(chrMutsContext, levels = pContexts),
            factor(chrMutsThatHit$Tumor_Seq_Allele2, levels = c("A", "C", "G", "T"))
        )
    }

    mutMat <- apply(simplify2array(chrMats), 1:2, sum)
    foldMat <- mutMat / (abundance * sum(mutMat))
    foldMatHash <- mat2hash(foldMat, w)

    attr(foldMatHash, "mutMat") <- mutMat
    attr(foldMatHash, "foldMat") <- foldMat
    attr(foldMatHash, "abundance") <- abundance
    attr(foldMatHash, "cohort") <- cohort
    attr(foldMatHash, "k") <- k
    attr(foldMatHash, "cvars") <- cvars
    return(foldMatHash)

}

mutmatApplyTr <- function(M, trGtf, trSNPfoldDt, genome, cgenomes) {

    trChr <- trGtf$Chromosome[1]
    trStrand <- trGtf$Strand[1]
    cvars <- attr(M, "cvars")
    k <- attr(M, "k")
    w <- (k - 1) / 2
    changeTo <- gsub("U", "T", trSNPfoldDt$Change_To)

    # get transcript kmers (gtf ranges appear in trancription order but still need to read them from end to start if the coding strand is -)
    if (trStrand == "+") {

        trSites <- unlist(mapply(':', trGtf$Start_Position, trGtf$End_Position, SIMPLIFY = FALSE))

    }
    if (trStrand == "-") {
        
        trSites <- unlist(mapply(':', trGtf$End_Position, trGtf$Start_Position, SIMPLIFY = FALSE))
    
    }
    trKmerRanges <- IRanges::IRanges(trSites - w, trSites + w)
    trKmers <- Biostrings::extractAt(genome[[trChr]], trKmerRanges)
    # if coding strand is "-", need to get reverse complement to get right kmer in strand
    countStrand <- rep("+", length(trKmerRanges))
    if (trStrand == "-") {
        trKmers <- Biostrings::reverseComplement(trKmers)
        countStrand <- rep("-", length(trKmerRanges))
    }

    # orient by pyrimidine
    trKmerIsPurine <- isPurine(trKmers, w)
    trKmers[trKmerIsPurine] <- Biostrings::reverseComplement(trKmers[trKmerIsPurine])
    countStrand[trKmerIsPurine] <- ifelse(countStrand[trKmerIsPurine] == "+", "-", "+")

    # get covariates of each kmer, if they are used
    if (cvars[1] != "simple") {

        chrSiteRanges <- IRanges::IRanges(trSites, trSites)
        cvarsKmers <- lapply(cgenomes[cvars], cvarKmer2, trChr, chrSiteRanges, countStrand)
        cvarsKmers <- apply(do.call(cbind, cvarsKmers), 1, stringi::stri_join, collapse = "_")
        trKmers <- stringi::stri_join(trKmers, cvarsKmers, sep = "_")

    }

    # find correct key for each possible substitution in each transcript site in mutation matrix
    changeKmers <- rep(as.character(trKmers), each = 3)
    changeIsPurine <- rep(trKmerIsPurine, each = 3)
    changeTo[changeIsPurine] <- chartr("ACGT", "TGCA", changeTo[changeIsPurine])
    trKeys <- stringi::stri_join(changeKmers, changeTo, sep = ">")

    # get probabilities, normalize and handle case when prob for all sites is zero
    trMutprob <- hash::values(M, trKeys)
    trMutprob <- trMutprob / sum(trMutprob)
    trMutprob <- ifelse(is.finite(trMutprob), trMutprob, rep(0, length(trMutprob)))
    return(trMutprob)

}

trCountMuts <- function(trGtf, trSNPfoldDt, mafs, cohort) {

    trChr <- trGtf$Chromosome[1]
    trStrand <- trGtf$Strand[1] 

    # find mutations overlapping the transcript
    chrCohortMaf <- mafs[[trChr]]
    chrCohortMaf <- chrCohortMaf[chrCohortMaf$Cohort == cohort]
    gtfRanges <- IRanges::IRanges(trGtf$Start_Position, trGtf$End_Position)
    mutRanges <- IRanges::IRanges(chrCohortMaf$Start_Position, chrCohortMaf$End_Position)
    mutOverlaps <- IRanges::countOverlaps(mutRanges, gtfRanges)
    trMaf <- chrCohortMaf[mutOverlaps > 0, ]
    trMaf$Reference_Allele <- gsub("T", "U", trMaf$Reference_Allele)
    trMaf$Tumor_Seq_Allele2 <- gsub("T", "U", trMaf$Tumor_Seq_Allele2)

    # find positions of mutations that hit inside the transcript
    if (trStrand == "+") {

        trSites <- unlist(mapply(':', trGtf$Start_Position, trGtf$End_Position, SIMPLIFY = FALSE))

    }
    if (trStrand == "-") {
        
        trSites <- unlist(mapply(':', trGtf$End_Position, trGtf$Start_Position, SIMPLIFY = FALSE))
    
    }
    trSites <- setNames(1:length(trSites), trSites)
    trMutsPos <- trSites[as.character(trMaf$Start_Position)] # works for SNPs only

    # according to coding strand of gene, get the RNA change
    s <- rep(trStrand, nrow(trMaf)) == "+"
    trWts <- ifelse(s, trMaf$Reference_Allele, chartr("ACGU", "UGCA", trMaf$Reference_Allele))
    trMuts <- ifelse(s, trMaf$Tumor_Seq_Allele2, chartr("ACGU", "UGCA", trMaf$Tumor_Seq_Allele2))
    trChanges <- stringi::stri_join(trWts, trMutsPos, trMuts)

    # count
    trMutCount <- table(factor(trChanges, levels = trSNPfoldDt$RNA_Change))

    return(as.integer(trMutCount))

}

selectionTest <- function(P, trSNPfoldDt) {

    nsims <- 1e+05

    obs <- rep(trSNPfoldDt$PCC, P$mutObs)
    nobs <- length(obs)
    if (length(obs) == 0) return(NA)

    sims <- rmultinom(nsims, nobs, P$transitProb)
    simSums <- apply(sims, 2, function(sim) sum(rep(trSNPfoldDt$PCC, sim)))

    return(sum(simSums <= sum(obs)) / nsims)

}

#' @export
ensemblSelectionTest <- function(rootDirs, gElms, snpfoldDir, rng, outFile) {

    set.seed(rng)
    nsims <- 1e+05

    # for each genomic element
    trSNPfoldPaths <- stringi::stri_join(snpfoldDir, gElms$Id, sep = "/")
    n <- nrow(gElms)
    pvalues <- numeric(n)
    nmuts <- integer(n)
    for (i in 1:n) {

        if (i %% 100 == 0) cat(stringi::stri_join(i, "/", n, "...\n", sep = ""))

        # load SNPfold results
        trSNPfoldDt <- data.table::fread(trSNPfoldPaths[i])

        # for each model in ensembl load the transition probabilities for each site
        pPaths <- stringi::stri_join(rootDirs, gElms[i], ".tsv")
        Ps <- lapply(pPaths, data.table::fread)

        # filter out models where no mutations are observed
        Ps <- Ps[sapply(Ps, function(P) sum(P$mutObs)) > 0]

        # for each model in ensembl get the PCC sum of nsims simulations
        simSums <- sapply(Ps, function(P) {
            obs <- rep(trSNPfoldDt$PCC, P$mutObs)
            nobs <- length(obs)
            sims <- rmultinom(nsims, nobs, P$transitProb)
            apply(sims, 2, function(sim) sum(rep(trSNPfoldDt$PCC, sim)))
        })
        simSums <- as.matrix(simSums)

        # sum of sums along ensembl
        simSS <- apply(simSums, 1, sum) 

        # get observed sum along ensembl
        obs <- as.matrix(sapply(Ps, function(P) P$mutObs))
        obsSum <- sum(apply(obs, 2, function(o) sum(rep(trSNPfoldDt$PCC, o))))
        
        # compute empirical pvalue of ensembl
        pvalues[i] <- sum(simSS <= obsSum) / nsims
        nmuts[i] <- sum(obs)

    }

    # save results
    R <- data.table::data.table(Id = gElms$Id, Name = gElms$Name, Mutations = nmuts, Pvalue = pvalues)
    data.table::fwrite(R, outFile, sep = "\t", na = "NA")
    cat("DONE\n")

}
