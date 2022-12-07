simplifyLm = function(cm) {

    cm$y = c()
    cm$model = c()

    cm$residuals = c()
    cm$fitted.values = c()
    cm$effects = c()
    cm$qr$qr = c()  
    cm$linear.predictors = c()
    cm$weights = c()
    cm$prior.weights = c()
    cm$data = c()


    cm$family$variance = c()
    cm$family$dev.resids = c()
    cm$family$aic = c()
    cm$family$validmu = c()
    cm$family$simulate = c()
    attr(cm$terms,".Environment") = c()
    attr(cm$formula,".Environment") = c()

    cm
}

mtypes <- function(kmer, nuc) {

    nucs <- c("A", "C", "G", "T")
    mutTypes <- paste0(kmer, ">", nucs[!(nucs %in% nuc)])

    return(mutTypes)

}
