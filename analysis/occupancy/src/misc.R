
logit <- function(x) {
    log(x/(1 - x))
}

expit <- function(x) {
    exp(x)/(1 + exp(x))
}


matchSpPersistFlex <- function(pattern,
                               sds,
                               means,
                               zscores,
                               drop.hyper.param,
                               model.input){
    ## match phi/psi values to interaction flexibility
    means <- means[grepl(pattern, names(means))]
    if(drop.hyper.param){
        to.drop <- names(means) %in% c(paste("mu", pattern,
                                             "0", sep="."),
                                       paste("sigma", pattern,
                                             "0", sep="."))
        means <- means[!to.drop]
        sds <- sds[grepl(pattern, names(sds))][!to.drop]
    }
    names(means) <-  names(sds) <-
        dimnames(model.input$data$X)$species

    zscores[pattern] <- means[match(zscores$GenusSpecies,
                                    names(means))]
    zscores[paste(pattern, "sd", sep=".")] <- sds[match(zscores$GenusSpecies,
                                                        names(sds))]

    return(zscores)
}


standardize <- function(x)
(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)



## function to make pollinator visitation matrices
make.mats <- function(pollinator.id,
                      null.mat,
                      pollinator,
                      var1,
                      var2) {
    make.mat <- function(P) {
        var1 <- var1[pollinator==P]
        var2 <- var2[pollinator==P]
        m <- tapply(rep(1, length(var1)),
                    list(sites=var1, dates=var2), sum)
        null.mat[rownames(m), colnames(m)][!is.na(m)] <- m[!is.na(m)]
        null.mat
    }
    mats <- lapply(pollinator.id, function(x) make.mat(x))
    names(mats) <- pollinator.id
    mats
}



## convert occurrence data into a site by species matrix
samp2site.spp <- function(site, yr, abund) {
    x <- tapply(abund, list(site=site,yr=yr), sum)
    x[is.na(x)] <- 0
    return(x)
}



## return sorted unique values
id <- function(x) as.character(unique(sort(x)))

## gets the means and sds from variables and round them 
calcMeanSd <- function(x){
  meanSd <- paste0(round(all.vars[,grep(x, colnames(all.vars))][1],3), " (", round(all.vars[,grep(x, colnames(all.vars))][2],3), ")")
  return(meanSd)
}
