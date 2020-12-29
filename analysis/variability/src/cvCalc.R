## regresses coefficent of variation against traits

corCv <- function(x,...){
    cv(x)*(1 + (1/4*length(x)))
}

calcCvTrait <- function(spec.dat,
                     byType,
                     trait1,
                     trait2,
                     time.col,
                     abund.col,
                     status.order=c("control", "maturing", "mature"),
                     cv.function=corCv,
                     zero2na =FALSE,
                     standard.cv=TRUE,
                     species.type="GenusSpecies",...){
    ##spec.dat: specimens data
    ## byType: table of species level network metrics
    ## trait1: first trait of interest (to regress cv
    ## against)
    ## trait2: second trait of interest
    ## time.col: name of the time column
    ## abund.col: network metric column name
    ## cv.function: function to use for calculating cv. cv.function
    ## can be corCv, cv, or sd

    ## zero2na: converts NAs to zeros, but NAs are not
    ## inclcued in the cv calculations whereas zeros are, so
    ## recommended to keep as FALSE
    ## standard cv logs the cv
    byStatus <- split(byType, byType$SiteStatus)
    bySite <- lapply(byStatus, function(x) {split(x, x$Site)})
    bySite <- unlist(bySite, recursive=FALSE)
    prep.cv <- lapply(bySite, function(y) {
        samp2site.spp(y[, time.col], y[,"GenusSpecies"], y[, abund.col])
    })
    ## to avoid zeros in calculation of cv
    if(zero2na){
        prep.cv <- lapply(prep.cv, function(x){
            x[x == 0] <- NA
            return(x)
        })
    }
    coeff.cv <- lapply(prep.cv, function(x){apply(x, 2, cv.function, ...)})
    dats <- data.frame(cv=unlist(coeff.cv))
    if(standard.cv){
        dats$cv[!is.na(dats$cv)] <- log(dats$cv[!is.na(dats$cv)])
    }
    dats$SiteStatus <-  gsub('\\..*', '', rownames(dats))
    dats$SiteStatus <- factor(dats$SiteStatus, levels=status.order)
    dats$GenusSpecies <- unlist(lapply(coeff.cv, names))
    dats$Site <-  sapply(strsplit(rownames(dats), "\\."),
                         function(x) x[2])
    rownames(dats) <- NULL

    dats <- cbind(dats, spec.dat[, c(trait1, trait2)][match(dats$GenusSpecies,
                                                            spec.dat[,
                                                                     species.type]),])
    lm.dats <- dats[!is.na(dats$cv) & !is.na(dats[,trait1]) & !is.na(dats[,trait2]),]
    return(list(data=dats,
                lm.data =lm.dats))
}

## variance inflation factor

vif.mer <- function (fit) {
    ## adapted from rms::vif

    v <- vcov(fit)
    nam <- names(fixef(fit))

    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
        v <- v[-(1:ns), -(1:ns), drop = FALSE]
        nam <- nam[-(1:ns)]
    }

    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    return(v)
}


formula.cv <- formula(cv ~ occ.date +  r.degree + (1|Site) +
                          (1|GenusSpecies))


formula.plant.cv <- formula(cv ~ occ.plant.date +  plant.r.degree +
                                (1|Site) + (1|GenusSpecies))
