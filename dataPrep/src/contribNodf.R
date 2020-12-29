## Modified from bipartite::nestedcontribution. I only commented out the calculation 
## for plants because it takes too long to run with 999

nestedcontributionHiguer <- function (web, nsimul = 999) 
{
    web <- ifelse(web > 0, 1, 0)
    if (is.null(rownames(web))) 
        rownames(web) <- paste0("L", seq.int(nrow(web)))
    if (is.null(colnames(web))) 
        colnames(web) <- paste0("H", seq.int(ncol(web)))
    lower.out <- data.frame(row.names = rownames(web))
    lower.out$nestedcontribution <- NA
    higher.out <- data.frame(row.names = colnames(web))
    higher.out$nestedcontribution <- NA
    if (any(dim(web) < 2)) {
        warning("Your web is too small for a meaningful computation of nestedcontrrank (and probably other indices)!")
    }
    else {
        nested.orig <- vegan::nestednodf(web)$statistic["NODF"]
        # for (i in rownames(web)) {
        #     message(i)
        #     probs <- (rowSums(web)[i]/ncol(web) + colSums(web)/nrow(web))/2
        #     nested.null <- sapply(1:nsimul, function(x) {
        #         web.null <- web
        #         web.null[i, ] <- rbinom(ncol(web), 1, probs)
        #         vegan::nestednodf(web.null)$statistic["NODF"]
        #     })
        #     lower.out[i, "nestedcontribution"] <- (nested.orig - 
        #                                                mean(nested.null))/sd(nested.null)
        # }
        for (i in colnames(web)) {
            message(i)
            probs <- (rowSums(web)/ncol(web) + colSums(web)[i]/nrow(web))/2
            nested.null <- sapply(1:nsimul, function(x) {
                web.null <- web
                web.null[, i] <- rbinom(nrow(web), 1, probs)
                vegan::nestednodf(web.null)$statistic["NODF"]
            })
            higher.out[i, "nestedcontribution"] <- (nested.orig - 
                                                        mean(nested.null))/sd(nested.null)
        }
    }
    #out <- list(`higher level` = higher.out, `lower level` = lower.out)
    out <- list(`higher level` = higher.out)
    return(out)
}

calcNofContrib <- function(networks){
    ## applies nestedcontribution from bipartite to networks
    ## calculates the species roles from a network and returns a dataframe
    ## with site status and ypr
    ## takes networks and specimen data
    nested.contr <- lapply(networks, nestedcontributionHiguer)
    all.nested.contr <- rapply(nested.contr, function(x) as.matrix(x),
                               how="replace")
    all.nested.contr <- lapply(nested.contr, function(x) {
        do.call(rbind, x)
    })
    all.nested.contr <- do.call(rbind, all.nested.contr)
    all.nested.contr$GenusSpecies <- sapply(strsplit(rownames(all.nested.contr), "[.]"),
                                            function(x)x[4])
    all.nested.contr$Site <- sapply(strsplit(rownames(all.nested.contr), "[.]"),
                                    function(x)x[1])
    all.nested.contr$Year <- sapply(strsplit(rownames(all.nested.contr), "[.]"),
                                    function(x)x[2])
    
    rownames(all.nested.contr) <- NULL
    return(all.nested.contr)
}
