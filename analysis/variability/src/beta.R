calcBetaStatus <- function(comm,
                           dis.method,
                           nulls,
                           occ= FALSE,
                           zscore=FALSE){ ## calculate zscores?
    ## computes dispersion of community matrices, returns output of
    ## vegan function
    ## create community dissimilarity matrix
    comm.dis <-  lapply(comm, function(x) {
        as.matrix(vegdist(x, method= dis.method, diag= TRUE, binary= occ))
    })
    ## null dissimilarity matrices: list of spp. Each 999 element is a null matrix of year x plant spp
    null.dis <- vector("list", length(nulls))
    for(i in 1:length(nulls)){
        this.null <- nulls[[i]]
        null.dis[[i]] <- lapply(this.null, function(x) {
            as.matrix(vegdist(x, method= dis.method, diag=TRUE, binary = occ))
        })
        null.dis[[i]][[length(nulls[[i]]) + 1]] <- comm.dis[[i]]
        
    }
    beta.disper.result <- vector("list", length(comm))
    for(i in 1:length(comm)){
        this.spp <- names(comm)[i]
        arr <- array(unlist(null.dis[[i]]), c(dim(comm.dis[[i]])[1],
                                              dim(comm.dis[[i]])[2],
                                              length(nulls[[i]]) + 1))
        ## standardize dissimilarities
        if(!zscore){
            less.than  <-   apply(arr, 1:2, function(x){
                sum(x[length(null.dis[[i]])] > x)
            })
            equal.2  <-   apply(arr, 1:2, function(x){
                sum(x[length(null.dis[[i]])] == x)
            })
            cor.dis <- as.dist((less.than + 0.5*equal.2)/
                                   length(null.dis[[i]]), diag= TRUE)
        }else{
            cor.dis  <-  (comm.dis[[i]] -
                              apply(arr , 1:2 , mean))/
                (apply(arr , 1:2 , sd) + 10^-10)
            cor.dis <-  try(as.dist(((cor.dis - min(cor.dis))/diff(range(cor.dis))),
                                    diag= TRUE), silent=TRUE)
            if(inherits(cor.dis, "try-error"))browser()
        }
        ## run model
        beta.disper.result[[i]] <- try(betadisper(cor.dis, rep("all", nrow(as.matrix(cor.dis))),
                                                  type="centroid", add = TRUE), silent = TRUE)
        if(inherits(beta.disper.result[[i]], "try-error"))browser()
    }
    return(beta.disper.result)
}

## wrangling the results
makeBetaDataPretty <- function(comm, dis){
    distances <- list()
    for(i in 1:length(dis)){
        this.site <- unlist(lapply(strsplit(names(dis[i]), "[.]"), function(x) x[1]))
        these.spp <- names(comm[[i]])
        names(dis[[i]]) <- these.spp
        distances[[i]] <- lapply(dis[[i]], function(y) y$distances)
    }
    names(distances) <- names(dis)
    dist <- unlist(distances)
    #browser()
    if(zscore){
        GenusSpecies <- unlist(lapply(strsplit(names(dist), "[.]"), function(x) x[2]))
        Site <-  unlist(lapply(strsplit(names(dist), "[.]"), function(x) x[1]))
        Year <-  unlist(lapply(strsplit(names(dist), "[.]"), function(x) x[3]))
        dats <- cbind(as.data.frame(dist), GenusSpecies, Site, Year)
    } else{
        GenusSpecies <- unlist(lapply(strsplit(names(dist), "[.]"), function(x) x[2]))
        GenusSpecies <- gsub('.{0,1}$', '', GenusSpecies)
        Site <-  unlist(lapply(strsplit(names(dist), "[.]"), function(x) x[1]))
        dats <- cbind(as.data.frame(dist), GenusSpecies, Site)
    }
    rownames(dats) <- NULL
    return(dats)
}
