prepOccModelInput <- function(nzero, ## if augmenting data
                              threshold, ## min times a pollinator is seen
                              spec, ## specimen data
                              sr.sched, ## sampling schedule
                              veg,
                              metrics,
                              traits){ ## floral availability matrix
    ## this function preps specimen data for the multi season multi
    ## species occupancy model, and returns a lsit of the inits, data,
    ## constants needed to run in jags or nimble
    ## create complete site x date x species matrix based on sample
    ## schedule]
    min.sets <- apply(sapply(metrics, dim), 1,
                      function(x) names(x)[x == min(x)][1])
    min.sites <- rownames(metrics[[min.sets[1]]])
    min.sp <- colnames(metrics[[min.sets[2]]])

    null.mat <- tapply(rep(0, nrow(sr.sched)),
                       list(sites=paste(sr.sched$Site),
                            dates=sr.sched$Date), sum)

    ## extract the relevant specimen, site, date data
    spec.data <- data.frame(pollinator=spec$GenusSpecies,
                            site=spec$Site,
                            date=spec$Date)

    flower.mat <- samp2site.ypr(site=veg$Site,
                                yr=veg$Year,
                                abund=veg[, "Div"])

    flower.mat <- t(apply(flower.mat, 1, function(x){
        nas <- which(is.na(x))
        if(length(nas) != 0){
            x[is.na(x)] <- mean(x, na.rm=TRUE)
        }
        return(x)
    }))


    ## drop sites without data

    min.sites <- min.sites[min.sites %in% rownames(flower.mat)]
    spec.data <- spec.data[spec.data$pollinator %in% min.sp,]
    spec.data <- spec.data[spec.data$site %in% min.sites,]
    spec.data$site <- as.character(spec.data$site)
    null.mat <- null.mat[rownames(null.mat) %in%
                         min.sites,]

    ## make data array
    pollinator.id <- id(spec.data$pollinator)
    occ.mat <- make.mats(pollinator.id,
                         null.mat,
                         pollinator=as.vector(spec.data$pollinator),
                         var1=as.vector(spec.data$site),
                         var2=spec.data$date)
    sites <- rownames(occ.mat[[1]])
    dates <- colnames(occ.mat[[1]])
    species <- names(occ.mat)
    mat <- array(unlist(occ.mat),
                 dim=c(dim(occ.mat[[1]]),
                       length(occ.mat)))
    dimnames(mat) <- list(site=sites, date=dates, species=species)
    mm <- make.mat(mat, threshold, nzero)
    mat <- mm$mat

    ## make 4D matrix
    occ.mat <- lapply(1:dim(mat)[2], function(x) mat[,x,])
    occ.mat.split <- split(occ.mat, dimnames(mat)$date)
    yr.table <- table(dimnames(mat)$date)
    X <- array(NA, dim=c(dim(mat)[1], length(yr.table),
                         max(yr.table), dim(mat)[3]))
    dimnames(X) <- list(site=dimnames(mat)$site,
                        year=unique(dimnames(mat)$date),
                        rep=1:max(yr.table),
                        species=dimnames(mat)$species)

    null.mat <- matrix(NA, dim(mat)[1], dim(mat)[3],
                       dimnames=dimnames(X)[c('site', 'species')])
    f <- function(i) {
        missing <- max(yr.table) - yr.table[i]
        if(missing == 0) return(occ.mat.split[[i]])
        c(occ.mat.split[[i]], lapply(1:missing, function(x) null.mat))
    }
    tmp <- lapply(seq_along(yr.table), f)

    for(i in 1:length(yr.table))
        for(j in 1:max(yr.table))
            X[,i,j,] <- tmp[[i]][[j]]

    ## only keep sites with some positive number of reps
    no.reps <- apply(mat, 1, function(x) sum(x >= 0, na.rm=TRUE))==0
    X <- X[!no.reps,,,,drop=FALSE]

    ## create an array of day of the year
    make.date.mat <- function(dates) {
        date.mat <- array(NA, dim=dim(X)[1:3],
                          dimnames=dimnames(X)[1:3])
        for(i in seq_along(dates)) {
            year <- as.numeric(format(as.Date(dates[[i]],
                                              format='%Y-%m-%d'),
                                      format = '%Y'))
            lengths <- rle(as.vector(year))$lengths
            ind <- cbind(rep(i, sum(lengths)),
                         match(year, dimnames(X)$year),
                         as.vector(unlist(sapply(lengths,
                                                 seq_len))))
            date.mat[ind] <- strptime(dates[[i]],
                                      '%Y-%m-%d')$yday+1
        }
        date.mat
    }
    date.occ.mat <- lapply(mm$dates, make.date.mat)
    date.mat <- array(unlist(date.occ.mat),
                      dim=c(dim(date.occ.mat[[1]]), length(date.occ.mat)))
    dimnames(date.mat) <- dimnames(X)

    ## function to re-arrange replicate dimension
    compress <- function(x) {
        if(!any(is.na(x))) return(x)
        return(c(x[!is.na(x)], x[is.na(x)]))
    }
    X <- aperm(apply(X, c(1,2,4), compress), c(2,3,1,4))
    names(dimnames(X)) <- c("site", "year", "rep", "species")
    date.mat <- aperm(apply(date.mat, c(1,2,4), compress), c(2,3,1,4))

    ## where missed data, take mean for species, and standardize
    fillNAsMean <- function(traits){
        apply(traits, 2, function(x){
            nas <- which(is.na(x))
            if(length(nas) != 0){
                x[is.na(x)] <- mean(x, na.rm=TRUE)
            }
            return(x)
        })
    }
  ## if the cell is NA, fill it with the mean of the spp across site
    metrics.filled <- lapply(metrics, fillNAsMean)
  ## x[min.sites, min.sp] is unnecessary beacuse the data is already all
  ## organized, but I'll leave it to be extra sure
  ## it is the raw data so it needs to be standardized
    sp.data <- lapply(metrics.filled, function(x){
      standardize(x[min.sites, min.sp])
    })
    ## rerganizing the data according to X (also unnecessary but 
    ## I'll leave it to be extra carefull)
    sp.data <- lapply(sp.data, function(x){
      x <- x[dimnames(X)$site, dimnames(X)$species]
    })
    ## drop last year of flower data
    flower.mat <- flower.mat[rownames(flower.mat) %in%
                             dimnames(X)$site,]
    flower.mat <- flower.mat[, as.numeric(colnames(flower.mat)) <
                               max(dimnames(X)$year)]
    flower.mat <- standardize(flower.mat)

    nsp <- dim(X)[4]
    inits <- getInits(nsp)
    constants <-  list(nrep=apply(X, c(1,2,4),
                                  function(x) sum(x >= 0,na.rm=TRUE)),
                       nsp=dim(X)[4],
                       nsite=dim(X)[1],
                       nyear=dim(X)[2],
                       max.nreps = dim(X)[3])

    ## create date polynomial terms
    day <- day.2 <- date.mat
    poly.terms <- poly(date.mat[!is.na(date.mat)], 2)
    day[!is.na(day)]  <- standardize(poly.terms[,1])
    day.2[!is.na(day.2)]  <- standardize(poly.terms[,2])

    ## fill is -1000 to avoid computation when using function for
    ## integrating
    X[is.na(X) ] <- -1000
    data <- c(list(X=X,
                   day=day,
                   day.2=day.2,
                   fra = flower.mat[dimnames(X)$site,]),
              sp.data)

    monitors <- getParams()
    my.inits <- inits
    model.input <- list(data=data,
                        monitors=monitors,
                        constants=constants,
                        inits=my.inits)
    return(model.input)
}

samp2site.ypr <- function(site, yr, abund) {
    x <- tapply(abund, list(site=site,yr=yr),
                sum, na.rm=TRUE)
    return(x)
}


getInits <- function(nsp){
    list(mu.p.0 = rnorm(1),
         sigma.p.0 = runif(1),
         mu.p.day.1 = rnorm(1),
         sigma.p.day.1 = runif(1),
         mu.p.day.2 = rnorm(1),
         sigma.p.day.2= runif(1),

         mu.phi.0 = rnorm(1),
         mu.gam.0 = rnorm(1),
         sigma.phi.0 = runif(1),
         sigma.gam.0 = runif(1),
         mu.phi.fra = rnorm(1),
         sigma.phi.fra = runif(1),
         mu.gam.fra = rnorm(1),
         sigma.gam.fra = runif(1),

         phi.var.partner=rnorm(1),
         gam.var.partner=rnorm(1),
         phi.var.role=rnorm(1), 
         gam.var.role=rnorm(1), 
         phi.var.cnodf=rnorm(1),
         gam.var.cnodf=rnorm(1), 

         p.0 = rnorm(nsp),
         p.day.1 = rnorm(nsp),
         p.day.2 = rnorm(nsp),
         phi.0 = rnorm(nsp),
         gam.0 = rnorm(nsp)
         )
}


getParams <- function(){
    c('mu.p.0',
      'sigma.p.0',
      'mu.p.day.1',
      'sigma.p.day.1',
      'mu.p.day.2',
      'sigma.p.day.2',

      'mu.phi.0',
      'sigma.phi.0',
      'mu.gam.0',
      'sigma.gam.0',

      "mu.phi.fra", "mu.gam.fra",
      "sigma.phi.fra", "sigma.gam.fra",
      "phi.var.partner",
      "gam.var.partner",
      "phi.var.role",
      "gam.var.role",
      "phi.var.cnodf",
      "gam.var.cnodf"
      )
}
