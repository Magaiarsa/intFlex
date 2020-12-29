
calcYearBeta <- function(x, species.type, spec, species.type.int,
                         year.cut.off = 3,
                         observation.cut.off = 2){
    ## creates a year by species community
    ## species.type referes to the column in the dataset with the species
    ## names to be aggregated over
    this.spec <- spec[spec$Site == x,]
    prep.comm <- aggregate(list(Abund=this.spec[, species.type]),
                           list(GenusSpecies=this.spec[, species.type],
                                InterGenusSpecies=this.spec[, species.type.int],
                                Year=this.spec$Year),
                           length)

    by.species <- split(prep.comm, prep.comm$GenusSpecies)

    num.observations <- sapply(by.species, function(x) sum(x$Abund))
    ## subset to species seen in > 3 years and at least 5 times
    ## (defaults)
    by.species <- by.species[num.observations >= observation.cut.off]
    num.years <- sapply(by.species, function(x) length(unique(x$Year)))
    by.species <- by.species[num.years >= year.cut.off]
    by.species <- by.species[!sapply(by.species, is.null)]

    ## year plant combinations
    empty.matrix <- matrix(0, nrow=length(unique(prep.comm$Year)),
                           ncol=length(unique(prep.comm$InterGenusSpecies)))
    rownames(empty.matrix) <- sort(unique(prep.comm$Year))
    colnames(empty.matrix) <-
        sort(unique(prep.comm$InterGenusSpecies))

    if(length(by.species) != 0){
        comm <- vector("list", length=length(by.species))
        for(i in 1:length(by.species)){
            comm[[i]] <- empty.matrix
            this.by.species <- by.species[[i]]
            for(j in 1:nrow(this.by.species)){
                this.row <- this.by.species[j,]
                comm[[i]][match(this.row["Year"], rownames(comm[[i]])),
                          match(this.row["InterGenusSpecies"],colnames(comm[[i]]))] <-
                    as.numeric(this.row[["Abund"]])
            }
            comm[[i]] <- comm[[i]][rowSums(comm[[i]]) > 0,]
            ## if(!is.matrix(comm[[i]])) comm[[i]] <- NA
        }
        names(comm) <- names(by.species)
        ## comm <- comm[!sapply(comm, function(x) all(is.na(x)))]
        return(list(comm=comm, site=x))
    }
}

makePretty <- function(comms, spec){
    sites <- sapply(comms, function(x) x$site)
    comms <- lapply(comms, function(x) x$comm)
    comms <- comms[!sapply(sites, is.null)]
    sites <- unlist(sites[!sapply(sites, is.null)])
    names(comms) <- sites
    years <- lapply(comms, function(x) sapply(x, rownames))

    comm.pp <- list(comm=comms,
                    years=years,
                    sites= rep(names(comms),
                               sapply(comms, length)))
    statuses <- spec$SiteStatus[match(comm.pp$sites,
                                      spec$Site)]
    comm.pp$status <- split(statuses,
                            comm.pp$sites)[names(comm.pp$comm)]
    return(comm.pp)
}
