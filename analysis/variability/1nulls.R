## this file creates the null models (reshuffling interactions within
## communities) necessary for the partner variability calculations

## setwd('~/Dropbox/speciesRoles/')
rm(list=ls())
setwd('analysis/variability')
source('src/initialize_nulls.R')

nnull <- 999
## ************************************************************
## year by species matrices pollinators
## ************************************************************
comms <- lapply(sites, calcYearBeta,
                species.type=species.type,
                spec=spec,
                species.type.int=species.type.int)

comm <- makePretty(comms, spec)

save(comm, file=file.path(save.dir.comm,
                          sprintf('%s-abund.Rdata', type)))

## ************************************************************
## alpha div nulls
## ************************************************************
method.alpha <- "r0_ind"

load(file=file.path(save.dir.comm,
                    sprintf('%s-abund.Rdata', type)))

occ.null <- function(web){
  simulate(vegan::nullmodel(web, method=method.alpha),1)[,,1]
}

rep.occ.null <- function(web, N){
  replicate(N, occ.null(web), simplify = FALSE)
}
nulls <- rapply(comm$comm, rep.occ.null, N=nnull, how="replace")

save(nulls, file=file.path(save.dir.nulls,
                           sprintf('%s-alpha.Rdata', type)))

## ************************************************************
## occurrence nulls
## ************************************************************
method.occ <- "r0"

load(file=file.path(save.dir.comm,
                    sprintf('%s-abund.Rdata', type)))

occ.null <- function(web){
    simulate(vegan::nullmodel(web, method=method.occ),1)[,,1]
}

rep.occ.null <- function(web, N){
  replicate(N, occ.null(web), simplify = FALSE)
}
nulls <- rapply(comm$comm, rep.occ.null, N=nnull, how="replace")

save(nulls, file=file.path(save.dir.nulls,
                           sprintf('%s-occ.Rdata', type)))
