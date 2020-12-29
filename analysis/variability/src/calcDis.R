

calcDisper <- function(arr, ind){
  pp.dis <- apply(arr, ind, function(y){
    vegdist(empty(t(y)), method="chao")
  })
  pp.dis <- pp.dis[sapply(pp.dis, function(q){
    all(dim(as.matrix(q)) > 3)
  })]
  if(length(pp.dis) == 0){
    return(NA)
  } else{
    pp.disper <- lapply(pp.dis, function(z){
      out <- try(betadisper(z,
                            group=rep("1",
                              nrow(as.matrix(z))))$distances,
                 silent=TRUE)
      if(inherits(out, "try-error")){
        out <- NA
        names(out) <- "NA.NA"
      }
      return(out)
    })
    return(pp.disper)
    
  }
}

calcDis <- function(site, ind, nets, spec){
  these.nets <- nets[sites == site]
  arr <- simplify2array(these.nets)
  pp <- calcDisper(arr, ind)
  if(all(is.na(pp))){
    return(NA)
  }else{
    out.dist <- try(data.frame(GenusSpecies=
                               rep(names(pp), sapply(pp, length)),
                               Dist=unlist(pp),
                               Site= site,
                               Year =  sapply(strsplit(
                                 unlist(sapply(pp, names)), "[.]"),
                                 function(x) x[2])))
    out.dist$SiteStatus <- spec$SiteStatus[match(paste(out.dist$Site,
                                                       out.dist$Year),
                                                 paste(spec$Site,
                                                       spec$Year))]
    rownames(out.dist) <- NULL
    return(out.dist)
  }
}


getDis <- function(sites, ind, nets, specs.agg, traits, spec){
  pp <- lapply(unique(sites), calcDis, ind=ind, nets, spec)
  pp <- pp[!sapply(pp, function(x) all(is.na(x)))]
  pp <- do.call(rbind, pp)
  pp$k <- specs.agg$k[match(pp$GenusSpecies,
                            specs.agg$GenusSpecies)]
  pp$closeness <- specs.agg$closeness[match(pp$GenusSpecies,
                                            specs.agg$GenusSpecies)]
  pp$d <- traits$d[match(pp$GenusSpecies,
                         traits$GenusSpecies)]
  pp$occ.date <- traits$occ.date[match(pp$GenusSpecies,
                                       traits$GenusSpecies)]
  pp$core <- "core"
  pp$core[pp$k <= 1] <- "peripheral"
  return(pp)
}
