
## This functions takes site-species-abundance data and creates a
## matrix where the sites are columns and the rows are species.

samp2site.spp <- function(site, spp, abund, FUN=sum) {
  x <- tapply(abund, list(site = site, spp = spp), FUN)
  x[is.na(x)] <- 0
  return(x)
}

## does the reverse of samp2site

comm.mat2sample <-  function (z) {
  temp <- data.frame(expand.grid(dimnames(z))[1:2],
                     as.vector(as.matrix(z)))
  temp <- temp[sort.list(temp[, 1]), ]
  data.frame(Site = temp[, 1], Samp = temp[, 3],
             Date = temp[, 2])
}


## function to clean up white-space in a column of data (replaces all
## instances of white-space with " " and empty cells with ""
fix.white.space <- function(d) {
  d <- as.character(d)
  remove.first <- function(s) substr(s, 2, nchar(s))
  d <- gsub("      ", " ", d, fixed=TRUE)
  d <- gsub("     ", " ", d, fixed=TRUE)
  d <- gsub("    ", " ", d, fixed=TRUE)
  d <- gsub("   ", " ", d, fixed=TRUE)
  d <- gsub("  ", " ", d, fixed=TRUE)

  tmp <- strsplit(as.character(d), " ")
  d <- sapply(tmp, function(x) paste(x, collapse=" "))

  first <- substr(d, 1, 1)
  d[first==" "] <- remove.first(d[first==" "])
  d
}


## return sorted unique values
id <- function(x) unique(sort(x))

## function to make pollinator visitation matrices
make.mats <- function(pollinator.id,
                      null.mat,
                      pollinator,
                      var1,
                      var2,
                      occ) {
  make.mat <- function(P) {
    var1 <- var1[pollinator==P]
    var2 <- var2[pollinator==P]
    occ <- occ[pollinator==P]
    null.mat[unique(var1),
             unique(as.character(var2))][!is.na(null.mat)] <- occ
    null.mat
  }
  mats <- lapply(pollinator.id, function(x) make.mat(x))
  names(mats) <- pollinator.id
  mats
}



## abundance of each species for all site-date combinations
make.by.species <- function(spec,
                            sr.sched,
                            site.date,
                            type="PlantGenusSpecies"){
  ## number of species at each site-date
  all.sp <- aggregate(list(Abund=spec[, type]),
                      list(Site=spec$Site,
                           Date=spec$Date,
                           GenusSpecies=spec[,type]),
                      length)
  ## all site-date species combinations
  sp <- expand.grid(
    SiteDate = unique(paste(sr.sched$Site,
      sr.sched$Date,
      sep=';')),
    GenusSpecies=unique(spec[,type]),
    ## Occ=0,
    Abund=0)
  ## fill in the matrix of all combinations
  match.dates <- match(paste(all.sp$Site,
                       all.sp$Date,
                             all.sp$GenusSpecies,
                             sep=';'),
                       paste(sp$SiteDate,
                             sp$GenusSpecies,
                             sep=';'))
  match.abund <- match(paste(sp$SiteDate,
                             sp$GenusSpecies,
                             sep=';'),
                       paste(all.sp$Site,
                       all.sp$Date,
                             all.sp$GenusSpecies,
                             sep=';'))

  sp$Abund <- all.sp$Abund[match.abund]

  ## create site, date, genus etc. columns
  sp$Site <- sapply(strsplit(as.character(sp$SiteDate), ";"),
                    function(x) x[1])
  sp$Date <- as.Date(sapply(strsplit(as.character(sp$SiteDate), ";"),
                            function(x) x[2]))

  ## make a site by date by species arrary

  null.mat <- site.date
  pollinator.id <- id(sp$GenusSpecies)
  mats <- make.mats(pollinator.id,
                    null.mat,
                    pollinator=as.vector(sp$GenusSpecies),
                    var1=as.vector(sp$Site),
                    var2=sp$Date,
                    occ=sp$Abund)

  sites <- rownames(mats[[1]])
  dates <- colnames(mats[[1]])
  species <- names(mats)
  mat <- array(unlist(mats),
               dim=c(dim(mats[[1]]),
                 length(mats)))
  dimnames(mat) <- list(site=sites,
                        date=dates,
                        species=species)
  return(mat)
}



findOcc <- function(x){
  out <- try(occ[x["GenusSpecies"], paste(x["Site"], x["SiteStatusBACI"],
                                          sep=":")], silent=TRUE)
  if(inherits(out, "try-error")) out <- NA
  return(out)
}

findOccPlant <- function(x){
  out <- try(occ.plant[x["PlantGenusSpecies"], x["Site"]], silent=TRUE)
  if(inherits(out, "try-error")) out <- NA
  return(out)
}

calcOccArray <-  function(x){
  x[x > 1] <- 1
  sum(x > 0, na.rm=TRUE)/sum(x >= 0, na.rm=TRUE)
}
