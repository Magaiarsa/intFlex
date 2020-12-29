load("../../data/speciesSet.Rdata")
species.to.analyze <- paste(species.to.analyze$GenusSpecies, species.to.analyze$Site, sep=".")

library(bmotif)
## renaming so it saves the sites, spp and simulation number
if(binary){
  for(i in 1: length(nulls.motif)){
    nome <- names(nulls.motif)[i]
    names(nulls.motif[[i]]) <- paste(nome, 1:nnull, sep=".")
  }
}else{
  for(i in 1: length(nulls.motif.q)){
    nome <- names(nulls.motif.q)[i]
    names(nulls.motif.q[[i]]) <- paste(nome, 1:nnull, sep=".")
  }
}


## function to calculate the motif position for each spp in each null community
## and wrangles the output
calcMotifPosNull <- function (networks, binary){
  calc.motif <- function(x) {
    if(binary){
      node_positions(x, level="columns", 
                     six_node = FALSE,
                     weights_method = "none", 
                     weights_combine = "none", 
                     normalisation = "none")
    } else {
      node_positions(x, level="columns", 
                     six_node = FALSE,
                     weights_method = "mean_motifweights", 
                     weights_combine = "sum", 
                     normalisation = "none")}
  }
  motifs.all <- lapply(networks, calc.motif)
  motifs.all <- do.call(rbind, motifs.all)
  
  motifs.all$site <- unlist(lapply(strsplit(rownames(motifs.all), "[.]"),
                                   function(x)
                                     x[1]))
  motifs.all$spp <- unlist(lapply(strsplit(rownames(motifs.all), "[.]"),
                                  function(x)
                                    x[4]))
  motifs.all$year <- unlist(lapply(strsplit(rownames(motifs.all), "[.]"),
                                   function(x)
                                     x[2]))
  motifs.all$sim <- unlist(lapply(strsplit(rownames(motifs.all), "[.]"),
                                  function(x)
                                    x[3]))
  motifs.all$sppSite <- paste(motifs.all$spp, motifs.all$site, sep = ".")
  return(motifs.all)
}

make.list.pretty <- function(lista){
  y <- plyr::ldply(lista)
  y$year <- sapply(strsplit(y$.id, "[.]"),
                   function(x)x[2])
  return(y)
}