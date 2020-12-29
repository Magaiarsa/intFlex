## modified version of the partner calculations
## considers the nulls to calculate the beta diversity
library(vegan)
library(dplyr)

calcBetaStatusMotif <- function(comm,
                                dis.method,
                                occ= FALSE){ ## binary?
  ## computes dispersion of community matrices, returns output of
  ## vegan function
  ## create community dissimilarity matrix
  comm.dis <-  lapply(comm, function(x) {
    vegdist(x, method= dis.method, add=TRUE, diag= TRUE, binary= occ)
  })
  beta.disper.result <- vector("list", length(comm))
  for(i in 1:length(comm)){
    this.spp <- names(comm)[i]
    cor.dis  <-  comm.dis[[i]]
  
    if(sum(colSums(as.matrix(cor.dis)))==0){
      beta.disper.result[[i]] <- NA
      names(beta.disper.result)[i] <- this.spp
    }else{
    beta.disper.result[[i]] <- try(betadisper(cor.dis, rep("all", nrow(as.matrix(cor.dis))),
                                              type="centroid", add=TRUE), silent = TRUE)
    if(inherits(beta.disper.result[[i]], "try-error"))browser()
    names(beta.disper.result)[i] <- this.spp
    }
  }
  return(beta.disper.result)
}

makeBetaDataPrettyMotifs <- function(res.motif){
  distances <- lapply(res.motif, function(x){
    lapply(x, function(y) y$distances)})
  dist <- unlist(distances)
  GenusSpecies <- unlist(lapply(strsplit(names(dist), "[.]"), function(x) x[2]))
  Site <-  unlist(lapply(strsplit(names(dist), "[.]"), function(x) x[1]))
  dats <- cbind(as.data.frame(dist), GenusSpecies, Site)
  rownames(dats) <- NULL
  return(dats)
}
