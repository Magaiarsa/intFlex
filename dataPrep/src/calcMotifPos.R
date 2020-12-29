# require(devtools)
# install.packages("Rcpp")
# devtools::install_github("SimmonsBI/bmotif", build_vignettes = TRUE) # install bmotif

library(bmotif)
## applies the bmotif::function node_positions  to bipartite to networks
## calculates the motif position of each spp at each site/year
## returns a dataframe with the spp.site.year as rows 

calcMotifPos <- function (networks, bin){
  calc.motif <- function(x) {
    if(bin){
    node_positions(x, level="columns", 
                   six_node = FALSE,
                   weights_method = "none", 
                   weights_combine = "none", 
                   normalisation = "sum")
    }else{
      node_positions(x, level="columns", 
                     six_node = FALSE,
                     weights_method = "mean_motifweights", 
                     weights_combine = "mean", 
                     normalisation = "none")}
      }
  motifs.all <- lapply(networks, calc.motif)

  motifs.all <- do.call(rbind, motifs.all)
  
  motifs.all$site <- unlist(lapply(strsplit(rownames(motifs.all), "[.]"),
                                   function(x)
                                     x[1]))
  motifs.all$spp <- unlist(lapply(strsplit(rownames(motifs.all), "[.]"),
                                  function(x)
                                    x[3]))
  motifs.all$year <- unlist(lapply(strsplit(rownames(motifs.all), "[.]"),
                                   function(x)
                                     x[2]))
  motifs.all$sppSite <- paste(motifs.all$spp, motifs.all$site, sep = ".")
  return(motifs.all)
}
