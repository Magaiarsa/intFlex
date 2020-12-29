## setwd('~/Dropbox/speciesRoles')
rm(list=ls())
setwd('analysis/occupancy')
source('src/initialize.R')
save.dir <- "../../../speciesRoles_saved/occupancy"
all.types <- c("partner", "role", "cnodf")#, "all", "degree")

all.samples <- list()

for(i in 1:length(all.types)){
  load(file=file.path(save.dir,
                      sprintf("ms-ms-samples-%s.Rdata",
                              all.types[i])))
  samples <- do.call(rbind, ms.ms.nimble$samples)
  samples.means <- apply(samples[-c(1:100),], 2, mean)
  all.samples[[i]] <- samples.means[grepl(paste(c("var", "mean"), collapse="|"), names(samples.means))]
  print(i)
}

names(all.samples) <- all.types
