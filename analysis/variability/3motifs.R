## setwd('~/Dropbox/speciesRoles')
rm(list=ls())
setwd('analysis/variability')

## loading the data calculated in the dataPrep file

load("../../data/motifPositionQ.Rdata") ##empirical
motif.position <- motif.position.Q
binary <- FALSE ## for the dissimilarity metric below
## dissimilarity metric
dis.method <- "altGower"

source("src/betaMotif.R")
source('src/misc.R')
s.path <- 'saved'
## this script calculates the "motif beta diversity" using Marti
## Anderson's multivariate dispersion method. It uses the null
## communities generated in the nulls.R script to control for alpha
## diversity in the beta diversity calculation

## organizing the empirical data
## making a list of sites with species as elements, and motifs position x years 
nps <- grep("np", colnames(motif.position))
comm.motif <- split(motif.position, f = motif.position$site)
comm.motif <- lapply(comm.motif, function(x) split(x[,nps], f = x$spp))

## this function applies vegdist to each year x motif position combination and 
## returns the betadisper outcome 
dis.motif <- mapply(function(a)
  calcBetaStatusMotif(comm= a, ## observed communities
                      dis.method, ## dissimilarity metric
                      occ=binary), ## binary or abundance weighted
  a= comm.motif,
  SIMPLIFY=FALSE)

dis.motif <- lapply(dis.motif, function(x) x[which(is.na(x)==FALSE)])
motifBeta <- makeBetaDataPrettyMotifs(dis.motif)

## calculate the sd per site and cv per species per site (across years)
pca.var <- motifBeta  %>%
  group_by(Site, GenusSpecies)  %>%
  summarize(pca.cv = cv(dist),
            pca.sd = sd(dist),
            pca.Mean = mean(dist))

save(pca.var, motifBeta,
       file=file.path(s.path,"roleQ.RData"))
