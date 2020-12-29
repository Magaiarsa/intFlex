library(piecewiseSEM)
library(dplyr)
library(nlme)
library(xtable)
library(car)
## all this for the plotting
library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)

source("src/misc.R")
load(file="saved/cnodf.RData")

if(binary){
  print("loading binaryRdata")
  load(file="../variability/saved/partner.RData")
  load(file="../variability/saved/role.RData")
}else{
  print("loading quantitative Rdata")
  load(file="../variability/saved/partnerQ.RData")
  load(file="../variability/saved/roleQ.RData")
}

load('../../data/speciesSet.Rdata')
load('../../data/traitsPath.Rdata')

BACI <- c("Barger", "Butler", "MullerB", "Sperandio", "Hrdy")

## merging the traits of interest with the variability measures
cv.path <- c("GenusSpecies", "Site", "cv", "median.abund", "median.days",
             "r.degree", "Lecty", "MeanITD")

mean.path <- c("GenusSpecies", "Site", "Mean", "median.abund", "median.days",
             "r.degree", "Lecty")

if(!use.mean){
  partner <- dplyr::inner_join(partner.var, traits.path) 
  partner <- partner[,grepl(paste(cv.path, collapse="|"),
                            colnames(partner))]
  
  role <- dplyr::inner_join(pca.var, traits.path) 
  role <- role[,grepl(paste(cv.path, collapse="|"),
                      colnames(role))]
  
  cnodf <- dplyr::inner_join(cnodf.var, traits.path) 
  cnodf <- cnodf[,grepl(paste(cv.path, collapse="|"),
                              colnames(cnodf))]
} else{
  partner <- dplyr::inner_join(partner.var, traits.path) 
  partner <- partner[,grepl(paste(mean.path, collapse="|"),
                            colnames(partner))]
  
  role <- dplyr::inner_join(pca.var, traits.path) 
  role <- role[,grepl(paste(mean.path, collapse="|"),
                      colnames(role))]
  
  cnodf <- dplyr::inner_join(cnodf.var, traits.path) 
  cnodf <- cnodf[,grepl(paste(mean.path, collapse="|"),
                                  colnames(cnodf))]
}
## subset to only assembling hedgerows (BACI sites, before, after
## control impact) for path analyses to avoid having the
## explanatory variables (abundance etc.) calculated from exactly the
## same data as the variability. We chose the assembling hedgerows
## because since the community context is shifting, species are most
## likely to be variable In addition subset to species seen more than
## three times at each BACI site

partner <- partner[partner$Site %in% BACI,]
role <- role[role$Site %in% BACI,]
cnodf <- cnodf[cnodf$Site %in% BACI,]

##############################
## creating a path analysis object
path.variables <- list()
## to make the 5piecewiseSEM code prettier, changing mean by cv here
if(use.mean){
  colnames(partner)[grep("Mean.partner", colnames(partner))] <- "cv.partner"
  colnames(role)[grep("pca.Mean", colnames(role))] <- "pca.cv"
  colnames(cnodf)[grep("cnodf.Mean", colnames(cnodf))] <- "cnodf.cv"
}

path.variables$partner <- partner
path.variables$role <- role[complete.cases(role),]
path.variables$cnodf <- as.data.frame(cnodf)


## standardizing everyone
for(i in 1: length(path.variables)){
  path.variables[[i]] <-
    path.variables[[i]] %>%
  mutate_if(is.numeric, standardize)}

median.spp <- path.variables
for(i in 1: length(median.spp)){
  median.spp[[i]] <- median.spp[[i]]   %>%
    group_by(Site, GenusSpecies)  %>%
    summarize_all(median)}