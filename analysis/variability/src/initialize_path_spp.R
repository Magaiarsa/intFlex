library(piecewiseSEM)
library(dplyr)
library(nlme)
library(xtable)
load('../../data/speciesSet.Rdata')

BACI <- c("Barger", "Butler", "MullerB", "Sperandio", "Hrdy")

load('../../data/traits.Rdata')
## cv of partner per species across the BACI sites ----------
load("saved/partnerBin.RData")
dis <- dis[dis$Site %in% BACI,]

partner.var.spp <- dis  %>%
  group_by(GenusSpecies)  %>%
  summarize(cv.partner = cv(dist),
            mean.partner = mean(dist))

partner.var.spp <- merge(partner.var.spp, traits, by = "GenusSpecies")

## mean area (role variability) per species across the landscape ----------
load("saved/role.RData")
pca.var <- pca.var[pca.var$Site %in% BACI,]

role.var.spp <- pca.var  %>%
  group_by(GenusSpecies)  %>%
  summarize(cv.role = cv(pca.cv),
            mean.role = mean(pca.cv))
## some spp only have one obs and thus their cv is NA, but should be the original value
spp.nas <- role.var.spp$GenusSpecies[which(is.na(role.var.spp$cv.role))]
role.var.spp[which(is.na(role.var.spp$cv.role)),"cv.role"] <- pca.var$pca.cv[match(spp.nas, pca.var$GenusSpecies)]
role.var.spp <- merge(role.var.spp, traits, by = "GenusSpecies")

## cv of cnodf per species across the landscape ----------
load('../../data/nestedContribution.Rdata')
contr.nodf  <- contr.nodf [contr.nodf$Site %in% BACI,]

## calculating mean and sd per species per site
cnodf.var.spp <- contr.nodf %>%
  group_by(GenusSpecies)  %>%
  dplyr::summarize(
    cnodf.cv = cv(nestedcontribution),
    cnodf.mean = mean(nestedcontribution,
                      na.rm = TRUE))

cnodf.var.spp <- merge(cnodf.var.spp, traits, by="GenusSpecies",
                       all.x=TRUE)

## getting only the species to analyze 
partner.var.spp <- partner.var.spp[partner.var.spp$GenusSpecies %in% species.to.analyze$GenusSpecies,]
role.var.spp <- role.var.spp[role.var.spp$GenusSpecies %in% species.to.analyze$GenusSpecies,]
cnodf.var.spp <- cnodf.var.spp[cnodf.var.spp$GenusSpecies %in% species.to.analyze$GenusSpecies,]


## getting only the variables of interest
vars.cv <- c("GenusSpecies", "mean.abund", "mean.days",
             "r.degree", "Lecty", "MeanITD", "cv")
vars.mean <- c("GenusSpecies", "mean.abund", "mean.days",
               "r.degree", "Lecty", "MeanITD", "mean")

##############################
## Getting only the necessary columns for path analysis
path.variables <- list()

#### partner ####
path.variables$partner <- partner.var.spp[,grepl(paste(vars.cv, collapse="|"),
                                             colnames(partner.var.spp))]
#### role motif  ######
path.variables$role.cv <- role.var.spp[,grepl(paste(vars.cv, collapse="|"),
                                         colnames(role.var.spp))]
#### cnodf ####
path.variables$cnodf.mean <- cnodf.var.spp[,grepl(paste(vars.mean, collapse="|"),
                                              colnames(cnodf.var.spp))]
path.variables$cnodf.cv <- cnodf.var.spp[,grepl(paste(vars.cv, collapse="|"),
                                            colnames(cnodf.var.spp))]

## standardizing everyone
for(i in 1: length(path.variables)){
  path.variables[[i]] <-
    path.variables[[i]] %>%
    mutate_if(is.numeric, standardize)}

