library(ggplot2)
library(dplyr)
library(ggfortify)
library(vegan)

f.path <- 'figures'
s.path <- 'saved'

source('src/misc.R')
load('../../data/networks/allSpecimens.Rdata')
load('../../data/networks/all_networks_years.Rdata')
#load('../../data/species_roles.Rdata')
load('../../data/traits.Rdata')
# 
# species.roles$SiteYr <- paste(species.roles$Site,
#                               species.roles$Year, sep='.')
# 
# ## combine GenusSpecies + site + year
# species.roles$sp.site.year <- paste(species.roles$GenusSpecies,
#                                     species.roles$SiteYr, sep='.')
# 
# species.roles$SiteStatus <-
#     spec$SiteStatus[match(species.roles$SiteYr,
#                           paste(spec$Site, spec$Year, sep="."))]
# 
# BACI <- c("Barger", "Butler", "MullerB", "Sperandio", "Hrdy")
# species.roles$SiteStatus[species.roles$Site %in% BACI] <- "hedgerow"
# 
# species.roles$SiteStatus[species.roles$SiteStatus == "mature" |
#                          species.roles$SiteStatus == "maturing" |
#                          species.roles$SiteStatus == "forb"] <- "hedgerow"
