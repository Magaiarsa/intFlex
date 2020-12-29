## this file loads and preps all of the specimen and trait data
## necessary for the analyses. It requires access to the raw data.
## rm(list=ls())
## setwd('~/Dropbox/speciesRoles/dataPrep')
load('~/Dropbox/hedgerow/data_sets/traditional/specimens-complete.RData')
source('src/misc.R')
source('src/prepNets.R')
source('src/contribNodf.R')
source('src/calcMotifPos.R')x
source('src/specialization.R')
library(bipartite)
library(fossil)
library(lubridate)
library(dplyr, plyr)
library(xtable)

## Assembling sites 
BACI <- c("Butler", "Hrdy", "MullerB", "Barger", "Sperandio")
## *******************************************************************
## load and prep data
## *******************************************************************
min.years <- 3 ## years to keep the sites

trait.dir <- '~/Dropbox/hedgerow/data_sets/traditional/functional_traits'
traits <- read.csv(file.path(trait.dir, 'bee.csv'), row.names=1)
spec <- dd
## subset to net specimens
spec <- spec[spec$NetPan == 'net',]
## create species column
spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                                 spec$PlantSpecies))
## drop pollinators and plants without identifications
spec <-  spec[spec$PlantGenusSpecies != '',]

## *******************************************************************
## quantities for manuscript - whole dataset
## *******************************************************************
## proportion bees for manuscript
tot.spec <- nrow(spec)
print("proprotion bee")
nrow(spec[spec$BeeNonbee == 'bee',])/tot.spec
spec <- spec[spec$BeeNonbee == 'bee',]
## drop species without an ID
spec <-  spec[spec$Species != '',]
## data used for analyses
## total specimens
print("total specimens")
nrow(spec)
print("total sites")
length(unique(spec$Site))
## total species
print("total bee species")
length(unique(spec$GenusSpecies))
print("total plant species")
length(unique(spec$PlantGenusSpecies))
## interactions
print("total interactions")
length(unique(paste(spec$GenusSpecies, spec$PlantGenusSpecies)))
## Total number of years each site was sampled
all.sites <- aggregate(list(n.samp.year=spec$Year),
                        list(Site=spec$Site, Status=spec$SiteStatus),
                       function(x) length(unique(x)))
## Unique combinations of site and year
print("total number of networks")
length(unique(paste(spec$Site, spec$Year)))

## *******************************************************************
## sampling table for manuscript - ALL sites
## *******************************************************************

site.table <- aggregate(list(Samples=spec$Date),
                        list(Year=spec$Year,
                             Site=spec$Site),
                        function(x) length(unique(x)))

ms.table <- samp2site.spp(site=site.table$Site,
                          spp=site.table$Year,
                          abund=site.table$Samples,
                          FUN=sum)
nYears <- apply(ms.table != 0, 1, sum)
print("mean times each site was sampled")
mean(rowSums(ms.table[,1:10])/nYears)
ms.table <- as.data.frame(cbind(ms.table, nYears))

write.csv(ms.table,"../data/samples.csv")

more.than.2 <- ms.table[which(ms.table$nYears>2),]
print("mean numer of years each site was sampled")
mean(more.than.2$nYears)
mean(ms.table$nYears)

## changing the tables site names and adding * to the baci ones
sites.table <- rownames(ms.table)
sites.S <- paste0("S", 1:length(sites.table))
sites.S[which(sites.table %in% BACI)] <- 
  paste0(sites.S[which(sites.table %in% BACI)], "*")
rownames(ms.table) <- sites.S
print(xtable(ms.table, type = "latex", digits=0), 
      file="../data/samples.txt")

## *******************************************************************
## subsetting the data for the spp and sites 
## to be included in the analysis
## *******************************************************************
## number of encounters per species
## *******************************************************************
all.spec <- spec ## saving before subset to calculate the traits from

## removing the sites that were samples less than 3 years
spec <- spec[spec$Site %in% rownames(more.than.2),]
## removing networks with any dim < 5 for the analysis 
num.dim <- 5
nets <- breakNet(spec.dat = spec, 'Site', 'Year')
site.year.nets <- names(nets)
spec$SiteYear <- paste(spec$Site, spec$Year, sep=".")
spec <- spec[spec$SiteYear %in% site.year.nets,]

## number of species per network
netws.size <- lapply(nets, dim)
netws.size <- lapply(netws.size, function(x) data.frame(row=x[1], col=x[2]))
netws.size <- as.data.frame(do.call(rbind, netws.size))
netws.size$site <- unlist(lapply(strsplit(rownames(netws.size), "[.]"), function(x) x[1]))
netws.size$year <- unlist(lapply(strsplit(rownames(netws.size), "[.]"), function(x) x[2]))
netws.size$total <- netws.size$row + netws.size$col

range(netws.size$total)
median(netws.size$total)
mean(netws.size$total)
mean(netws.size$row)
range(netws.size$col)

## save networks for each site
f.path <- '../data/networks'
save(nets, file=file.path(f.path, 'all_networks_years.Rdata'))

## number of encounters of each spp per site & year
spp.site.year <- spec %>%
  count(GenusSpecies, Site, Year, sort = TRUE) %>%
  arrange(GenusSpecies, Site, Year) 

## how many years in a site a spp was seen
n.spp.year <- spp.site.year %>%
  group_by(GenusSpecies, Site) %>%
  dplyr::summarise(nyears = n_distinct(Year))

## selecting spp & sites that were seen at least in 3 different years
spp.sites.to.keep <- n.spp.year[n.spp.year$nyears>2,]

## number of years each site was sampled
n.year.site <- spp.site.year %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(nyears = n_distinct(Year))

## selecting only sites that were samples at least twice  
n.year.site <- n.year.site[n.year.site$nyears>=min.years,]

## subset the spp in the sites that were sampled at least twice 
spp.sites.to.keep <- spp.sites.to.keep[spp.sites.to.keep$Site %in% n.year.site$Site,]

species.to.analyze <- spp.sites.to.keep ## it is required for the occupancy models

print("total number of species seen at least three years across all sites")
length(unique(spp.sites.to.keep$GenusSpecies))
print("total number of sites sampled at least three years")
length(unique(spp.sites.to.keep$Site))

## average number of of sampling rounds in the kept sites
n.sample.sites.keep <- site.table[site.table$Site %in% spp.sites.to.keep$Site,]
print("mean number of samples per year in the kept sites")
mean(n.sample.sites.keep$Samples)
print("mean number of years sampled")
mean(n.year.site$nyears)

save(species.to.analyze, file="../data/speciesSet.Rdata")

## number of encounters of each spp per site & year - Table S2
## spp, site, NYears, Nobs
table.S2 <- spp.site.year %>%
  dplyr::group_by(GenusSpecies, Site) %>%
  dplyr::summarise(NYears = n_distinct(Year),
                   NObs = sum(n))

## subseting to the included ones
table.S2.sppSite <- paste(table.S2$GenusSpecies, table.S2$Site, sep=".")
sppSite.keep <- paste(species.to.analyze$GenusSpecies, species.to.analyze$Site, sep=".")

table.S2 <- table.S2[table.S2.sppSite %in% sppSite.keep,]

## transforming the site names into S
sites.table.S1 <- as.data.frame(cbind("Site"= sites.table, "Sites.S"= sites.S))

table.S2 <- dplyr::inner_join(table.S2, sites.table.S1) %>%
  arrange(GenusSpecies)
table.S2 <- table.S2[,c("GenusSpecies", "Sites.S", "NYears", "NObs")]
colnames(table.S2) <- c("Species", "Site", "NYears", "NObs")
table.S2 <- as.data.frame(table.S2)
print(xtable(table.S2, type = "latex"),
      include.rownames=FALSE,
      file="../data/s2.txt")

## only at the BACI sites for the Structural Eq. Models
spp.baci <- species.to.analyze[species.to.analyze$Site %in% BACI,]
samples.baci <- n.sample.sites.keep[n.sample.sites.keep$Site %in% BACI,]
n.years.baci <- n.year.site[n.year.site$Site %in% BACI,]

print("Total number of spp in BACI sites")
length(unique(spp.baci$GenusSpecies))
print("mean number of samples per year in the kept sites")
mean(samples.baci$Samples)
print("range of years sampled")
range(n.years.baci$nyears)
print("mean number of years sampled")
mean(n.years.baci$nyears)

print("Total number of networks in BACI sites")
nets.baci <- spec[spec$Site %in% BACI,]
length(unique(paste(nets.baci$Site, nets.baci$Year, sep = ".")))
## *******************************************************************
## traits using the whole dataset
## *******************************************************************
## trait 1 & 2: pollinator lecty and body size 
## *******************************************************************
## Getting the lecty data
traits$Lecty[traits$Lecty == "parasite"] <- "polylectic"
traits$Lecty <- as.character(traits$Lecty)
traits$Lecty <- as.factor(traits$Lecty)
traits$Lecty <- as.numeric(traits$Lecty)
traits$Lecty[traits$Lecty == 2] <- 0 ## for path analysis
## subseting Lecty and mean body size
traits <- traits[,c("GenusSpecies", "Lecty", "MeanITD")]
traits.path <- traits
## *******************************************************************
## trait 3: phenological range
## *******************************************************************
## transforming date to date format
spp.dates <- as.Date(all.spec$Date,format="%Y-%m-%d")
day.observed <- day(spp.dates)
month.observed <- month(spp.dates)
day.month <- as.Date(paste(month.observed,day.observed, sep = "."),
                     format = "%m.%d")
## calculate the min and max date which spp was observed, then drop
## the year and calculate the period)
all.spec.dates <- cbind(all.spec, spp.dates, day.month)

## calculating the min and max date observed, per year
spp.max.min.days <- all.spec.dates %>%
    group_by(GenusSpecies, Year) %>%
    dplyr::summarise(min.day=min(day.month), max.day=max(day.month))%>%
    mutate(days.year=max.day - min.day)
spp.max.min.days$days.year <- as.numeric(spp.max.min.days$days.year)
spp.max.min.days$days.year[which(spp.max.min.days$days.year == 0)] <- 1

## calculating the min and max date observed, ignoring the year
spp.pheno.pol <- spp.max.min.days %>%
    group_by(GenusSpecies) %>%
    dplyr::summarise(min.day=min(min.day),
              max.day=max(max.day),
              total.days=max.day-min.day,
              max.days=max(days.year),
              min.days=min(days.year),
              mean.days=mean(days.year),
              median.days=median(days.year))

## merge all.specimen and phenology data
traits <- merge(traits, spp.pheno.pol, all=TRUE)

## *******************************************************************
## trait 4: abundance
## *******************************************************************
## abundance of each all.species at each site/sampling round
pol.SR <- aggregate(list(abund=all.spec$GenusSpecies),
                    list(GenusSpecies=
                             all.spec$GenusSpecies,
                                    site=all.spec$Site,
                                    SR=all.spec$Date),
                      length)
## create a all.species by site-date matrix.
pol.samp.mat <- samp2site.spp(spp=paste(pol.SR$site,
                                        pol.SR$SR, sep="."),
                              site=pol.SR$GenusSpecies,
                              abund=pol.SR$abund,
                              FUN= function(x) x)
## rarefy abundaunce using chao method
r.abund <- apply(pol.samp.mat, 1, chao1)

pol.year <- aggregate(list(abund=all.spec$GenusSpecies),
                      list(GenusSpecies=all.spec$GenusSpecies,
                                    site=all.spec$Site,
                                    year=all.spec$Year),
                      length)

spp.abund.pol <- pol.year %>%
    group_by(GenusSpecies) %>%
    dplyr::summarise(min.abund.year=min(abund),
              max.abund.year=max(abund),
              total.abund.across.yr=sum(abund),
              mean.abund=mean(abund),
              median.abund=median(abund))

spp.abund.pol$r.abund <- r.abund[match(spp.abund.pol$GenusSpecies,
                                       names(r.abund))]
spp.abund.pol$all.speciesType <- "pollinator"

## merge abundance and all.specimen data
traits <- merge(traits, spp.abund.pol, all=TRUE)

## *******************************************************************
## trait 5: diet breadth
## *******************************************************************
## pollinator visitation at each site across the entire dataset of all
## the sites in the region
agg.all.spec <- aggregate(list(abund=all.spec$GenusSpecies),
                      list(GenusSpecies=all.spec$GenusSpecies,
                           PlantGenusSpecies=all.spec$PlantGenusSpecies),
                      length)

## create a network where columns are pollinators and rows are plants
nets.all <- samp2site.spp(agg.all.spec$PlantGenusSpecies,
                          agg.all.spec$GenusSpecies,
                          agg.all.spec$abund,
                          FUN=sum)

all.all.specializations <- specieslevel(nets.all,
                                    index=c("degree"))

## calculate rarefied plant.pol degree
rare.pols.degree <- apply(nets.all, 2, chao1)
pol.degree <- data.frame(GenusSpecies= unlist(sapply(all.all.specializations,
                                                 rownames)),
                     do.call(rbind, all.all.specializations))
pol.degree$r.degree <-  rare.pols.degree[match(pol.degree$GenusSpecies,
                                           names(rare.pols.degree))]

rownames(pol.degree) <- NULL

traits <- merge(traits, pol.degree, all=TRUE)

## *******************************************************************
## traits for the PATH ANALYSIS
## *******************************************************************
## trait 1 & 2: pollinator lecty and body size - cant exclude BACI sites
## *******************************************************************

## *******************************************************************
## trait 3: phenological range
## *******************************************************************
## transforming date to date format
all.spec.path <- all.spec[!all.spec$Site%in%BACI,]
spp.dates <- as.Date(all.spec.path$Date,format="%Y-%m-%d")
day.observed <- day(spp.dates)
month.observed <- month(spp.dates)
day.month <- as.Date(paste(month.observed,day.observed, sep = "."),
                     format = "%m.%d")
## calculate the min and max date which spp was observed, then drop
## the year and calculate the period)
all.spec.path.dates <- cbind(all.spec.path, spp.dates, day.month)

## calculating the min and max date observed, per year
spp.max.min.days <- all.spec.path.dates %>%
  group_by(GenusSpecies, Year) %>%
  dplyr::summarise(min.day=min(day.month), max.day=max(day.month))%>%
  mutate(days.year=max.day - min.day)
spp.max.min.days$days.year <- as.numeric(spp.max.min.days$days.year)
spp.max.min.days$days.year[which(spp.max.min.days$days.year == 0)] <- 1

## calculating the min and max date observed, ignoring the year
spp.pheno.pol <- spp.max.min.days %>%
  group_by(GenusSpecies) %>%
  dplyr::summarise(min.day=min(min.day),
                   max.day=max(max.day),
                   total.days=max.day-min.day,
                   max.days=max(days.year),
                   min.days=min(days.year),
                   mean.days=mean(days.year),
                   median.days=median(days.year))

## merge all.spec.path and phenology data
traits.path <- inner_join(spp.pheno.pol, traits.path)

## *******************************************************************
## trait 4: abundance
## *******************************************************************
## abundance of each all.spec.pathies at each site/sampling round
pol.SR <- aggregate(list(abund=all.spec.path$GenusSpecies),
                    list(GenusSpecies=
                           all.spec.path$GenusSpecies,
                         site=all.spec.path$Site,
                         SR=all.spec.path$Date),
                    length)
## create a all.spec.pathies by site-date matrix.
pol.samp.mat <- samp2site.spp(spp=paste(pol.SR$site,
                                        pol.SR$SR, sep="."),
                              site=pol.SR$GenusSpecies,
                              abund=pol.SR$abund,
                              FUN= function(x) x)
## rarefy abundaunce using chao method
r.abund <- apply(pol.samp.mat, 1, chao1)

pol.year <- aggregate(list(abund=all.spec.path$GenusSpecies),
                      list(GenusSpecies=all.spec.path$GenusSpecies,
                           site=all.spec.path$Site,
                           year=all.spec.path$Year),
                      length)

spp.abund.pol <- pol.year %>%
  group_by(GenusSpecies) %>%
  dplyr::summarise(min.abund.year=min(abund),
                   max.abund.year=max(abund),
                   total.abund.across.yr=sum(abund),
                   mean.abund=mean(abund),
                   median.abund=median(abund))

spp.abund.pol$r.abund <- r.abund[match(spp.abund.pol$GenusSpecies,
                                       names(r.abund))]
spp.abund.pol$all.spec.pathiesType <- "pollinator"

## merge abundance and all.spec.pathimen data
traits.path <- inner_join(spp.abund.pol, traits.path)

## *******************************************************************
## trait 5: diet breadth
## *******************************************************************
## pollinator visitation at each site across the entire dataset of all
## the sites in the region
agg.all.spec.path <- aggregate(list(abund=all.spec.path$GenusSpecies),
                          list(GenusSpecies=all.spec.path$GenusSpecies,
                               PlantGenusSpecies=all.spec.path$PlantGenusSpecies),
                          length)

## create a network where columns are pollinators and rows are plants
nets.all <- samp2site.spp(agg.all.spec.path$PlantGenusSpecies,
                          agg.all.spec.path$GenusSpecies,
                          agg.all.spec.path$abund,
                          FUN=sum)

all.all.spec.pathializations <- specieslevel(nets.all,
                                        index=c("degree"))

## calculate rarefied plant.pol degree
rare.pols.degree <- apply(nets.all, 2, chao1)
pol.degree <- data.frame(GenusSpecies= unlist(sapply(all.all.spec.pathializations,
                                                     rownames)),
                         do.call(rbind, all.all.spec.pathializations))
pol.degree$r.degree <-  rare.pols.degree[match(pol.degree$GenusSpecies,
                                               names(rare.pols.degree))]

rownames(pol.degree) <- NULL
traits.path <- inner_join(traits.path, pol.degree)

## *******************************************************************
## Analyses of nestedness contribution & Motif
## *******************************************************************
## nestedness contribution
contr.nodf <- calcNofContrib(nets)
contr.nodf$sppSite <- paste(contr.nodf$GenusSpecies, contr.nodf$Site, sep=".")

## motif position
motif.position.Q <- calcMotifPos(nets, bin=FALSE)

## subsetting to the spp x site combinations
motif.position.Q <- motif.position.Q[motif.position.Q$sppSite %in% sppSite.keep,]
motif.position.Q <- motif.position.Q %>% replace(is.na(.), 0) ## replacing NAs

contr.nodf <- contr.nodf[contr.nodf$sppSite %in% sppSite.keep,]

save(contr.nodf, file='../data/nestedContribution.Rdata')
save(motif.position.Q, file='../data/motifPositionQ.Rdata')

## save the final specimen data
spec$siteYear <- paste(spec$Site, spec$Year, sep=".")
save(spec, file='../data/networks/allSpecimens.Rdata') ## used in the partner nulls
save(traits,  file='../data/traits.Rdata')
save(traits.path,  file='../data/traitsPath.Rdata')
