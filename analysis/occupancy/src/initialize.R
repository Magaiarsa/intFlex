library(nimble)
library(igraph)
library(parallel)
library(dplyr)
library(reshape2)
library(xtable)
source("src/dynamicOcc.R")
source("src/checkChains.R")
source("src/plot_posteriors_all.R")
source('src/setup_data_all.R')
source("src/misc.R")
source("src/setup.R")
source('src/model.R')
source("src/prep.R")
source("src/make-matrix.R")
source("src/prepOccupancyData.R")
source("src/testingDimModelData.R")

save.dir <- "../../../speciesRoles_saved/occupancy"
f.path <- 'figures'
hedgerow.dir <- "../../data"
## spatial data
geo <- read.csv(file.path(hedgerow.dir, 'tables/geography.csv'),
                as.is=TRUE)

## sampling schedule
sr.sched <- read.csv(file.path(hedgerow.dir,
                               'tables/conditions.csv'),
                     as.is=TRUE)
sr.sched$Date <- as.Date(sr.sched$Date)

sr.sched$Site <- geo$Site[match(sr.sched$GeographyFK,
                                geo$GeographyPK)]

## veg data
load('../../data/tables/veg.Rdata')
load('../../data/networks/allSpecimens.Rdata')
load('../../data/traits.Rdata')

## for the prepOccupancy
load(file="../variability/saved/cnodf.RData") ##only binary

if(binary){
  print("loading binary Rdata")
  load(file="../variability/saved/partner.RData")
  load(file="../variability/saved/role.RData")
  tipo <- "B"
} else {
  print("loading quantitative Rdata")
  load(file="../variability/saved/partnerQ.RData")
  load(file="../variability/saved/roleQ.RData")
  tipo <- "Q"
}

load('../../data/speciesSet.Rdata')

var.model <- prepOccupancyData()
