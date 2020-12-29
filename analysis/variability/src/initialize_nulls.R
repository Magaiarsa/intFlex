load('../../data/networks/allSpecimens.Rdata')
source('src/commPrep.R')
save.dir.comm <- "saved/communities"
save.dir.nulls <- "saved/nulls"
source('src/misc.R')
source('src/vaznull2.R')
load("../../data/networks/all_networks_years.Rdata")

sites <- unique(spec$Site)
spec$Int <- paste(spec$GenusSpecies,
                  spec$PlantGenusSpecies)

args <- commandArgs(trailingOnly=TRUE)
type <- args[1]
nnull <- args[2]

if(length(args) == 0){
    type <- "pols"
    nnull <- nnull
}

if(type=="pols"){
    species.type="GenusSpecies"
    species.type.int="PlantGenusSpecies"
}

if(type=="plants"){
    species.type="PlantGenusSpecies"
    species.type.int="GenusSpecies"
}
