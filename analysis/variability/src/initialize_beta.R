library(vegan)
library(dplyr)
source('src/misc.R')
source('src/beta.R')

if(length(args) == 0){
    type <- "pols"
} else{
    type <- args[1]
}

if(type == "pols"){
    speciesType <- "pollinator"
} else{
    speciesType <- "plants"
}

if(!binary & alpha){
  print(paste(binary, alpha))
  occ <- "abund"
  print(occ)
  dis.method <- dis.method
  load(file=file.path('saved/communities',
         sprintf('%s-abund.Rdata', type)))
  load(file=file.path('saved/nulls',
         sprintf('%s-alpha.Rdata', type)))
}

if(!binary & !alpha){
  print(paste(binary, alpha))
  occ <- "indiv"
  print(occ)
  dis.method <- dis.method
  load(file=file.path('saved/communities',
         sprintf('%s-abund.Rdata', type)))
  load(file=file.path('saved/nulls',
         sprintf('%s-indiv.Rdata', type)))
}

if(binary){
  print(paste(binary, alpha))
  occ <- "occ"
  print(occ)
  dis.method <- dis.method
  load(file=file.path('saved/communities',
         sprintf('%s-abund.Rdata', type)))
  load(file=file.path('saved/nulls',
         sprintf('%s-occ.Rdata', type)))
}

if(type=="pols"){
  ylabel <- "Pollinator species turnover"
}
if(type=="ints"){
  ylabel <- "Interaction turnover"
}
if(type=="plants"){
  ylabel <- "Plant species turnover"
}