## rm(list=ls())
## setwd('~/Dropbox/speciesRoles')
## load("../../data/motifPosition.Rdata")
## setwd('analysis/variability')
library(ggridges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)

## plot for each spp, one ridge per site of the values of motif postions
motif.position.ggplot <- reshape2::melt(motif.position[,1:49], 
                                        id.vars = c("spp", "site", "year"))
motif.position.ggplot.bin <- reshape2::melt(motif.position.bin[,1:151], 
                                        id.vars = c("spp", "site", "year"))
##-------------------------------------------------------------------
## plot the motif position of the spp of interest - quantitative
motif.spp.ggplot <- function(dados, spp, site){
  #browser()
  spp.interst <- dados[dados$spp==spp,]
  if(!is.na(site)){
    spp.interst <- spp.interst[spp.interst$site==site,]
  }
  spp.interst <- spp.interst[!is.na(spp.interst$value),]

  plote <- ggplot(spp.interst, aes(y=value, x=year, fill=year)) + 
    geom_col() + facet_grid(site~variable)+
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                     vjust=0.5, size=12, family="Times"),
          axis.text.y = element_text(size=12, family="Times"),
          legend.position="none",
          strip.text.x = element_text(size = 12, family="Times")) +
    xlab("") + ylab(paste("Frequency of motif position", 
                          spp, sep=" "))
  return(print(plote))
}


## plot the motif position of the spp of interest - binary
motif.spp.ggplot.bin <- function(dados, spp, site){
  #browser()
  spp.interst <- dados[dados$spp==spp,]
  if(!is.na(site)){
    spp.interst <- spp.interst[spp.interst$site==site,]
  }
  spp.interst <- spp.interst[spp.interst$value!=0,]
  
  plote <- ggplot(spp.interst, aes(y=value, x=year, fill=year)) + 
    geom_col() + facet_grid(site~variable)+
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                     vjust=0.5, size=12, family="Times"),
          axis.text.y = element_text(size=12, family="Times"),
          legend.position="none",
          strip.text.x = element_text(size = 12, family="Times")) +
    xlab("") + ylab(paste("Frequency of motif position", 
                          spp, sep=" "))
  return(print(plote))
}

motif.spp.ggplot.years <- function(dados, spp, site, QB){
  spp.interst <- dados[dados$spp==spp,]
  if(!is.na(site)){
    spp.interst <- spp.interst[spp.interst$site==site,]
  }
  spp.interst <- spp.interst[!is.na(spp.interst$value),]
  if(QB=="B"){
    spp.interst <- spp.interst[spp.interst$value>0,]
  }
  plote <- ggplot(spp.interst, aes(y=value, x=site)) + 
    geom_col() + facet_grid(year~variable)+
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                     vjust=0.5, size=12, family="Times"),
          axis.text.y = element_text(size=12, family="Times"),
          legend.position="none",
          strip.text.x = element_text(size = 12, family="Times")) +
    xlab("") + ylab(paste("Frequency of motif position", 
                          spp, sep=" "))
  return(print(plote))
}

# areas.na.bin
# g5 <- motif.spp.ggplot.bin(motif.position.ggplot.bin, "Andrena candida", "Rominger")
# g6 <- motif.spp.ggplot.bin(motif.position.ggplot.bin, "Bombus melanopygus", "Barger")
# g7 <- motif.spp.ggplot.bin(motif.position.ggplot.bin, "Lasioglossum (Dialictus) impavidum", "Barger")
# g8 <- motif.spp.ggplot.bin(motif.position.ggplot.bin, "Megachile apicalis", "Gregory")
# 
# g5/g6/g7/g8
# ggsave(filename = "figures/motifs/motif_bin_problems.pdf",
#        width = 10, height = 15)

# areas.na
# g1 <- motif.spp.ggplot(motif.position.ggplot, "Andrena candida", "Rominger")
# g2 <- motif.spp.ggplot(motif.position.ggplot, "Bombus melanopygus", "Barger")
# g3 <- motif.spp.ggplot(motif.position.ggplot, "Megachile apicalis", "H16")
# g4 <- motif.spp.ggplot(motif.position.ggplot, "Svastra obliqua expurgata", "MullerM")
# g1/g2/g3/g4
# ggsave(filename = "figures/motifs/motif_quant_problems.pdf",
#        width = 10, height = 15)
# 

# 
motif.spp.ggplot.years(motif.position.ggplot.bin, "Lasioglossum (Dialictus) impavidum", "Barger", "B")
ggsave(filename = "figures/motifs/lassImpav_Barguer.pdf",
       width = 10, height = 7)

motif.spp.ggplot.years(motif.position.ggplot.bin, "Megachile apicalis", "Gregory", "B")
ggsave(filename = "figures/motifs/megaAcp_Gregory.pdf",
       width = 10, height = 7)

motif.spp.ggplot.years(motif.position.ggplot.bin, "Halictus tripartitus", "MullerB", "B")
ggsave(filename = "figures/motifs/halicTrip_MullerB.pdf",
       width = 20, height = 10)
