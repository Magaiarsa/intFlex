## setwd('~/Dropbox/speciesRoles')
rm(list=ls())
library(randomcoloR); library(ggplot2); library(dplyr); library(RColorBrewer)
setwd('analysis/variability')
source('src/misc.R')
f.path <- 'figures'

######## Plotting the cvs used in the occupancy model
load(file="saved/partnerQ.RData") #Partner variability
load(file="saved/roleQ.RData") #role variability
load(file="saved/cnodf.RData") #contribution to nodf variability
load('../../data/nestedContribution.Rdata') ## raw data of cnodf
###############################################
## Plotting all the cvs per species
###############################################
partner.gg.cv <- partner.var[,c("GenusSpecies", "Site", "cv.partner")]

role.gg.cv <- pca.var[,c("GenusSpecies", "Site", "pca.cv")]

str.gg.cv <- cnodf.var[,c("GenusSpecies", "Site", "cnodf.cv")]

cv.ggplot <- plyr::join_all(list(partner.gg.cv, role.gg.cv, str.gg.cv), 
                            by=c("GenusSpecies", "Site"), type='left')
cv.ggplott <- cv.ggplot[complete.cases(cv.ggplot),]

## organizing data for ggplot
cv.ggplott <- plyr::rename(cv.ggplott, c("cv.partner"="Partner", 
                                         "pca.cv"="Role",
                                         "cnodf.cv"="Structural"))
cv.ggplott <- reshape2::melt(cv.ggplott, id.vars = c("GenusSpecies", "Site"))

cv.ggplott <- cv.ggplott[!duplicated(cv.ggplott),]

## getting color pallete
palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(cv.ggplott$Site)))

ggplot(cv.ggplott, aes(x=factor(GenusSpecies, levels = sort(unique(GenusSpecies), decreasing = TRUE)), 
                       y=value, fill=Site)) + 
  geom_jitter(width = 0.01, size=2.5, shape=21, color="black") + 
  scale_fill_manual(values=c(palette))+
  theme_bw() + 
  facet_grid(cols=vars(variable))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                   vjust=0.5, size=12, family="Times"),
        axis.text.y = element_text(size=10, family="Times", 
                                   face="italic"),
        legend.position="none",
        strip.text.x = element_text(size = 12, family="Times")) + 
  coord_flip() + xlab("") + ylab("") #+ ylim(-3, 3)

ggsave(filename = file.path(f.path, sprintf("%s.pdf", "cvValues")),
       width=12, height=8)

###############################################
#### Plotting the raw values per species ----
###############################################
raw.partner <- dis.beta[,c("GenusSpecies", "Site", "dist")]
colnames(raw.partner)[3] <- "Partner"
#partner.gg.cv$cv.partner <- standardize(partner.gg.cv$cv.partner)
raw.role <- motifBeta[,c("GenusSpecies", "Site", "dist")]
colnames(raw.role)[3] <- "Role"
raw.cnodf <- contr.nodf[,c("GenusSpecies", "Site", "nestedcontribution")]
colnames(raw.cnodf)[3] <- "Structural"
raw.ggplot <- plyr::join_all(list(raw.partner, raw.role, raw.cnodf), 
                            by=c("GenusSpecies", "Site"), type='left')

raw.ggplot <- reshape2::melt(raw.ggplot, id.vars = c("GenusSpecies", "Site"))
raw.ggplot <- raw.ggplot[!duplicated(raw.ggplot),]

#cv.ggplott$value <- as.numeric(as.factor(cv.ggplott$value))
## getting color pallete
palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(raw.ggplot$Site)))

ggplot(raw.ggplot, aes(x=GenusSpecies, y=value, fill=Site)) + 
  geom_jitter(width = 0.01, size=2.5, shape=21, color="black") + 
  scale_fill_manual(values=c(palette))+
  theme_bw() + 
  facet_grid(cols=vars(variable), scales="free_x") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                   vjust=0.5, size=12, family="Times"),
        axis.text.y = element_text(size=10, family="Times", 
                                   face="italic"),
        legend.position="none",
        strip.text.x = element_text(size = 12, family="Times")) + 
  coord_flip() + xlab("") + ylab("") #+ ylim(-3, 3)

ggsave(filename = file.path(f.path, sprintf("%s.pdf", "rawData")),
       width=12, height=8)


