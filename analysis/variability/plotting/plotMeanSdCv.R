## setwd('~/Dropbox/speciesRoles')
## rm(list=ls())
setwd('analysis/variability')
f.path <- 'figures'

library(ggplot2)
load(file="saved/partnerQ.RData") #Partner varialbility
load(file="saved/roleQ.RData") #role variability
load(file="saved/cnodf.RData") #contribution to nodf variabiity
load('../../data/nestedContribution.Rdata')
load('../../data/speciesSet.Rdata')

## selecting sites that were sampled at least 3x 
## and prepping the data for plotting
partner.var <- base::merge(partner.var, species.to.analyze, all.x=FALSE)
partner.plot <- partner.var[,c("GenusSpecies", "Site", "cv.partner", "nyears")]
partner.plot$var <- "partner"
colnames(partner.plot)[3] <- "value"

pca.var <- merge(pca.var, species.to.analyze, all.x=FALSE)
role.plot <- pca.var[,c("GenusSpecies", "Site", "pca.cv", "nyears")]
role.plot$var <- "role"
colnames(role.plot)[3] <- "value"

cnodf.var <- merge(cnodf.var, species.to.analyze, all.x=FALSE)
cnodf.plot <- cnodf.var[,c("GenusSpecies", "Site", "cnodf.cv", "nyears")]
cnodf.plot$var <- "cnodf"
colnames(cnodf.plot)[3] <- "value"

all.vars <- rbind(partner.plot, role.plot, cnodf.plot)

###############################################
## Ploting
###############################################
all.vars$var <- as.factor(all.vars$var)
levels(all.vars$var) <- c(partner = "Partner variability", 
                             role = "Role variability", 
                             cnodf = "Structural variability")

all.vars$nyears <- log(all.vars$nyears)

ggplot(all.vars, aes(x=nyears, y=value, group=var)) + 
  geom_point() + 
  theme_bw() + 
    facet_grid(.~ var) + 
  ylab("Total observations across years (log)") + xlab("Coefficient of variation") 

ggsave(filename = file.path(f.path, sprintf("%s.pdf", "CvAcrossYears20")),
       width=10, height=6)
