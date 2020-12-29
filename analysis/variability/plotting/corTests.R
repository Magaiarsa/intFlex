## rm(list=ls())
setwd('~/Dropbox/speciesRoles')
## load the data 
setwd('analysis/variability')
load(file="../variability/saved/cnodf.RData") 
load(file="../variability/saved/partnerQ.RData")
load(file="../variability/saved/roleQ.RData")

## ordering everyone for testing the correlation
cnodf.var$SppSite <- paste(cnodf.var$GenusSpecies, cnodf.var$Site, sep=".")

partner.var$SppSite <- paste(partner.var$GenusSpecies, partner.var$Site, sep=".")
partner.var <- partner.var[partner.var$SppSite %in% cnodf.var$SppSite, ]
partner.var <- partner.var[order(match(partner.var$SppSite, cnodf.var$SppSite)),]

pca.var$SppSite <- paste(pca.var$GenusSpecies, pca.var$Site, sep=".")
pca.var <- pca.var[pca.var$SppSite %in% cnodf.var$SppSite, ]
pca.var <- pca.var[order(match(pca.var$SppSite, cnodf.var$SppSite)),]

cor.test(partner.var$cv.partner, pca.var$pca.cv)
cor.test(partner.var$cv.partner, cnodf.var$cnodf.cv)
cor.test(cnodf.var$cnodf.cv, pca.var$pca.cv)

all.cvs <- cbind(partner.var$cv.partner, pca.var$pca.cv, cnodf.var$cnodf.cv)
cor(all.cvs)
cor.cvs <- cor(all.cvs)

print(xtable(cor.cvs, type = "latex"), include.rownames=TRUE, 
      file = "saved/corCvsTable.txt")