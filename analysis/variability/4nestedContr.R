## setwd('~/Dropbox/speciesRoles')
rm(list=ls())
setwd('analysis/variability')
source('src/initialize.R')
load('../../data/nestedContribution.Rdata')

## This script calculates the  variability of nestedness contribution
## (calculated in the dataPrep.R file using the method from the
## function nestedcontribution (in the bipartite package)

## Adding the minimum to evryone so the cv makes sense
smallest <- abs(min(contr.nodf$nestedcontribution))
contr.nodf$nestedcontribution <- contr.nodf$nestedcontribution + smallest

## calculating mean and cv per species per site
cnodf.var <- contr.nodf %>%
  group_by(GenusSpecies, Site)  %>%
  dplyr::summarize(
    cnodf.cv = cv(nestedcontribution),
    cnodf.sd = sd(nestedcontribution),
    cnodf.Mean = mean(nestedcontribution,
                      na.rm = TRUE))

save(cnodf.var, contr.nodf,
     file=file.path(s.path,"cnodf.RData"))
