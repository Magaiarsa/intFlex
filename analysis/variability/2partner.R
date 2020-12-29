## setwd('~/Dropbox/speciesRoles')
setwd('analysis/variability')
source('src/initialize.R')

args <- commandArgs(trailingOnly=TRUE)

## this script calculates the partner beta diversity using Marti
## Anderson's miultivariate dispersion method. It uses the null
## communities generated in the nulls.R script to control for alpha
## diversity in the beta diversity calculation

binary <- FALSE ## calculate occurrence-based beta-div metric?
alpha <- TRUE ## use null models that control for alpha diversity?

if(binary){
  dis.method <- "jaccard" ## distance method to use in the distance matrix 
}else{
  dis.method <- "chao" ## distance method to use in the distance matrix 
}

zscore <- TRUE
source('src/initialize_beta.R')

dis <- mapply(function(a, c, d)
    calcBetaStatus(comm= a, ## observed communities
                   dis.method, ## dissimilarity metric
                   nulls=c, ## null communities
                   occ=binary, ## binary or abundance weighted?
                   zscore=zscore), ## use Chase method not zscores
    a=comm$comm,
    c= nulls,
    d= comm$comm,
    SIMPLIFY=FALSE)


dis.beta <- makeBetaDataPretty(comm = comm$comm, dis=dis)

## calculate the sd per site and cv per species per site (across years)
partner.var <- dis.beta  %>%
  group_by(Site, GenusSpecies)  %>%
    summarize(cv.partner = cv(dist),
              sd.partner = sd(dist),
              Mean.partner = mean(dist))

if(binary){
save(partner.var, dis.beta,
     file=file.path(s.path,"partner.RData"))
}else{
  save(partner.var, dis.beta,
       file=file.path(s.path,"partnerQ.RData"))
}
