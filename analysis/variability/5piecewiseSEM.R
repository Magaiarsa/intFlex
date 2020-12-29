## setwd('~/Dropbox/speciesRoles')
rm(list=ls())
setwd('analysis/variability')
s.path <- "saved/pathAn"
s.path.fig <- "figures"

binary <- FALSE ## binary or quantitative
use.mean <- FALSE ## cv or mean
source("src/initialize_path.R")

## *******************************************************************
### Partner variability
## *******************************************************************
## fit piecewise model with random effect of species and site
## indirect effect through rdegree, test the effect of cv
partner.var.cv = psem(
  r.degree = lme(r.degree ~ median.abund + median.days + Lecty + MeanITD,
                 random = ~ 1 | Site/GenusSpecies,
                 data = path.variables$partner),
  cv.dist = lme(cv.partner ~ median.abund + median.days+ r.degree  + Lecty + MeanITD,
                random = ~ 1 | Site/GenusSpecies,
                data = path.variables$partner)) 

# Evaluate path significance
summary(partner.var.cv)

## removing Bombus and re-running = still significant
## path.variables$partner <- path.variables$partner[path.variables$partner$GenusSpecies!="Bombus melanopygus",]

## *******************************************************************
### Role variability
## *******************************************************************
## effect of cv
role.var.cv = psem(
  r.degree = lme(r.degree ~ median.abund + median.days + Lecty + MeanITD,
                 random = ~ 1 | Site/GenusSpecies,
                 data = path.variables$role),
  pca.cv = lme(pca.cv ~ median.abund + median.days+ r.degree  +  Lecty + MeanITD,
               random = ~ 1 | Site/GenusSpecies,
               data = path.variables$role))

summary(role.var.cv)

## *******************************************************************
### Nodf contribution ----
## *******************************************************************
## effect of cv
cnodf.var.cv = psem(
  r.degree = lme(r.degree ~ median.abund + median.days + Lecty + MeanITD,
                 random = ~ 1 | Site/GenusSpecies,
                 data = path.variables$cnodf),
  cnodf.cv = lme(cnodf.cv ~ median.abund + median.days + r.degree + Lecty + MeanITD,
               random = ~ 1 | Site/GenusSpecies,
               data = path.variables$cnodf))

summary(cnodf.var.cv)

## *******************************************************************
## saving the results and creating a table
## *******************************************************************
## Plotting 
if(use.mean){
  save(partner.var.cv, role.var.cv, cnodf.var.cv,
       file = file.path(s.path, "pathResultsMean.Rdata"))
  export_graph(plot(partner.var.cv, return = TRUE), 
               file_name = file.path(s.path.fig,"path_MeanPartner.pdf"),file_type = "pdf")
  export_graph(plot(role.var.cv, return = TRUE), 
               file_name = file.path(s.path.fig,"path_MeanRole.pdf"),file_type = "pdf")
  export_graph(plot(cnodf.var.cv, return = TRUE), 
               file_name = file.path(s.path.fig,"path_MeanStr.pdf"),file_type = "pdf")
}else{
  save(partner.var.cv, role.var.cv, cnodf.var.cv,
       file = file.path(s.path, "pathResultsCV.Rdata"))
  export_graph(plot(partner.var.cv, return = TRUE), 
               file_name = file.path(s.path.fig,"path_CVpartner.pdf"),file_type = "pdf")
  export_graph(plot(role.var.cv, return = TRUE), 
               file_name = file.path(s.path.fig,"path_CVrole.pdf"),file_type = "pdf")
  export_graph(plot(cnodf.var.cv, return = TRUE), 
               file_name = file.path(s.path.fig,"path_CVstr.pdf"),file_type = "pdf")
}

## Creating the table
## partner
path.results.all <- cbind("partner.var.cv",
                          summary(partner.var.cv)$coefficients)
colnames(path.results.all)[1] <- "var"

## role
role.cv <- cbind("role.var.cv", summary(role.var.cv)$coefficients)
colnames(role.cv)[1] <- "var"

# role.mean <- cbind("role.var.mean", summary(role.var.mean)$coefficients)
# colnames(role.mean)[1] <- "var"
##cnodf
cnodf.cv <- cbind("cnodf.var.cv", summary(cnodf.var.cv)$coefficients)
colnames(cnodf.cv)[1] <- "var"

path.results.all <- rbind(path.results.all,
                          role.cv,
                          cnodf.cv)

rownames(path.results.all) <- NULL 

## renaming stuff to make my life easier 
## variables
if(use.mean){
  path.results.all$var[path.results.all$var== "partner.var.cv"] <- "Mean Partner"
  path.results.all$var[path.results.all$var== "role.var.cv"] <- "Mean Role"
  path.results.all$var[path.results.all$var== "cnodf.var.cv"] <- "Mean Struct."
  
}else{
  path.results.all$var[path.results.all$var== "partner.var.cv"] <- "Partner var."
  path.results.all$var[path.results.all$var== "role.var.cv"] <- "Role var."
  path.results.all$var[path.results.all$var== "cnodf.var.cv"] <- "Struct var."
}

## responses 
path.results.all$Response[path.results.all$Response== "r.degree"] <- "diet breadth"
path.results.all$Response[path.results.all$Response== "cv.partner"] <- "Partner var."
path.results.all$Response[path.results.all$Response== "pca.cv"] <- "Role var."
path.results.all$Response[path.results.all$Response== "cnodf.cv"] <- "Structural var."

## responses 
path.results.all$Predictor[path.results.all$Predictor== "median.abund"] <- "abund."
path.results.all$Predictor[path.results.all$Predictor== "median.days"] <- "phenology"
path.results.all$Predictor[path.results.all$Predictor== "Lecty"] <- "speciallization"
path.results.all$Predictor[path.results.all$Predictor== "MeanITD"] <- "body size"
path.results.all$Predictor[path.results.all$Predictor== "r.degree"] <- "diet breadth"

if(use.mean){
  print(xtable(path.results.all, type = "latex"), include.rownames=FALSE, 
        file = "saved/pathResultsTableMean.txt")
}else{
  print(xtable(path.results.all, type = "latex"), include.rownames=FALSE, 
        file = "saved/pathResultsTableCV.txt")
}

