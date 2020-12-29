##This function creates a list with the variables that will go into 
##the occupancy model
prepOccupancyData <- function(){
  #browser()
  var.model <- list()
  ## *******************************************************************
  ## Creating the matrix with the values that will be used in the
  ## occupancy model
  
  partner.var <- base::merge(partner.var, species.to.analyze, all.x=FALSE)
  partner.var <- partner.var[partner.var$nyears >= 3,]
  
  pca.var <- merge(pca.var, species.to.analyze, all.x=FALSE)
  pca.var <- pca.var[pca.var$nyears >= 3,]
  
  cnodf.var <- merge(cnodf.var, species.to.analyze, all.x=FALSE)
  cnodf.var <- cnodf.var[cnodf.var$nyears >= 3,]
  
  ##combination of spp and site 
  min.set <- unique(paste(pca.var$Site, pca.var$GenusSpecies))
  
  partner.var <- partner.var[paste(partner.var$Site,
                                   partner.var$GenusSpecies) %in%
                               min.set,]
  cnodf.var <- cnodf.var[paste(cnodf.var$Site,
                               cnodf.var$GenusSpecies) %in% min.set,]
  
  ## *******************************************************************
  ## organizing data for the model and standardizing
  ## *******************************************************************
  #### partner ####
  partner.model <- partner.var[,c("GenusSpecies", "Site", "cv.partner")]
  partner.model$cv.partner <- standardize(partner.model$cv.partner)
  ## creating a species x site matrix
  var.model$var.partner <- acast(partner.model,
                                 Site~GenusSpecies, value.var="cv.partner")
  
  ## *******************************************************************
  ## role cv
  role.cv.model <- pca.var[,c("GenusSpecies", "Site", "pca.cv")]
  role.cv.model$pca.cv <- standardize(role.cv.model$pca.cv)
  var.model$var.role <- acast(role.cv.model,
                              Site~GenusSpecies, value.var="pca.cv")
  
  ## *******************************************************************
  ## cnodf cv
  cnodf.cv.model <- cnodf.var[,c("GenusSpecies", "Site", "cnodf.cv")]
  cnodf.cv.model$cnodf.cv <- standardize(cnodf.cv.model$cnodf.cv)
  var.model$var.cnodf <- acast(cnodf.cv.model,
                               Site~GenusSpecies, value.var="cnodf.cv")
  
  ## *******************************************************************
  traits.occupancy <- traits[,c("GenusSpecies",
                                "r.degree",
                                "median.abund",
                                "median.days",
                                "Lecty")]
  
  # standardizing everyone
  traits.occupancy <- traits.occupancy %>%
    mutate_if(is.numeric, standardize)
  
  traits.occupancy <- traits.occupancy[traits.occupancy$GenusSpecies%in%unique(cnodf.var$GenusSpecies),]
  
  save(var.model, traits.occupancy,
       file="../../../speciesRoles_saved/occupancy/varModelData.RData")  
  return(var.model)
}
