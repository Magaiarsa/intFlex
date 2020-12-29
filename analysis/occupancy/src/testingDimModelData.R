
testOccData <- function(){
  ## ************************************************************
  ## check all the dimensions of the data and the names line up 
  ## for all 5 variability measures
  ## ************************************************************
  print(paste("all site names across variability data matches",
  all(all(dimnames(model.input$data$var.partner)[[1]] ==
        dimnames(model.input$data$mean.role)[[1]])==TRUE) &
  (all(dimnames(model.input$data$mean.role)[[1]] ==
        dimnames(model.input$data$var.cnodf)[[1]])) &
  (all(dimnames(model.input$data$var.cnodf)[[1]] ==
        dimnames(model.input$data$mean.cnodf)[[1]]))))

  print(paste("all species names across variability data matches",
              all(all(dimnames(model.input$data$var.partner)[[2]] ==
                        dimnames(model.input$data$mean.role)[[2]])==TRUE) &
                (all(dimnames(model.input$data$mean.role)[[2]] ==
                       dimnames(model.input$data$var.cnodf)[[2]])) &
                (all(dimnames(model.input$data$var.cnodf)[[2]] ==
                       dimnames(model.input$data$mean.cnodf)[[2]]))))

  print(paste("data sites match X matrix",
              all(dimnames(model.input$data$X)$site ==
                    dimnames(model.input$data$var.partner)[[1]])))
  
  print(paste("data species match X matrix",
              all(dimnames(model.input$data$X)$species ==
                    dimnames(model.input$data$var.partner)[[2]])))
  
  

  ## ************************************************************
  ## check x matrix against raw data
  ## ************************************************************
  
  raw.check1 <- unique(spec$GenusSpecies[spec$Site == dimnames(model.input$data$X)$site[1]&
                                           spec$Year == dimnames(model.input$data$X)$year[1] &
                                           spec$Date == min(sr.sched$Date[sr.sched$Site == dimnames(model.input$data$X)$site[1]])])

  model.check1 <- names(model.input$data$X[1,1,1,][ model.input$data$X[1,1,1,] == 1])
  
  print(paste("random raw, model data match 1", all(raw.check1 %in%
                                                      model.check1)))
  
  raw.check2 <- unique(spec$GenusSpecies[spec$Site == dimnames(model.input$data$X)$site[5]&
                                           spec$Year == dimnames(model.input$data$X)$year[2] &
                                           spec$Date == min(sr.sched$Date[sr.sched$Site == dimnames(model.input$data$X)$site[5]])])
  
  model.check2 <- names(model.input$data$X[5,2,1,][ model.input$data$X[5,2,1,] == 1])
  
  print(paste("random raw, model data match 2", all(raw.check2 %in%
                                                      model.check2)))
  
  raw.check3 <- unique(spec$GenusSpecies[spec$Site == dimnames(model.input$data$X)$site[15]&
                                           spec$Year == dimnames(model.input$data$X)$year[5] &
                                           spec$Date == min(sr.sched$Date[sr.sched$Site == dimnames(model.input$data$X)$site[15]])])
  
  model.check3 <- names(model.input$data$X[15,5,1,][ model.input$data$X[15,5,1,] == 1])
  
  print(paste("random raw, model data match 3", all(raw.check3 %in%
                                                      model.check3)))
}
