
prepData4Plotting <- function(ms.ms.nimble, all.model.input){
    ## function to pull out coefficent estimates for plotting
    if(is.list(ms.ms.nimble$samples)){
        samples <- do.call(rbind, ms.ms.nimble$samples)
    } else {
        samples <- ms.ms.nimble$samples
    }

    means <- apply(samples, 2, mean, na.rm=TRUE)
    sds <- apply(samples, 2, sd,  na.rm=TRUE)
    all.vars <- rbind(means, sds)

    all.model.input <- list()
    all.all.vars <- list()

    all.model.input[[1]] <- unlist(model.input$data[grep("var.partner", names(model.input$data))])
    all.model.input[[2]] <- unlist(model.input$data[grep("var.role", names(model.input$data))])
    all.model.input[[3]] <- unlist(model.input$data[grep("var.cnodf", names(model.input$data))])
    
    var.types <- c("var.partner","var.role","var.cnodf")

    all.all.vars[[1]] <- all.vars[,grep("var.partner", dimnames(all.vars)[[2]])]
    all.all.vars[[2]] <- all.vars[,grep("var.role", dimnames(all.vars)[[2]])]
    all.all.vars[[3]] <- all.vars[,grep("var.cnodf", dimnames(all.vars)[[2]])]
    all.all.vars[[4]] <- cbind(all.vars[,grep("mu.gam", dimnames(all.vars)[[2]])],
                               all.vars[,grep("gam.var", dimnames(all.vars)[[2]])],
                               all.vars[,grep("mu.phi", dimnames(all.vars)[[2]])],
                               all.vars[,grep("phi.var", dimnames(all.vars)[[2]])],
                               all.vars[,grep("gam.mean", dimnames(all.vars)[[2]])],
                               all.vars[,grep("phi.mean", dimnames(all.vars)[[2]])])

    names(all.all.vars) <- c(var.types, "mu.gam.phi")
    names(all.model.input) <- var.types
    return(list(all.all.vars=all.all.vars,
                all.model.input=all.model.input))

}
