## setwd('~/Dropbox/speciesRoles')
rm(list = ls())
setwd('analysis/occupancy')

source('src/initialize.R')
type <- c("all")

scale <- 50
burnin <- 1e1 * scale
niter <- (1e3) * scale
nthin <- 2
nchain <- 3

model.input <- prepOccModelInput(nzero=0,
                                 threshold=2,
                                 spec=spec,
                                 sr.sched=sr.sched,
                                 veg=by.site,
                                 metrics=var.model,
                                 traits=traits.occupancy)

## check names of dimensions of model and data
testOccData()

print(type)
ms.ms.occ <- setModel(type)
model.input <- setup(type, model.input)
ms.ms.model <- nimbleModel(
    code = ms.ms.occ,
    constants = model.input$constants,
    data = model.input$data,
    inits = model.input$inits,
    check = FALSE,
    calculate = FALSE
)
C.model <- compileNimble(ms.ms.model)
## configure and build mcmc
mcmc.spec <- configureMCMC(ms.ms.model,
                           print = FALSE,
                           monitors =
                               C.model$getNodeNames(stochOnly = TRUE,
                                                    includeData = FALSE),
                           enableWAIC = TRUE)
mcmc <- buildMCMC(mcmc.spec,
                  enableWAIC = TRUE)
C.mcmc <- compileNimble(mcmc, project = ms.ms.model)
## run model
ms.ms.nimble <- runMCMC(C.mcmc,
                        niter = niter,
                        nchains = nchain,
                        nburnin = burnin,
                        WAIC=TRUE)

save(ms.ms.nimble, model.input,
     file = file.path(save.dir,
                      sprintf("ms-ms-samples-%s.Rdata",
                              type)))
runMCMCcheckChains(ms.ms.nimble$samples, type = type)
all.out <- prepData4Plotting(ms.ms.nimble, model.input)
xlabs.cv <- c("Partner variability", "Role variability",
              "Structural variability")
pdf.f(f.plot.all.var,
      file=sprintf("figures/phiGamPsiMeanVar%s.pdf", type),
      height=5, width=13.5)
