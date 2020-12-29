## setwd('~/Dropbox/speciesRoles')
rm(list=ls())
setwd('analysis/occupancy')
library(xtable)

save.dir <- "../../../speciesRoles_saved/occupancy"

### getting the data
load(file = file.path(save.dir, "ms-ms-samples-all.Rdata"))
samples <- do.call(rbind, ms.ms.nimble$samples)
means <- apply(samples, 2, mean, na.rm = TRUE)
sds <- apply(samples, 2, sd,  na.rm = TRUE)
all.vars <- rbind(means, sds)
quants <- apply(samples, 2, function(x) quantile (x, probs = c(0.025, 0.975))) 
samples.vars <- c("var.partner", "var.role", "var.cnodf")

## occupancy is a result of gam and phi {gam/(1-phi _ gam)}
bySample <- function(samp, vars){
  out <- (samp[grep(paste0("gam.", vars), names(samp))]) /
    (1 - samp[grep(paste0("phi.", vars), names(samp))] +
       samp[grep(paste0("gam.", vars), names(samp))])
  return(out)
}

## calculating occupancy
calcOcc <- function(vars, samples){
  x <- apply(samples, 1, bySample, vars)
  return(x)
}

#takes forever to run
samples.occ <- sapply(samples.vars, calcOcc, samples)
colnames(samples.occ) <- paste0("occ.", colnames(samples.occ))

## getting the occurence CI
quants.occ <- apply(samples.occ, 2, function(x) quantile (x, probs = c(0.025, 0.975))) 
quants.occ <- round(quants.occ, 3)

samples <- cbind(samples, samples.occ)
## Calculating the posteriors
## how may are greater than zero
h1s <- round(apply(samples, 2, function(x)
  sum(x > 0) / length(x)),3)

post.table.h1 <- data.frame(
  Partner = h1s[grep("partner", names(h1s))],
  RoleVar = h1s[grep("var.role", names(h1s))],
  CnodfVar = h1s[grep("var.cnodf", names(h1s))])

post.table.h2 <- 1-post.table.h1

## calcMeanSd uses the all.vars object to grep 
## the mean and sd of each variable of interest
calcMeanSd <- function(x){
  meanSd <- paste0(round(all.vars[,grep(x, colnames(all.vars))][1],3), " (", round(all.vars[,grep(x, colnames(all.vars))][2],3), ")")
  return(meanSd)
}

phi.mean.sd <- data.frame(
  Partner = calcMeanSd("phi.var.partner"),
  RoleVar = calcMeanSd("phi.var.role"),
  CnodfVar = calcMeanSd("phi.var.cnodf"))

gam.mean.sd <- data.frame(
  Partner = calcMeanSd("gam.var.partner"),
  RoleVar = calcMeanSd("gam.var.role"),
  CnodfVar = calcMeanSd("gam.var.cnodf"))

## credible intervals
quants <- round(quants, 3)
phi.quants <- data.frame(
  Partner = paste0("[", quants[colnames(quants) == "phi.var.partner"][1], ", ", 
                   quants[colnames(quants) == "phi.var.partner"][2], "]"),
  RoleVar = paste0("[", quants[colnames(quants) == "phi.var.role"][1], ", ", 
                   quants[colnames(quants) == "phi.var.role"][2], "]"),
  CnodfVar = paste0("[", quants[colnames(quants) == "phi.var.cnodf"][1], ", ", 
                    quants[colnames(quants) == "phi.var.cnodf"][2], "]"))

gam.quants <- data.frame(
  Partner = paste0("[", quants[colnames(quants) == "gam.var.partner"][1], ", ", 
                   quants[colnames(quants) == "gam.var.partner"][2], "]"),
  RoleVar = paste0("[", quants[colnames(quants) == "gam.var.role"][1], ", ", 
                   quants[colnames(quants) == "gam.var.role"][2], "]"),
  CnodfVar = paste0("[", quants[colnames(quants) == "gam.var.cnodf"][1], ", ", 
                    quants[colnames(quants) == "gam.var.cnodf"][2], "]"))

occ.quants <- data.frame(
  Partner = paste0("[", quants.occ[colnames(quants.occ) == "occ.var.partner"][1], ", ", 
                   quants.occ[colnames(quants.occ) == "occ.var.partner"][2], "]"),
  RoleVar = paste0("[", quants.occ[colnames(quants.occ) == "occ.var.role"][1], ", ", 
                   quants.occ[colnames(quants.occ) == "occ.var.role"][2], "]"),
  CnodfVar = paste0("[", quants.occ[colnames(quants.occ) == "occ.var.cnodf"][1], ", ", 
                    quants.occ[colnames(quants.occ) == "occ.var.cnodf"][2], "]"))

post.table.all <- rbind(post.table.h1, post.table.h2,
                        gam.mean.sd, phi.mean.sd, 
                        gam.quants, phi.quants)
rownames(post.table.all) <- c("Gam>1", "Phi>1", "Occ>1",
                              "Gam<1", "Phi<1", "Occ<1",
                              "GamMeanSd", "PhiMeanSd",
                              "GamCI", "PhiCI")

post.table.all <- as.data.frame(t(post.table.all))
col.order <- c("GamMeanSd", "GamCI", "Gam>1", "Gam<1",
               "PhiMeanSd", "PhiCI", "Phi>1", "Phi<1")
post.table.ms <- post.table.all[, col.order]
post.table.occ <- post.table.all[, c("Occ>1", "Occ<1")]
post.table.occ$occCI <- t(occ.quants)
post.table.occ <- post.table.occ[,c("occCI", "Occ>1", "Occ<1")]

print(xtable(post.table.ms, type = "latex"),
      file = file.path(save.dir, "table/postProbTableMs.txt"))

print(xtable(post.table.occ, type = "latex"),
      file = file.path(save.dir, "table/postProbTableOcc.txt"))