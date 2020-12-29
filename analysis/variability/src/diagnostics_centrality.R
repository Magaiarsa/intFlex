resid.plot <- function(){
    layout(matrix(1:4))
        plot(fitted(pol.mod.pca), residuals(pol.mod.pca),
         xlab = "Fitted Values", ylab = "Residuals", main="abund + degree")
    abline(h=0, lty=2)
    lines(smooth.spline(fitted(pol.mod.pca),
                        residuals(pol.mod.pca)))

    plot(fitted(abund.pol.mod.pca), residuals(abund.pol.mod.pca),
         xlab = "Fitted Values", ylab = "Residuals", main="just abund")
    abline(h=0, lty=2)
    lines(smooth.spline(fitted(abund.pol.mod.pca),
                        residuals(abund.pol.mod.pca)))

    plot(fitted(degree.pol.mod.pca), residuals(degree.pol.mod.pca),
         xlab = "Fitted Values", ylab = "Residuals", main="just degree")
    abline(h=0, lty=2)
    lines(smooth.spline(fitted(degree.pol.mod.pca),
                        residuals(degree.pol.mod.pca)))

        plot(fitted(plant.mod.pca), residuals(plant.mod.pca),
         xlab = "Fitted Values", ylab = "Residuals",
         main="abund + degree plants")
    abline(h=0, lty=2)
    lines(smooth.spline(fitted(plant.mod.pca),
                        residuals(plant.mod.pca)))

}

pol.pca.scores$pca.var <-
    pol.pca.scores$pca.var[!is.na(pol.pca.scores$pca.var[,"var.pca1"]),]
plant.pca.scores$pca.var <-
    plant.pca.scores$pca.var[!is.na(plant.pca.scores$pca.var[,"var.pca1"]),]


box.site <- function(){
    layout(matrix(1:4, nrow=2))

    boxplot(residuals(pol.mod.pca) ~ Site,
            data = pol.pca.scores$pca.var, main = "Site",
            ylab = "Residuals")
    abline(h=0, lty=2)

    boxplot(residuals(pol.mod.pca) ~ GenusSpecies,
            data = pol.pca.scores$pca.var, main = "Species",
            ylab = "Residuals")
    abline(h=0, lty=2)

    boxplot(residuals(plant.mod.pca) ~ Site,
            data = plant.pca.scores$pca.var, main = "Site",
            ylab = "Residuals")
    abline(h=0, lty=2)
    boxplot(residuals(plant.mod.pca) ~ GenusSpecies,
            data = plant.pca.scores$pca.var, main = "Species",
            ylab = "Residuals")
    abline(h=0, lty=2)
}

