resid.plot <- function(){
    layout(matrix(1:3))
    plot(fitted(type.mod.beta), residuals(type.mod.beta),
         xlab = "Fitted Values", ylab = "Residuals", main="abund + degree")
    abline(h=0, lty=2)
    lines(smooth.spline(fitted(type.mod.beta),
                        residuals(type.mod.beta)))

    plot(fitted(abund.type.mod.beta), residuals(abund.type.mod.beta),
         xlab = "Fitted Values", ylab = "Residuals", main="just abund")
    abline(h=0, lty=2)
    lines(smooth.spline(fitted(abund.type.mod.beta),
                        residuals(abund.type.mod.beta)))

    plot(fitted(degree.type.mod.beta), residuals(degree.type.mod.beta),
         xlab = "Fitted Values", ylab = "Residuals", main="just degree")
    abline(h=0, lty=2)
    lines(smooth.spline(fitted(degree.type.mod.beta),
                        residuals(degree.type.mod.beta)))

}


box.site <- function(){
    layout(matrix(1:3, nrow=3))

    boxplot(residuals(type.mod.beta) ~ site,
            data = dats.type, main = "Site",
            ylab = "Residuals")
    abline(h=0, lty=2)

    boxplot(residuals(type.mod.beta) ~ species,
            data = dats.type, main = "Species",
            ylab = "Residuals")
    abline(h=0, lty=2)

    boxplot(residuals(type.mod.beta) ~ year,
            data = dats.type, main = "year",
            ylab = "Residuals")
    abline(h=0, lty=2)
}

