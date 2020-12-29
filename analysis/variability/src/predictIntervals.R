
predict.int <- function(mod, xlabel, dats, method){
  pdf.f <- function(f, file, ...) {
    cat(sprintf("Writing %s\n", file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
  }
  dd <- expand.grid(traits=seq(from= min(dats, na.rm=TRUE),
                      to= max(dats, na.rm=TRUE), length=10),
                    status=c('control', 'mature', 'maturing'),
                    cv=0)
  mm <- model.matrix(terms(mod), dd)
  dd$cv <- mm %*% fixef(mod) 
  pvar1 <- diag(mm %*% vcov(mod) %*% t(mm))
  tvar1 <- pvar1 + VarCorr(mod)$site[1] + VarCorr(mod)$sp[1] 
  new.dd <- data.frame(dd,
                       plo = dd$cv-2*sqrt(pvar1),
                       phi = dd$cv+2*sqrt(pvar1),
                       tlo = dd$cv-2*sqrt(tvar1),
                       thi = dd$cv+2*sqrt(tvar1))

  plot.ci <- function(trait, new.dd){
    cols <- hsv(c(0.75, 0.35, 0.1), alpha=0.2)
    cols.lines <- hsv(c(0.75, 0.35, 0.1), alpha=1)
    plot(NA, xlim=range(dats, na.rm=TRUE), ylim=c(0,510),
         xlab=trait, ylab='Coefficient of variation')
    legend("bottomright", legend= c("Mature", "Maturing",
                            "Unrestored"), col=cols.lines,
           pch=16, bty="n", cex=1.5)
    statuses <- c("mature", "maturing", "control")
    for(i in 1:3){
      lines(x=new.dd$traits[new.dd$status == statuses[i]],
            y=new.dd$cv[new.dd$status == statuses[i]], col=
            cols.lines[i])
      lines(x=new.dd$traits[new.dd$status == statuses[i]],
            y=new.dd$plo[new.dd$status == statuses[i]], lty=2,
            col= cols.lines[i])
      lines(x=new.dd$traits[new.dd$status == statuses[i]],
            y=new.dd$phi[new.dd$status == statuses[i]], lty=2,
            col=cols.lines[i])
      polygon(c(new.dd$traits[new.dd$status == statuses[i]],
                rev(new.dd$traits[new.dd$status == statuses[i]])),
              c(new.dd$phi[new.dd$status == statuses[i]],
                rev(new.dd$plo[new.dd$status == statuses[i]])),
              col = cols[i], border = NA)
    }
  }

  f <- function(){
    par(mfrow=c(1,1), oma=c(2.5, 2.5, 0.25, 1),
        mar=c(5, 5, 2, 0.5), cex.axis=1.3, cex.lab=1.3)
    plot.ci(xlabel, new.dd)
  }

  path <- 'figures' 
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste(xlabel, method, sep="_"))),
        width=6, height=6)

}
