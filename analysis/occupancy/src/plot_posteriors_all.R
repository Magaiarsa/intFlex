f.plot.all.var <- function(){
    plotAll(mean.var="var",
            all.model.input=all.out$all.model.input,
            all.all.vars=all.out$all.all.vars,
            xlabs.mean, xlabs.cv)
}


f.plot.all.mean <- function(){
    plotAll(mean.var="mean",
            all.model.input=all.out$all.model.input,
            all.all.vars=all.out$all.all.vars,
            xlabs.mean, xlabs.cv)
}

plotAll <- function(mean.var, all.model.input, all.all.vars,
                    xlabs.mean, xlabs.cv){
  if(mean.var=="mean"){
    layout(matrix(c(1:2), ncol=2))
  }else{
  layout(matrix(c(1:3), ncol=3))}
  par(oma=c(5, 5, 2, 1),
      mar=c(4, 2, 1, 1.5), cex.axis=1.5)
  ## creating the lists of data
  if(mean.var=="var"){
    this.input.var <-
      all.model.input[grep(mean.var, names(all.model.input))]
    this.vars <- all.all.vars[grep(mean.var, names(all.all.vars))]
    xlabs.this <- xlabs.cv
    }
  if(mean.var=="mean"){
        this.input.var <-
          all.model.input[grep(mean.var, names(all.model.input))]
        this.vars <- all.all.vars[grep(mean.var, names(all.all.vars))]
        xlabs.this <- xlabs.mean
      }
  mu.gam.phi.var <- all.all.vars[grep("mu.gam.phi", names(all.all.vars))]
  for(t in 1:length(this.input.var)){
    print(t)
    makeOccPLotAll(t,
                mean.var,
                xlabs=xlabs.this[t],
                this.input=this.input.var[t],
                this.var=this.vars[t],
                mu.gam.phi.var=mu.gam.phi.var)
    }
}

makeOccPLotAll <- function(this.type,
                        mean.var, ## var or mean
                        #mean.var.model.input, ## mean or cv
                        xlabs,
                        this.input,
                        this.var,
                        mu.gam.phi.var){
      plot(NA, ylim=c(0, 1),
         xlim=range(this.input),
         las=1, ylab="", xlab="", yaxt="n")
    mtext(xlabs, 1, line=4, cex=1.3)

    if(this.type==1){
        axis(2, at=pretty(c(0,1)))
        mtext("Population parameter", 2, line=4, cex=1.3)
        if(mean.var=="var"){
            #legend("topleft", inset=0.15,
              legend(x=-1.5, y=1,
                   legend=c("Persistence", "Colonization",
                                   "Occupancy"),
               lty=c("solid", "dashed", "dotted"), bty="n", cex=1.2)
        }else{
            #legend("bottomright",
            legend(x=-1, y=0.9,
                   legend=c("Persistence", "Colonization",
                                            "Occupancy"),
                        lty=c("solid", "dashed", "dotted"), bty="n",
                        cex=0.8)
               }
    }

    curve(expit(mu.gam.phi.var$mu.gam.phi['means','mu.phi.0'] +
                mu.gam.phi.var$mu.gam.phi['means',
                (paste0('phi.', names(this.var)))] * x),
          from=range(this.input)[1],
          to=range(this.input)[2],
          col="black",
          lwd=2,
          add=TRUE)

    curve(expit(mu.gam.phi.var$mu.gam.phi['means','mu.gam.0'] +
                mu.gam.phi.var$mu.gam.phi['means',
                (paste0('gam.', names(this.var)))] * x),
          from=range(this.input)[1],
          to=range(this.input)[2],
          col="black",
          lty="dashed",
          lwd=2,
          add=TRUE)

    curve((expit(mu.gam.phi.var$mu.gam.phi['means','mu.gam.0'] +
                 mu.gam.phi.var$mu.gam.phi['means',
                 (paste0('gam.', names(this.var)))] * x))/
          (1 + expit(mu.gam.phi.var$mu.gam.phi['means','mu.gam.0'] +
                     mu.gam.phi.var$mu.gam.phi['means',
                     (paste0('gam.', names(this.var)))] * x) -
           expit(mu.gam.phi.var$mu.gam.phi['means','mu.phi.0'] +
                 mu.gam.phi.var$mu.gam.phi['means',
                 (paste0('phi.', names(this.var)))] * x)),
          from=range(this.input)[1],
          to=range(this.input)[2],
          col="grey51",
          lty="dotted",
          lwd=2,
          add=TRUE)
}

