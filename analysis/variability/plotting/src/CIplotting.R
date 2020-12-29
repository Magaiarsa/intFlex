plot.panel <- function(dats,
                       new.dd,
                       y1,
                       xs,
                       treatments,
                       col.lines,
                       col.fill,
                       col.points=col.fill,
                       ylabel="",
                       agg.col="GenusSpecies",
                       plot.x=TRUE,
                       scaled=FALSE,
                       plot.y=TRUE,
                       pchs=c(16),
                       dec=1){
    plotting.loop <- function(){
        ## subset data to specific treatment
        if(!is.null(agg.col)){
            ys <- aggregate(list(y=dats[,y1]),
                            list(sp=dats[, agg.col],
                                 x=dats[,xs]),
                            mean, na.rm=TRUE)
        }else{
            ys <- data.frame(y=dats[,y1],
                             sp=dats[, agg.col],
                             x=dats[,xs])
        }
        ## plots means
        points(x=jitter(ys$x, factor=0.25),
               y=ys$y,
               col=col.points,
               cex=1,
               pch=pchs)
        ## plots CI
        lines(x=new.dd[,xs],
              y=new.dd[,y1],
              col=col.lines,
              lwd=2)
        lines(x=new.dd[,xs],
              y=new.dd$plo,
              col=col.lines,
              lty="dashed")
        lines(x=new.dd[,xs],
              y=new.dd$phi,
              col=col.lines,
              lty="dashed")
        ## add fill from ci.up to ci.lb
        polygon(c(new.dd[,xs],
                  rev(new.dd[,xs])),
                c(new.dd$phi,
                  rev(new.dd$plo)),
                col=col.fill, border=NA)

    }
    plot(NA,
         xlim=range(new.dd[,xs], na.rm=TRUE),
         ylim=range(c(new.dd$phi,  new.dd$plo), na.rm=TRUE),
                      ## dats[,y1]), na.rm=TRUE),
         xlab="",
         ylab="",
         xaxt="n",
         yaxt="n",
         las=1)
    if(plot.y){
        axis(2,
             pretty(range(c(new.dd$phi,  new.dd$plo),
                            ## dats[,y1]),
                          na.rm=TRUE),
                    min.n=3),
             las=1)
        mtext(ylabel, 2, line=3, cex=1.5)
    }
    if(plot.x){
           axis(1, pretty(new.dd[,xs], 5))
    }
    plotting.loop()
}


plot.predict.ypr <- function(new.dd,
                             ylabel,
                             dats,
                             y1,
                             #y2=NA,
                             xs='ypr',
                             xlabel,
                             path = 'figures',
                             extinction.method,
                             agg.col="Site"){
    plot.ci <- function(){
        col.lines <-  brewer.pal(4, "Greys")[3]
        col.fill <- add.alpha(col.lines, alpha=0.2)
        layout(matrix(1, ncol=1))
        par(oma=c(6, 6, 2, 1),
            mar=c(0.5, 0, 0.5, 1))
        plot.panel(dats, new.dd, y1=y1, y2=y2,  xs=xs,
                   col.lines=col.lines, col.fill=col.fill,
                   agg.col=agg.col)
        axis(1, pretty(dats[,xs]), labels=pretty(dats[,xs]))
        mtext(ylabel, 2, line=3, cex=1.5)
        mtext(xlabel, 1, line=3.5, cex=1.5)
    }
    pdf.f(plot.ci, file=file.path(path,
                                  sprintf('%s.pdf', paste(
                                                        gsub('[[:space:]]', '_', xlabel),
                                                        gsub('[[:space:]]', '_', ylabel),
                                                        extinction.method,  sep='_'))),
          width=4, height=4)

}


makePlots <- function(pp, xvar, ys, dd, mods, ylabs, all.specs, xs,
    xlab, agg.col="GenusSpecies"){
    for(j in pp){
        for(i in 1:length(ys)){
            dd1 <- cbind(dd, 0)
            colnames(dd1) <- c(colnames(dd), ys[i])
            dd.pi <- predict.int(mod= mods[[j]][[i]],
                                 dd=dd1,
                                 y=ys[i],
                                 family="gaussian")
            plot.predict.ypr(new.dd=dd.pi,
                             ylabel=ylabs[i],
                             dats=all.specs[[j]],
                             y1=ys[i],
                             extinction.method=j,
                             agg.col=agg.col,
                             xs=xs,
                             xlabel= xlab)
        }
    }
}
