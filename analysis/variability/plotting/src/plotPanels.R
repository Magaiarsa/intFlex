plot.panels <- function(){
  f <- function(){
    col.white <- add.alpha("white", alpha=0)
    col.lines <- brewer.pal(4, "Greys")[3]
    col.fill <- add.alpha(col.lines, alpha=0.5)
    cols.points <- add.alpha(brewer.pal(5, "Set1"), alpha=0.6)
    layout(matrix(1:9, ncol=3, byrow = TRUE))
    par(oma=c(6, 9, 2, 1),
        mar=c(1, 1, 0.5, 1), cex.axis=1.5)
    ## ************************************************************
    ## contribution to nodf
    ## ************************************************************
    ys <- c("contr.nodf")
    agg.col <- "GenusSpecies"
    y2 <- aggregate(list(y=traits.contr.nodf[,ys]),
                    list(sp=traits.contr.nodf[, agg.col]),
                    mean, na.rm=TRUE)
    ## phenology
    plot.panel(new.dd=nodf.days.pi,
               dats=traits.contr.nodf,
               y1="contr.nodf",
               xs="median.days",
               col.fill=col.fill,
               col.lines=col.lines,
               col.points= cols.points,
               agg.col = "GenusSpecies",
               plot.x=FALSE,
               plot.y=TRUE)#,
    mtext("Structure contribution's \n variability", 2, line=5, cex=1.5)
    
    ## degree
    plot.panel(new.dd=nodf.rdegree.pi,
               dats=traits.contr.nodf,
               y1="contr.nodf",
               xs="r.degree",
               col.fill=col.fill,
               col.lines=col.lines,
               col.points= cols.points,
               agg.col = "GenusSpecies",
               plot.x=FALSE,
               plot.y=FALSE)#,
    
    ## abundance
    plot.panel(new.dd=nodf.abund.pi,
               dats=traits.contr.nodf,
               y1="contr.nodf",
               xs="median.abund",
               col.fill=col.white,
               col.lines=col.white,
               col.points= cols.points,
               agg.col = "GenusSpecies",
               plot.x=FALSE,
               plot.y=FALSE)#,
    
    ## ************************************************************
    ## role variability
    ## ************************************************************
    ys <- c("var.pca1")
    agg.col <- "GenusSpecies"
    y2 <- aggregate(list(y=pol.pca.scores$pca.var[,ys]),
                    list(sp=pol.pca.scores$pca.var[, agg.col]),
                    mean, na.rm=TRUE)
    ## phenology 
    plot.panel(new.dd=role.days.pi,
               dats=pol.pca.scores$pca.var,
               y1="var.pca1",
               xs="median.days",
               col.fill=col.fill,
               col.lines=col.lines,
               col.points= cols.points,
               agg.col = "GenusSpecies",
               plot.x=FALSE,
               plot.y=TRUE)#,
    mtext("Role \n variability", 2, line=5, cex=1.5)
    
     ## degree
    plot.panel(new.dd=role.rdegree.pi,
               dats=pol.pca.scores$pca.var,
               y1="var.pca1",
               xs="r.degree",
               col.fill=col.fill,
               col.lines=col.lines,
               col.points= cols.points,
               agg.col = "GenusSpecies",
               plot.x=FALSE,
               plot.y=FALSE)#,
    
    ## abundance
    plot.panel(new.dd=role.abund.pi,
               dats=pol.pca.scores$pca.var,
               y1="var.pca1",
               xs="median.abund",
               col.fill=col.white,
               col.lines=col.white,
               col.points= cols.points,
               agg.col = "GenusSpecies",
               plot.x=FALSE,
               plot.y=FALSE)#,
    
    ## ************************************************************
    ## partner variability
    ## ************************************************************
    ## phenology
    plot.panel(new.dd=partner.days.pi,
               dats=dats.type,
               y1="dist",
               xs="median.days",
               col.fill=col.white,
               col.lines=col.white,
               col.points= cols.points,
               agg.col = "species",
               plot.x=TRUE,
               plot.y=TRUE)#,
    mtext("Phenological breadth", 1, line=5, cex=1.5)
    mtext("Partner \n variability", 2, line=5, cex=1.5)
    
    ## degree 
    plot.panel(new.dd=partner.rdegree.pi,
               dats=dats.type,
               y1="dist",
               xs="r.degree",
               col.fill=col.fill,
               col.lines=col.lines,
               col.points= cols.points,
               agg.col = "species",
               plot.x=TRUE,
               plot.y=FALSE)
    mtext("Degree", 1, line=5, cex=1.5)
    
    ## abundance 
    plot.panel(new.dd=partner.abund.pi,
               dats=dats.type,
               y1="dist",
               xs="median.abund",
               col.fill=col.fill,
               col.lines=col.lines,
               col.points= cols.points,
               agg.col = "species",
               plot.x=TRUE,
               plot.y=FALSE)
    mtext("Abundance", 1, line=5, cex=1.5)
    
  
   }
  path <- 'figures'
  pdf.f(f, file=file.path(path,
                          sprintf("%s.pdf", "baci")),
        width=12, height=10)
  
}
