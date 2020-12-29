library(RColorBrewer)
library(coda)

pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

runMCMCcheckChains <- function(all.chains.samps, f.path=save.dir,
                               params=NULL, num.samps=1000,
                               type){
    ## browser()
    ## function to plot values of parameters as a function of
    ## iterations. Input is the output of runMCMC
    if(class(all.chains.samps) != "list"){
        all.chains.samps <- list(all.chains.samps)
    }
    niter <- length(all.chains.samps[[1]])
    if(is.null(params)){
        params <- colnames(all.chains.samps[[1]])
    }

    f <- function(){
        cols <- c("dodgerblue", "darkolivegreen","goldenrod")
        layout(matrix(1:4, ncol=2))
        for(j in params){
           all.pars <- do.call(rbind, all.chains.samps)
            plot(NA, xlim=c(0,num.samps),
                ylim=range(all.pars[,j]),
                xlab = 'iteration', ylab="", main=j)
            for(i in 1:length(all.chains.samps)){
                this.chain <- all.chains.samps[[i]][,j]
                points(this.chain[seq(from=1, to=length(this.chain),
                                      length.out=num.samps)],
                       type="l",
                       col=cols[i]
                       )


            } # end i loop
        } # end j loop
    }

    pdf.f(f,
          file= file.path(f.path,
                          sprintf("chains_%s.pdf", type)),
          height=11, width=8.5)
}

