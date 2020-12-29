setModel <- function(type){
  ##*************************************
  ## model all
  ##*************************************
  if(type == "all"){
    print("all model")
    ms.ms.occ <- nimbleCode({
      ## multi-species priors
      ## detectablility
      mu.p.0     ~ dnorm(0,0.001)
      mu.p.day.1 ~ dnorm(0,0.001)
      mu.p.day.2 ~ dnorm(0,0.001)
      sigma.p.0     ~ dunif(0,100)
      sigma.p.day.1 ~ dunif(0,100)
      sigma.p.day.2 ~ dunif(0,100)
      
      ## phi/gam random intercepts
      mu.phi.0  ~ dnorm(0,0.001)
      mu.gam.0  ~ dnorm(0,0.001)
      sigma.phi.0 ~ dunif(0,100)
      sigma.gam.0 ~ dunif(0,100)
      mu.phi.fra  ~ dnorm(0,0.001)
      mu.gam.fra  ~ dnorm(0,0.001)
      sigma.phi.fra ~ dunif(0,100)
      sigma.gam.fra ~ dunif(0,100)
      
      phi.var.partner  ~ dnorm(0,0.001)
      gam.var.partner  ~ dnorm(0,0.001)
      
      phi.var.role  ~ dnorm(0,0.001)
      gam.var.role  ~ dnorm(0,0.001)
      
      phi.var.cnodf  ~ dnorm(0,0.001)
      gam.var.cnodf  ~ dnorm(0,0.001)
      
      ## species-specific  parameters
      for(sp in 1:nsp) {
        ## day
        p.0[sp]     ~ dnorm(mu.p.0,     sd=sigma.p.0)
        p.day.1[sp] ~ dnorm(mu.p.day.1, sd=sigma.p.day.1)
        p.day.2[sp] ~ dnorm(mu.p.day.2, sd=sigma.p.day.2)
        
        ## species specific intercept
        phi.0[sp] ~ dnorm(mu.phi.0, sd=sigma.phi.0)
        gam.0[sp] ~ dnorm(mu.gam.0, sd=sigma.gam.0)
        phi.fra[sp] ~ dnorm(mu.phi.fra, sd=sigma.phi.fra)
        gam.fra[sp] ~ dnorm(mu.gam.fra, sd=sigma.gam.fra)
        
      }
      for(sp in 1:nsp) {
        for(site in 1:nsite) {
          for(yr in 1:nyear) {
            for(rep in 1:nrep[site,yr,sp]) {
              logit(p[site,yr,rep,sp]) <-
                p.0[sp] +
                p.day.1[sp]*day[site,yr,rep,sp] +
                p.day.2[sp]*day.2[site,yr,rep,sp]
            }
          }
          ## start off at the average for each species, site
          ## across years
          logit(phi.site.sp.mean[site,sp]) <-
            mean(phi[site, 1:(nyear-1),sp])
          logit(gam.site.sp.mean[site,sp]) <-
            mean(gam[site,1:(nyear-1),sp])
          
          psi.1[site,sp] <- gam.site.sp.mean[site,sp]/
            (1 - phi.site.sp.mean[site,sp] +
               gam.site.sp.mean[site,sp])
          
          ## occupancy in year 1
          ## psi is on a logit scale
          psi[site,1,sp] <- psi.1[site,sp]
          
          ## occupancy in subsequent years
          ## phi and gam are on a linear scale
          for(yr in 1:(nyear-1)) {
            phi[site,yr,sp] <-
              phi.0[sp] +
              phi.fra[sp]*fra[site,yr] +
              phi.var.partner*var.partner[site,sp] +
              phi.var.role*var.role[site,sp] +
              phi.var.cnodf*var.cnodf[site,sp]
            
            gam[site,yr,sp] <-
              gam.0[sp] +
              gam.fra[sp]*fra[site,yr] +
              gam.var.partner*var.partner[site,sp] +
              gam.var.role*var.role[site,sp] +
              gam.var.cnodf*var.cnodf[site,sp] 
          }
        }
      }
      
      
      for(site in 1:nsite) {
        for(sp in 1:nsp) {
          X[site, 1:nyear, 1:max.nreps, sp] ~
            dDynamicOccupancy(nrep=nrep[site, 1:nyear, sp],
                              psi1=psi[site,1,sp],
                              phi=expit(phi[site,1:(nyear-1),sp]),
                              gamma=expit(gam[site,1:(nyear-1),sp]),
                              p=p[site, 1:nyear, 1:max.nreps, sp])
          
        }
      }
      
      
    })
  }
  return(ms.ms.occ)
}