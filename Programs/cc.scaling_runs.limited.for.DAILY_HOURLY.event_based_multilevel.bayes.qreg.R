rm(list=ls())
# jags
library(runjags)
library(rjags)
library(R2jags)
library(MCMCpack)
library(parallel)
library(mcmcplots)
library(superdiag)
library(coda)
library(MCMCvis)

setwd("D:/Projects/CA_CCscaling")



#------------------------------------------------------#
#------------------------------------------------------#
# Multilevel Model # for Ta and Td
#------------------------------------------------------#
#------------------------------------------------------#
# for 2 states of ARs: non-AR, ARx formulation

#temporal_selec <- 'hourly'

#lst.precip.ta.td.ar.nevents.all <- readRDS(file="./Data/processed.data/lst.CA.HOURLY.combo.ERA5.850hPa.precip.ta.td.event_based.AR.states.nevents.all.rds")
#lst.precip.ta.td.ar.nevents.all <- readRDS(file="./Data/processed.data/lst.CA.HOURLY.precip.ta.td.event_based.AR.states.nevents.all.rds")

temporal_selec <- 'daily'
#lst.precip.ta.td.ar.nevents.all <- readRDS(file="./Data/processed.data/lst.CA.DAILY.combo.ERA5.850hPa.precip.ta.td.event_based.AR.states.nevents.all.rds")
lst.precip.ta.td.ar.nevents.all <- readRDS(file="./Data/processed.data/lst.CA.DAILY.precip.ta.td.event_based.AR.states.nevents.all.rds")

# removing those sites with small number of events:
lst.precip.ta.td.ar.nevents.all.dummy <- lst.precip.ta.td.ar.nevents.all
min.num.events <- 150 # e.g., at least 100, 150, 200, etc events
nevents <- lst.precip.ta.td.ar.nevents.all$nevents
idx.keep <- which(nevents>min.num.events)
nsites <- length(idx.keep)
lst.precip.ta.td.ar.nevents.all$precip <- lst.precip.ta.td.ar.nevents.all.dummy$precip[idx.keep,]
lst.precip.ta.td.ar.nevents.all$totalprecip <- lst.precip.ta.td.ar.nevents.all.dummy$totalprecip[idx.keep,]
lst.precip.ta.td.ar.nevents.all$ta <- lst.precip.ta.td.ar.nevents.all.dummy$ta[idx.keep,]
lst.precip.ta.td.ar.nevents.all$priorta <- lst.precip.ta.td.ar.nevents.all.dummy$priorta[idx.keep,]
lst.precip.ta.td.ar.nevents.all$td <- lst.precip.ta.td.ar.nevents.all.dummy$td[idx.keep,]
lst.precip.ta.td.ar.nevents.all$priortd <- lst.precip.ta.td.ar.nevents.all.dummy$priortd[idx.keep,]
lst.precip.ta.td.ar.nevents.all$ar <- lst.precip.ta.td.ar.nevents.all.dummy$ar[idx.keep,,]
lst.precip.ta.td.ar.nevents.all$fdates <- lst.precip.ta.td.ar.nevents.all.dummy$fdates[idx.keep,]
lst.precip.ta.td.ar.nevents.all$my.AR.info <- lst.precip.ta.td.ar.nevents.all.dummy$my.AR.info[idx.keep,]
lst.precip.ta.td.ar.nevents.all$nevents <- lst.precip.ta.td.ar.nevents.all.dummy$nevents[idx.keep]
nevents <- lst.precip.ta.td.ar.nevents.all$nevents


bayes_cc_qreg_ar_multilevel <- function() {
  for(i in 1:nsites)
  {
    for(t in 1:nevents[i])
    {
      # Level 1 #
      precip[i,t] ~ dnorm(me[i,t],pe[i,t])
      me[i,t] <- (1-2*p)/(p*(1-p))*w[i,t] + mu[i,t]
      mu[i,t] <- alpha.ars[i,t] + beta.ars[i,t]*ta_td[i,t]
      pe[i,t] <- (tau.q[i]*p*(1-p))/(2*w[i,t])
      w[i,t] ~ dexp(tau.q[i])
      # Level 2 #
      beta.ars[i,t] <- b1[i]*ar1[i,t]+b2[i]*ar2[i,t]
      alpha.ars[i,t] <- a1[i]*ar1[i,t]+a2[i]*ar2[i,t]
    }
    # Prior #
    lsigma[i] ~ dunif(-10,10)
    sigma[i] <- exp(lsigma[i]/2)
    tau.q[i] <- pow(sigma[i],-2)
    a1[i] ~ dnorm(0,0.01)
    a2[i] ~ dnorm(0,0.01)
    b1[i] ~ dnorm(mu.b1,tau.b1)
    b2[i] ~ dnorm(mu.b2,tau.b2)
  }
  # Prior #
  mu.b1 ~ dnorm(0,1)
  mu.b2 ~ dnorm(0,1)
  tau.b1 ~ dgamma(1,1)
  tau.b2 ~ dgamma(1,1)
  sigma.b1 <- 1/sqrt(tau.b1)
  sigma.b2 <- 1/sqrt(tau.b2)
}

# Run/Save for all precipitations (2), percentiles (3), all temperatures (4) #
# 2*3*4 = 24
AR.all.states <- c('non-AR', 'ARx')
num.AR.states <- length(AR.all.states)

ta_td_selec.opts <- c("Ta","Td","PreTa","PreTd")
precip_selec.opts <- c("Precip","TotalPrecip")

p.opts=c(0.5,0.9,0.99)

# pr = 2 # 1, 2
# te = 3 # 1,2,3,4
# 
# 
# per= 3 # 1,2,3

for (pr in 1:length(precip_selec.opts)){
  precip_selec <- precip_selec.opts[pr]
  if (pr==1){
    precip <- log(lst.precip.ta.td.ar.nevents.all$precip)
  } else{
    precip <- log(lst.precip.ta.td.ar.nevents.all$totalprecip)
  }
  
  for (te in 1:length(ta_td_selec.opts)){
    ta_td_selec <- ta_td_selec.opts[te]
    if (te==1){
      ta_td <- lst.precip.ta.td.ar.nevents.all$ta
    }
    if (te==2){
      ta_td <- lst.precip.ta.td.ar.nevents.all$td
    }
    if (te==3){
      ta_td <- lst.precip.ta.td.ar.nevents.all$priorta
    }
    if (te==4){
      ta_td <- lst.precip.ta.td.ar.nevents.all$priortd
    }
    
    for (per in 1:length(p.opts)){
      p <- p.opts[per]
      
      print(paste('Started:::: Percentile:',p,'Temperature:',ta_td_selec,
                  'Prcp:',precip_selec))
      
      
      ar1 <- lst.precip.ta.td.ar.nevents.all$ar[,,1]
      ar2 <- lst.precip.ta.td.ar.nevents.all$ar[,,2]
      
      nevents <- lst.precip.ta.td.ar.nevents.all$nevents
      nsites <- dim(precip)[1]
      
      # (extra--only for saving variables later)
      fdates <- lst.precip.ta.td.ar.nevents.all$fdates
      AR.info <- lst.precip.ta.td.ar.nevents.all$my.AR.info
      
      jags.data <- list('precip' = precip,
                        'ta_td' = ta_td,
                        'p' = p,
                        'ar1' = ar1,
                        'ar2' = ar2,
                        'nevents' = nevents,
                        'nsites' = nsites)
      
      jags.params <- c('b1','b2',
                       'me',
                       'beta.ars','alpha.ars',
                       'sigma.b1','sigma.b2',
                       'a1','a2',
                       'mu.b1','mu.b2')
      
      jags.inits <- NULL
      n.iter <- 8000
      #memory.limit(size=1e+13)
      gc()
      #=============#
      # RUN jags  #
      #=============#
      # Fit the model using run.jags -- parallel #
      start_time <- Sys.time()
      print(paste("started at:",start_time))
      cc.scaling.mod.fit  <- jags.parallel(model.file = bayes_cc_qreg_ar_multilevel, 
                                           data = jags.data, 
                                           parameters.to.save = jags.params,
                                           inits = jags.inits,
                                           n.iter=8000,
                                           n.chains = 4)
      end_time <- Sys.time()
      run.time.jags <- end_time - start_time
      print(paste("ended at:",end_time))
      print(run.time.jags)
      
      gc()
      
      
      attach.jags(cc.scaling.mod.fit)
      
      #plot(cc.scaling.mod.fit)
      
      
      beta.ars.prct <- (exp(beta.ars)-1)*100 # n.iter -by- nsites -by- maximum number of events
      beta.ars.prct.est <- apply(beta.ars.prct,FUN=median,c(2,3)) # nsites -by- maximum number of events
      
      beta.ars.prct.est.final <- array(NA,c(num.AR.states,nsites))
      for (k in 1:num.AR.states){
        for (s in 1:nsites){
          idx.extract <- which(lst.precip.ta.td.ar.nevents.all$ar[s,,k]>0)
          if (length(idx.extract)>0){
            beta.ars.prct.est.final[k,s] <- median(beta.ars.prct.est[s,idx.extract])
          }
        }
      }
      
      alpha.ars.est <- apply(alpha.ars,FUN=median,c(2,3)) # nsites -by- maximum number of events
      
      alpha.ars.est.final <- array(NA,c(num.AR.states,nsites))
      for (k in 1:num.AR.states){
        for (s in 1:nsites){
          idx.extract <- which(lst.precip.ta.td.ar.nevents.all$ar[s,,k]>0)
          if (length(idx.extract)>0){
            alpha.ars.est.final[k,s] <- median(alpha.ars.est[s,idx.extract])
          }
        }
      }
      
      colm <- function(x){return(apply(x,FUN=median,2))}
      
      b.coefs <- rbind(colm(b1),colm(b2))
      b.coefs.prct <- (exp(b.coefs)-1)*100
      a.coefs <- rbind(colm(a1),colm(a2))
      
      sigma.b.coefs <- cbind(sigma.b1,sigma.b2)
      
      sigma.b.coefs.prct <- (exp(sigma.b.coefs)-1)*100
      
      mu.b.coefs <- cbind(mu.b1,mu.b2)
      
      mu.b.coefs.prct <- (exp(mu.b.coefs)-1)*100
      
      # calculating Rhats and convergence checks
      # NEW
      mod_fit_summary = cc.scaling.mod.fit$BUGSoutput$summary
      mod_all_convergence = cc.scaling.mod.fit$BUGSoutput$summary[,8:9]
      
      mod_rhat_summary = c(min(mod_all_convergence[,1]),quantile(mod_all_convergence[,1],c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),max(mod_all_convergence[,1]),mean(mod_all_convergence[,1]))
      mod_neff_summary = c(min(mod_all_convergence[,2]),quantile(mod_all_convergence[,2],c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),max(mod_all_convergence[,2]),mean(mod_all_convergence[,2]))
      
      
      # saving the parameters #
      lst.names.order <- c(tolower(precip_selec.opts[pr]),
                           tolower(ta_td_selec[te]),
                           'p','nevents','nsites','n.iter',
                           'sigma.b.coefs ','sigma.b.coefs.prct',
                           'b.coefs.prct','a.coefs','beta.ars.prct.est',
                           'beta.ars.prct.est.final','alpha.ars.est','alpha.ars.est.final',
                           'mu.b.coefs.prct','mu.b.coefs',
                           'fdates','AR.info',
                           'rhat_summary',
                           'neff_summary',
                           'idx.sites.keep')
      
      lst.cc.scaling.mod.fit <- list(precip,
                                     ta_td,
                                     p,nevents,nsites,n.iter,
                                     sigma.b.coefs, sigma.b.coefs.prct,
                                     b.coefs.prct,a.coefs,beta.ars.prct.est,
                                     beta.ars.prct.est.final,alpha.ars.est,alpha.ars.est.final,
                                     mu.b.coefs.prct,mu.b.coefs,
                                     fdates,AR.info,
                                     mod_rhat_summary,
                                     mod_neff_summary,
                                     idx.keep)
      names(lst.cc.scaling.mod.fit) <- lst.names.order
      
      saveRDS(lst.cc.scaling.mod.fit,
              #file=paste0("./Data/model.mcmc.run.limited/cc.scaling.mod.multi.bayes.ARs.event_based.fit.850Combo.",
                          file=paste0("./Data/model.mcmc.run.limited/cc.scaling.mod.multi.bayes.ARs.event_based.fit.",
                          temporal_selec,".",
                          p,".prct.",ta_td_selec,
                          ".",precip_selec,
                          ".",n.iter,"iter.rds"))
      print(paste('Finished:::: Percentile:',p,'Temperature:',ta_td_selec))
      rm(lst.cc.scaling.mod.fit,p,cc.scaling.mod.fit,
         beta.ars.prct,alpha.ars.est,mu.b.coefs.prct,mu.b.coefs)
      gc()
      #
    }
  }
}
