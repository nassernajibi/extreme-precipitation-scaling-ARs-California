rm(list=ls())
# jags
library(dplyr)
library(runjags)
library(rjags)
library(R2jags)
library(MCMCpack)
library(parallel)
library(mcmcplots)
library(superdiag)
library(coda)
library(MCMCvis)
#library(jagsUI)
library(superdiag)
library(parallel)
cl <- makeCluster(4)
library(stringr)

setwd("D:/Projects/CA_CCscaling")

#****************************** ##
##/--  Hierarchical Bayesian Trend Scaling ----  \##
#****************************** ##

#------------------------------------------------------#
#------------------------------------------------------#
# Multilevel Model # for Ta and Td
#------------------------------------------------------#
#------------------------------------------------------#
# for 2 states of ARs: non-AR, ARx formulation
# processed P and T matrices, ARs tags, threshold, lambda

#@ surface 2 m
temporal_selec <- 'daily'
lst.precip.ta.td.ar.nevents.all <- readRDS(file='./Data/processed.data/lst.CA.trend.GPD_DAILY.precip.ta.td.event_based.AR.states.nevents.all.rds')

# temporal_selec <- 'hourly'
# lst.precip.ta.td.ar.nevents.all <- readRDS(file='./Data/processed.data/lst.CA.trend.GPD_HOURLY.precip.ta.td.event_based.AR.states.nevents.all.rds')

## @ 850-hPa
# temporal_selec <- 'daily'
# lst.precip.ta.td.ar.nevents.all <- readRDS(file='./Data/processed.data/lst.CA.trend.GPD_Combo.850.DAILY.precip.ta.td.event_based.AR.states.nevents.all.rds')

# temporal_selec <- 'hourly'
# lst.precip.ta.td.ar.nevents.all <- readRDS(file='./Data/processed.data/lst.CA.trend.GPD_Combo.850.HOURLY.precip.ta.td.event_based.AR.states.nevents.all.rds')


lst.precip.ta.td.ar.nevents.all0 <- lst.precip.ta.td.ar.nevents.all


# removing those sites with small number of events:
lst.precip.ta.td.ar.nevents.all.dummy <- lst.precip.ta.td.ar.nevents.all0
min.num.events <- 30 # e.g., at least 30, etc events
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
lst.precip.ta.td.ar.nevents.all$threshold_mm <- lst.precip.ta.td.ar.nevents.all.dummy$threshold_mm[idx.keep]
lst.precip.ta.td.ar.nevents.all$lambda_mean <- lst.precip.ta.td.ar.nevents.all.dummy$lambda_mean[idx.keep]
nevents <- lst.precip.ta.td.ar.nevents.all$nevents


bayes_trend_scaling_ar_multilevel <- "model {
  for(i in 1:nsites)
  {
    for(t in 1:nevents[i])
    {
      precip[i,t] ~ dgenpar(sigma[i,t],threshold_mm[i],xi[i])
      sigma[i,t] <- xi[i]*max(0.00001,(exp(alpha_d[i,t]+beta_d[i,t]*ta_td[i,t])-threshold_mm[i]))/(pow(return.period*lambda_mean[i],xi[i])-1)
      beta_d[i,t] <- db1[i]*ar1[i,t]+db2[i]*ar2[i,t] # slope
      alpha_d[i,t] <- da1[i]*ar1[i,t]+da2[i]*ar2[i,t] # intercept
    }
    # Prior #
    xi[i] ~ dnorm(0,0.01) T(0,) #truncated normal distribution
    #xi[i] ~ dnorm(0,0.01);I(0,) #truncated normal distribution
    da1[i] ~ dnorm(0,0.01)
    da2[i] ~ dnorm(0,0.01)
    db1[i] ~ dnorm(mu.db1,tau.db1)
    db2[i] ~ dnorm(mu.db2,tau.db2)
  }
  # Prior #
  mu.db1 ~ dnorm(0,1)# then, fit competing time-varying models based on PRCP vs. TEMP
  mu.db2 ~ dnorm(0,1)# then, fit competing time-varying models based on PRCP vs. TEMP
  tau.db1 ~ dgamma(1,1)
  tau.db2 ~ dgamma(1,1)
  sigma.db1 <- 1/sqrt(tau.db1)
  sigma.db2 <- 1/sqrt(tau.db2)
}"

# Run/Save for all precipitations (2), all temperatures (4) #
# 2*4 = 8
AR.all.states <- c('non-AR', 'ARx')
num.AR.states <- length(AR.all.states)

ta_td_selec.opts <- c("Ta","Td","PreTa","PreTd")
precip_selec.opts <- c("Precip","TotalPrecip")

# pr = 2 # 1, 2
# te = 2 # 1,2,3,4

return.period <- 20; qq <- 0.05

# (extra--only for saving variables later)
fdates <- lst.precip.ta.td.ar.nevents.all$fdates
AR.info <- lst.precip.ta.td.ar.nevents.all$my.AR.info


for (pr in 1:length(precip_selec.opts)){
  
  precip_selec <- precip_selec.opts[pr]
  if (pr==1){
    precip <- lst.precip.ta.td.ar.nevents.all$precip
  } else{
    precip <- lst.precip.ta.td.ar.nevents.all$totalprecip
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
    
    print(paste('Started:::: Temperature:',ta_td_selec,'Prcp:',precip_selec))
    
    ar1 <- lst.precip.ta.td.ar.nevents.all$ar[,,1]
    ar2 <- lst.precip.ta.td.ar.nevents.all$ar[,,2]
    
    nevents <- lst.precip.ta.td.ar.nevents.all$nevents
    
    threshold_mm <- lst.precip.ta.td.ar.nevents.all$threshold_mm
    lambda_mean <- lst.precip.ta.td.ar.nevents.all$lambda_mean
    
    jags.data <- list('precip' = precip,
                      'ta_td' = ta_td,
                      'ar1' = ar1,
                      'ar2' = ar2,
                      'nevents' = nevents,
                      'nsites' = nsites,
                      'threshold_mm' = threshold_mm,
                      'return.period' = return.period,
                      'lambda_mean' = lambda_mean)
    
    jags.params <- c('mu.db1','mu.db2','sigma.db1','sigma.db2',
                     'xi','sigma',
                     'db1','db2','da1','da2')
    
    jags.inits <- NULL
    n.iter <- 1000
    n.iter2 <- 4*n.iter # 4 chains combined
    # memory.limit(size=1e+13)
    gc()
    
    #Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.1")
    #Sys.setenv(JAGS_HOME="C:/Users/nn289/Documents/R/win-library/3.6/JAGS/JAGS-4.3.1")
    #Sys.setenv(JAGS_ROOT="C:/Program Files/JAGS/JAGS-4.3.1")
    #load.runjagsmodule()
    #load.module("runjags")
    #list.modules()
    #=============#
    # RUN jags  #
    #=============#
    # Fit the model using run.jags -- parallel #
    start_time <- Sys.time()
    print(paste("started at:",start_time))
    ## trend.scaling.mod.fit  <- jags.parallel(model.file = bayes_trend_scaling_ar_multilevel, 
    ##                                         data = jags.data, 
    ##                                         parameters.to.save = jags.params,
    ##                                         inits = jags.inits,
    ##                                         n.iter=50000,
    ##                                         n.chains = 4,
    ##                                         jags.module='runjags')
    
    trend.scaling.mod.fit.run.jags <- run.jags(model=bayes_trend_scaling_ar_multilevel,
                                               data=jags.data,
                                               n.chains=4,
                                               burnin=4000,
                                               sample=n.iter,
                                               monitor=jags.params,
                                               method="rjparallel", cl=cl,
                                               modules="runjags")
    
    
    end_time <- Sys.time()
    
    run.time.jags <- end_time - start_time
    print(paste("ended at:",end_time))
    print(run.time.jags)
    
    gc()
    
    ##trend.scaling.mod.fit <- as.mcmc.list(trend.scaling.mod.fit.run.jags)
    trend.scaling.mod.fit.mcmc <- as.mcmc(trend.scaling.mod.fit.run.jags)
    df.trend.scaling.mod.fit.mcmc <- as.data.frame(trend.scaling.mod.fit.mcmc)
    ##attach.jags(trend.scaling.mod.fit)
    
    ## diagnostics plots
    ## plot(trend.scaling.mod.fit)
    #plot(trend.scaling.mod.fit.mcmc)
    # 
    # trend.scaling.mod.fit.mcmc <- as.mcmc(trend.scaling.mod.fit)
    # denplot(trend.scaling.mod.fit.mcmc,parms = c('mu.db1','mu.db2'))
    # traplot(trend.scaling.mod.fit.mcmc,parms = c('mu.db1','mu.db2'))
    # caterplot(trend.scaling.mod.fit.mcmc,parms = c('mu.db1','mu.db2'))
    # ##1: non-ARs; 2: ARx
    # summary(cbind((exp(mu.db1)-1)*100,(exp(mu.db2)-1)*100))
    # 
    # par(mfcol=c(2,2))
    # boxplot(cbind((exp(mu.db1)-1)*100,(exp(mu.db2)-1)*100),col=c('gray50','red'),frame=F,main='mu: non-AR (gray), ARx (red)')
    # boxplot(cbind((exp(sigma.db1)-1)*100,(exp(sigma.db2)-1)*100),col=c('gray50','red'),main='sigma')
    # boxplot((exp(db1)-1)*100,outline=F,col='gray50',main='non-AR')
    # boxplot((exp(db2)-1)*100,outline=F,col='red',main='ARx')
    # 
    #
    mu.db1 <- df.trend.scaling.mod.fit.mcmc$mu.db1
    mu.db2 <- df.trend.scaling.mod.fit.mcmc$mu.db2
    
    sigma.db1 <- df.trend.scaling.mod.fit.mcmc$sigma.db1
    sigma.db2 <- df.trend.scaling.mod.fit.mcmc$sigma.db2
    
    # names1 <- paste0('xi[',seq(1,length(nevents)),']')
    # names2 <- array(NA,c(length(nevents),max(nevents)))
    # for (n2 in 1:max(nevents)){
    #   for (n1 in 1:length(nevents)){
    #     names2[n1,n2] <- paste0('sigma[',n1,',',n2,']')
    #   }
    # }
    # names22 <- as.vector(names2)
    # xi <- df.trend.scaling.mod.fit.mcmc[,names(df.trend.scaling.mod.fit.mcmc)%in%names1]
    # sigma <- df.trend.scaling.mod.fit.mcmc[,names(df.trend.scaling.mod.fit.mcmc)%in%names22]
    # 
    # #fit competing time-varying models based on PRCP vs. TEMP
    # #u + sigma[t]/xi*((return.period*lambda)^xi-1)
    # my.precip <- precip
    # my.temp <- ta_td
    # my.thresh <- threshold_mm
    # my.lambda <- lambda_mean
    # for(i in 1:nsites)
    # {
    #   for(t in 1:nevents[i])
    #   {
    #     my.precip[i,t]<- threshold_mm[i]+sigma[i,t]/xi[i]*((return.period*lambda_mean[i])^xi[i]-1)
    #   }
    #   my.lm.d <- lm(log(my.precip[i,]) ~ my.temp[i,])
    #   site.scaling <- c(100*(exp(my.lm.d$coefficients[2])-1))
    # }
    
    mu.db.coefs <- cbind(mu.db1,mu.db2)
    mu.db.coefs.prct <- (exp(mu.db.coefs)-1)*100
    
    sigma.db.coefs <- cbind(sigma.db1,sigma.db2)
    sigma.db.coefs.prct <- (exp(sigma.db.coefs)-1)*100
    
    # db1.coefs.prct <- (exp(db1)-1)*100
    # db2.coefs.prct <- (exp(db2)-1)*100
    # 
    # da1.coefs <- da1
    # da2.coefs <- da2
    
    #sigma.GPD.param <- apply(sigma,FUN=median,c(2,3))
    # calculating Rhats and convergence checks
    # mod_fit_summary = trend.scaling.mod.fit$BUGSoutput$summary
    # mod_all_convergence = trend.scaling.mod.fit$BUGSoutput$summary[,8:9]
    # 
    # mod_rhat_summary = c(min(mod_all_convergence[,1]),quantile(mod_all_convergence[,1],c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),max(mod_all_convergence[,1]),mean(mod_all_convergence[,1]))
    # mod_neff_summary = c(min(mod_all_convergence[,2]),quantile(mod_all_convergence[,2],c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),max(mod_all_convergence[,2]),mean(mod_all_convergence[,2]))
    #print(rbind(mod_rhat_summary,mod_neff_summary))
    
    # saving the parameters #
    lst.names.order <- c(tolower(precip_selec.opts[pr]),
                         tolower(ta_td_selec[te]),
                         'nevents','nsites','n.iter',
                         'mu.db.coefs.prct','sigma.db.coefs.prct',
                         #'db1.coefs.prct','db2.coefs.prct',
                         #'da1.coefs','da2.coefs',
                         #'xi','sigma.GPD.param',
                         'fdates','AR.info',
                         #'rhat_summary',
                         #'neff_summary',
                         'idx.keep',
                         'df.trend.scaling.mod.fit.mcmc')
    
    lst.trend.scaling.mod.fit <- list(precip,
                                      ta_td,
                                      nevents,nsites,n.iter,
                                      mu.db.coefs.prct,sigma.db.coefs.prct,
                                      #db1.coefs.prct,db2.coefs.prct,
                                      #da1.coefs,da2.coefs,
                                      #xi,sigma.GPD.param,
                                      fdates,AR.info,
                                      #mod_rhat_summary,
                                      #mod_neff_summary,
                                      idx.keep,
                                      df.trend.scaling.mod.fit.mcmc)
    
    names(lst.trend.scaling.mod.fit) <- lst.names.order
    
    saveRDS(lst.trend.scaling.mod.fit,
            file=paste0("./Data/model.mcmc.run.limited/trend.scaling.mod.multi.bayes.ARs.event_based.fit.850Combo.",
                        #file=paste0("./Data/model.mcmc.run.limited/trend.scaling.mod.multi.bayes.ARs.event_based.fit.",
                        temporal_selec,".",
                        ta_td_selec,
                        ".",precip_selec,
                        ".",n.iter2,"iter.rds"))
    print(paste('Finished:::: Temperature:',ta_td_selec,
                '; Precipitation:',precip_selec))
    rm(df.trend.scaling.mod.fit.mcmc,trend.scaling.mod.fit.run.jags,
       jags.data)
    gc()
  }
}

#


