# Tadoussac Bird Observatory (OOT), located in Tadoussac, Québec, Canada 
library (nimble)
library (nimbleHMC)
library (MCMCvis)
library (coda)
library(parallel)
set.seed(5757575)
load("data//data_s1.Rdata")
source("R//plot_diag_function.R")
datl$y <- datl$y_juv
# create a list of models for 
# Goodness of fit checks 
# using Poisson, negative binomial, and zero-inflated Poisson
m <- list()
########################
#*******************
#* Poisson model
#*******************
########################
run <- function(seed, datl, constl){
  library(nimble)
  library(coda)
  library(nimbleHMC)
  code <- nimbleCode(
    {
      # priors
      r ~ dexp(1)  
      mu ~ dnorm(0, sd=10) 
      for (i in 1:8){
        delta[i] ~ dnorm(0,sd=10)
      }
      sigma ~ dunif(0,10)
      # likelihood
      for (t in 1:ntime){
        y[t] ~ dpois(lambda[t])
        log(lambda[t]) <- mu + log(effort[t]) + 
          delta[1]*time[t] + 
          delta[2]*time[t]^2 +
          delta[3]*cos(2*3.141593*time2[t]) + 
          delta[4]*sin(2*3.141593*time2[t]) +
          delta[5]*cos(2*3.141593*time2[t])*time[t] + 
          delta[6]*sin(2*3.141593*time2[t])*time[t] +
          delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
          delta[8]*sin(2*3.141593*time2[t])*time[t]^2 +
          eps[t]
        
        eps[t] ~ dnorm(0, sd=sigma)
      } # t
      ###################
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      ###################
      # GOF for model
      for (t in 1:ntime){
          c.obs[t] <- y[t] # observed counts
          c.exp[t] <- lambda[t] # expected counts adult breeder
          c.rep[t] ~ dpois(lambda[t]) # expected counts
          # Compute fit statistics, Mean absolute error
          dssm.obs[t] <- abs( ( (c.obs[t]) - (c.exp[t]) ) / (c.obs[t]+0.001) )
          dssm.rep[t] <- abs( ( (c.rep[t]) - (c.exp[t]) ) / (c.rep[t]+0.001) )
        } # t
      dmape.obs <- sum(dssm.obs[1:ntime])
      dmape.rep <- sum(dssm.rep[1:ntime])
      # variance-mean ratio 
      tvm.rep <- sd(c.rep[1:ntime])^2/mean(c.rep[1:ntime])
      tvm.obs <- sd(y[1:ntime])^2/mean(y[1:ntime])
    }
  ) # end model
  
  params <- c( "mu", "delta", "sigma",
               "dssm.obs", "dssm.rep",
               "dmape.obs", "dmape.rep",
               "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             delta = runif(8, -2, 2),
                             sigma = runif(1)
  )}
  
  #n.chains=1; n.thin=1; n.iter=2000; n.burnin=1000
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  #n.chains=1; n.thin=1000; n.iter=2000000; n.burnin=1000000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod )
  
  confhmc <- configureHMC(mod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  post <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = T)
  
  return(post)
} # run model function end

this_cluster <- makeCluster(4)
m[[1]] <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)
params_p <-c( "mu", "delta", "sigma")

pois <- lapply(m[[1]], as.mcmc)

MCMCtrace(pois, params_p, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
MCMCsummary(pois, params_p, 
            HPD=TRUE,  hpd_prob = 0.95,
            round=3, pg0 = TRUE, func = median)
par(mfrow=c(1,1))
MCMCplot(object = pois, params = 'mu', HPD=TRUE, 
         ci=c(85, 95))
MCMCplot(object = pois, params = 'delta', HPD=TRUE, 
         ci=c(85, 95))

########################
#***********************
#* Negative binomial model
#***********************
########################
run <- function(seed, datl, constl){
  library(nimble)
  library(coda)
  library(nimbleHMC)
  code <- nimbleCode(
    {
      # priors
      r ~ dexp(1)  
      mu ~ dnorm(0, sd=10) 
      for (i in 1:8){
      delta[i] ~ dnorm(0, sd=10)
      }
      sigma ~ dunif(0, 10)
      # likelihood
      for (t in 1:ntime){
          p[t] <- r/(r+lambda[t])
          y[t] ~ dnegbin(p[t], r)
          log(lambda[t]) <- mu + log(effort[t]) + 
            delta[1]*time[t] + 
            delta[2]*time[t]^2 +
            delta[3]*cos(2*3.141593*time2[t]) + 
            delta[4]*sin(2*3.141593*time2[t]) +
            delta[5]*cos(2*3.141593*time2[t])*time[t] + 
            delta[6]*sin(2*3.141593*time2[t])*time[t] +
            delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
            delta[8]*sin(2*3.141593*time2[t])*time[t]^2 +
            eps[t]
          
          eps[t] ~ dnorm(0, sd=sigma)
      } # t
      ###################
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      ###################
      # GOF for model
      for (t in 1:ntime){
          c.obs[t] <- y[t] # observed counts
          c.exp[t] <- lambda[t] # expected counts adult breeder
          c.rep[t] ~ dnegbin(p[t],r) # expected counts
          # Compute fit statistics, Mean absolute error
          dssm.obs[t] <- abs( ( (c.obs[t]) - (c.exp[t]) ) / (c.obs[t]+0.001) )
          dssm.rep[t] <- abs( ( (c.rep[t]) - (c.exp[t]) ) / (c.rep[t]+0.001) )
        } # t
      dmape.obs <- sum(dssm.obs[1:ntime])
      dmape.rep <- sum(dssm.rep[1:ntime])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime])^2/mean(c.rep[1:ntime])
      tvm.obs <- sd(y[1:ntime])^2/mean(y[1:ntime])
    }
  ) # end model
  
  params <- c( "mu", "delta", "sigma", 
               "r", "p", "lambda",
               "dssm.obs", "dssm.rep",
               "dmape.obs", "dmape.rep",
               "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             delta = runif(8, -2, 2),
                             r = runif(1),
                             sigma = runif(1)
  )}
  
  #n.chains=1; n.thin=1; n.iter=2000; n.burnin=1000
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  #n.chains=1; n.thin=1000; n.iter=2000000; n.burnin=1000000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod )
  
  confhmc <- configureHMC(mod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  post <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = T)
  
  return(post)
} # run model function end

this_cluster <- makeCluster(4)
m[[2]] <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)
params_nb <-c( "mu", "delta", "r", "sigma")

nb <- lapply(m[[2]], as.mcmc)

MCMCtrace(nb, params_nb, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
MCMCsummary(nb, params_nb, 
            HPD=TRUE,  hpd_prob = 0.95,
            round=3, pg0 = TRUE, func = median)
par(mfrow=c(1,1))
MCMCplot(object = nb, params = 'mu', HPD=TRUE, ci=c(85, 95))
MCMCplot(object = nb, params = 'r', HPD=TRUE, ci=c(85, 95))
MCMCplot(object = nb, params = c('delta'), HPD=TRUE, ci=c(85, 95))

########################
#***********************
#* Zero-inflated Poisson
#***********************
########################
run <- function(seed, datl, constl){
  library(nimble, warn.conflicts = F)
  library(coda)
  library(nimbleHMC)
  
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=10) 
      for (i in 1:8){
        delta[i] ~ dnorm(0, sd=10)
      }
      psi ~ dunif(0,1)
      sigma ~ dunif(0,10)
      # likelihood
      for (t in 1:ntime){
        y[t] ~ dpois (lambda[t]*z[t])
        z[t] ~ dbern(psi)
        log(lambda[t]) <- mu + log(effort[t]) + 
          delta[1]*time[t] + 
          delta[2]*time[t]^2 +
          delta[3]*cos(2*3.141593*time2[t]) + 
          delta[4]*sin(2*3.141593*time2[t]) +
          delta[5]*cos(2*3.141593*time2[t])*time[t] + 
          delta[6]*sin(2*3.141593*time2[t])*time[t] +
          delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
          delta[8]*sin(2*3.141593*time2[t])*time[t]^2 +
          eps[t]
        
        eps[t] ~ dnorm(0, sd=sigma)
      } # t
      ###################
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      ###################
      # GOF for model
      for (t in 1:ntime){
        c.obs[t] <- y[t] # observed counts
        c.exp[t] <- z[t]*lambda[t] # expected counts 
        c.rep[t] ~ dpois(lambda[t]*psi) # expected counts
        # Compute fit statistics, Mean absolute error
        dssm.obs[t] <- abs( ( (c.obs[t]) - (c.exp[t]) ) / (c.obs[t]+0.001) )
        dssm.rep[t] <- abs( ( (c.rep[t]) - (c.exp[t]) ) / (c.rep[t]+0.001) )
      } # t
      dmape.obs <- sum(dssm.obs[1:ntime])
      dmape.rep <- sum(dssm.rep[1:ntime])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime])^2/mean(c.rep[1:ntime])
      tvm.obs <- sd(y[1:ntime])^2/mean(y[1:ntime])
    }
  ) # end model
  
  params <- c( "mu", "delta", "sigma", 
               "psi", "lambda", #"z",
               "dssm.obs", "dssm.rep",
               "dmape.obs", "dmape.rep",
               "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             delta = runif(8, -1, 1),
                             psi = runif(1, 0.9, 1), 
                             lambda = rpois(constl$ntime, mean(datl$y)),
                             sigma = dunif(1)
  )}
  
  n.chains=1; n.thin=50; n.iter=100000; n.burnin=50000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  cmod <- compileNimble(mod )
  confhmc <- configureMCMC(mod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  post <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = T)
  
  return(post)
} # run model function end
this_cluster <- makeCluster(4)
m[[3]] <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)

params_zip <-c( "mu", "delta", "psi", "sigma")
zip <- lapply(m[[3]], as.mcmc)

MCMCtrace(zip, params_zip, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
MCMCsummary(zip, params_zip, 
            HPD=TRUE,  hpd_prob = 0.95,
            round=3, pg0 = TRUE, func = median)
par(mfrow=c(1,1))
MCMCplot(object = zip, params = 'mu', HPD=TRUE, ci=c(85, 95))
MCMCplot(object = zip, params = 'psi', HPD=TRUE, ci=c(85, 95))
MCMCplot(object = zip, params = c('delta'), HPD=TRUE, ci=c(85, 95))

save(m=m,
     file="outputs\\site1_gof-juvs.Rdata")

#***************************
# compare goodness-of-fit
#***************************
plot.diag(pois)
plot.diag(pois, ratio=TRUE)
plot.diag(nb)
plot.diag(zip)

########################
#***********************
#* Poisson model
#* with longer run
#***********************
########################
##################
# M1- Drop quadratic interaction btw cycles and time
##################
run <- function(seed, datl, constl){
  library(nimble)
  library(coda)
  library(nimbleHMC)
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5)
      for (n in 1:8){
        delta[n] ~ dnorm(0,sd=5)
      }
      sigma ~ dunif(0,10)
      
      # likelihood
      for (t in 1:ntime){
        y[t] ~ dpois(lambda[t])
        log(lambda[t]) <- mu + log(effort[t]) + 
          delta[1]*time[t] + 
          #delta[2]*time[t]^2 +
          delta[3]*cos(2*3.141593*time2[t]) + 
          delta[4]*sin(2*3.141593*time2[t]) +
          delta[5]*cos(2*3.141593*time2[t])*time[t] + 
          delta[6]*sin(2*3.141593*time2[t])*time[t] +
          #delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
          #delta[8]*sin(2*3.141593*time2[t])*time[t]^2 +
          eps[t]
        
        eps[t] ~ dnorm(0, sd=sigma)
      } # t
    }
  ) # end model
  
  params <- c( "mu", "delta", "sigma", "eps", "lambda", "logProb_y")
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             lambda = rpois(constl$ntime, mean(datl$y)), 
                             delta = runif(8, -2, 2), 
                             sigma = runif(1) 
  )}
  
  n.chains=1; n.thin=100; n.iter=200000; n.burnin=100000
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  mod$calculate()
  cmod <- compileNimble(mod )
  confhmc <- configureHMC(mod, enableWAIC = TRUE)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  post <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = T)
  
  return(post)
} # run function end

  this_cluster <- makeCluster(4)
  m.tj <- parLapply(cl = this_cluster, 
                    X = 1:4, 
                    fun = run, 
                    dat=datl, 
                    const=constl)
  stopCluster(this_cluster)

save(m.tj=m.tj,  
     file="outputs\\full-model-site1-juvs.Rdata")
