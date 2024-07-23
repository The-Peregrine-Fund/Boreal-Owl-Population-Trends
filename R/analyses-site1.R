dat <- read.csv("data\\data_T.csv")
datl <- list(y = dat$caps)
constl <- list(time = (dat$year - median(dat$year))/ max(dat$year - median(dat$year)),
               effort = dat$effort/500,
               ntime = length(dat$year),
               time2 = 1:length(1996:2023)/4)
d1 <- list(dat=dat, datl=datl, constl=constl)
save(d1, file="data//data_s1.Rdata")
# caps are number of captures
# effort are net hours
library (nimble)
library (nimbleHMC)
library (MCMCvis)
library (coda)
library(parallel)
set.seed(5757575)

########################
#*******************
#* Check goodness-of-fit
#* for Pois, NB, and ZIP
#*******************
########################

###############
# function to plot posterior predictive checks
###############
plot.diag <- function(out, ratio=FALSE, lab=""){
  par(mfrow=c(1,1))
  # plot mean absolute percentage error
  samps <- MCMCpstr(out, "all", type="chains")
  mx <- max(c(samps$dmape.rep[1,], samps$dmape.obs[1,]))
  mn <- min(c(samps$dmape.rep[1,], samps$dmape.obs[1,]))
  plot(jitter(samps$dmape.obs[1,]), 
       jitter(samps$dmape.rep[1,]),
       main=paste0("Mean absolute percentage error\nmodel\n",lab),
       ylab="Discrepancy replicate values",
       xlab="Discrepancy observed values", 
       xlim=c(mn,mx), ylim=c(mn,mx), 
       pch=16, cex=0.5, col="gray10")
  curve(1*x, from=mn, to=mx, add=T, lty=2, lwd=2, col="blue")
  bp1 <- round(mean(samps$dmape.rep[1,] > samps$dmape.obs[1,]),2)
  loc <- ifelse(bp1 < 0.5, "topleft", "bottomright")
  legend(loc, legend=bquote(p[B]==.(bp1)), bty="n", cex=2)
  
  if (ratio==TRUE){
    # plot variance/mean ratio
    hist(samps$tvm.rep[1,], nclass=50,
         xlab="variance/mean ", main=NA, axes=FALSE)
    abline(v=samps$tvm.obs[1,1], col="red")
    axis(1); axis(2)
  }
  return(list('Bayesian p-value'=bp1))
}
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
        #p[t] <- r/(r+lambda[t])
        #y[t] ~ dnegbin(p[t], r)
        y[t] ~ dpois(lambda[t])
        # y[t] ~ dZIP(lambda[t], psi)
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
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)
params_p <-c( "mu", "delta", "sigma")

pois <- list(as.mcmc(post[[1]]), 
           as.mcmc(post[[2]]), 
           as.mcmc(post[[3]]),
           as.mcmc(post[[4]]))

MCMCtrace(pois, params_p, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
MCMCsummary(pois, params_p, 
            HPD=TRUE,  hpd_prob = 0.95,
            round=3, pg0 = TRUE, func = median)
MCMCsummary(pois, params_p, 
            HPD=TRUE,  hpd_prob = 0.85,
            round=3, pg0 = TRUE, func = median)
par(mfrow=c(1,1))
MCMCplot(object = pois, params = 'mu', HPD=TRUE, 
         ci=c(85, 95))
MCMCplot(object = pois, params = 'delta', HPD=TRUE, 
         ci=c(85, 95))

plot.diag(pois)
plot.diag(post, ratio=TRUE)

save(out=pois, post=post,
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\poisson.Rdata")

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
          # y[t] ~ dpois(lambda[t])
          # y[t] ~ dZIP(lambda[t], psi)
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
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)
params_nb <-c( "mu", "delta", "r", "sigma")

nb <- list(as.mcmc(post[[1]]), 
           as.mcmc(post[[2]]), 
           as.mcmc(post[[3]]),
           as.mcmc(post[[4]]))

MCMCtrace(nb, params_nb, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
MCMCsummary(nb, params_nb, 
            HPD=TRUE,  hpd_prob = 0.95,
            round=3, pg0 = TRUE, func = median)
MCMCsummary(nb, params_nb, 
            HPD=TRUE,  hpd_prob = 0.85,
            round=3, pg0 = TRUE, func = median)
par(mfrow=c(1,1))
MCMCplot(object = nb, params = 'mu', HPD=TRUE, ci=c(85, 95))
MCMCplot(object = nb, params = 'r', HPD=TRUE, ci=c(85, 95))
MCMCplot(object = nb, params = c('delta'), HPD=TRUE, ci=c(85, 95))

plot.diag(nb)

save(out=nb, post=post,
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\nb.Rdata")

########################
#***********************
#* Zero-inflated Poisson
#***********************
########################
run <- function(seed, datl, constl){
  library(nimble, warn.conflicts = F)
  library(coda)
  library(nimbleHMC)
  
  dZIP <- nimbleFunction(
    run = function(x = integer(), lambda = double(), 
                   zeroProb = double(), log = logical(0, default = 0)) {
      returnType(double())
      ## First handle non-zero data
      if (x != 0) {
        ## return the log probability if log = TRUE
        if (log) return(dpois(x, lambda, log = TRUE) + log(1 - zeroProb))
        ## or the probability if log = FALSE
        else return((1 - zeroProb) * dpois(x, lambda, log = FALSE))
      }
      ## From here down we know x is 0
      totalProbZero <- zeroProb + (1 - zeroProb) * dpois(0, lambda, log = FALSE)
      if (log) return(log(totalProbZero))
      return(totalProbZero)
    })
  assign('dZIP', dZIP, envir = .GlobalEnv)
  
  rZIP <- nimbleFunction(
    run = function(n = integer(), lambda = double(), zeroProb = double()) {
      returnType(integer())
      isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
      if (isStructuralZero) return(0)
      return(rpois(1, lambda))
    })
  assign('rZIP', rZIP, envir = .GlobalEnv)
  
  registerDistributions(list(
    dZIP = list(
      BUGSdist = "dZIP(lambda, zeroProb)",
      discrete = TRUE,
      range = c(0, Inf),
      types = c('value = integer()', 'lambda = double()', 'zeroProb = double()')
    )))
  
  code <- nimbleCode(
    {
      # priors
      #r ~ dexp(1)  
      mu ~ dnorm(0, sd=10) 
      for (i in 1:8){
        delta[i] ~ dnorm(0, sd=10)
      }
      psi ~ dunif(0,1)
      sigma ~ dunif(0,10)
      # likelihood
      for (t in 1:ntime){
        # p[t] <- r/(r+lambda[t])
        # y[t] ~ dnegbin(p[t], r)
        # y[t] ~ dpois(lambda[t])
        #lambda.star[t] <- lambda[t]*z[t]
        #y[t] ~ dZIP (lambda[t], zeroProb = psi)
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
        #c.rep[t] ~ dZIP(lambda[t], zeroProb= psi) # expected counts
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
  
  #n.chains=1; n.thin=1; n.iter=2000; n.burnin=1000
  n.chains=1; n.thin=50; n.iter=100000; n.burnin=50000
  #n.chains=1; n.thin=1000; n.iter=2000000; n.burnin=1000000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod )
  
  #confhmc <- configureHMC(mod)
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
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run, 
                  dat=datl, 
                  const=constl)

stopCluster(this_cluster)

params_zip <-c( "mu", "delta", "psi", "sigma")
zip <- list(as.mcmc(post[[1]]), 
           as.mcmc(post[[2]]), 
           as.mcmc(post[[3]]),
           as.mcmc(post[[4]]))

MCMCtrace(zip, params_zip, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
MCMCsummary(zip, params_zip, 
            HPD=TRUE,  hpd_prob = 0.95,
            round=3, pg0 = TRUE, func = median)
MCMCsummary(zip, params_zip, 
            HPD=TRUE,  hpd_prob = 0.85,
            round=3, pg0 = TRUE, func = median)
par(mfrow=c(1,1))
MCMCplot(object = zip, params = 'mu', HPD=TRUE, ci=c(85, 95))
MCMCplot(object = zip, params = 'psi', HPD=TRUE, ci=c(85, 95))
MCMCplot(object = zip, params = c('delta'), HPD=TRUE, ci=c(85, 95))

plot.diag(zip)

save(out=zip, post=post,
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\zip.Rdata")




########################
#***********************
#* Negative binomial model
#* with population cycles on lambda
#* and WAIC
#***********************
########################
m <- list()
##################
# M1- GLOBAL MODEL
##################
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
    }
  ) # end model
  
  params <- c( "mu", "delta", "sigma", "eps", "lambda", "logProb_y")
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             lambda = rpois(constl$ntime, mean(datl$y)), 
                             delta = runif(8, -2, 2), 
                             sigma = runif(1) 
  )}
  
  n.chains=4; n.thin=100; n.iter=200000; n.burnin=100000
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  mod$calculate()
  cmod <- compileNimble(mod )
  confhmc <- configureMCMC(mod, enableWAIC = TRUE)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  m[[1]] <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = TRUE,
                  WAIC=TRUE)

##################
# M2 - no quad interaction
##################
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5)
      for (n in 1:6){
        delta[n] ~ dnorm(0,sd=5)
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
          #delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
          #delta[8]*sin(2*3.141593*time2[t])*time[t]^2
          eps[t]
        
        eps[t] ~ dnorm(0, sd=sigma)
      } # t
    }
  ) # end model

  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             lambda = rpois(constl$ntime, mean(datl$y)), 
                             delta = runif(6, -2, 2), 
                             sigma = runif(1)
  )}
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  mod$calculate()
  cmod <- compileNimble(mod )
  confhmc <- configureMCMC(mod, enableWAIC = TRUE)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  m[[2]] <- runMCMC(chmc,
                niter = n.iter, 
                nburnin = n.burnin,
                nchains = n.chains,
                thin = n.thin,
                samplesAsCodaMCMC = TRUE,
                WAIC=TRUE)
##################
# M3 - no interaction
##################  
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5)
      for (n in 1:4){
        delta[n] ~ dnorm(0,sd=5)
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
          eps[t]
        
        eps[t] ~ dnorm(0, sd=sigma)
        # delta[5]*cos(2*3.141593*time2[t])*time[t] + 
        #delta[6]*sin(2*3.141593*time2[t])*time[t] 
        #delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
        #delta[8]*sin(2*3.141593*time2[t])*time[t]^2
      } # t
    }
  ) # end model
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             lambda = rpois(constl$ntime, mean(datl$y)), 
                             delta = runif(4, -2, 2), 
                             sigma = runif(1)
  )}
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  mod$calculate()
  cmod <- compileNimble(mod )
  confhmc <- configureMCMC(mod, enableWAIC = TRUE)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  m[[3]] <- runMCMC(chmc,
                niter = n.iter, 
                nburnin = n.burnin,
                nchains = n.chains,
                thin = n.thin,
                samplesAsCodaMCMC = TRUE,
                WAIC=TRUE)

##################
# M4 - no quad year
##################  
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5)
      for (n in 1:3){
        delta[n] ~ dnorm(0,sd=5)
      }
      sigma ~ dunif(0,10)
      
      # likelihood
      for (t in 1:ntime){
        y[t] ~ dpois(lambda[t])
        log(lambda[t]) <- mu + log(effort[t]) + 
          delta[1]*time[t] + 
        #  delta[2]*time[t]^2 +
          delta[2]*cos(2*3.141593*time2[t]) + 
          delta[3]*sin(2*3.141593*time2[t]) +
          eps[t]
        
        eps[t] ~ dnorm(0, sd=sigma)
        # delta[5]*cos(2*3.141593*time2[t])*time[t] + 
        #delta[6]*sin(2*3.141593*time2[t])*time[t] 
        #delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
        #delta[8]*sin(2*3.141593*time2[t])*time[t]^2
      } # t
    }
  ) # end model
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             lambda = rpois(constl$ntime, mean(datl$y)), 
                             delta = runif(3, -2, 2), 
                             sigma = runif(1)
  )}
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  mod$calculate()
  cmod <- compileNimble(mod )
  confhmc <- configureMCMC(mod, enableWAIC = TRUE)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  m[[4]] <- runMCMC(chmc,
                niter = n.iter, 
                nburnin = n.burnin,
                nchains = n.chains,
                thin = n.thin,
                samplesAsCodaMCMC = TRUE,
                WAIC=TRUE)
  
##################
# M5 - no year
##################  
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5)
      for (n in 1:2){
        delta[n] ~ dnorm(0,sd=5)
      }
      sigma ~ dunif(0,10)
      
      # likelihood
      for (t in 1:ntime){
        y[t] ~ dpois(lambda[t])
        log(lambda[t]) <- mu + log(effort[t]) + 
          # delta[1]*time[t] + 
          #  delta[2]*time[t]^2 +
          delta[1]*cos(2*3.141593*time2[t]) + 
          delta[2]*sin(2*3.141593*time2[t]) +
          eps[t]
        
        eps[t] ~ dnorm(0, sd=sigma)
        # delta[5]*cos(2*3.141593*time2[t])*time[t] + 
        #delta[6]*sin(2*3.141593*time2[t])*time[t] 
        #delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
        #delta[8]*sin(2*3.141593*time2[t])*time[t]^2
      } # t
    }
  ) # end model
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             lambda = rpois(constl$ntime, mean(datl$y)), 
                             delta = runif(2, -2, 2), 
                             sigma = runif(1) 
  )}
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  mod$calculate()
  cmod <- compileNimble(mod )
  confhmc <- configureMCMC(mod, enableWAIC = TRUE)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  m[[5]] <- runMCMC(chmc,
                niter = n.iter, 
                nburnin = n.burnin,
                nchains = n.chains,
                thin = n.thin,
                samplesAsCodaMCMC = TRUE,
                WAIC=TRUE)
  
##################
# M6 - no season
##################  
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5)
      # for (n in 1:2){
      #   delta[n] ~ dnorm(0,sd=5)
      # }
      sigma ~ dunif(0,10)
      
      # likelihood
      for (t in 1:ntime){
        y[t] ~ dpois(lambda[t])
        log(lambda[t]) <- mu + log(effort[t]) +
          eps[t]
        
        eps[t] ~ dnorm(0, sd=sigma)  
          # delta[1]*time[t] + 
          #  delta[2]*time[t]^2 +
          #delta[1]*cos(2*3.141593*time2[t]) + 
          #delta[2]*sin(2*3.141593*time2[t]) 
        # delta[5]*cos(2*3.141593*time2[t])*time[t] + 
        #delta[6]*sin(2*3.141593*time2[t])*time[t] 
        #delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
        #delta[8]*sin(2*3.141593*time2[t])*time[t]^2
      } # t
    }
  ) # end model
  
  params <- c( "mu", "sigma", "eps", "lambda", "logProb_y")
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             lambda = rpois(constl$ntime, mean(datl$y)),
                             sigma = runif(1) 
                             #delta = runif(2, -2, 2)
  )}
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  mod$calculate()
  cmod <- compileNimble(mod )
  confhmc <- configureMCMC(mod, enableWAIC = TRUE)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  m[[6]] <- runMCMC(chmc,
                niter = n.iter, 
                nburnin = n.burnin,
                nchains = n.chains,
                thin = n.thin,
                samplesAsCodaMCMC = TRUE,
                WAIC=TRUE)  
###############
# M7 TIME + YEAR + TIME:YEAR
##############
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5)
      for (n in 1:5){
        delta[n] ~ dnorm(0,sd=5)
      }
      sigma ~ dunif(0,10)
      
      # likelihood
      for (t in 1:ntime){
        y[t] ~ dpois(lambda[t])
        log(lambda[t]) <- mu + log(effort[t]) + 
          delta[1]*time[t] + 
          #  delta[2]*time[t]^2 +
          delta[2]*cos(2*3.141593*time2[t]) + 
          delta[3]*sin(2*3.141593*time2[t]) +
          delta[4]*cos(2*3.141593*time2[t])*time[t]+
          delta[5]*sin(2*3.141593*time2[t])*time[t] +
          eps[t]
        
        eps[t] ~ dnorm(0, sd=sigma)
        # delta[5]*cos(2*3.141593*time2[t])*time[t] + 
        #delta[6]*sin(2*3.141593*time2[t])*time[t] 
        #delta[7]*cos(2*3.141593*time2[t])*time[t]^2 + 
        #delta[8]*sin(2*3.141593*time2[t])*time[t]^2
      } # t
    }
  ) # end model
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(1, -2, 2),
                             lambda = rpois(constl$ntime, mean(datl$y)), 
                             delta = runif(5, -2, 2), 
                             sigma = runif(1)
  )}
  params <- c( "mu", "sigma", "delta", "eps", "lambda", "logProb_y")
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  mod$calculate()
  cmod <- compileNimble(mod )
  confhmc <- configureMCMC(mod, enableWAIC = TRUE)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
  
  m[[7]] <- runMCMC(chmc,
                    niter = n.iter, 
                    nburnin = n.burnin,
                    nchains = n.chains,
                    thin = n.thin,
                    samplesAsCodaMCMC = TRUE,
                    WAIC=TRUE)

save(m=m, 
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-set_population-cycles.Rdata")
