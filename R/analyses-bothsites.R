library (nimble)
library (nimbleHMC)
library (MCMCvis)
library (coda)
library(parallel)
set.seed(5757575)
dat <- read.csv("data\\data.csv")

dat2 <- data.frame(caps=c(dat$TadoussacCapture , dat$MCFallCapture),
                   effort= c(dat$TadoussacEffort , dat$MCFallEffort),
                   season= c(rep("fall", length(dat$TadoussacEffort)),
                             rep("fall", length(dat$MCFallEffort)) ),
                   year=c(dat$Year, dat$Year), 
                   site = as.numeric( factor(c(rep("T", length(dat$TadoussacEffort )),
                            rep("M", length(dat$MCFallEffort)))) )
)
dat2 <- dat2[ !is.na(dat2$effort) & dat2$effort>0,]

datl <- list(y = dat2$caps)
constl <- list(time = (dat2$year - median(1995:2023))/ (max(dat2$year - median(1995:2023))-1),
               effort = dat2$effort/500,
               season = dat2$season,
               ntime = length(dat2$year2),
               site = dat2$site)

constl$cycle <- 1:length(1995:2023)/4
timemerge <- data.frame(year = seq(1,length(1995:2023), by=1)+1994,  
                        year1995 = seq(1,length(1995:2023), by=1),
                        year4 = seq(1,length(1995:2023), by=1)/4
)
dat3 <- merge(dat2, timemerge, by="year")
constl$time2 <- dat3$year4
cbind(constl$time, constl$time2)

plot(dat2$year2, datl$y/constl$effort, type="l")
acf(datl$y/constl$effort)
acf(datl$y/constl$effort, type="p")

run <- function(seed, datl, constl){
  library(nimble)
  library(coda)
  library(nimbleHMC)
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=10) 
      sd1 ~ dexp(1)
      sd2 ~ dexp(1)
      # likelihood
      for (i in 1:nobs){
        y[t] ~ dpois(lambda[t])
        log(lambda[t]) <- mu + log(effort[t]) + 
          eps[year[i]] + eta[site[i], year[i]]
      } # t
      for (t in 1:nyear){
        eps[t] ~ dnorm(0, sd=sd1)
        for (s in 1:nsite){
          eta[s,t] ~ dnorm(0, sd=sd2)  
      }}
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
  
  params <- c( "mu", "sd1", "sd2",
               "eps", "eta",
               "dssm.obs", "dssm.rep",
               "dmape.obs", "dmape.rep",
               "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  inits <- function(){ list (mu = runif(2, -2, 2),
                             sd1 = runif(1),
                             sd2 = runif(1)
  )}
  
  #n.chains=1; n.thin=1; n.iter=2000; n.burnin=1000
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  #n.chains=1; n.thin=1000; n.iter=2000000; n.burnin=1000000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE, 
                     calculate = TRUE)
  
  #mod$calculate()
  
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
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\poisson-site2.Rdata")
