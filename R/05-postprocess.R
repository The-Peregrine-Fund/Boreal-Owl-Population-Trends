library(MCMCvis)
library (coda)
library (ggplot2)
library (reshape2)
library (tidybayes)
library(tidyr)
library(cubelyr)
library(HDInterval)
library(cowplot)
library (bayestestR)
options(scipen=999)
source("R//plot_diag_function.R")

load("data//data_s1.Rdata") # data for site 1
d1 <- list(dat=dat, datl=datl, constl=constl)
load("data//data_s2.Rdata") # data for site 2
d2 <- list(dat=dat3, datl=datl, constl=constl)

# Look at global models 
# to assess covariates for significance 
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-site1.Rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\site1_gof.Rdata")
m.site1 <- list(m1=lapply(m[[1]], as.mcmc),  
                m2=lapply(m[[2]], as.mcmc),
                m3=lapply(m[[3]], as.mcmc),
                mr=lapply(m.t, as.mcmc))
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-site2.Rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\site2_gof.Rdata")
m.site2 <- list(m1=lapply(m[[1]], as.mcmc),  
                m2=lapply(m[[2]], as.mcmc),
                m3=lapply(m[[3]], as.mcmc),
                mr=lapply(m.wp, as.mcmc))
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-site1-juvs.Rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\site1_gof-juvs.Rdata")
m.site1j <- list(m1=lapply(m[[1]], as.mcmc),  
                 m2=lapply(m[[2]], as.mcmc),
                 m3=lapply(m[[3]], as.mcmc),
                 mr=lapply(m.tj, as.mcmc)) 

#*****************
#* Check goodness of fit 
#*****************
#*Tadoussac
plot.diag(m.site1$m1)
plot.diag(m.site1$m1, ratio=TRUE)
plot.diag(m.site1$m2)
plot.diag(m.site1$m3)

plot.diag(m.site1j$m1)
plot.diag(m.site1j$m1, ratio=TRUE)
plot.diag(m.site1j$m2)
plot.diag(m.site1j$m3)

plot.diag(m.site2$m1)
plot.diag(m.site2$m1, ratio=TRUE)
plot.diag(m.site2$m2)
plot.diag(m.site2$m3)

# Inspect posteriors from global models for convergence
pars <- c( "mu", "delta", "sigma")
MCMCtrace(m.site1[[2]], pars, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
ps1 <- MCMCpstr(m.site1[[2]], pars, type="chains")

MCMCtrace(m.site1j[[2]], pars, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
ps1j <- MCMCpstr(m.site1j[[2]], pars, type="chains")

MCMCtrace(m.site2[[2]], pars, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
ps2 <- MCMCpstr(m.site2[[2]], pars, type="chains") 

#*********************
#* Model estimates
#*********************
#################
# Check probability of direction from global models
##################
pars <-c( "mu", "delta", "sigma")
mp <- function(x) {cbind( median(x), pd(x))}
# Global model site 1 Tadoussac
coeftab.s1 <- MCMCsummary(m.site1[[1]], pars, 
                          HPD=TRUE,  hpd_prob = 0.95, round=2,  
                          func = mp, func_name = c("median", "pd"))
coeftab.s1$site <- "Tadoussac"
coeftab.s1 <- coeftab.s1[, c(9,7,3,4,5,8)]

# Global model site 1 Tadoussac juveniles
coeftab.s1j <- MCMCsummary(m.site1j[[1]], pars, 
                           HPD=TRUE,  hpd_prob = 0.95, round=2,  
                           func = mp, func_name = c("median", "pd"))
coeftab.s1j$site <- "Tadoussac juveniles"
coeftab.s1j <- coeftab.s1j[, c(9,7,3,4,5,8)]

# Global model site 2 Whitefish Point
coeftab.s2 <- MCMCsummary(m.site2[[1]], pars, 
                          HPD=TRUE,  hpd_prob = 0.95, round=3,  
                          func = mp, func_name = c("median", "pd"))
coeftab.s2$site <- "Whitefish Point"
coeftab.s2 <- coeftab.s2[, c(9,7,3,4,5,8)]


rbind(coeftab.s1, coeftab.s1j, coeftab.s2) |>
  write.csv( 'docs\\coeftab_global models.csv')

# Get estimates from models that
# retained covariates with significant effects 
# Determined from probability of direction
# Use these estimates for inference

# Subset of global model site 1 Tadoussac
coeftab.s1.r <- MCMCsummary(m.site1[[2]], pars, 
                            HPD=TRUE,  hpd_prob = 0.95, round=2, 
                            func = mp, func_name = c("median", "pd"))
coeftab.s1.r$site <- "Tadoussac"
coeftab.s1.r <- coeftab.s1.r[c(1:7,10), c(9,7,3,4,5,8)]

# Subset of global model Tadoussac juveniles
coeftab.s1j.r <- MCMCsummary(m.site1j[[2]], pars, 
                             HPD=TRUE,  hpd_prob = 0.95, round=2, 
                             func = mp, func_name = c("median", "pd"))
coeftab.s1j.r$site <- "Tadoussac juveniles"
coeftab.s1j.r <- coeftab.s1j.r[c(1,2, 4:7, 10), c(9,7,3,4,5,8)]

# subset of model site 2 Whitefish Point
coeftab.s2.r <- MCMCsummary(m.site2[[2]], pars, 
                          HPD=TRUE,  hpd_prob = 0.95, round=3, 
                          func = mp, func_name = c("median", "pd"))
coeftab.s2.r$site <- "Whitefish Point"
coeftab.s2.r <- coeftab.s2.r[, c(9,7,3,4,5,8)]

rbind(coeftab.s1.r, coeftab.s1j.r, coeftab.s2.r) |>
  write.csv('docs\\coeftab_reduced models.csv')

#**********************
#* Plot autocorrelation
#**********************
tiff("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\autocorrelation.tiff", 
     width=8, height=12, 
     units= "in", res=300)
par(mfrow=c(3,1))
acf(d1$datl$y/d1$constl$effort, main="Tadoussac Bird Observatory")
acf(d1$datl$y_juv/d1$constl$effort, main="Tadoussac Bird Observatory\nJuveniles")
acf(d2$datl$y/d2$constl$effort, main="Whitefish Point", na.action = na.pass)
dev.off()

#**********************
#* plot population cycles over time
#**********************
# Tadoussac
ps1 <- MCMCpstr(m.site1[[2]], c(pars, "eps"), type="chains")
yr <- seq(1996, 2023, by=0.25)
cycle.sc2 <- (yr-1995)/4
yr.sc <- (yr-median(d1$dat$year)) / max(d1$dat$year - median(d1$dat$year))
  
l.ab2 <- l.ab3 <- array(NA, dim=c(length(yr), 4000), 
             dimnames=list(year=yr, iter=c(1:4000)))
for (t in 1:length(yr)){
  l.ab2[t, ] <- ps1$mu + 
                ps1$delta[1,]*yr.sc[t] + 
                ps1$delta[2,]*yr.sc[t]^2 +
                ps1$delta[3,]*cos(2*pi*cycle.sc2[t]) + 
                ps1$delta[4,]*sin(2*pi*cycle.sc2[t]) +
                ps1$delta[5,]*cos(2*pi*cycle.sc2[t])*yr.sc[t] + 
                ps1$delta[6,]*sin(2*pi*cycle.sc2[t])*yr.sc[t] +
                ps1$eps[round(yr[t])-1995,]
}

mab2 <- l.ab2 %>%
  as.tbl_cube(met_name = "l.ab") %>%
  as_tibble

df.ab2 <- data.frame( year = yr,
                     md = apply(l.ab2, 1, median, na.rm=T) |> exp(),
                     lci95 = apply(l.ab2, 1, HDInterval::hdi)[1,]  |> exp(),
                     uci95 = apply(l.ab2, 1, HDInterval::hdi)[2,] |> exp(),
                     lci80 = apply(l.ab2, 1, HDInterval::hdi, credMass=.8)[1,] |> exp(),
                     uci80 = apply(l.ab2, 1, HDInterval::hdi, credMass=.8)[2,] |> exp()
)

df.dat <- data.frame(count=d1$datl$y/d1$constl$effort,
                     year=1996:2023)
  
col <- "black"
cols <- c("Data"="blue", "Median"="black", 
          "85% HDIs"= "black", "95% HDIs"="black")
szs <- c("Data"= 1, "Median"=1, 
         "80% HDIs"= 0.5, "95% HDIs"=0.25)
p2 <- ggplot() + theme_minimal() + 
  geom_line(data=mab2, aes(x=year, y=exp(l.ab), group=iter), 
            color="gray40", linewidth=0.25, alpha=0.05 ) +
  #geom_line(data=df.ab3, aes(x=year, y=md), color=col, linewidth=1, linetype="dashed" ) +
  xlab("Year") + ylab("Captures\nper 500 mist net hours") +
  geom_line(data=df.ab2, aes(x=year, y=lci95), color=col, linewidth=0.25 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci95), color=col, linewidth=0.25) +
  geom_line(data=df.ab2, aes(x=year, y=lci80), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci80), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab2, aes(x=year, y=md), color=col, linewidth=1 ) +
  geom_line(data=df.dat, aes(x=year, y=count), color="blue", linewidth=1, alpha=0.75) +
  #geom_line(data=df.ab3, aes(x=year, y=md), color="black", linewidth=1, linetype=2) +
  coord_cartesian(xlim=c(1995, 2023), ylim=c(0, 175)) +
  ggtitle("Tadoussac Bird Observatory")

# summaries
# which year had max median catches
df.ab2$year[which.max(df.ab2$md)]
df.ab2[df.ab2$year>=2004 & df.ab2$year<=2005,]
df.ab2[df.ab2$year>=2022 & df.ab2$year<=2023,]

#####################
# Tadoussac JUVENILES
ps1j <- MCMCpstr(m.site1j[[2]], c(pars, "eps"), type="chains")
yr <- seq(1996, 2023, by=0.25)
cycle.sc2 <- (yr-1995)/4
yr.sc <- (yr-median(d1$dat$year)) / max(d1$dat$year - median(d1$dat$year))

l.ab2 <- l.ab3 <- array(NA, dim=c(length(yr), 4000), 
                        dimnames=list(year=yr, iter=c(1:4000)))
for (t in 1:length(yr)){
  l.ab2[t, ] <- ps1j$mu + 
    ps1j$delta[1,]*yr.sc[t] + 
    ps1j$delta[3,]*cos(2*pi*cycle.sc2[t]) + 
    ps1j$delta[4,]*sin(2*pi*cycle.sc2[t]) +
    ps1j$delta[5,]*cos(2*pi*cycle.sc2[t])*yr.sc[t] + 
    ps1j$delta[6,]*sin(2*pi*cycle.sc2[t])*yr.sc[t] +
    ps1j$eps[round(yr[t])-1995,]
}

mab2 <- l.ab2 %>%
  as.tbl_cube(met_name = "l.ab") %>%
  as_tibble

df.ab2 <- data.frame( year = yr,
                      md = apply(l.ab2, 1, median, na.rm=T) |> exp(),
                      lci95 = apply(l.ab2, 1, HDInterval::hdi)[1,]  |> exp(),
                      uci95 = apply(l.ab2, 1, HDInterval::hdi)[2,] |> exp(),
                      lci80 = apply(l.ab2, 1, HDInterval::hdi, credMass=.8)[1,] |> exp(),
                      uci80 = apply(l.ab2, 1, HDInterval::hdi, credMass=.8)[2,] |> exp()
)

df.dat <- data.frame(count=d1$datl$y_juv /d1$constl$effort,
                     year=1996:2023)

col <- "black"
cols <- c("Data"="blue", "Median"="black", 
          "85% HDIs"= "black", "95% HDIs"="black")
szs <- c("Data"= 1, "Median"=1, 
         "80% HDIs"= 0.5, "95% HDIs"=0.25)
p3 <- ggplot() + theme_minimal() + 
  geom_line(data=mab2, aes(x=year, y=exp(l.ab), group=iter), 
            color="gray40", linewidth=0.25, alpha=0.05 ) +
  #geom_line(data=df.ab3, aes(x=year, y=md), color=col, linewidth=1, linetype="dashed" ) +
  xlab("Year") + ylab("Captures\nper 500 mist net hours") +
  geom_line(data=df.ab2, aes(x=year, y=lci95), color=col, linewidth=0.25 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci95), color=col, linewidth=0.25) +
  geom_line(data=df.ab2, aes(x=year, y=lci80), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci80), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab2, aes(x=year, y=md), color=col, linewidth=1 ) +
  geom_line(data=df.dat, aes(x=year, y=count), color="blue", linewidth=1, alpha=0.75) +
  #geom_line(data=df.ab3, aes(x=year, y=md), color="black", linewidth=1, linetype=2) +
  coord_cartesian(xlim=c(1995, 2023), ylim=c(0, 75)) +
  ggtitle("Tadoussac Bird Observatory, Juveniles")

# summaries
# which year had max median catches
df.ab2$year[which.max(df.ab2$md)]
df.ab2[df.ab2$year>=2004 & df.ab2$year<=2005,]
df.ab2[df.ab2$year>=2022 & df.ab2$year<=2023,]

#################
# Plot captures per mist nest hours 
# for Whitefish Point site
ps2 <- MCMCpstr(m.site2[[2]], c(pars, "eps"), type="chains")
range(d2$dat$year)
yr <- seq(1995, 2022, by=0.25)
cycle.sc2 <- (yr-1994)/4
yr.sc <- (yr-2008.5) / ((2022-1995)/2)

l.ab2 <- l.ab3 <- array(NA, dim=c(length(yr), 4000), 
                        dimnames=list(year=yr, iter=c(1:4000)))
for (t in 1:length(yr)){
  l.ab2[t, ] <- ps2$mu + 
    ps2$delta[1,]*yr.sc[t] + 
    ps2$delta[2,]*yr.sc[t]^2 +
    ps2$delta[3,]*cos(2*pi*cycle.sc2[t]) + 
    ps2$delta[4,]*sin(2*pi*cycle.sc2[t]) +
    ps2$delta[5,]*cos(2*pi*cycle.sc2[t])*yr.sc[t] + 
    ps2$delta[6,]*sin(2*pi*cycle.sc2[t])*yr.sc[t] +
    ps2$delta[7,]*cos(2*pi*cycle.sc2[t])*yr.sc[t]^2 + 
    ps2$delta[8,]*sin(2*pi*cycle.sc2[t])*yr.sc[t]^2 +
    ps2$eps[round(yr[t])-1994,]
}

mab2 <- l.ab2 %>%
  as.tbl_cube(met_name = "l.ab") %>%
  as_tibble

df.ab2 <- data.frame( year = yr,
                      md = apply(l.ab2, 1, median, na.rm=T) |> exp(),
                      lci95 = apply(l.ab2, 1, HDInterval::hdi)[1,]  |> exp(),
                      uci95 = apply(l.ab2, 1, HDInterval::hdi)[2,] |> exp(),
                      lci80 = apply(l.ab2, 1, HDInterval::hdi, credMass=.8)[1,] |> exp(),
                      uci80 = apply(l.ab2, 1, HDInterval::hdi, credMass=.8)[2,] |> exp()
)

df.dat <- data.frame(count_per_mnh=d2$datl$y/d2$constl$effort,
                     year=d2$dat$year)

col <- "black"
cols <- c("Data"="blue", "Median"="black", 
          "85% HDIs"= "black", "95% HDIs"="black")
szs <- c("Data"= 1, "Median"=1, 
         "80% HDIs"= 0.5, "95% HDIs"=0.25)
p4 <- ggplot() + theme_minimal() + 
  geom_line(data=mab2, aes(x=year, y=exp(l.ab), group=iter), 
            color="gray40", linewidth=0.25, alpha=0.05 ) +
  #geom_line(data=df.ab3, aes(x=year, y=md), color=col, linewidth=1, linetype="dashed" ) +
  xlab("Year") + ylab("Captures\nper 500 mist net hours") +
  geom_line(data=df.ab2, aes(x=year, y=lci95), color=col, linewidth=0.25 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci95), color=col, linewidth=0.25) +
  geom_line(data=df.ab2, aes(x=year, y=lci80), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci80), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab2, aes(x=year, y=md), color=col, linewidth=1 ) +
  geom_line(data=df.dat, aes(x=year, y=count_per_mnh), color="blue", linewidth=1, alpha=0.75) +
  geom_point(data=data.frame(year=c(1994+16, 1994+19), y=c(0,0)), 
             aes(x=year, y=y), col="red", shape=4, size=2, stroke=2) +
  coord_cartesian(xlim=c(1995, 2023), ylim=c(0, 75)) +
  ggtitle("Whitefish Point") 
#+  scale_color_manual(values = cols)

ap234 <- align_plots(p2, p3, p4, align="v", axis="l")
p234 <- plot_grid(ap234[[1]], ap234[[2]], ap234[[3]], 
                 nrow = 3, align="v", 
                 label_y = "Captures")
## ---- lambdaplot --------
p234

## ---- saveplot12 --------
ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\populationtrends.tiff",
       plot=p234, device="tiff", width=6, height=12, dpi=300)

# summaries
# which year had max median catches
df.ab2$year[which.max(df.ab2$md)]
df.ab2[df.ab2$year==1996,]
df.ab2[df.ab2$year>=2022,]
df.ab2[df.ab2$year %in% c(2014:2022),]
