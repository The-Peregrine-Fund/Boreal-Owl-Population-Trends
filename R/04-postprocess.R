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

load("data//data_s1.Rdata") # data for site 1
d1 <- list(dat=dat, datl=datl, constl=constl)
load("data//data_s2.Rdata") # data for site 2
d2 <- list(dat=dat3, datl=datl, constl=constl)

# Look at global models 
# to assess covariates for significance 
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-site1.Rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\site1_gof.Rdata")
m.site1 <- list(m1=lapply(m[[1]], as.mcmc),  mr=lapply(m.t, as.mcmc))
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-site2.Rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\site2_gof.Rdata")
m.site2 <- list(m1=lapply(m[[1]], as.mcmc),  mr=lapply(m.wp, as.mcmc))
 

# Inspect estimates from global models
# Use probability of direction to evaluate 
# whether an effect exists
mp <- function(x) {cbind( median(x), pd(x))}
pars <- c( "mu", "delta", "sigma")
MCMCsummary(m.site1[[1]], pars, 
            HPD=TRUE,  hpd_prob = 0.95, round=2, 
            func = mp, func_name = c("median", "pd"))
MCMCtrace(m.site1[[1]], pars, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
ps1 <- MCMCpstr(m.site1[[2]], pars, type="chains")
# Interaction with time is significant using p>0.975
# use reduced model dropping delta 7 and 8

MCMCsummary(m.site2[[1]], pars, 
            HPD=TRUE,  hpd_prob = 0.80, round=2, 
            func = mp, func_name = c("median", "pd"))
MCMCtrace(m.site2[[1]], pars, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
ps2 <- MCMCpstr(m.site2[[2]], pars, type="chains") 
# Interaction with time^2 is sig
# use global model

#################
# Retain covariates with the effects that exist
# Determined from probability of direction
# Use these estimates for inference
pars <-c( "mu", "delta", "sigma")
# Global model site 1
coeftab.s1 <- MCMCsummary(m.site1[[1]], pars, 
            HPD=TRUE,  hpd_prob = 0.95, round=2,  
            func = mp, func_name = c("median", "pd"))
coeftab.s1 <- coeftab.s1[, c(7,3,4,5,8)]
write.csv(coeftab.s1, 'docs\\coeftab_s1.csv')
# Subset of global model
coeftab.s1.r <- MCMCsummary(m.site1[[2]], pars, 
                            HPD=TRUE,  hpd_prob = 0.95, round=2, 
                            func = mp, func_name = c("median", "pd"))
coeftab.s1.r <- coeftab.s1.r[c(1:7,10), c(7,3,4,5,8)]
write.csv(coeftab.s1.r, 'docs\\coeftab_s1_retained.csv')

# Global model site 2
coeftab.s2 <- MCMCsummary(m.site2[[1]], pars, 
            HPD=TRUE,  hpd_prob = 0.95, round=3,  
            func = mp, func_name = c("median", "pd"))
coeftab.s2 <- coeftab.s2[, c(7,3,4,5,8)]
write.csv(coeftab.s2, 'docs\\coeftab_s2.csv')
# subset of model site 2
coeftab.s2.r <- MCMCsummary(m.site2[[2]], pars, 
                          HPD=TRUE,  hpd_prob = 0.95, round=3, 
                          func = mp, func_name = c("median", "pd"))
coeftab.s2.r <- coeftab.s2.r[, c(7,3,4,5,8)]
write.csv(coeftab.s2.r, 'docs\\coeftab_s2_retained.csv')

#**********************
#* Plot autocorrelation
#**********************
tiff("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\autocorrelation.tiff", 
     width=8, height=12, 
     units= "in", res=300)
par(mfrow=c(2,1))
acf(d1$datl$y/d1$constl$effort, main="Tadoussac Bird Observatory")
acf(d2$datl$y/d2$constl$effort, main="Whitefish Point", na.action = na.pass)
dev.off()


#**********************
#* plot population cycles over time
#**********************
# Tadoussac
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
                ps1$delta[6,]*sin(2*pi*cycle.sc2[t])*yr.sc[t] 
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
  coord_cartesian(xlim=c(1994, 2023), ylim=c(0, 350)) +
  ggtitle("Tadoussac Bird Observatory")

# summaries
# which year had max median catches
df.ab2$year[which.max(df.ab2$md)]
df.ab2[df.ab2$year>=2004 & df.ab2$year<=2005,]
df.ab2[df.ab2$year>=2022 & df.ab2$year<=2023,]

#################
# Plot captures per mist nest hours 
# for Whitefish Point site
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
    ps2$delta[8,]*sin(2*pi*cycle.sc2[t])*yr.sc[t]^2
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
  geom_line(data=df.dat, aes(x=year, y=count_per_mnh), color="blue", linewidth=1, alpha=0.75) +
  geom_point(data=data.frame(year=c(1994+16, 1994+19), y=c(0,0)), 
             aes(x=year, y=y), col="red", shape=4, size=2, stroke=2) +
  coord_cartesian(xlim=c(1994, 2023), ylim=c(0, 150)) +
  ggtitle("Whitefish Point") 
#+  scale_color_manual(values = cols)

ap23 <- align_plots(p2, p3, align="v", axis="l")
p23 <- plot_grid(ap23[[1]], ap23[[2]], nrow = 2, align="v", 
                 label_y = "Captures")
## ---- lambdaplot --------
p23

## ---- saveplot12 --------
# ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\populationtrends.tiff",
#        plot=p23, device="tiff", width=6, height=8, dpi=300)

# summaries
# which year had max median catches
df.ab2$year[which.max(df.ab2$md)]
