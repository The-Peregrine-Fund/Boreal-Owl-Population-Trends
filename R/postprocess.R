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
load("data//data_s2.Rdata") # data for site 2

# Look at global models 
# to assess covariates for significance 
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-set_population-cycles.Rdata")
m.site1 <- m
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-set_population-cycles-site2.Rdata")
m.site2 <- m
 

# Inspect estimates from global models
# Use probability of direction to evaluate 
# whether an effect exists
# log(lambda[t]) <- mu + log(effort[t]) + 
#   delta[1]*year[t] + 
#   delta[2]*year[t]^2 +
#   delta[3]*cos(2*3.141593*cycle[t]) + 
#   delta[4]*sin(2*3.141593*cycle[t]) +
#   delta[5]*cos(2*3.141593*cycle[t])*year[t] + 
#   delta[6]*sin(2*3.141593*cycle[t])*year[t] +
#   delta[7]*cos(2*3.141593*cycle[t])*year[t]^2 + 
#   delta[8]*sin(2*3.141593*cycle[t])*year[t]^2 +
#   eps[t]
pars <- c( "mu", "delta", "sigma")
MCMCsummary(m.site1[[1]]$samples, pars, 
            HPD=TRUE,  hpd_prob = 0.80,
            round=2, pg0 = TRUE, 
            func = median, func_name = "median")

ps1 <- MCMCpstr(m.site1[[2]]$samples, pars, type="chains")
# Interaction with time is significant using p<0.1 and p>0.9
# use model m[[2]]

MCMCsummary(m.site2[[8]]$samples, pars, 
            HPD=TRUE,  hpd_prob = 0.80,
            round=2, pg0 = TRUE, 
            func = median, func_name = "median")

ps2 <- MCMCpstr(m.site2[[8]]$samples, pars, type="chains") 
# calculate probability of direction for seasonal differences. 
pmu.diff <- ps2$mu[2,]-ps2$mu[1,]
mean(pmu.diff>0) # retain season

# Interaction with time^2 is sig
# use m[[8]]

#################
# Retain covariates with the effects that exist
# Determined from probability of direction
# Use these estimates for inference
pars <-c( "mu", "delta", "sigma")
# Global model site 1
coeftab.s1 <- MCMCsummary(m.site1[[1]]$samples, pars, 
            HPD=TRUE,  hpd_prob = 0.95,
            round=2, pg0 = TRUE, 
            func = median, func_name = "median")
coeftab.s1 <- coeftab.s1[, c(8,3,4,5,7)]
write.csv(coeftab.s1, 'docs\\coeftab_s1.csv')
# Subset of global model
coeftab.s1.r <- MCMCsummary(m.site1[[2]]$samples, pars, 
                            HPD=TRUE,  hpd_prob = 0.95,
                            round=2, pg0 = TRUE, 
                            func = median, func_name = "median")
coeftab.s1.r <- coeftab.s1.r[, c(8,3,4,5,7)]
write.csv(coeftab.s1.r, 'docs\\coeftab_s1_retained.csv')

# Global model site 2
coeftab.s2 <- MCMCsummary(m.site2[[8]]$samples, pars, 
            HPD=TRUE,  hpd_prob = 0.85,
            round=3, pg0 = TRUE, 
            func = median, func_name = "median")
coeftab.s2 <- coeftab.s2[, c(8,3,4,5,7)]
write.csv(coeftab.s2, 'docs\\coeftab_s2.csv')
# subset of model site 2
coeftab.s2.r <- MCMCsummary(m.site2[[8]]$samples, pars, 
                          HPD=TRUE,  hpd_prob = 0.85,
                          round=3, pg0 = TRUE, 
                          func = median, func_name = "median")
coeftab.s2.r <- coeftab.s2.r[, c(8,3,4,5,7)]
write.csv(coeftab.s2.r, 'docs\\coeftab_s2_retained.csv')

#**********************
#* Plot autocorrelation
#**********************
tiff("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\autocorrelation.tiff", 
     width=8, height=12, 
     units= "in", res=300)
par(mfrow=c(2,1))
acf(d1$datl$y/d1$constl$effort, main="Tadoussac Bird Observatory")
acf(d2$datl$y/d2$constl$effort, main="Whitefish Point")
#acf(d1$datl$y/d1$constl$effort, main="Tadoussac Bird Observatory", type="p")
#acf(d2$datl$y/d2$constl$effort, main="MC", type="p")
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
# calculate a four year moving average
for (t in 17:nrow(l.ab2)){
for (i in 1:ncol(l.ab2)){
l.ab3[t,i] <- mean(l.ab2[(t-16):t, i])
}}

df.ab3 <- data.frame( year = yr,
                      md = apply(l.ab3, 1, median, na.rm=T) |> exp(),
                      lci95 = apply(l.ab3, 1, HDInterval::hdi)[1,]  |> exp(),
                      uci95 = apply(l.ab3, 1, HDInterval::hdi)[2,] |> exp(),
                      lci80 = apply(l.ab3, 1, HDInterval::hdi, credMass=.8)[1,] |> exp(),
                      uci80 = apply(l.ab3, 1, HDInterval::hdi, credMass=.8)[2,] |> exp()
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
#+  scale_color_manual(values = cols)

# ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\abundance_cycles.tiff", 
#        plot=ptab2, 
#        device="tiff",
#        width=6, 
#        height=4,
#        dpi=300)

# summaries
# which year had max median catches
df.ab2$year[which.max(df.ab2$md)]
df.ab2[df.ab2$year>=2004 & df.ab2$year<=2005,]
df.ab2[df.ab2$year>=2022 & df.ab2$year<=2023,]

#################
# Plot captures per mist nest hours 
# for Whitefish Point site
range(d2$dat$year)
yr <- seq(1994, 2023, by=0.5)
cycle.sc2 <- 1:length(yr)/8
yr.sc <- (yr-2008.5) / ((2023-1994)/2)
sf <- ifelse(yr %% 1 ==0, 2, 1) # 2=spring, 1=fall

l.ab2 <- l.ab3 <- array(NA, dim=c(length(yr), 4000), 
                        dimnames=list(year=yr, iter=c(1:4000)))
for (t in 1:length(yr)){
  l.ab2[t, ] <- ps2$mu[sf[t],] + 
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
# calculate a four year moving average
for (t in 1:nrow(l.ab2)){
  for (i in 1:ncol(l.ab2)){
    ind1 <- ifelse( (t-16)<1, 1, t-16)
    ind2 <- ifelse( (t+16)>nrow(l.ab2), nrow(l.ab2), t+16)
    l.ab3[t,i] <- mean(l.ab2[ind1:ind2, i])
  }}

df.ab3 <- data.frame( year = yr,
                      md = apply(l.ab3, 1, median, na.rm=T) |> exp(),
                      lci95 = apply(l.ab3, 1, HDInterval::hdi)[1,]  |> exp(),
                      uci95 = apply(l.ab3, 1, HDInterval::hdi)[2,] |> exp(),
                      lci80 = apply(l.ab3, 1, HDInterval::hdi, credMass=.8)[1,] |> exp(),
                      uci80 = apply(l.ab3, 1, HDInterval::hdi, credMass=.8)[2,] |> exp()
)

df.dat <- data.frame(count=d2$datl$y/d2$constl$effort,
                     year=d2$dat$year2)

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
  coord_cartesian(xlim=c(1994, 2023), ylim=c(0, 150)) +
  ggtitle("Whitefish Point") +
#+  scale_color_manual(values = cols)

ap23 <- align_plots(p2, p3, align="v", axis="l")
p23 <- plot_grid(ap23[[1]], ap23[[2]], nrow = 2, align="v", 
                 label_y = "Captures")
## ---- lambdaplot --------
p23

## ---- saveplot12 --------
ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\populationtrends.tiff",
       plot=p23, device="tiff", width=6, height=8, dpi=300)


# summaries
# which year had max median catches
df.ab2$year[which.max(df.ab2$md)]
#**********************
#* Plot population growth over time
#**********************
plam <- array(NA, dim=c(d1$constl$ntime-1, 4000), 
              dimnames=list(year=c(1997:2023), iter=c(1:4000)))
lst2 <- ifelse(lst$lambda.star==0, 0.1, lst$lambda.star )
for (t in 1:(d1$constl$ntime-1)){
  plam[t,] <- lst$lambda.star[t+1,]/lst2[t,]
}

mlam <- plam %>%
  as.tbl_cube(met_name = "lam") %>%
  as_tibble

df.lam <- data.frame( year = 1997:2023,
                      md = apply(plam, 1, median, na.rm=T),
                      lci95 = apply(plam, 1, hdi)[1,],
                      uci95 = apply(plam, 1, hdi)[2,],
                      lci80 = apply(plam, 1, hdi, credMass=.8)[1,],
                      uci80 = apply(plam, 1, hdi, credMass=.8)[2,]
)

col <- "black"

plam <- ggplot() + theme_minimal() + 
  geom_line(data=mlam, aes(x=year, y=lam, group=iter), 
            color="gray40", linewidth=0.5, alpha=0.05 ) +
  xlim(1997, 2023) + ylim(0,2) +
  geom_hline(yintercept=1, lwd=2, color="black", linetype="dashed") +
  geom_line(data=df.lam, aes(x=year, y=lci95), color=col, linewidth=0.5 ) +
  geom_line(data=df.lam, aes(x=year, y=uci95), color=col, linewidth=0.5) +
  geom_line(data=df.lam, aes(x=year, y=lci80), color=col, linewidth=1 ) +
  geom_line(data=df.lam, aes(x=year, y=uci80), color=col, linewidth=1 ) +
  geom_line(data=df.lam, aes(x=year, y=md), color=col, linewidth=2 ) +
  xlab("Year") +
  ylab( expression(paste("Population growth (", lambda,")")) )+
  annotate(geom = "text", x = 1990, y = 1.2, label = "B", size=6)

ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\lambda.tiff", 
       plot=plam, 
       device="tiff",
       width=6, 
       height=4,
       dpi=300)