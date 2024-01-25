library(MCMCvis)
library (coda)
library (ggplot2)
library (reshape2)
library (tidybayes)
library(tidyr)
library(cubelyr)
library(HDInterval)
library(cowplot)
options(scipen=999)
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\nb_lambdacycles.Rdata")
dat <- read.csv("data\\data.csv")
datl <- list(y = dat$caps)
constl <- list(time = (dat$year - median(dat$year))/ max(dat$year - median(dat$year)),
               effort = dat$effort/500,
               ntime = length(dat$year),
               time2 = 1:length(1996:2023)/4)

pars <-c( "mu", "beta", "r", "delta")
MCMCtrace(nbzip, pars, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
MCMCsummary(nbzip, pars, 
            HPD=TRUE,  hpd_prob = 0.95,
            round=3, pg0 = TRUE, func = median)
MCMCsummary(nbzip, pars, 
            HPD=TRUE,  hpd_prob = 0.85,
            round=3, pg0 = TRUE, func = median)

p <- MCMCpstr(nbzip, pars, type="chains") 
lst <- MCMCpstr(nbzip, pars, type="chains") 

#**********************
#* plot abundance over time
#**********************
yr <- 1996:2023
yr.sc <- (yr-2009.5)/13.5

pab <- array(NA, dim=c(constl$ntime, 4000), 
             dimnames=list(year=c(1996:2023), iter=c(1:4000)))
for (t in 1:(constl$ntime)){
  pab[t, ] <- lst$mu + lst$beta[1,]*yr.sc[t] #+ lst$beta[2,]*yr.sc[t]^2
}
mab <- pab %>%
  as.tbl_cube(met_name = "l.ab") %>%
  as_tibble

mab$ab <- exp(mab$l.ab)

df.ab <- data.frame( year = yr,
                     md = apply(pab, 1, median, na.rm=T) |> exp(),
                     lci95 = apply(pab, 1, hdi)[1,]  |> exp(),
                     uci95 = apply(pab, 1, hdi)[2,] |> exp(),
                     lci80 = apply(pab, 1, hdi, credMass=.8)[1,] |> exp(),
                     uci80 = apply(pab, 1, hdi, credMass=.8)[2,] |> exp()
)

col <- "black"

ptab <- ggplot() + theme_minimal() + 
  geom_line(data=mab, aes(x=year, y=ab, group=iter), 
            color="gray40", linewidth=0.5, alpha=0.05 ) +
  xlim(1996, 2023) + 
  xlab("Year") + ylab("Abundance") +
  geom_line(data=df.ab, aes(x=year, y=lci95), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab, aes(x=year, y=uci95), color=col, linewidth=0.5) +
  geom_line(data=df.ab, aes(x=year, y=lci80), color=col, linewidth=1 ) +
  geom_line(data=df.ab, aes(x=year, y=uci80), color=col, linewidth=1 ) +
  geom_line(data=df.ab, aes(x=year, y=md), color=col, linewidth=2 ) +
  coord_cartesian(xlim=c(1996, 2023), ylim=c(0, 200)) 

ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\abundance_year.tiff", 
       plot=ptab, 
       device="tiff",
       width=6, 
       height=4,
       dpi=300)

#**********************
#* Plot population growth over time
#**********************
plam <- array(NA, dim=c(constl$ntime-1, 4000), 
                    dimnames=list(year=c(1997:2023), iter=c(1:4000)))
lst2 <- ifelse(lst$lambda.star==0, 0.1, lst$lambda.star )
for (t in 1:(constl$ntime-1)){
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

ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\popgrowth_year.tiff", 
       plot=plam, 
       device="tiff",
       width=6, 
       height=4,
       dpi=300)

#**********************
#* plot population cycles over time
#**********************
yr <- seq(1996, 2023, by=0.25)
yr.sc2 <- (yr-1995)/4
yr.sc <- (yr-median(dat$year)) / max(dat$year - median(dat$year))
  
l.ab2 <- l.ab3 <- array(NA, dim=c(length(yr), 4000), 
             dimnames=list(year=yr, iter=c(1:4000)))
for (t in 1:length(yr)){
  l.ab2[t, ] <- lst$mu + 
                lst$beta[1,]*yr.sc[t] + lst$beta[2,]*yr.sc[t]^2 +
                lst$delta[1,]*cos(2*pi*yr.sc2[t]) + lst$delta[2,]*sin(2*pi*yr.sc2[t]) +
                lst$delta[3,]*cos(2*pi*yr.sc2[t])*yr.sc[t] + lst$delta[4,]*sin(2*pi*yr.sc2[t])*yr.sc[t] +
                lst$delta[5,]*cos(2*pi*yr.sc2[t])*yr.sc[t]^2 + lst$delta[6,]*sin(2*pi*yr.sc2[t])*yr.sc[t]^2
  
  l.ab3[t, ] <- lst$mu + 
    lst$beta[1,]*yr.sc[t] + lst$beta[2,]*yr.sc[t]^2 +
    lst$delta[1,]*cos(2*pi*median(yr.sc2)) + lst$delta[2,]*sin(2*pi*median(yr.sc2)) +
    lst$delta[3,]*cos(2*pi*median(yr.sc2))*median(yr.sc[t]) + lst$delta[4,]*sin(2*pi*median(yr.sc2))*median(yr.sc[t]) +
    lst$delta[5,]*cos(2*pi*median(yr.sc2))*yr.sc[t]^2 + lst$delta[6,]*sin(2*pi*median(yr.sc2))*yr.sc[t]^2
}

mab2 <- l.ab2 %>%
  as.tbl_cube(met_name = "l.ab") %>%
  as_tibble

mab3 <- l.ab3 %>%
  as.tbl_cube(met_name = "l.ab") %>%
  as_tibble

df.ab2 <- data.frame( year = yr,
                     md = apply(l.ab2, 1, median, na.rm=T) |> exp(),
                     lci95 = apply(l.ab2, 1, hdi)[1,]  |> exp(),
                     uci95 = apply(l.ab2, 1, hdi)[2,] |> exp(),
                     lci80 = apply(l.ab2, 1, hdi, credMass=.8)[1,] |> exp(),
                     uci80 = apply(l.ab2, 1, hdi, credMass=.8)[2,] |> exp()
)

df.ab3 <- data.frame( year = yr,
                      md = apply(l.ab3, 1, median, na.rm=T) |> exp(),
                      lci95 = apply(l.ab3, 1, hdi)[1,]  |> exp(),
                      uci95 = apply(l.ab3, 1, hdi)[2,] |> exp(),
                      lci80 = apply(l.ab3, 1, hdi, credMass=.8)[1,] |> exp(),
                      uci80 = apply(l.ab3, 1, hdi, credMass=.8)[2,] |> exp()
)

df.dat <- data.frame(count=datl$y/constl$effort,
                     year=1996:2023)
  
col <- "black"


acf(datl$y/constl$effort, 
    main="Autocorrelation of counts per 500 net hours")

ptab2 <- ggplot() + theme_minimal() + 
  geom_line(data=mab2, aes(x=year, y=exp(l.ab), group=iter), 
            color="gray40", linewidth=0.5, alpha=0.05 ) +
  xlab("Year") + ylab("Abundance per 500 net hours") +
  geom_line(data=df.ab2, aes(x=year, y=lci95), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci95), color=col, linewidth=0.5) +
  geom_line(data=df.ab2, aes(x=year, y=lci80), color=col, linewidth=1 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci80), color=col, linewidth=1 ) +
  geom_line(data=df.ab2, aes(x=year, y=md), color=col, linewidth=2 ) +
  geom_line(data=df.dat, aes(x=year, y=count), color="blue", linewidth=1, alpha=0.75) +
  coord_cartesian(xlim=c(1996, 2023), ylim=c(0, 350)) 

ptab3 <- ggplot() + theme_minimal() + 
  geom_line(data=mab3, aes(x=year, y=exp(l.ab), group=iter), 
            color="gray40", linewidth=0.5, alpha=0.05 ) +
  xlab("Year") + ylab("Abundance per 500 net hours") +
  geom_line(data=df.ab3, aes(x=year, y=lci95), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab3, aes(x=year, y=uci95), color=col, linewidth=0.5) +
  geom_line(data=df.ab3, aes(x=year, y=lci80), color=col, linewidth=1 ) +
  geom_line(data=df.ab3, aes(x=year, y=uci80), color=col, linewidth=1 ) +
  geom_line(data=df.ab3, aes(x=year, y=md), color=col, linewidth=2 ) +
  #geom_line(data=df.dat, aes(x=year, y=count), color="blue", linewidth=1, alpha=0.75) +
  coord_cartesian(xlim=c(1996, 2023), ylim=c(0, 30))

#align_plots(ptab2, ptab3, align="v", axis="l")
both_aligned <- plot_grid(ptab2, ptab3, nrow = 2, align="v")

ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\abundance_cycles.tiff", 
       plot=both_aligned, 
       device="tiff",
       width=6, 
       height=8,
       dpi=300)
