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

dat <- read.csv("data\\data_T.csv")
datl1 <- list(y = dat$caps)
constl1 <- list(time = (dat$year - median(dat$year))/ max(dat$year - median(dat$year)),
               effort = dat$effort/500,
               ntime = length(dat$year),
               time2 = 1:length(1996:2023)/4)

dat <- read.csv("data\\data_MC.csv")

dat2 <- data.frame(caps=c(dat$MCCaptureSpring, dat$MCFallCapture),
                   effort= c(dat$MCSpringEffort, dat$MCFallEffort),
                   season= c(rep("spring", length(dat$MCSpringEffort)),
                             rep("fall", length(dat$MCFallEffort)) ),
                   year=c(dat$Year, dat$Year)
)
dat2$year2 <- ifelse(dat2$season=="spring", dat2$year, dat2$year+0.5)
dat2 <- dat2[ !is.na(dat2$effort) & dat2$effort>0,]
dat2 <- dat2[order(dat2$year2),]

datl2 <- list(y = dat2$caps)
constl2 <- list(time = (dat2$year2 - median(1994:2023))/ max(dat2$year2 - median(1994:2023)),
               effort = dat2$effort/500,
               season = as.numeric(factor(dat2$season, 
                                          levels=c("fall", "spring"), 
                                          ordered=TRUE)),
               ntime = length(dat2$year2) )

constl2$time2 <- 1:length(1994:2023)/4
timemerge <- data.frame(year = seq(1,length(1994:2023), by=0.5)+1993,  
                        year1994 = seq(1,length(1994:2023), by=0.5),
                        year4 = seq(1,length(1994:2023), by=0.5)/4
)
dat3 <- merge(dat2, timemerge, by="year")
constl2$time2 <- dat3$year4
cbind(constl2$time, constl2$time2)

plot(dat2$year2, datl2$y/constl$effort, type="l")
acf(datl2$y/constl2$effort)
acf(datl2$y/constl2$effort, type="p")

# Look at global models for significance 
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-set_population-cycles.Rdata")
m1 <- m
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\outputs\\full-model-set_population-cycles-site2.Rdata")
m2 <- m
# inspect estimates from global models
pars <-c( "mu", "delta", "sigma")
MCMCsummary(m1[[1]]$samples, pars, 
            HPD=TRUE,  hpd_prob = 0.85,
            round=3, pg0 = TRUE, func = median)

MCMCsummary(m2[[8]]$samples, pars, 
            HPD=TRUE,  hpd_prob = 0.85,
            round=3, pg0 = TRUE, func = median)




#################
# Older stuff very useful though
pars <-c( "mu", "delta", "sigma")
MCMCtrace(m[[2]]$samples, pars, pdf=F, 
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)

coeftab <- MCMCsummary(m[[2]]$samples, pars, 
            HPD=TRUE,  hpd_prob = 0.95,
            round=2, pg0 = TRUE, func = median, 
            func_name = "Median")
coeftab <- coeftab[, c(8,3,4,5,7)]
write.csv(coeftab, 'docs\\coeftab.csv')

MCMCsummary(m[[2]]$samples, pars, 
            HPD=TRUE,  hpd_prob = 0.85,
            round=3, pg0 = TRUE, func = median)

p <- MCMCpstr(m[[2]]$samples, pars, type="chains") 
lst <- MCMCpstr(m[[2]]$samples, pars, type="chains") 

acf(datl$y/constl$effort)
acf(datl$y/constl$effort, type="p")


#**********************
#* PLot autocorrelation
#**********************
tiff("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\autocorrelation.tiff", 
     width=8, height=4, 
     units= "in", res=300)
par(mfrow=c(1,2))
acf(datl$y/constl$effort, main="")
acf(datl$y/constl$effort, main="", type="p")
dev.off()

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
                lst$delta[1,]*yr.sc[t] + 
                lst$delta[2,]*yr.sc[t]^2 +
                lst$delta[3,]*cos(2*pi*yr.sc2[t]) + 
                lst$delta[4,]*sin(2*pi*yr.sc2[t]) +
                lst$delta[5,]*cos(2*pi*yr.sc2[t])*yr.sc[t] + 
                lst$delta[6,]*sin(2*pi*yr.sc2[t])*yr.sc[t] 

}

mab2 <- l.ab2 %>%
  as.tbl_cube(met_name = "l.ab") %>%
  as_tibble

df.ab2 <- data.frame( year = yr,
                     md = apply(l.ab2, 1, median, na.rm=T) |> exp(),
                     lci95 = apply(l.ab2, 1, hdi)[1,]  |> exp(),
                     uci95 = apply(l.ab2, 1, hdi)[2,] |> exp(),
                     lci80 = apply(l.ab2, 1, hdi, credMass=.8)[1,] |> exp(),
                     uci80 = apply(l.ab2, 1, hdi, credMass=.8)[2,] |> exp()
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
                      lci95 = apply(l.ab3, 1, hdi)[1,]  |> exp(),
                      uci95 = apply(l.ab3, 1, hdi)[2,] |> exp(),
                      lci80 = apply(l.ab3, 1, hdi, credMass=.8)[1,] |> exp(),
                      uci80 = apply(l.ab3, 1, hdi, credMass=.8)[2,] |> exp()
)

df.dat <- data.frame(count=datl$y/constl$effort,
                     year=1996:2023)
  
col <- "black"
cols <- c("Data"="blue", "Median"="black", 
          "85% HDIs"= "black", "95% HDIs"="black")
szs <- c("Data"= 1, "Median"=1, 
         "85% HDIs"= 0.5, "95% HDIs"=0.25)
ptab2 <- ggplot() + theme_minimal() + 
  geom_line(data=mab2, aes(x=year, y=exp(l.ab), group=iter), 
            color="gray40", linewidth=0.25, alpha=0.05 ) +
  #geom_line(data=df.ab3, aes(x=year, y=md), color=col, linewidth=1, linetype="dashed" ) +
  xlab("Year") + ylab("Captures of Boreal Owl\nper 500 mist net hours") +
  geom_line(data=df.ab2, aes(x=year, y=lci95), color=col, linewidth=0.25 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci95), color=col, linewidth=0.25) +
  geom_line(data=df.ab2, aes(x=year, y=lci80), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab2, aes(x=year, y=uci80), color=col, linewidth=0.5 ) +
  geom_line(data=df.ab2, aes(x=year, y=md), color=col, linewidth=1 ) +
  geom_line(data=df.dat, aes(x=year, y=count), color="blue", linewidth=1, alpha=0.75) +
  coord_cartesian(xlim=c(1996, 2023), ylim=c(0, 350)) 
#+  scale_color_manual(values = cols)

ggsave(filename="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Boreal Owl pop trends\\docs\\figs\\abundance_cycles.tiff", 
       plot=ptab2, 
       device="tiff",
       width=6, 
       height=4,
       dpi=300)

# summaries
# which year had max median catches
df.ab2$year[which.max(df.ab2$md)]

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