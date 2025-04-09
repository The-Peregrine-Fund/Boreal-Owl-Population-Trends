
# Tadoussac Bird Observatory (OOT), located in Tadoussac, Québec, Canada 
dat <- read.csv("data\\data_T.csv")
datl <- list(y = dat$caps)
constl <- list(time = (dat$year - median(dat$year))/ max(dat$year - median(dat$year)),
               effort = dat$effort/500,
               ntime = length(dat$year),
               time2 = 1:length(1996:2023)/4)
save(dat=dat, datl=datl, constl=constl, file="data//data_s1.Rdata")

# Whitefish Point Bird Observatory (WPBO), located at Whitefish Point, Chippewa County, MI
dat <- read.csv("data\\data_WP.csv")

dat2 <- data.frame(caps=dat$MCFallCapture,
                   effort= dat$MCFallEffort,
                   season= rep("fall", length(dat$MCFallEffort)),
                   year=dat$Year)
dat2 <- dat2[!is.na(dat2$effort), ]
dat2 <- dat2[order(dat2$year),]

datl <- list(y = dat2$caps)
constl <- list(time = (dat2$year - median(1995:2022))/ max(dat2$year - median(1995:2022)),
               effort = dat2$effort/500,
               ntime = length(dat2$year) )

constl$time2 <- 1:length(1995:2022)/4
timemerge <- data.frame(year = seq(1,length(1995:2022), by=1)+1994,  
                        year1995 = seq(1,length(1995:2022), by=1),
                        year4 = seq(1,length(1995:2022), by=1)/4
)
dat3 <- merge(dat2, timemerge, by="year")
constl$time2 <- dat3$year4
cbind(constl$time, constl$time2)
# Insert NAs to years without surveys
datl$y[constl$effort==0] <- NA # no effort so missing value
constl$effort[constl$effort==0] <- mean(constl$effort[constl$effort!=0])

save(dat=dat3, datl=datl, constl=constl, file="data//data_s2.Rdata")