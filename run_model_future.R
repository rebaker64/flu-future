require("deSolve")
require("rgdal")
require("pals")
require("readxl")
#set parameters for model
Lset = 35
R0maxset = 2.2
R0minset = 1.2

#Load shape file
# county shape file from census bureau
usashp <- readOGR("cb_2016_us_county_5m.shp")
shpdt <- usashp@data

# load historic climate data
load("historic_us_average.RData")

# load disease model
source("flu_model.R", encoding = 'UTF-8')

list_mods <- c("ACC-CM2","ACC-ESM","AWI-CM","BCC-CSM","CanESMO","CanESM5","CAS-ESM2","CESM2","CIESM","CMC-CM2","CMC-ESM",
               "CNR-CMH","CNR-CM6","CNRM-ES","E3SM-0","E3SM-E","E3SM-1","EC-CC","EC-LR","EC-Veg","FGOAL-f3",
               "FGOAL-g3","GFDL-CM4","GFDL-ESM","GISS-G","GISS-H","GISS-2","Had-LL","Had-MM","IITM-ESM","INM-CM4",
               "INM-CM5","IPSL-CM6","KACE-1","MCM-UA","MIROC-E","MIROC6","MPI-ESM1","MRI-ESM2","NorE-LM","NorE-MM","TaiESM1","UKESM1")

# pop projections
# load population projections from CIESIN
pop_proj <- read_excel("hauer_county_totpop_SSPs.xlsx")
pop_proj$COUNTYFP10[pop_proj$STATEFP10=="02" & pop_proj$COUNTYFP10 == "270"] <- "158" # county fips changed
pop_proj$COUNTYFP10[pop_proj$STATEFP10=="46" & pop_proj$COUNTYFP10 == "113"] <- "102" # county fips changed

pop_proj_locations <- unique(pop_proj$STATEFP10)

# pop historic
# load historic population from the US Census Bureau
pop_obs <- read.csv("co-est2019-alldata.csv")
pop_obs <- pop_obs[c("STATE","COUNTY","POPESTIMATE2010")]
pop_locations <- unique(pop_obs$STATE)

for(i in 1:43){ # model index
setwd("CMIP6_output")
# load specific humidity data from CMIP6 model is list_mods
load(paste0(list_mods[i],".RData"))

projected_us_average <- historic_us_average + mod_out

# set up model
xstart = c(S = 0.7, I = 0.1, R = 0.2)
times = seq(1, 5000, by = 1)
yearst <- seq(1,1000,by = 1/52)[1:length(times)]

#metrics to fill in for each county
shpdt <- usashp@data
shpdt$statefp_num <- as.numeric(as.character(shpdt$STATEFP))
shpdt$countyfp_num <- as.numeric(as.character(shpdt$COUNTYFP))
shpdt$peak_size_change <- NA
shpdt$intensity_change <- NA
shpdt$R0_mean_change <- NA
shpdt$R0_range_change <- NA
shpdt$q_mean_change <- NA
shpdt$q_range_change <- NA

for(j in 1:3233){ # county index
  if(shpdt$statefp_num[j]  %in% pop_locations){
    
# get pop data
pop_proj_sub <- pop_proj[pop_proj$STATEFP10==shpdt$STATEFP[j] & pop_proj$COUNTYFP10 == shpdt$COUNTYFP[j],]
pop_proj_select <- pop_proj_sub$ssp52090/(as.numeric(as.character(shpdt$ALAND[j]))/1000000)

if(pop_proj_select != 0){ # some population projections are zero - remove these as can't run model in this case
pop_obs_sub <-  pop_obs$POPESTIMATE2010[pop_obs$STATE == shpdt$statefp_num[j] & 
                                          pop_obs$COUNTY == shpdt$countyfp_num[j]]/ (as.numeric(as.character(shpdt$ALAND[j]))/1000000)

# calculate "a" parameter using population data
ause_hist <- -286.0764942 + 77.6571653*log(pop_obs_sub)  + -13.5141550*(log(pop_obs_sub))^2 + 0.7758193*(log(pop_obs_sub))^3 
ause_proj <- -286.0764942 + 77.6571653*log(pop_proj_select)  + -13.5141550*(log(pop_proj_select))^2 + 0.7758193*(log(pop_proj_select))^3 

# run model for historic data
qhist <- rep(historic_us_average[j,], length = length(times))
paras <- list(D = 4/7, L = Lset, R0max = R0maxset, R0min = R0minset, qList = qhist, ause = ause_hist)
out_hist <- as.data.frame(ode(xstart, times, sirs_ode, paras))
out_hist <- out_hist[yearst >= 50 & yearst < 55,] # 50 year burn-in period
plot(yearst[1:260],out_hist$I,type="l", main = j, ylim = c(0, 0.05))

# run model for future projections, baseline pop
qproj <- rep(projected_us_average[j,], length = length(times))
paras <- list(D = 4/7, L = Lset, R0max = R0maxset, R0min = R0minset, qList = qproj, ause = ause_hist)
out_proj_int <- as.data.frame(ode(xstart, times, sirs_ode, paras))
out_proj_int <- out_proj_int[yearst >= 50 & yearst < 55,] # 50 year burn-in period
lines(yearst[1:260],out_proj_int$I,type="l", col="orange")

# run model for future projections
qproj <- rep(projected_us_average[j,], length = length(times))
paras <- list(D = 4/7, L = Lset, R0max = R0maxset, R0min = R0minset, qList = qproj, ause = ause_proj)
out_proj <- as.data.frame(ode(xstart, times, sirs_ode, paras))
out_proj <- out_proj[yearst >= 50 & yearst < 55,] # 50 year burn-in period
lines(yearst[1:260],out_proj$I,type="l", col="red")

# calculate change in peak size
peak_size_change <- max(out_proj$I) - max(out_hist$I)

# calculate change in entropy/intensity
ps <- out_proj$I[1:52]/sum(out_proj$I[1:52])
ps <- ps[ps != 0]
entropy_proj <- 1- (-1*sum(ps*log(ps)))/(-52*sum(1/52*log(1/52)))

ps <- out_hist$I[1:52]/sum(out_hist$I[1:52])
ps <- ps[ps != 0]
entropy_hist <- 1- (-1*sum(ps*log(ps)))/(-52*sum(1/52*log(1/52)))

intensity_change <- entropy_proj - entropy_hist

# R0 calculation
R0_hist = exp(ause_hist*historic_us_average[j,] + log(R0maxset - R0minset)) + R0minset
R0_proj_int = exp(ause_hist*projected_us_average[j,] + log(R0maxset - R0minset)) + R0minset
R0_proj = exp(ause_proj*projected_us_average[j,] + log(R0maxset - R0minset)) + R0minset
R0_mean_change = mean(R0_proj) - mean(R0_hist)
R0_range_change = diff(range(R0_proj)) -diff(range(R0_hist))

# q change calculation
q_mean_change = mean(projected_us_average[j,]) - mean(historic_us_average[j,])
q_range_change = diff(range(projected_us_average[j,])) -diff(range(historic_us_average[j,]))

# update values for county
shpdt$peak_size_historic[j] <- max(out_hist$I)
shpdt$peak_size_projint[j] <-  max(out_proj_int$I)
shpdt$peak_size_proj[j] <-  max(out_proj$I)
shpdt$peak_size_change[j] <- peak_size_change

shpdt$intensity_historic[j] <- entropy_hist
shpdt$intensity_proj_int[j] <-entropy_proj_int
shpdt$intensity_proj[j] <- entropy_proj
shpdt$intensity_change[j] <- intensity_change

shpdt$R0_historic[j] <- R0_hist
shpdt$R0_proj_int[j] <- R0_proj_int
shpdt$R0_proj[j] <- R0_proj
shpdt$R0_mean_change[j] <- R0_mean_change
shpdt$R0_range_change[j] <- R0_range_change
shpdt$q_hist[j] <-mean(historic_us_average[j,])
shpdt$q_proj[j] <- mean(projected_us_average[j,])
shpdt$q_mean_change[j] <- q_mean_change
shpdt$q_range_change[j] <- q_range_change
print(paste0("county = ",j," and ","model = ",i))
}
}
}

setwd("~/output_metrics_pop_projections")
save(shpdt, file = paste0(list_mods[i],".RData"))
}


# example plotting 
shpdt <- as.data.frame(shpdt)
usashp <- merge(usashp, shpdt, by.x=c(names(shpdt)[1:9]), by.y = c(names(shpdt)[1:9]))
usashp <- usashp[!usashp$STATEFP %in% c("02","15","60","66","69","78","72"),]
# make plots
spplot(usashp, "peak_size_change",  sub = " ", 
       col = "transparent", col.regions=rev(brewer.spectral(31)), bty="n", 
       par.settings = list(axis.line = list(col =  'transparent')), main = list_mods[i])

