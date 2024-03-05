require("raster")
require("rgdal")
require("ggplot2")
theme_set(theme_bw())
library("sf")
require("pals")

#Load shape file
# get county shapefile from census bureau
usashp <- readOGR("cb_2016_us_county_5m.shp")

list_mods <- c("ACC-CM2","ACC-ESM","AWI-CM","BCC-CSM","CanESMO","CanESM5","CAS-ESM2","CESM2","CIESM","CMC-CM2","CMC-ESM",
               "CNR-CMH","CNR-CM6","CNRM-ES","E3SM-0","E3SM-E","E3SM-1","EC-CC","EC-LR","EC-Veg","FGOAL-f3",
               "FGOAL-g3","GFDL-CM4","GFDL-ESM","GISS-G","GISS-H","GISS-2","Had-LL","Had-MM","IITM-ESM","INM-CM4",
               "INM-CM5","IPSL-CM6","KACE-1","MCM-UA","MIROC-E","MIROC6","MPI-ESM1","MRI-ESM2","NorE-LM","NorE-MM","TaiESM1","UKESM1")

# get distribution of 43 projected peak incidence, intensity, q range and r0 range values
gen_pk_dist <- NULL
for(i in 1:length(list_mods)){
setwd("output_metrics_pop_projections")
load(file = paste0(list_mods[i],".RData"))

get_var <- shpdt$peak_size_change 
gen_pk_dist <- cbind(gen_pk_dist, get_var)

}

shpdt$peak_size_50 <-apply(gen_pk_dist, 1, function(x) quantile(x,probs=0.5, na.rm = T)) 
shpdt$peak_size_20 <-apply(gen_pk_dist, 1, function(x) quantile(x,probs=0.2, na.rm = T)) 
shpdt$peak_size_80 <-apply(gen_pk_dist, 1, function(x) quantile(x,probs=0.8, na.rm = T)) 

# example plotting for later
shpdt <- as.data.frame(shpdt)
usashp <- merge(usashp, shpdt , by.x=c(names(shpdt)[1:9]), by.y = c(names(shpdt)[1:9]))
usashp <- usashp[!usashp$STATEFP %in% c("02","15","60","66","69","78","72"),]


# make plots

usashp$peak_size_80[usashp$peak_size_80 >= 0.34] <- 0.34 # limit projection range for visualization

pdf("peak_size_80.pdf")
spplot(usashp, "peak_size_80",  sub = " ", 
       col = "transparent", col.regions=rev(brewer.rdbu(36)),at = seq(-0.034,0.034,0.002), main = "80th percentile change", 
       par.settings = list(axis.line = list(col =  'transparent')))
dev.off()

pdf("peak_size_50.pdf")
spplot(usashp, "peak_size_50",  sub = " ", 
       col = "transparent", col.regions=rev(brewer.rdbu(36)), at = seq(-0.034,0.034,0.002), main = "50th percentile change", 
       par.settings = list(axis.line = list(col =  'transparent')))
dev.off()

usashp$peak_size_80[usashp$peak_size_20 <= -0.34] <- -0.34 # limit projection range for visualization
pdf("peak_size_20.pdf")
spplot(usashp, "peak_size_20",  sub = " ", 
       col = "transparent", col.regions=rev(brewer.rdbu(36)),at = seq(-0.034,0.034,0.002), main = "20th percentile change", 
       par.settings = list(axis.line = list(col =  'transparent')))
dev.off()

