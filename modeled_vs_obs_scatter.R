require("pals")

load(file ="observed_intensity.RData") 
load(file = "modeled_intensity.RData")

mg <- merge(mg, shpdt, by=c("STATEFP","COUNTYFP"))

pdf("obs_mod_fitting_poly.pdf", height = 2.75, width = 3)
par(mar=c(3,3,1,1))
coluse <- brewer.spectral(5)[5]
mg$col <- brewer.spectral(9)[as.numeric(cut(mg$average_historic_humidity, breaks = c(seq(0,0.016,0.002), 0.45)))]
mg$col <- brewer.spectral(9)[as.numeric(cut(mg$flu.Longitude, breaks = 9))]
plot(mg$intensity, mg$intensity_obs, pch = 16, ylim = c(0,0.31), col = mg$col, bty = "n",xlab = "",ylab = "")
title(xlab = "Modeled intensity", line = 2)
title(ylab = "Observed intensity", line = 2)

statefipsuse <- "36"
countyfipsuse <- "061"
points(mg$intensity[mg$STATEFP==statefipsuse & mg$COUNTYFP== countyfipsuse],
       mg$intensity_obs[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], 
       bg = mg$col[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], pch = 21, cex = 1.5)
text(mg$intensity[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse],
     mg$intensity_obs[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], 
     label = "NY", pos = 1, cex = 0.7)

statefipsuse <- "06"
countyfipsuse <- "075"
points(mg$intensity[mg$STATEFP==statefipsuse & mg$COUNTYFP== countyfipsuse],
       mg$intensity_obs[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], 
       bg = mg$col[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], pch = 21, cex = 1.5)
text(mg$intensity[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse],
     mg$intensity_obs[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], 
     label = "SF", pos = 1, cex = 0.7)

statefipsuse <- "56"
countyfipsuse <- "021"
points(mg$intensity[mg$STATEFP==statefipsuse & mg$COUNTYFP== countyfipsuse],
       mg$intensity_obs[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], 
       bg = mg$col[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], pch = 21, cex = 1.5)
text(mg$intensity[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse],
     mg$intensity_obs[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], 
     label = "Cheyenne", pos = 3, cex = 0.7)

statefipsuse <- "44"
countyfipsuse <- "005"
points(mg$intensity[mg$STATEFP==statefipsuse & mg$COUNTYFP== countyfipsuse],
       mg$intensity_obs[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], 
       bg = mg$col[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], pch = 21, cex = 1.5)
text(mg$intensity[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse],
     mg$intensity_obs[mg$STATEFP==statefipsuse & mg$COUNTYFP==countyfipsuse], 
     label = "Newport", pos = 3, cex = 0.7)


lm.out <- lm(intensity_obs ~ intensity, data= mg)
intensity = seq(0,0.5,by = 0.05)
conf_interval <- predict(lm.out, newdata=data.frame(intensity=intensity), interval="confidence")
lines(intensity, conf_interval[,2], col=coluse, lty=2, lwd = 2)
lines(intensity, conf_interval[,3], col=coluse, lty=2, lwd = 2)
#polygon(c(intensity, rev(intensity)), c( conf_interval[,2], rev(conf_interval[,3])), border = NA, col = "lightblue")
abline(lm.out, col=coluse, lwd= 2, lty = 1)
dev.off()
