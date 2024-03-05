library("AICcmodavg")
require("pals")


load(file = "fitted_a_values.RData")
mg$pop_density <- mg$POPESTIMATE2010/(as.numeric(as.character(mg$ALAND))/1000000)
mg$lgpop <- log(mg$pop_density)
mg$lgpop2 <- log(mg$pop_density)^2
mg$lgpop3 <- log(mg$pop_density)^3
mg$lgpop4 <- log(mg$pop_density)^4

#linear model
lm1 <- lm(ause ~ log(pop_density), data = mg)
summary(lm1)

#quadratic model
lm2 <- lm(ause ~ lgpop + lgpop2, data = mg)
summary(lm2)

#3rd order polynomial model
lm3 <- lm(ause ~ lgpop + lgpop2 + lgpop3 , data = mg)
summary(lm3)

#4th order polynomial model
lm4 <- lm(ause ~  lgpop + lgpop2 + lgpop3 + lgpop4, data = mg)
summary(lm4)

models <- list(lm1, lm2,lm3, lm4)

#specify model names
mod.names <- c('log1', 'log2','log3', 'log4')

#calculate AIC of each model
aictab(cand.set = models, modnames = mod.names)
