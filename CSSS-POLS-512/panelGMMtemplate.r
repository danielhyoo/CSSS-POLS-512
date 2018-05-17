# Dynamic Panel Data Models (Arellano-Bond, etc.)
# for fixed effects in short-T panels
# Model Fitting and Interpretation via Simulation
#
# Christopher Adolph   faculty.washington.edu/cadolph
# 20 May 2015

# Clear memory
rm(list=ls())

# Load libraries
library(plm)            # Econometrics package for linear panel models
library(tseries)        # For ADF unit root test
library(simcf)          # For panel functions and simulators
library(tile)           # For visualization of model inference
library(RColorBrewer)   # For nice colors
library(MASS)           # For mvrnorm()
source("helperCigs.R")  # For graphics functions

# Load cigarette consumption data (Jonathan Gruber, MIT)
# Variables (see codebook):
# state	year	cpi	pop	packpc	income	tax	avgprs	taxs
data <- read.csv("cigarette.csv")

# Quick inflation adjustment to 1995 dollars
inflAdjust <- function(x,cpi,year,target) {
    unique(cpi[year==target])*x/cpi
}
    
data$income95 <- with(data, inflAdjust(income, cpi, year, 1995))
data$tax95 <- with(data, inflAdjust(tax, cpi, year, 1995))
data$avgprs95 <- with(data, inflAdjust(avgprs, cpi, year, 1995))
data$taxs95 <- with(data, inflAdjust(taxs, cpi, year, 1995))

# Create per capita income (in k)
data$income95pc <- data$income95/data$pop

# Create pretax price, 1995 dollars
data$pretax95 <- data$avgprs95 - data$taxs95


# Look at some time series
stateList <- unique(as.character(data$state))
nstate <- length(stateList)
yearList <- unique(data$year)
nyear <- length(yearList)
packpcMat <- taxs95Mat <- popMat <- NULL
for (istate in 1:nstate) {
    curstate <- stateList[istate]
    packpcMat <- cbind(packpcMat, data$packpc[as.character(data$state)==curstate])
    taxs95Mat <- cbind(taxs95Mat, data$taxs95[as.character(data$state)==curstate])
    popMat <- cbind(popMat, data$pop[as.character(data$state)==curstate])
}
packpcWAvg <- apply(packpcMat*popMat, 1, sum)/apply(popMat, 1, sum)
taxs95WAvg <- apply(taxs95Mat*popMat, 1, sum)/apply(popMat, 1, sum)


packpcTrace <- polylinesTile(y=packpcMat,
                             x=matrix(yearList, nrow=nyear, ncol=nstate),
                             id=matrix(1:nstate, nrow=nyear, ncol=nstate, byrow=TRUE),
                             col="gray80",
                             plot=1)

taxs95Trace <- polylinesTile(y=taxs95Mat/100,
                             x=matrix(yearList, nrow=nyear, ncol=nstate),
                             id=matrix(1:nstate, nrow=nyear, ncol=nstate, byrow=TRUE),
                             col="gray80",
                             plot=2)

packpcLabels <- textTile(x=rep(1995.75, nstate),
                         y=packpcMat[11,],
                         labels=labelTopBottomK(stateList, packpcMat[11,], 3),
                         cex=0.8,
                         col="gray80",
                         plot=1)

taxs95Labels <- textTile(x=rep(1995.75, nstate),
                         y=taxs95Mat[11,]/100,
                         labels=labelTopBottomK(stateList, taxs95Mat[11,], 3),
                         cex=0.8,
                         col="gray80",
                         plot=2)

packpcWAvgTrace <- linesTile(y=packpcWAvg, x=yearList,
                             col="red",
                             plot=1)

taxs95WAvgTrace <- linesTile(y=taxs95WAvg/100, x=yearList,
                             col="red",
                             plot=2)

packpcWAvgLabel <- textTile(y=packpcWAvg[11], x=1995.75, labels="US", cex=0.8, col="red", plot=1)

taxs95WAvgLabel <- textTile(y=taxs95WAvg[11]/100, x=1995.75, labels="US", cex=0.8, col="red", plot=2)


tile(packpcTrace, taxs95Trace,
     packpcLabels, taxs95Labels,
     packpcWAvgTrace, taxs95WAvgTrace,
     packpcWAvgLabel, taxs95WAvgLabel,
     limits=rbind(c(1984.5, 1996.5, 45, 210),
                  c(1984.5, 1996.5, 0.25, 1.20)),
     xaxis=list(at=c(1986, 1988, 1990, 1992, 1994)),
     yaxis=list(log=TRUE, major=FALSE),
     plottitle=list(labels1="Packs consumed per capita",
         labels2="Avg taxes per pack, $1995", y=1),
     height=list(plottile=3),
     output=list(file="cigByState", width=8)
     )
     


# Are there trends or lags? *to do

# Look at some unit root tests *to do



# Alternative model specifications
model1 <- packpc ~ income95pc + avgprs95 
model2 <- packpc ~ income95pc + pretax95 + taxs95
model3 <- log(packpc) ~ log(income95pc) + log(avgprs95)

# Simple linear models
lm.res1 <- lm(model1, data)
lm.res2 <- lm(model2, data)
lm.res3 <- lm(model3, data)


# Do plm model *to do

# Do RE model *to do


# Panel based diagnostics available in the plm library
# (This package recently expanded to contain many many panel data tests
#  for serial correlation, fixed effects, and unit roots)

# First, create a plm data frame (special data frame that "knows" the
# unit variable and time variable
pdata <- pdata.frame(data, index=c("state", "year"))


# Do an panel unit root test on the undifferenced cigarette data;
# there are many options; see ?purtest

# Note:  for some reason this isn't working
##purtest(packpc~1, data=pdata, test="ips")


#####
# Estimate Arellano-Bond GMM for fixed effects with lagged DV
#
# pgmm needs formulas in a specific format:
# 1. in the first part of the RHS, include lags of DV and covariates, as shown
# 2. in the second part, include the panel data instruments (99 here means use
#    up to the 99th lag of the difference as an instrument)
# 3. in an optional (not shown) third part of the RHS, include any other instruments
#
# note that pgmm formulas construct lag() properly for panel data,
# though lag() usually doesn't
pgmmformula.1a <- packpc ~ lag(packpc, 1) + income95pc + avgprs95 | lag(packpc, 2:99)

# We'll run GMM with only unit fixed effects,
# but we could include period fixed effects as well by setting effect to "two-way"
# (often a good practice in short T panels)
pgmm.res1a <- pgmm(pgmmformula.1a,
                   data = pdata,
                   effect = "individual",   # should consider two-way for small T
                   transformation = "d")    # should do ld if T=3

summary(pgmm.res1a)
# Good Sargan test, Good AR(2) test

# Let's consider alternative sets of instruments; concern: distant lags are weak instruments
pgmmformula.1b <- packpc ~ lag(packpc, 1) + income95pc + avgprs95 | lag(packpc, 2:5)
pgmm.res1b <- pgmm(pgmmformula.1b,
                  data = pdata,
                  effect = "individual",   # should consider two-way for small T
                  transformation = "d")    # should do ld if T=3

summary(pgmm.res1b)
# Poor Sargan test, Good AR(2) test

# Keeping just the most recent two instruments makes no substantive difference
pgmmformula.1c <- packpc ~ lag(packpc, 1) + income95pc + avgprs95 | lag(packpc, 2:3)
pgmm.res1c <- pgmm(pgmmformula.1c,
                  data = pdata,
                  effect = "individual",   # should consider two-way for small T
                  transformation = "d")    # should do ld if T=3

summary(pgmm.res1c)
# Poor Sargan test, Good AR(2) test


# Slight difference with one instrument, but not substatively noteworthy?
pgmmformula.1d <- packpc ~ lag(packpc, 1) + income95pc + avgprs95 | lag(packpc, 2)
pgmm.res1d <- pgmm(pgmmformula.1d,
                  data = pdata,
                  effect = "individual",   # should consider two-way for small T
                  transformation = "d")    # should do ld if T=3

summary(pgmm.res1d)
# Poor Sargan test, Good AR(2) test


# Try system GMM with all lags
pgmm.res1e <- pgmm(pgmmformula.1a,
                  data = pdata,
                  effect = "individual",   # should consider two-way for small T
                  transformation = "ld")    # should do ld if T=3

summary(pgmm.res1e)
# Good Sargan test, Good AR(2) test


# Try system GMM with only recent lag
pgmm.res1f <- pgmm(pgmmformula.1d,
                  data = pdata,
                  effect = "individual",   # should consider two-way for small T
                  transformation = "ld")    # should do ld if T=3

summary(pgmm.res1f)
# Poor Sargan test, Good AR(2) test


# Try difference GMM with two way effects
pgmm.res1g <- pgmm(pgmmformula.1a,
                  data = pdata,
                  effect = "twoways",   # should consider two-way for small T
                  transformation = "d")    # should do ld if T=3

summary(pgmm.res1g)
# Good Sargan test, Good AR(2) test, Wald supports 2-way


# Try system GMM with two way effects
pgmm.res1h <- pgmm(pgmmformula.1a,
                  data = pdata,
                  effect = "twoways",   # should consider two-way for small T
                  transformation = "ld")    # should do ld if T=3

summary(pgmm.res1h)
# Good Sargan test, Good AR(2) test, Wald supports 2-way


# Now consider last two models with alternative specifications
pgmmformula.2a <- packpc ~ lag(packpc, 1) + income95pc + pretax95 + taxs95 | lag(packpc, 2:99)
pgmmformula.3a <- log(packpc) ~ lag(log(packpc), 1) + log(income95pc) + log(avgprs95) | lag(log(packpc), 2:99)
pgmmformula.4a <- log(packpc) ~ lag(log(packpc), 1) + log(income95pc) + log(pretax95) + log(taxs95) | lag(log(packpc), 2:99)


# Model 2: Unique tax effects

# Try difference GMM with two way effects
pgmm.res2g <- pgmm(pgmmformula.2a,
                   data = pdata,
                   effect = "twoways",   # should consider two-way for small T
                   transformation = "d")    # should do ld if T=3

summary(pgmm.res2g)
# Good Sargan test, Good AR(2) test, Wald supports 2-way

# Try system GMM with two way effects
pgmm.res2h <- pgmm(pgmmformula.2a,
                   data = pdata,
                   effect = "twoways",   # should consider two-way for small T
                   transformation = "ld")    # should do ld if T=3

summary(pgmm.res2h)
# Good Sargan test, Good AR(2) test, Wald supports 2-way


# Model 3: Elasticity specification

# Try difference GMM with only unit fixed effects
pgmm.res3a <- pgmm(pgmmformula.3a,
                   data = pdata,
                   effect = "individual",   # should consider two-way for small T
                   transformation = "d")    # should do ld if T=3

summary(pgmm.res3a)
# Good Sargan test, Good AR(2) test


# Try system GMM with all lags
pgmm.res3e <- pgmm(pgmmformula.3a,
                  data = pdata,
                  effect = "individual",   # should consider two-way for small T
                  transformation = "ld")    # should do ld if T=3

summary(pgmm.res3e)
# Good Sargan test, Good AR(2) test


# Try difference GMM with two way effects
pgmm.res3g <- pgmm(pgmmformula.3a,
                   data = pdata,
                   effect = "twoways",   # should consider two-way for small T
                   transformation = "d")    # should do ld if T=3

summary(pgmm.res3g)
# Good Sargan test, Good AR(2) test, Wald supports 2-way

# Try system GMM with two way effects
pgmm.res3h <- pgmm(pgmmformula.3a,
                   data = pdata,
                   effect = "twoways",   # should consider two-way for small T
                   transformation = "ld")    # should do ld if T=3

summary(pgmm.res3h)
# Good Sargan test, Good AR(2) test, Wald supports 2-way


# Model 4: Elasticity specification, components of price

# Try difference GMM with two way effects
pgmm.res4g <- pgmm(pgmmformula.4a,
                   data = pdata,
                   effect = "twoways",   # should consider two-way for small T
                   transformation = "d")    # should do ld if T=3

summary(pgmm.res4g)
# Good Sargan test, Good AR(2) test, Wald supports 2-way

# Try system GMM with two way effects
pgmm.res4h <- pgmm(pgmmformula.4a,
                   data = pdata,
                   effect = "twoways",   # should consider two-way for small T
                   transformation = "ld" # should do ld if T=3
                   )    

summary(pgmm.res4h)
# Good Sargan test, Good AR(2) test, Wald supports 2-way





#######
# Simulate conditional forecasts from 8 different models
# 1a, 1e, 1g, 1h, 3a, 3e, 3g, 3h


# Forecast for 3 years from 1996 to 1998
periods.out <- 3
sims <- 1000

# How big a change in price to simulate?
# How about "double" the average tax in the most recent year?
summary(pdata$taxs95[pdata$year==1995])
# The average (and median) tax is about 60 cents/pack
sd(pdata$taxs95[pdata$year==1995])
# A 60 cent increase would also be about 3 sd's,
# and raise the tax to a bit more than the max observed

# Other possibilities:
# (2) A 10 cent increase
# (3) Raise every sate to the max observed for any state in 1995 (112.60 cents)


# Construct the year dummies
yearfe <- makeFEdummies(pdata$year)
yearfe <- yearfe[,3:ncol(yearfe)]  # Why drop first 2 col's?
yearlist <- unique(pdata$year)
yearlist <- yearlist[3:length(yearlist)]
colnames(yearfe) <- paste0("y",yearlist)

# Construct formulas -- without year dummies (1a)
formula.1a <- packpc ~ income95pc + avgprs95 -1

# Construct formulas -- without year dummies but with intercept (1e)
formula.1e <- packpc ~ income95pc + avgprs95

# Construct formulas -- with year dummies (1g)
formula <- "packpc ~ income95pc + avgprs95 -1"
datayearfe <- cbind(pdata,yearfe)
yearfenames <- NULL
for (i in 1:ncol(yearfe)) {
  formula <- paste0(formula,"+ y",yearlist[i]," ")
  yearfenames <- c(yearfenames,paste0("y",yearlist[i]))
}
names(datayearfe) <- c(names(data),yearfenames)
formula.1g <- as.formula(formula)

# Construct formulas -- with year dummies and intercept (1h)
formula <- "packpc ~ income95pc + avgprs95"
datayearfe <- cbind(pdata,yearfe)
yearfenames <- NULL
for (i in 1:ncol(yearfe)) {
  formula <- paste0(formula,"+ y",yearlist[i]," ")
  yearfenames <- c(yearfenames,paste0("y",yearlist[i]))
}
names(datayearfe) <- c(names(data),yearfenames)
formula.1h <- as.formula(formula)

# Population in 1995 in average state
avgpop1995 <- mean(pdata$pop[pdata$year==1995])


# Forecast: Model 1a, +60

# Simulate parameters
simparam.1a <- mvrnorm(sims, coefficients(pgmm.res1a), vcovHC(pgmm.res1a))
simphis.1a <- simparam.1a[,1]
simbetas.1a <- simparam.1a[,2:ncol(simparam.1a)]


# Make matrix of hypothetical x's:
# Assume an average state raised taxes 60 cents starting 1996
#

# Make matrix of hypothetical x's: covariates
xhyp.1a <- cfMake(formula.1a,datayearfe, periods.out)

# pgmm uses covariates in differenced form
# so we want most of them to be 0 (no change)
# exceptions:
# (1) changes in covariates of interest
# (2) differenced time dummies require special care
xhyp.1a$x <- xhyp.1a$xpre <- 0*xhyp.1a$x
xhyp.1a <- cfChange(xhyp.1a, "avgprs95", x=60, scen=1)

# We can "ignore" the state fixed effects for now and add them later
# because model is total linear

# Create baseline scenario
xbase.1a <- xhyp.1a
xbase.1a$x <- xbase.1a$xpre

# We need a lag of the price per pack
lagY.1a <- NULL # Hypothetical previous change in Y for simulation
for (i in 1:length(pgmm.res1a$model))
    lagY.1a <- c(lagY.1a, as.data.frame(pgmm.res1a$model[[i]])["1995",]$packpc)
lagY.1a <- mean(lagY.1a, na.rm=TRUE)

# Hypothetical initial level of Y for simulation
initialY <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)


# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev1a <- ldvsimev(xhyp.1a,               # The matrix of hypothetical x's
                     simbetas.1a,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # NA indicates no constant!
                     phi=simphis.1a,            # estimated AR parameters; length must match lagY 
                     lagY=lagY.1a,          # lags of y, most recent last
                     transform="diff",   # "log" to undo log transformation,
                                         # "diff" to under first differencing
                                         # "difflog" to do both
                     initialY=initialY   # for differenced models, the lag of the level of y
                     )

# Simulate expected values of Y given no change in covariates
sim.base1a <- ldvsimev(xbase.1a,               # The matrix of hypothetical x's
                       simbetas.1a,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=NA,        # NA indicates no constant!
                       phi=simphis.1a,            # estimated AR parameters; length must match lagY 
                       lagY=lagY.1a,          # lags of y, most recent last
                       transform="diff",   # "log" to undo log transformation,
                                         # "diff" to under first differencing
                                         # "difflog" to do both
                       initialY=initialY   # for differenced models, the lag of the level of y
                       )    

# Simulate first differences in y 
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.fd1a <- ldvsimfd(xhyp.1a,            # The matrix of hypothetical x's
                     simbetas.1a,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.1a,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.1a,       # lags of y, most recent last
                     transform="diff",   # Model is differenced
                     initialY=initialY   # Redundant in this case (fd of linear differenced Y)
                     )

# Simulate relative risks in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.rr1a <- ldvsimrr(xhyp.1a,            # The matrix of hypothetical x's
                     simbetas.1a,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.1a,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.1a,       # lags of y, most recent last
                     transform="diff",   # Model is differenced
                     initialY=initialY   # Required for differenced Y in ldvsimrr
                     )


# Forecast: Model 1e, +60

# Simulate parameters
simparam.1e <- mvrnorm(sims, coefficients(pgmm.res1e), vcovHC(pgmm.res1e))
simphis.1e <- simparam.1e[,1]
simbetas.1e <- simparam.1e[,2:ncol(simparam.1e)]

# System GMM does NOT difference the covariates

# Make matrix of hypothetical x's:
# Assume an average state raised taxes 60 cents starting 1996
#

# Make matrix of hypothetical x's: covariates
xhyp.1e <- cfMake(formula.1a, datayearfe, periods.out)

# system pgmm uses covariates in *level* form
# -> back to our usual use of simcf; note apply to all 3 periods!
xhyp.1e <- cfChange(xhyp.1e, "avgprs95", x=60 + mean(pdata$avgprs95), scen=1:3)

# State fixed effects are not removed from the covariates,
# but from the instruments (so we can ignore them here)

# Create baseline scenario
xbase.1e <- xhyp.1e
xbase.1e$x <- xbase.1e$xpre

# We need a lag of the price per pack, now in levels
# But the code above to extract it from the pgmm object won't work!
lagY.1e <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)

# Hypothetical initial level of Y for simulation
initialY <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)


# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev1e <- ldvsimev(xhyp.1e,               # The matrix of hypothetical x's
                     simbetas.1e,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # NA indicates no constant!
                     phi=simphis.1e,            # estimated AR parameters; length must match lagY 
                     lagY=lagY.1e,          # lags of y, most recent last
                     transform="none",   # NOTE: System GMM is not differenced!
                     initialY=initialY
                     )

# Simulate expected values of Y given no change in covariates
sim.base1e <- ldvsimev(xbase.1e,               # The matrix of hypothetical x's
                       simbetas.1e,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=NA,        # NA indicates no constant!
                       phi=simphis.1e,            # estimated AR parameters; length must match lagY 
                       lagY=lagY.1e,          # lags of y, most recent last
                       transform="none"   # NOTE: System GMM is not differenced!
                       )    

# Simulate first differences in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.fd1e <- ldvsimfd(xhyp.1e,            # The matrix of hypothetical x's
                     simbetas.1e,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.1e,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.1e,       # lags of y, most recent last
                     transform="none"    # NOTE: System GMM is not differenced!
                     )

# Simulate relative risks in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.rr1e <- ldvsimrr(xhyp.1e,            # The matrix of hypothetical x's
                     simbetas.1e,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.1e,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.1e,       # lags of y, most recent last
                     transform="none"    # NOTE: System GMM is not differenced!
                     )



# Forecast: Model 1g, +60

# Simulate parameters
simparam.1g <- mvrnorm(sims, coefficients(pgmm.res1g), vcovHC(pgmm.res1g))
simphis.1g <- simparam.1g[,1]
simbetas.1g <- simparam.1g[,2:ncol(simparam.1g)]

# Make matrix of hypothetical x's:
# Assume an average state raised taxes 60 cents starting 1996
#
# Issues -- we need to somehow include the state and year FEs:
#           Let's set the state to be an "average" state in 1995,
#           and year to be like the last year (1995)

# Make matrix of hypothetical x's: covariates
xhyp.1g <- cfMake(formula.1g, datayearfe, periods.out)

# pgmm uses covariates in differenced form
# so we want most of them to be 0 (no change)
# exceptions:
# (1) changes in covariates of interest
# (2) differenced time dummies require special care
xhyp.1g$x <- xhyp.1g$xpre <- 0*xhyp.1g$x
xhyp.1g <- cfChange(xhyp.1g, "avgprs95", x=60, scen=1)

# The lines below would "repeat" the 1995 year effect in each
# subsequent period -- but this is an aberrant year effect
#xhyp.1g <- cfChange(xhyp.1g, "y1994", x=-1, xpre=-1, scen=1:3)
#xhyp.1g <- cfChange(xhyp.1g, "y1995", x=1, xpre=1, scen=1:3)

# The current trend seems to start in 1993; we will average over the
# the last three years of year effects:
#   = 1/3*( (tau93 - tau92) - (tau94 - tau93) - (tau95 - tau94) )
#   = 1/3*tau95 - 1/3*tau92
xhyp.1g <- cfChange(xhyp.1g, "y1992", x=-1/3, xpre=-1/3, scen=1:3)
xhyp.1g <- cfChange(xhyp.1g, "y1995", x=1/3, xpre=1/3, scen=1:3)

# We can "ignore" the state fixed effects for now and add them later
# because model is total linear

# Create baseline scenario
xbase.1g <- xhyp.1g
xbase.1g$x <- xbase.1g$xpre

# We need a lag of the price per pack
lagY.1g <- NULL # Hypothetical previous change in Y for simulation
for (i in 1:length(pgmm.res1g$model))
    lagY.1g <- c(lagY.1g, as.data.frame(pgmm.res1g$model[[i]])["1995",]$packpc)
lagY.1g <- mean(lagY.1g, na.rm=TRUE)

# Hypothetical initial level of Y for simulation
initialY <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)

# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev1g <- ldvsimev(xhyp.1g,               # The matrix of hypothetical x's
                     simbetas.1g,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # NA indicates no constant!
                     phi=simphis.1g,            # estimated AR parameters; length must match lagY 
                     lagY=lagY.1g,          # lags of y, most recent last
                     transform="diff",   # "log" to undo log transformation,
                                         # "diff" to under first differencing
                                         # "difflog" to do both
                     initialY=initialY   # for differenced models, the lag of the level of y
                     )

# Simulate expected values of Y given no change in covariates
sim.base1g <- ldvsimev(xbase.1g,               # The matrix of hypothetical x's
                       simbetas.1g,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=NA,        # NA indicates no constant!
                       phi=simphis.1g,            # estimated AR parameters; length must match lagY 
                       lagY=lagY.1g,          # lags of y, most recent last
                       transform="diff",   # "log" to undo log transformation,
                                         # "diff" to under first differencing
                                         # "difflog" to do both
                       initialY=initialY   # for differenced models, the lag of the level of y
                       )   

# Simulate first differences in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.fd1g <- ldvsimfd(xhyp.1g,            # The matrix of hypothetical x's
                     simbetas.1g,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.1g,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.1g,       # lags of y, most recent last
                     transform="diff",   # Model is differenced
                     initialY=initialY   # Redundant in this case (fd of linear differenced Y)
                     )

# Simulate relative risks in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.rr1g <- ldvsimrr(xhyp.1g,            # The matrix of hypothetical x's
                     simbetas.1g,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.1g,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.1g,       # lags of y, most recent last
                     transform="diff",   # Model is differenced
                     initialY=initialY   # Required for differenced Y in ldvsimrr
                     )




# Forecast: Model 1h, +60

# Simulate parameters
simparam.1h <- mvrnorm(sims, coefficients(pgmm.res1h), vcovHC(pgmm.res1h))
simphis.1h <- simparam.1h[,1]
simbetas.1h <- simparam.1h[,2:ncol(simparam.1h)]

# System GMM does NOT difference the covariates
# -> with 2-way effects, the model has a constant,
# which pgmm() puts in an odd place
# Move the constant to the front of the matrix!
simbetas.1h <- cbind(simbetas.1h[,3], simbetas.1h[,-3])

# Make matrix of hypothetical x's:
# Assume an average state raised taxes 60 cents starting 1996
#
# Issues -- we need to somehow include the state and year FEs:
#           Let's set the state to be an "average" state in 1995,
#           and year to be like the last year (1995)

# Make matrix of hypothetical x's: covariates
xhyp.1h <- cfMake(formula.1h, datayearfe, periods.out)

# system pgmm uses covariates in *level* form
# -> back to our usual use of simcf; note apply to all 3 periods!
xhyp.1h <- cfChange(xhyp.1h, "avgprs95", x=60 + mean(pdata$avgprs95), scen=1:3)

# The current trend seems to start in 1993; we will average over the
# the last three years of year effects:
xhyp.1h <- cfChange(xhyp.1h, "y1987", x=0, xpre=0, scen=1:3)
xhyp.1h <- cfChange(xhyp.1h, "y1988", x=0, xpre=0, scen=1:3)
xhyp.1h <- cfChange(xhyp.1h, "y1989", x=0, xpre=0, scen=1:3)
xhyp.1h <- cfChange(xhyp.1h, "y1990", x=0, xpre=0, scen=1:3)
xhyp.1h <- cfChange(xhyp.1h, "y1991", x=0, xpre=0, scen=1:3)
xhyp.1h <- cfChange(xhyp.1h, "y1992", x=0, xpre=0, scen=1:3)
xhyp.1h <- cfChange(xhyp.1h, "y1993", x=1/3, xpre=1/3, scen=1:3)
xhyp.1h <- cfChange(xhyp.1h, "y1994", x=1/3, xpre=1/3, scen=1:3)
xhyp.1h <- cfChange(xhyp.1h, "y1995", x=1/3, xpre=1/3, scen=1:3)

# State fixed effects are not removed from the covariates,
# but from the instruments (so we can ignore them here)

# Create baseline scenario
xbase.1h <- xhyp.1h
xbase.1h$x <- xbase.1h$xpre

# We need a lag of the price per pack, now in levels
# But the code above to extract it from the pgmm object won't work!
lagY.1h <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)

# Hypothetical initial level of Y for simulation
initialY <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)

# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev1h <- ldvsimev(xhyp.1h,               # The matrix of hypothetical x's
                     simbetas.1h,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=1,        # NOTE: System GMM has a constant!
                                         # You will need to note the column of the constant in simbetas
                     phi=simphis.1h,            # estimated AR parameters; length must match lagY 
                     lagY=lagY.1h,          # lags of y, most recent last
                     transform="none"   # NOTE: System GMM is not differenced!
                     )

# Simulate expected values of Y given no change in covariates
sim.base1h <- ldvsimev(xbase.1h,               # The matrix of hypothetical x's
                       simbetas.1h,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=1,        # NOTE: System GMM has a constant!
                                        # You will need to note the column of the constant in simbetas
                       phi=simphis.1h,            # estimated AR parameters; length must match lagY 
                       lagY=lagY.1h,          # lags of y, most recent last
                       transform="none"   # NOTE: System GMM is not differenced!
                       ) 

# Simulate first differences in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.fd1h <- ldvsimfd(xhyp.1h,            # The matrix of hypothetical x's
                     simbetas.1h,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=1,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.1h,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.1h,       # lags of y, most recent last
                     transform="none"   # NOTE: System GMM is not differenced!
                     )

# Simulate relative risks in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.rr1h <- ldvsimrr(xhyp.1h,            # The matrix of hypothetical x's
                     simbetas.1h,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=1,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.1h,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.1h,       # lags of y, most recent last
                     transform="none"   # NOTE: System GMM is not differenced!
                     )


#######
# Try log specification

# Forecast: Model 3a, +60

# Simulate parameters
simparam.3a <- mvrnorm(sims, coefficients(pgmm.res3a), vcovHC(pgmm.res3a))
simphis.3a <- simparam.3a[,1]
simbetas.3a <- simparam.3a[,2:ncol(simparam.3a)]

# Make matrix of hypothetical x's:
# Assume an average state raised taxes 60 cents starting 1996
#

# Make matrix of hypothetical x's: covariates
xhyp.3a <- cfMake(formula.1a, datayearfe, periods.out)

# pgmm uses covariates in differenced form
# so we want most of them to be 0 (no change)
# exceptions:
# (1) changes in covariates of interest
# (2) differenced time dummies require special care
xhyp.3a$x <- xhyp.3a$xpre <- 0*xhyp.3a$x

# Need log version of differenced key covariate (doubling tax in avg state)
meanPrice95 <- mean(pdata$avgprs95[pdata$year==1995], na.rm=TRUE)
meanTaxs95 <- mean(pdata$taxs95[pdata$year==1995], na.rm=TRUE)

xhyp.3a <- cfChange(xhyp.3a, "avgprs95",
                    x=log(meanPrice95+meanTaxs95) - log(meanPrice95),
                    scen=1)

# We can "ignore" the state fixed effects for now and add them later
# because model is total linear

# Create baseline scenario
xbase.3a <- xhyp.3a
xbase.3a$x <- xbase.3a$xpre

# We need a lag of the price per pack
lagY.3a <- NULL # Hypothetical previous change in Y for simulation
for (i in 1:length(pgmm.res3a$model))
    lagY.3a <- c(lagY.3a, as.data.frame(pgmm.res3a$model[[i]])["1995",1])
lagY.3a <- mean(lagY.3a, na.rm=TRUE)

# Hypothetical initial level of Y for simulation
initialY <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)


# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev3a <- ldvsimev(xhyp.3a,               # The matrix of hypothetical x's
                     simbetas.3a,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # NA indicates no constant!
                     phi=simphis.3a,            # estimated AR parameters; length must match lagY 
                     lagY=lagY.3a,          # lags of y, most recent last
                     transform="difflog",   # "log" to undo log transformation,
                                         # "diff" to under first differencing
                                         # "difflog" to do both
                     initialY=initialY   # for differenced models, the lag of the level of y
                     )

# Simulate expected values of Y given no change in covariates
sim.base3a <- ldvsimev(xbase.3a,               # The matrix of hypothetical x's
                       simbetas.3a,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=NA,        # NA indicates no constant!
                       phi=simphis.3a,            # estimated AR parameters; length must match lagY 
                       lagY=lagY.3a,          # lags of y, most recent last
                       transform="difflog",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                         # "difflog" to do both
                       initialY=initialY   # for differenced models, the lag of the level of y
                       )

# Simulate first differences in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.fd3a <- ldvsimfd(xhyp.3a,            # The matrix of hypothetical x's
                     simbetas.3a,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.3a,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.3a,       # lags of y, most recent last
                     transform="difflog",# Model is differenced
                     initialY=initialY   # Required in this case (fd of differenced log Y)
                     )

# Simulate relative risks in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.rr3a <- ldvsimrr(xhyp.3a,            # The matrix of hypothetical x's
                     simbetas.3a,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.3a,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.3a,       # lags of y, most recent last
                     transform="difflog",# Model is differenced
                     initialY=initialY   # Required for differenced Y in ldvsimrr
                     )



# Forecast: Model 3e, +60

# Because system GMM is in levels, it is convenient to
# handle logging through the formula combined with simcf
formula.3e <- log(packpc) ~ log(income95pc) + log(avgprs95) -1

# Simulate parameters
simparam.3e <- mvrnorm(sims, coefficients(pgmm.res3e), vcovHC(pgmm.res3e))
simphis.3e <- simparam.3e[,1]
simbetas.3e <- simparam.3e[,2:ncol(simparam.3e)]


# System GMM does NOT difference the covariates

# Make matrix of hypothetical x's:
# Assume an average state raised taxes 60 cents starting 1996
#
# Make matrix of hypothetical x's: covariates
xhyp.3e <- cfMake(formula.3e, datayearfe, periods.out)

# system pgmm uses covariates in *level* form
# -> back to our usual use of simcf; note apply to all 3 periods!
xhyp.3e <- cfChange(xhyp.3e, "avgprs95", x=60 + mean(pdata$avgprs95), scen=1:3)

# State fixed effects are not removed from the covariates,
# but from the instruments (so we can ignore them here)

# Create baseline scenario
xbase.3e <- xhyp.3e
xbase.3e$x <- xbase.3e$xpre

# We need a lag of the price per pack, now in logged levels
# But the code above to extract it from the pgmm object won't work!
# Getting this right is crucial
lagY.3e <- log(mean(pdata$packpc[pdata$year==1995], na.rm=TRUE))

# Hypothetical initial level of Y for simulation
# Still in linear levels
initialY <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)


# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev3e <- ldvsimev(xhyp.3e,               # The matrix of hypothetical x's
                     simbetas.3e,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # NA indicates no constant!
                     phi=simphis.3e,            # estimated AR parameters; length must match lagY 
                     lagY=lagY.3e,          # lags of y, most recent last
                     transform="log"   # NOTE: System GMM is not differenced!
                     )

# Simulate expected values of Y given no change in covariates
sim.base3e <- ldvsimev(xbase.3e,               # The matrix of hypothetical x's
                     simbetas.3e,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # NA indicates no constant!
                     phi=simphis.3e,            # estimated AR parameters; length must match lagY 
                     lagY=lagY.3e,          # lags of y, most recent last
                     transform="log",   # NOTE: System GMM is not differenced!
                     )    

# Simulate first differences in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.fd3e <- ldvsimfd(xhyp.3e,            # The matrix of hypothetical x's
                     simbetas.3e,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.3e,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.3e,       # lags of y, most recent last
                     transform="log"     # NOTE: System GMM is not differenced!
                     )

# Simulate relative risks in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.rr3e <- ldvsimrr(xhyp.3e,            # The matrix of hypothetical x's
                     simbetas.3e,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.3e,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.3e,       # lags of y, most recent last
                     transform="log"     # NOTE: System GMM is not differenced!
                     )


# Forecast: Model 3g, +60

simparam.3g <- mvrnorm(sims, coefficients(pgmm.res3g), vcovHC(pgmm.res3g))
simphis.3g <- simparam.3g[,1]
simbetas.3g <- simparam.3g[,2:ncol(simparam.3g)]

# Make matrix of hypothetical x's:
# Assume an average state raised taxes 60 cents starting 1996
#
# Issues -- we need to somehow include the state and year FEs:
#           Let's set the state to be an "average" state in 1995,
#           and year to be like the last year (1995)

# Make matrix of hypothetical x's: covariates
# Still use the 1g formula (no logs) -- we will handle logging manually
#  to get the differences of logs right
xhyp.3g <- cfMake(formula.1g, datayearfe, periods.out)

# pgmm uses covariates in differenced form
# so we want most of them to be 0 (no change)
# exceptions:
# (1) changes in covariates of interest
# (2) differenced time dummies need special care
xhyp.3g$x <- xhyp.3g$xpre <- 0*xhyp.3g$x

# Need log version of differenced key covariate (doubling tax in avg state)
meanPrice95 <- mean(pdata$avgprs95[pdata$year==1995], na.rm=TRUE)
meanTaxs95 <- mean(pdata$taxs95[pdata$year==1995], na.rm=TRUE)

xhyp.3g <- cfChange(xhyp.3g, "avgprs95",
                    x=log(meanPrice95+meanTaxs95) - log(meanPrice95),
                    scen=1)

# The lines below would "repeat" the 1995 year effect in each
# subsequent period -- but this is an aberrant year effect
#xhyp.1g <- cfChange(xhyp.1g, "y1994", x=-1, xpre=-1, scen=1:3)
#xhyp.1g <- cfChange(xhyp.1g, "y1995", x=1, xpre=1, scen=1:3)

# The current trend seems to start in 1993; we will average over the
# the last three years of year effects:
#   = 1/3*( (tau93 - tau92) - (tau94 - tau93) - (tau95 - tau94) )
#   = 1/3*tau95 - 1/3*tau92
xhyp.3g <- cfChange(xhyp.3g, "y1992", x=-1/3, xpre=-1/3, scen=1:3)
xhyp.3g <- cfChange(xhyp.3g, "y1995", x=1/3, xpre=1/3, scen=1:3)

# We can "ignore" the state fixed effects for now and add them later
# because model is total linear

# Create baseline scenario
xbase.3g <- xhyp.3g
xbase.3g$x <- xbase.3g$xpre

# We need a lag of the price per pack
lagY.3g <- NULL # Hypothetical previous change in Y for simulation
for (i in 1:length(pgmm.res3g$model))
    lagY.3g <- c(lagY.3g, as.data.frame(pgmm.res3g$model[[i]])["1995",1])
lagY.3g <- mean(lagY.3g, na.rm=TRUE)


initialY <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)
               # Hypothetical initial level of Y for simulation

# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev3g <- ldvsimev(xhyp.3g,               # The matrix of hypothetical x's
                     simbetas.3g,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # NA indicates no constant!
                     phi=simphis.3g,            # estimated AR parameters; length must match lagY 
                     lagY=lagY.3g,          # lags of y, most recent last
                     transform="difflog",   # "log" to undo log transformation,
                                         # "diff" to under first differencing
                                         # "difflog" to do both
                     initialY=initialY   # for differenced models, the lag of the level of y
                     )

# Simulate expected values of Y given no change in covariates
sim.base3g <- ldvsimev(xbase.3g,               # The matrix of hypothetical x's
                       simbetas.3g,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=NA,        # NA indicates no constant!
                       phi=simphis.3g,            # estimated AR parameters; length must match lagY 
                       lagY=lagY.3g,          # lags of y, most recent last
                       transform="difflog",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                       initialY=initialY   # for differenced models, the lag of the level of y
                       ) 

# Simulate first differences in y (as cumulated changes)
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.fd3g <- ldvsimfd(xhyp.3g,            # The matrix of hypothetical x's
                     simbetas.3g,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.3g,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.3g,       # lags of y, most recent last
                     transform="difflog",# Model is differenced logs
                     initialY=initialY   # Required in this case (fd of differenced log Y)
                     )

# Simulate relative risks in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.rr3g <- ldvsimrr(xhyp.3g,            # The matrix of hypothetical x's
                     simbetas.3g,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.3g,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.3g,       # lags of y, most recent last
                     transform="difflog",# Model is differenced logs
                     initialY=initialY   # Required for differenced Y in ldvsimrr
                     )



# Forecast: Model 3h, +60

# Because system GMM is in levels, it is convenient to
# handle logging through the formula combined with simcf
formula <- "log(packpc) ~ log(income95pc) + log(avgprs95)"
datayearfe <- cbind(pdata,yearfe)
yearfenames <- NULL
for (i in 1:ncol(yearfe)) {
  formula <- paste0(formula,"+ y",yearlist[i]," ")
  yearfenames <- c(yearfenames,paste0("y",yearlist[i]))
}
names(datayearfe) <- c(names(data),yearfenames)
formula.3h <- as.formula(formula)

# Simulate parameters
simparam.3h <- mvrnorm(sims, coefficients(pgmm.res3h), vcovHC(pgmm.res3h))
simphis.3h <- simparam.3h[,1]
simbetas.3h <- simparam.3h[,2:ncol(simparam.3h)]

# System GMM does NOT difference the covariates
# -> the model has a constant, which pgmm() puts in an odd place
# Move the constant to the front of the matrix!
simbetas.3h <- cbind(simbetas.3h[,3], simbetas.3h[,-3])

# Make matrix of hypothetical x's:
# Assume an average state raised taxes 60 cents starting 1996
#
# Issues -- we need to somehow include the state and year FEs:
#           Let's set the state to be an "average" state in 1995,
#           and year to be like the last year (1995)

# Make matrix of hypothetical x's: covariates
xhyp.3h <- cfMake(formula.3h, datayearfe, periods.out)

# system pgmm uses covariates in *level* form
# -> back to our usual use of simcf; let simcf handle logging here.
# The current trend seems to start in 1993; we will average over the
# the last three years of year effects:
xhyp.3h <- cfChange(xhyp.3h, "avgprs95", x=60 + mean(pdata$avgprs95), scen=1:3)
xhyp.3h <- cfChange(xhyp.3h, "y1987", x=0, xpre=0, scen=1:3)
xhyp.3h <- cfChange(xhyp.3h, "y1988", x=0, xpre=0, scen=1:3)
xhyp.3h <- cfChange(xhyp.3h, "y1989", x=0, xpre=0, scen=1:3)
xhyp.3h <- cfChange(xhyp.3h, "y1990", x=0, xpre=0, scen=1:3)
xhyp.3h <- cfChange(xhyp.3h, "y1991", x=0, xpre=0, scen=1:3)
xhyp.3h <- cfChange(xhyp.3h, "y1992", x=0, xpre=0, scen=1:3)
xhyp.3h <- cfChange(xhyp.3h, "y1993", x=1/3, xpre=1/3, scen=1:3)
xhyp.3h <- cfChange(xhyp.3h, "y1994", x=1/3, xpre=1/3, scen=1:3)
xhyp.3h <- cfChange(xhyp.3h, "y1995", x=1/3, xpre=1/3, scen=1:3)

# State fixed effects are not removed from the covariates,
# but from the instruments (so we can ignore them here)

# Create baseline scenario
xbase.3h <- xhyp.3h
xbase.3h$x <- xbase.3h$xpre

# We need a lag of the price per pack, now in logged levels
# But the code above to extract it from the pgmm object won't work!
# Getting this right is crucial
lagY.3h <- log(mean(pdata$packpc[pdata$year==1995], na.rm=TRUE))

# Hypothetical initial level of Y for simulation
# Still in linear levels
initialY <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)

# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev3h <- ldvsimev(xhyp.3h,               # The matrix of hypothetical x's
                     simbetas.3h,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=1,        # NOTE: System GMM with two-way effects has a constant!
                                         # You will need to note the column of the constant in simbetas
                     phi=simphis.3h,            # estimated AR parameters; length must match lagY 
                     lagY=lagY.3h,          # lags of y, most recent last
                     transform="log"   # NOTE: System GMM is not differenced!
                     )

# Simulate expected values of Y given no change in covariates
sim.base3h <- ldvsimev(xbase.3h,               # The matrix of hypothetical x's
                       simbetas.3h,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=1,        # NOTE: System GMM with two-way effects has a constant!
                                        # You will need to note the column of the constant in simbetas
                       phi=simphis.3h,            # estimated AR parameters; length must match lagY 
                       lagY=lagY.3h,          # lags of y, most recent last
                       transform="log"   # NOTE: System GMM is not differenced!
                       )

# Simulate first differences in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.fd3h <- ldvsimfd(xhyp.3h,            # The matrix of hypothetical x's
                     simbetas.3h,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=1,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.3h,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.3h,       # lags of y, most recent last
                     transform="log"   # NOTE: System GMM is not differenced!
                     )

# Simulate relative risks in y
# out to periods.out given hypothetical future values of x, xpre,
# and initial lags of the change in y
sim.rr3h <- ldvsimrr(xhyp.3h,            # The matrix of hypothetical x's
                     simbetas.3h,        # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=1,        # Column containing the constant
                                         # set to NA for no constant
                     phi=simphis.3h,     # estimated AR parameters; length must match lagY 
                     lagY=lagY.3h,       # lags of y, most recent last
                     transform="log"   # NOTE: System GMM is not differenced!
                     )



# Make plots of expected values, first differences, and percent changes
# using custom tile code in helperCigs.R

# Hypothetical initial level of Y for simulation
initialY <- mean(pdata$packpc[pdata$year==1995], na.rm=TRUE)

# axis limits
limitsEV <- c(50, 105)
limitsFD <- c(-40, 5)
limitsRR <- c(-40, 5)

# Model 1a Lineplot: Packs
cigLineplots(sim.ev1a, sim.base1a, sim.fd1a, sim.rr1a,
             limitsEV, limitsFD, limitsRR, initialY,
             "Cigarette Taxes & Consumption: 1a. Linear Difference GMM, Individual Effects",
             "m1aEVFDRR"
             )

# Model 1e Lineplot: Packs
cigLineplots(sim.ev1e, sim.base1e, sim.fd1e, sim.rr1e,
             limitsEV, limitsFD, limitsRR, initialY,
             "Cigarette Taxes & Consumption: 1e. Linear System GMM, Individual Effects",
             "m1eEVFDRR"
             )

# Model 1g Lineplot: Packs
cigLineplots(sim.ev1g, sim.base1g, sim.fd1g, sim.rr1g,
             limitsEV, limitsFD, limitsRR, initialY,
             "Cigarette Taxes & Consumption: 1g. Linear Difference GMM, Two-Way Effects",
             "m1gEVFDRR"
             )

# Model 1h Lineplot: Packs
cigLineplots(sim.ev1h, sim.base1h, sim.fd1h, sim.rr1h,
             limitsEV, limitsFD, limitsRR, initialY,
             "Cigarette Taxes & Consumption: 1h. Linear System GMM, Two-Way Effects",
             "m1hEVFDRR"
             )

# Model 3a Lineplot: Packs
cigLineplots(sim.ev3a, sim.base3a, sim.fd3a, sim.rr3a,
             limitsEV, limitsFD, limitsRR, initialY,
             "Cigarette Taxes & Consumption: 3a. Log-Log Difference GMM, Individual Effects",
             "m3aEVFDRR"
             )

# Model 3e Lineplot: Packs
cigLineplots(sim.ev3e, sim.base3e, sim.fd3e, sim.rr3e,
             limitsEV, limitsFD, limitsRR, initialY,
             "Cigarette Taxes & Consumption: 3e. Log-Log System GMM, Individual Effects",
             "m3eEVFDRR"
             )

# Model 3g Lineplot: Packs
cigLineplots(sim.ev3g, sim.base3g, sim.fd3g, sim.rr3g,
             limitsEV, limitsFD, limitsRR, initialY,
             "Cigarette Taxes & Consumption: 3g. Log-Log Difference GMM, Two-Way Effects",
             "m3gEVFDRR"
             )

# Model 3h Lineplot: Packs
cigLineplots(sim.ev3h, sim.base3h, sim.fd3h, sim.rr3h,
             limitsEV, limitsFD, limitsRR, initialY,
             "Cigarette Taxes & Consumption: 3h. Log-Log System GMM, Two-Way Effects",
             "m3hEVFDRR"
             )


#####
# Ropeladder to show model robustness (with CIs) : FD
k <- 1

traceDiffGMM1yearFD <- ropeladder(x=c(sim.fd1a$pe[k], sim.fd1g$pe[k],
                                      sim.fd3a$pe[k], sim.fd3g$pe[k]),
                                  lower=c(sim.fd1a$lower[k], sim.fd1g$lower[k],
                                      sim.fd3a$lower[k], sim.fd3g$lower[k]),
                                  upper=c(sim.fd1a$upper[k], sim.fd1g$upper[k],
                                      sim.fd3a$upper[k], sim.fd3g$upper[k]),
                                  labels=c("Linear / no trend", "Linear / year effects",
                                      "Log-Log / no trend", "Log-Log / year effects  "),
                                  sublabels="diff GMM",
                                  entryheight=0.25,
                                  subentryheight=0.8,
                                  sublabelsX=.14,
                                  plot=1
                                  )

traceSysGMM1yearFD <- ropeladder(x=c(sim.fd1e$pe[k], sim.fd1h$pe[k],
                                     sim.fd3e$pe[k], sim.fd3h$pe[k]),
                                 lower=c(sim.fd1e$lower[k], sim.fd1h$lower[k],
                                     sim.fd3e$lower[k], sim.fd3h$lower[k]),
                                 upper=c(sim.fd1e$upper[k], sim.fd1h$upper[k],
                                     sim.fd3e$upper[k], sim.fd3h$upper[k]),
                                 sublabels="sys GMM",
                                 sublabelsX=.14,
                                 plot=1
                                 )

k <- 3

traceDiffGMM3yearFD <- ropeladder(x=c(sim.fd1a$pe[k], sim.fd1g$pe[k],
                                      sim.fd3a$pe[k], sim.fd3g$pe[k]),
                                  lower=c(sim.fd1a$lower[k], sim.fd1g$lower[k],
                                      sim.fd3a$lower[k], sim.fd3g$lower[k]),
                                  upper=c(sim.fd1a$upper[k], sim.fd1g$upper[k],
                                      sim.fd3a$upper[k], sim.fd3g$upper[k]),
                                  sublabelsX=-45,
                                  plot=2
                                  )

traceSysGMM3yearFD <- ropeladder(x=c(sim.fd1e$pe[k], sim.fd1h$pe[k],
                                     sim.fd3e$pe[k], sim.fd3h$pe[k]),
                                 lower=c(sim.fd1e$lower[k], sim.fd1h$lower[k],
                                     sim.fd3e$lower[k], sim.fd3h$lower[k]),
                                 upper=c(sim.fd1e$upper[k], sim.fd1h$upper[k],
                                     sim.fd3e$upper[k], sim.fd3h$upper[k]),
                                 sublabelsX=-45,
                                 plot=2
                                 )

vertmark <- linesTile(x=c(0,0), y=c(0,1), plot=1:2)

tile(traceDiffGMM1yearFD, traceSysGMM1yearFD, traceDiffGMM3yearFD, traceSysGMM3yearFD,
     vertmark,
     limits=rbind(c(-50,10), c(-50,10)),
     xaxis=list(at1=seq(-30,0,10), at2=seq(-40,0,10)),
     topaxis=list(at1=seq(-30,0,10), at2=seq(-40,0,10), add=TRUE),
     topaxistitle=list(labels1="Change in Packs pc, +1 years",
         labels2="Change in Packs pc, +3 years"),
     gridlines=list(type="xt"),
     output=list(file="FDrobustnessCigsHC", width=8.5)      
     )



# Make a ropeladder to show model robustness (with CIs) : RR

# helper function to convert risk ratios to percent change
rr2pc <- function(x) 100*(x-1)

k <- 1

traceDiffGMM1yearRR <- ropeladder(x=rr2pc(c(sim.rr1a$pe[k], sim.rr1g$pe[k],
                                      sim.rr3a$pe[k], sim.rr3g$pe[k])),
                                  lower=rr2pc(c(sim.rr1a$lower[k], sim.rr1g$lower[k],
                                      sim.rr3a$lower[k], sim.rr3g$lower[k])),
                                  upper=rr2pc(c(sim.rr1a$upper[k], sim.rr1g$upper[k],
                                      sim.rr3a$upper[k], sim.rr3g$upper[k])),
                                  labels=c("Linear / no trend", "Linear / year effects",
                                      "Log-Log / no trend", "Log-Log / year effects  "),
                                  sublabels="diff GMM",
                                  entryheight=0.25,
                                  subentryheight=0.8,
                                  sublabelsX=.14,
                                  plot=1
                                  )

traceSysGMM1yearRR <- ropeladder(x=rr2pc(c(sim.rr1e$pe[k], sim.rr1h$pe[k],
                                     sim.rr3e$pe[k], sim.rr3h$pe[k])),
                                 lower=rr2pc(c(sim.rr1e$lower[k], sim.rr1h$lower[k],
                                     sim.rr3e$lower[k], sim.rr3h$lower[k])),
                                 upper=rr2pc(c(sim.rr1e$upper[k], sim.rr1h$upper[k],
                                     sim.rr3e$upper[k], sim.rr3h$upper[k])),
                                 sublabels="sys GMM",
                                 sublabelsX=.14,
                                 plot=1
                                 )

k <- 3

traceDiffGMM3yearRR <- ropeladder(x=rr2pc(c(sim.rr1a$pe[k], sim.rr1g$pe[k],
                                      sim.rr3a$pe[k], sim.rr3g$pe[k])),
                                  lower=rr2pc(c(sim.rr1a$lower[k], sim.rr1g$lower[k],
                                      sim.rr3a$lower[k], sim.rr3g$lower[k])),
                                  upper=rr2pc(c(sim.rr1a$upper[k], sim.rr1g$upper[k],
                                      sim.rr3a$upper[k], sim.rr3g$upper[k])),
                                  sublabelsX=-45,
                                  plot=2
                                  )

traceSysGMM3yearRR <- ropeladder(x=rr2pc(c(sim.rr1e$pe[k], sim.rr1h$pe[k],
                                     sim.rr3e$pe[k], sim.rr3h$pe[k])),
                                 lower=rr2pc(c(sim.rr1e$lower[k], sim.rr1h$lower[k],
                                     sim.rr3e$lower[k], sim.rr3h$lower[k])),
                                 upper=rr2pc(c(sim.rr1e$upper[k], sim.rr1h$upper[k],
                                     sim.rr3e$upper[k], sim.rr3h$upper[k])),
                                 sublabelsX=-45,
                                 plot=2
                                 )

vertmark <- linesTile(x=c(0,0), y=c(0,1), plot=1:2)

tile(traceDiffGMM1yearRR, traceSysGMM1yearRR, traceDiffGMM3yearRR, traceSysGMM3yearRR,
     vertmark,
     limits=rbind(c(-50,10), c(-50,10)),
     xaxis=list(at1=seq(-30,0,10), at2=seq(-40,0,10)),
     topaxis=list(at1=seq(-30,0,10), at2=seq(-40,0,10), add=TRUE),
     topaxistitle=list(labels1="% Change in Packs pc, +1 years",
         labels2="% Change in Packs pc, +3 years"),
     gridlines=list(type="xt"),
     output=list(file="RRrobustnessCigsHC", width=8.5)      
     )



##### Supplement
##### Simulate year effects, all else equal, model 3g

allPeriods <- 1987:1995

# All else equal Forecast: Model 3g, all periods

# Make matrix of hypothetical x's: covariates
# Still use the 1g formula (no logs) -- we will handle logging manually
#  to get the differences of logs right
xhyp0.3g <- cfMake(formula.1g, datayearfe, length(allPeriods))

# Set everything to zero for a "static" model with no changes in covariates
xhyp0.3g$x <- xhyp0.3g.xpre <- 0*xhyp0.3g$x

# Extract the actual averages across states from the analyzed data;
xhypMat <- Reduce("+", pgmm.res3g$model)/length(pgmm.res3g$model)

# Pull the year effects from the above and insert as counterfactuals;
# check the indexes manually!
xhyp0.3g$x[,4:ncol(xhyp0.3g$x)] <- xhypMat[,5:ncol(xhypMat)]
xhyp0.3g$xpre[,4:ncol(xhyp0.3g$x)] <- xhypMat[,5:ncol(xhypMat)]

# We need a lag of the price per pack
lagY0.3g <- NULL # Hypothetical previous change in Y for simulation
for (i in 1:length(pgmm.res3g$model))
    lagY0.3g <- c(lagY0.3g, as.data.frame(pgmm.res3g$model[[i]])["1987",1])
lagY0.3g <- mean(lagY0.3g, na.rm=TRUE)

initialY0 <- mean(pdata$packpc[pdata$year==1986], na.rm=TRUE)
               # Hypothetical initial level of Y for simulation

# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev3g <- ldvsimev(xhyp0.3g,               # The matrix of hypothetical x's
                     simbetas.3g,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # NA indicates no constant!
                     phi=simphis.3g,            # estimated AR parameters; length must match lagY 
                     lagY=lagY0.3g,          # lags of y, most recent last
                     transform="difflog",   # "log" to undo log transformation,
                                         # "diff" to under first differencing
                                         # "difflog" to do both
                     initialY=initialY0   # for differenced models, the lag of the level of y
                     )

# Plot results using tile
colFits <- brewer.pal(3,"Set1")[2:3]

yrsTrace <- scatter(x=1987:1995,
                    y=sim.ev3g$pe,
                    ylower=sim.ev3g$lower,
                    yupper=sim.ev3g$upper,
                    pch=1, size=.8,
                    plot=1
                    )

yrsTraceWLS <- scatter(x=1987:1995,
                    y=sim.ev3g$pe,
                    ylower=sim.ev3g$lower,
                    yupper=sim.ev3g$upper,
                    fit=list(method="wls",
                        weights=1/sim.ev3g$se^2,
                        col=colFits[1]),
                    pch=1, size=.8,
                    plot=1
                    )

yrsTraceLOESS <- scatter(x=1987:1995,
                    y=sim.ev3g$pe,
                    ylower=sim.ev3g$lower,
                    yupper=sim.ev3g$upper,
                    fit=list(method="loess",
                        weights=1/sim.ev3g$se^2,
                        col=colFits[2]),
                    pch=1, size=.8,
                    plot=1
                    )

yrsLabel <- textTile(labels="EV and 95% CI",
                     x=1988.75, y=105,
                     plot=1)

wlsLabel <- textTile(labels="WLS trend and 95% CI",
                     x=1988.75, y=100, col=colFits[1],
                     plot=1)

loessLabel <- textTile(labels="Loess trend and 95% CI",
                     x=1988.75, y=100, col=colFits[2],
                     plot=1)

tile(yrsTrace, yrsLabel,
     output = list(file="yearEffectsModel3g", width=5.5),
     limits = c(1986.5,1995.5,95,140),
     yaxis=list(major=FALSE),
     xaxistitle = list(labels="Year"),
     yaxistitle = list(labels="Trend due to year effects only", x=0),
     height=list(plot="golden"),
     width=list(null=5, # widen plot area for visibility
         yaxis.labelspace=0)
     )


tile(yrsTraceWLS, yrsLabel, wlsLabel,
     output = list(file="yearEffectsModel3gWLS", width=5.5),
     limits = c(1986.5,1995.5,95,140),
     yaxis=list(major=FALSE),
     xaxistitle = list(labels="Year"),
     yaxistitle = list(labels="Trend due to year effects only", x=0),
     height=list(plot="golden"),
     width=list(null=5, # widen plot area for visibility
         yaxis.labelspace=0)
     )


tile(yrsTraceLOESS, yrsLabel, loessLabel,
     output = list(file="yearEffectsModel3gLOESS", width=5.5),
     limits = c(1986.5,1995.5,95,140),
     yaxis=list(major=FALSE),
     xaxistitle = list(labels="Year"),
     yaxistitle = list(labels="Trend due to year effects only", x=0),
     height=list(plot="golden"),
     width=list(null=5, # widen plot area for visibility
         yaxis.labelspace=0)
     )


