###################################################################
# CSSS 512
# Lab Session 7 - In Sample Simulation for Panel Data Models
# 5/25/18

# Christopher Adolph   faculty.washington.edu/cadolph
###################################################################
# Clear memory
rm(list=ls())
setwd("~/desktop")

# Load libraries
library(plm)            # Econometrics package for linear panel models
library(simcf)          # For panel functions and simulators
library(tile)           # For visualization of model inference
library(RColorBrewer)   # For nice colors
library(MASS)           # For mvrnorm()
library(abind)          # For combining matrices into arrays
source("helperCigs.R")  # For processing and graphics functions

# Load cigarette consumption data (Jonathan Gruber, MIT)
# Variables (see codebook):
# state	year	cpi	pop	packpc	income	tax	avgprs	taxs
data <- read.csv("cigarette.csv")

library(Ecdat)
help(Cigarette)

# Quick inflation adjustment to 1995 dollars
inflAdjust <- function(x,cpi,year,target) {
    unique(cpi[year==target])*x/cpi		#Multiply x with cpi in target year then divide by cpi in observed year
}
    
data$income95 <- with(data, inflAdjust(income, cpi, year, 1995))	#Make adjustments to state personal income
data$tax95 <- with(data, inflAdjust(tax, cpi, year, 1995))			#Average state, federal, and average local excise taxes
data$avgprs95 <- with(data, inflAdjust(avgprs, cpi, year, 1995))	#Average price, including sales taxes
data$taxs95 <- with(data, inflAdjust(taxs, cpi, year, 1995))		#Average excise taxes, including sales taxes

# Create per capita income (in k)
data$income95pc <- data$income95/data$pop

# Create pretax price, 1995 dollars
data$pretax95 <- data$avgprs95 - data$taxs95

data[1:15,]
attach(data)

########################################################################################################
# Construct the weighted avg of consumption over states

stateList <- unique(as.character(data$state))	#Store a vector of the state names
nstate <- length(stateList)						#Store the number of states
yearList <- unique(data$year)					#Store the years
nyear <- length(yearList)						#Store the number of years
packpcMat <- taxs95Mat <- popMat <- NULL		#Create empty matrices for key variables

stateList
nstate
yearList
nyear

for (istate in 1:nstate) {						#Start a loop over the number of states
    curstate <- stateList[istate]				#Store the current state
    packpcMat <- cbind(packpcMat, data$packpc[as.character(data$state)==curstate])	#Store the packpc values of the current state to packpcMat
    taxs95Mat <- cbind(taxs95Mat, data$taxs95[as.character(data$state)==curstate])	#Do the same for taxs95
    popMat <- cbind(popMat, data$pop[as.character(data$state)==curstate])			#Do the same for population
}

packpcMat
taxs95Mat
popMat

packpcWAvg <- apply(packpcMat*popMat, 1, sum)/apply(popMat, 1, sum)	#Find the cumulative sum of packs per capita*population across states, then divide by the cumulative population for each year

packpcWAvg

taxs95WAvg <- apply(taxs95Mat*popMat, 1, sum)/apply(popMat, 1, sum) #Find the cumulative sum of taxes*population across states, then divide by the cumulative population for each year

taxs95WAvg

# Make a list of complete state names using built in objects
stateFull <- sapply(stateList, function(x) {state.name[grep(x, state.abb)]})	#Match the state's full name with stateList

# First, create a plm data frame (special data frame that "knows" the
# unit variable and time variable
pdata <- pdata.frame(data, index=c("state", "year"))
########################################################################################################
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
pgmmformula.3a <- log(packpc) ~ lag(log(packpc), 1) + log(income95pc) + log(avgprs95) | lag(log(packpc), 2:99)


# Model 3: Elasticity specification

# Try difference GMM with two way effects
pgmm.res3g <- pgmm(pgmmformula.3a,
                   data = pdata,
                   effect = "twoways",   # should consider two-way for small T
                   transformation = "d")    # should do ld if T=3

summary(pgmm.res3g)
# Good Sargan test, Good AR(2) test, Wald supports 2-way



########################################################################################################
# Simulate conditional forecasts from model 3g

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
yearfe <- makeFEdummies(pdata$year)								# Construct the dummies for each year
yearfe <- yearfe[,3:ncol(yearfe)] 								# Why drop first 2 col's?	
yearlist <- unique(pdata$year)									# List all the years
yearlist <- yearlist[3:length(yearlist)]						# List the years less the first two
colnames(yearfe) <- paste0("y",yearlist)						# Create names for the year dummies

# Construct formulas -- with year dummies (1g)
formula <- "packpc ~ income95pc + avgprs95 -1"					#with Income and Price as covariates
datayearfe <- cbind(pdata,yearfe)								#Combine pdata variables with the year dummies 
colnames(datayearfe)

yearfenames <- NULL
for (i in 1:ncol(yearfe)) {
  formula <- paste0(formula,"+ y",yearlist[i]," ")				#Add the year dummies to the initial formula
  yearfenames <- c(yearfenames,paste0("y",yearlist[i]))			#Make a vector of names for the years
}
names(datayearfe) <- c(names(data),yearfenames)
formula.1g <- as.formula(formula)
formula.1g

########################################################################################################
# Forecast: Model 3g, +60

# Recall model 3g: log(packpc) ~ lag(log(packpc), 1) + log(income95pc) + log(avgprs95) | lag(log(packpc), 2:99)
# log log Difference GMM with state and year fixed effects

simparam.3g <- mvrnorm(sims, coefficients(pgmm.res3g), vcovHC(pgmm.res3g))		#Sample parameters from multivariate normal distribution
simphis.3g <- simparam.3g[,1]													#Extract the phis
simbetas.3g <- simparam.3g[,2:ncol(simparam.3g)]								#Extract the betas

# Make matrix of hypothetical x's:
# Assume an average state raised taxes 60 cents starting 1996
#
# Issues -- we need to somehow include the state and year FEs:
#           Let's set the state to be an "average" state in 1995,
#           and year to be like the last year (1995)

# Make matrix of hypothetical x's: covariates
# Still use the 1g formula (no logs) -- we will handle logging manually
#  to get the differences of logs right
xhyp.3g <- cfMake(formula.1g, datayearfe, periods.out)			#Create matrix of hypothetical x's at their means
xhyp.3g

# pgmm uses covariates in differenced form
# so we want most of them to be 0 (no change)
# exceptions:
# (1) changes in covariates of interest
# (2) time dummies aren't differenced


changeTax <- c(60,0,0)											#Difference GMM requires hypothetical x's as changes

# Need log version of differenced key covariate (doubling tax in avg state)
meanPrice <- mean(pdata$avgprs, na.rm=TRUE)						#Find the mean average price
meanIncome <- mean(pdata$income95pc, na.rm=TRUE)				#Find the mean income per capita

xhyp.3g <- cfChange(xhyp.3g, "avgprs95",						#Change avg price to be the log of the change in mean price for 1996
                    x=log(meanPrice+changeTax[1]) - log(meanPrice),
                    xpre=log(meanPrice) - log(meanPrice),
                    scen=1)

xhyp.3g <- cfChange(xhyp.3g, "avgprs95",						#Change avg price to be the log of the change in mean price for 1997
                    x=log(meanPrice+changeTax[2]) - log(meanPrice),
                    xpre=log(meanPrice) - log(meanPrice),
                    scen=2)

xhyp.3g <- cfChange(xhyp.3g, "avgprs95",						#Change avg price to be the log of the change in mean price for 1998
                    x=log(meanPrice+changeTax[3]) - log(meanPrice),
                    xpre=log(meanPrice) - log(meanPrice),
                    scen=3)

xhyp.3g <- cfChange(xhyp.3g, "income95pc",						#Change income per capita to be the log of the change in mean income for 1996-1998
                    x=log(meanIncome) - log(meanIncome),
                    xpre=log(meanIncome) - log(meanIncome),
                    scen=1:3)
xhyp.3g

for (iyear in 1985:1991)
    xhyp.3g <- cfChange(xhyp.3g, paste0("y",iyear), x=0, xpre=0, scen=1:3)	#Set the change in year fixed effects for 1985-1991 to be zero

xhyp.3g <- cfChange(xhyp.3g, "y1992", x=-1, xpre=-1, scen=1)	#Set the change in the 1992 dummy to be -1
xhyp.3g <- cfChange(xhyp.3g, "y1992", x=0, xpre=0, scen=2)
xhyp.3g <- cfChange(xhyp.3g, "y1992", x=0, xpre=0, scen=3)

xhyp.3g <- cfChange(xhyp.3g, "y1993", x=1, xpre=1, scen=1)		#Set the change in 1993 dummy to be 1
xhyp.3g <- cfChange(xhyp.3g, "y1993", x=-1, xpre=-1, scen=2)	#Set the change in 1993 dummy to be -1
xhyp.3g <- cfChange(xhyp.3g, "y1993", x=0, xpre=0, scen=3)

xhyp.3g <- cfChange(xhyp.3g, "y1994", x=0, xpre=0, scen=1)
xhyp.3g <- cfChange(xhyp.3g, "y1994", x=1, xpre=1, scen=2)		#Set the change in 1994 dummy to be 1
xhyp.3g <- cfChange(xhyp.3g, "y1994", x=-1, xpre=-1, scen=3)	#Set the change in 1994 dummy to be -1

xhyp.3g <- cfChange(xhyp.3g, "y1995", x=0, xpre=0, scen=1)		
xhyp.3g <- cfChange(xhyp.3g, "y1995", x=0, xpre=0, scen=2)
xhyp.3g <- cfChange(xhyp.3g, "y1995", x=1, xpre=1, scen=3)		#Set the change in the 1995 dummy to be 1

# We can "ignore" the state fixed effects for now and add them later
# because model is total linear

# Create baseline scenario
xbase.3g <- xhyp.3g
xbase.3g$x <- xbase.3g$xpre

# We need a lag of the price per pack
lagY.3g <- NULL # Hypothetical previous change in Y for simulation

for (iunit in 1:length(pgmm.res3g$model))
    lagY.3g <- c(lagY.3g, as.data.frame(pgmm.res3g$model[[iunit]])["1995",1])	#Store the changes in log(packpc) in 1995 for all states
pgmm.res3g$model
lagY.3g <- mean(lagY.3g, na.rm=TRUE)											#Find the mean


# Hypothetical initial level of Y for simulation
initialY <- mean(pdata$packpc[pdata$year==1993], na.rm=TRUE)					#Store the mean packpc for 1993
               

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


########################################################################################################
# Forecasts: avg state, by state, and weighted avg of all states
#         for +.6 and raise to max scenarios
#         outcomes of interest: packs FD, packs RR/%change,
#                               revenue pc FD, revenue mils FD, revenue %GSP FD


# Visuals: line plots.  1x2 and 1x3 for two diff scenarios (x2)

# List for results
simEV.3g.scen1 <- simBASE.3g.scen1 <- simFD.3g.scen1 <- simRR.3g.scen1 <- list()		#Empty lists for results

# All states modeled
allstates <- names(pgmm.res3g$model)		#Names of all states considered in the model

# Storage vector for weights
popweight <- rep(NA, length(allstates))		#Temporary empty vector that will store populations

# Loop over states in model (we want to generate forecasts for each state separately)
for (i in 1:length(pgmm.res3g$model)) {			#loop from 1 to the number of states

    curstate <- allstates[i]		#refer to the current state being considered

    # Grab population weight for later
    popweight[i] <- pdata$pop[(pdata$year==1995)&(pdata$state==curstate)]	#Find the population in 1995 for the current state
    
    # Make matrix of hypothetical x's: covariates
    # Still use the 1g formula (no logs) -- we will handle logging manually
    #  to get the differences of logs right
    xhyp.3g <- cfMake(formula.1g, datayearfe, periods.out)					#Create the matrix of hypothetical x values at their mean

    # pgmm uses covariates in differenced form
    # so we want most of them to be 0 (no change)
    # exceptions:
    # (1) changes in covariates of interest
    # (2) time dummies aren't differenced
    #xhyp.3g$x <- xhyp.3g$xpre <- 0*xhyp.3g$x

    # Need log version of differenced key covariate (doubling tax in avg state)
    curPrice92 <- pdata$avgprs95[(pdata$state==curstate)&(pdata$year==1992)]		#Find the average price of the current state for 1992		
    curPrice93 <- pdata$avgprs95[(pdata$state==curstate)&(pdata$year==1993)]		#For 1993, etc.
    curPrice94 <- pdata$avgprs95[(pdata$state==curstate)&(pdata$year==1994)]
    curPrice95 <- pdata$avgprs95[(pdata$state==curstate)&(pdata$year==1995)]
    curIncome92 <-  pdata$income95pc[(pdata$state==curstate)&(pdata$year==1992)]		#Find the per capita income for the current state for 1992
    curIncome93 <-  pdata$income95pc[(pdata$state==curstate)&(pdata$year==1993)]		#For 1993, etc.
    curIncome94 <-  pdata$income95pc[(pdata$state==curstate)&(pdata$year==1994)]
    curIncome95 <-  pdata$income95pc[(pdata$state==curstate)&(pdata$year==1995)]

    changeTax <- c(60, 0, 0)		#Set the hypothetical change in tax we want to observe
    
    #Make changes to the initial counterfactual matrix by finding the log differences in prices and income
    xhyp.3g <- cfChange(xhyp.3g, "avgprs95",
                        x=log(curPrice93+changeTax[1]) - log(curPrice92),
                        xpre=log(curPrice93) - log(curPrice92),
                        scen=1)

    xhyp.3g <- cfChange(xhyp.3g, "avgprs95",
                        x=log(curPrice94+changeTax[2]) - log(curPrice93),
                        xpre=log(curPrice94) - log(curPrice93),
                        scen=2)
    
    xhyp.3g <- cfChange(xhyp.3g, "avgprs95",
                        x=log(curPrice95+changeTax[3]) - log(curPrice94),
                        xpre=log(curPrice95) - log(curPrice94),
                        scen=3)

    xhyp.3g <- cfChange(xhyp.3g, "income95pc",
                        x=log(curIncome93) - log(curIncome92),
                        xpre=log(curIncome93) - log(curIncome92),
                        scen=1)
    
    xhyp.3g <- cfChange(xhyp.3g, "income95pc",
                        x=log(curIncome94) - log(curIncome93),
                        xpre=log(curIncome94) - log(curIncome93),
                        scen=2)
    
    xhyp.3g <- cfChange(xhyp.3g, "income95pc",
                        x=log(curIncome95) - log(curIncome94),
                        xpre=log(curIncome95) - log(curIncome94),
                        scen=3)
    
    
    #Make the same changes to the year fixed effects as earlier
    for (iyear in 1985:1991)
        xhyp.3g <- cfChange(xhyp.3g, paste0("y",iyear), x=0, xpre=0, scen=1:3)
    
    xhyp.3g <- cfChange(xhyp.3g, "y1992", x=-1, xpre=-1, scen=1)
    xhyp.3g <- cfChange(xhyp.3g, "y1992", x=0, xpre=0, scen=2)
    xhyp.3g <- cfChange(xhyp.3g, "y1992", x=0, xpre=0, scen=3)
    
    xhyp.3g <- cfChange(xhyp.3g, "y1993", x=1, xpre=1, scen=1)
    xhyp.3g <- cfChange(xhyp.3g, "y1993", x=-1, xpre=-1, scen=2)
    xhyp.3g <- cfChange(xhyp.3g, "y1993", x=0, xpre=0, scen=3)
    
    xhyp.3g <- cfChange(xhyp.3g, "y1994", x=0, xpre=0, scen=1)
    xhyp.3g <- cfChange(xhyp.3g, "y1994", x=1, xpre=1, scen=2)
    xhyp.3g <- cfChange(xhyp.3g, "y1994", x=-1, xpre=-1, scen=3)
    
    xhyp.3g <- cfChange(xhyp.3g, "y1995", x=0, xpre=0, scen=1)
    xhyp.3g <- cfChange(xhyp.3g, "y1995", x=0, xpre=0, scen=2)
    xhyp.3g <- cfChange(xhyp.3g, "y1995", x=1, xpre=1, scen=3)

    # Baseline scenario
    xbase.3g <- xhyp.3g
    xbase.3g$x <- xbase.3g$xpre
    
    # We need a lag of the price per pack
    lagY.3g <- as.data.frame(pgmm.res3g$model[[i]])["1992",1]

    # This restores the fixed effect by state
    initialY <- pdata$packpc[(pdata$state==curstate)&(pdata$year==1992)]

    
    # Simulate expected values of Y (on original level scale)
    # out to periods.out given hypothetical future values of X,
    # initial lags of the change in Y, and an initial level of Y
    simState.ev3g <- ldvsimev(xhyp.3g,               # The matrix of hypothetical x's
                              simbetas.3g,           # The matrix of simulated betas
                              ci=0.95,            # Desired confidence interval
                              constant=NA,        # NA indicates no constant!
                              phi=simphis.3g,            # estimated AR parameters; length must match lagY 
                              lagY=lagY.3g,          # lags of y, most recent last
                              transform="difflog",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                              initialY=initialY,   # for differenced models, the lag of the level of y
                              simulates=TRUE
                              )

    # Simulate expected values of Y given no change in covariates
    simState.base3g <- ldvsimev(xbase.3g,               # The matrix of hypothetical x's
                                simbetas.3g,           # The matrix of simulated betas
                                ci=0.95,            # Desired confidence interval
                                constant=NA,        # NA indicates no constant!
                                phi=simphis.3g,            # estimated AR parameters; length must match lagY 
                                lagY=lagY.3g,          # lags of y, most recent last
                                transform="difflog",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                                initialY=initialY,   # for differenced models, the lag of the level of y
                                simulates=TRUE
                                ) 
    
    # Simulate first differences in y (as cumulated changes)
    # out to periods.out given hypothetical future values of x, xpre,
    # and initial lags of the change in y
    simState.fd3g <- ldvsimfd(xhyp.3g,            # The matrix of hypothetical x's
                              simbetas.3g,        # The matrix of simulated betas
                              ci=0.95,            # Desired confidence interval
                              constant=NA,        # Column containing the constant
                                        # set to NA for no constant
                              phi=simphis.3g,     # estimated AR parameters; length must match lagY 
                              lagY=lagY.3g,       # lags of y, most recent last
                               transform="difflog",# Model is differenced logs
                              initialY=initialY,  # Required in this case (fd of differenced log Y)
                              simulates=TRUE
                              )

    # Simulate relative risks in y
    # out to periods.out given hypothetical future values of x, xpre,
    # and initial lags of the change in y
    simState.rr3g <- ldvsimrr(xhyp.3g,            # The matrix of hypothetical x's
                              simbetas.3g,        # The matrix of simulated betas
                              ci=0.95,            # Desired confidence interval
                              constant=NA,        # Column containing the constant
                                        # set to NA for no constant
                              phi=simphis.3g,     # estimated AR parameters; length must match lagY 
                              lagY=lagY.3g,       # lags of y, most recent last
                              transform="difflog",# Model is differenced logs
                              initialY=initialY,  # Required for differenced Y in ldvsimrr
                              simulates=TRUE
                              )


    # Collect results (only works for 1 CI)
    simEV.3g.scen1 <- collectUnitSims(simEV.3g.scen1, simState.ev3g)		#Store the simulated expected values for each state
    simBASE.3g.scen1 <- collectUnitSims(simBASE.3g.scen1, simState.base3g)
    simFD.3g.scen1 <- collectUnitSims(simFD.3g.scen1, simState.fd3g)
    simRR.3g.scen1 <- collectUnitSims(simRR.3g.scen1, simState.rr3g)
}

# Compute weighted average across states
simEV.3g.scen1 <- aggregateUnitSims(simEV.3g.scen1, weighted.mean, w=popweight)	#Store the weighted averages of the simulated 
simBASE.3g.scen1 <- aggregateUnitSims(simBASE.3g.scen1, weighted.mean, w=popweight)
simFD.3g.scen1 <- aggregateUnitSims(simFD.3g.scen1, weighted.mean, w=popweight)
simRR.3g.scen1 <- aggregateUnitSims(simRR.3g.scen1, weighted.mean, w=popweight)


########################################################################################################
# Plot the FD and RR results using a custom tile-graphics function

# Nice colors
col <- brewer.pal(4,"Dark2")[3:4]

# Other settings common to all plots
avg <- "weighted.mean"
periods <- 1993:1995
limits <- c(1991.9, 1995.25, -43, 2)
xat <- c(1992, 1993, 1994, 1995)
plottitles <- c("Change, Packs pc after +$.60 tax",
                "Percent Change, Packs pc after +$.60 tax")
xtitles <- c("Forecast Year", "Forecast Year")
ytitles <- c("","")

# Plot the results by state
cigUnitLineplots(FD=simFD.3g.scen1, RR=simRR.3g.scen1,
                 periods=periods, units=stateList,
                 showUnits=TRUE, showMean=FALSE,
                 avg=avg, avgLab="US",
                 unitLabX=1995.125, labSet="tb", labK=1,
                 initialZero=TRUE, CI=TRUE,
                 limitsFD=limits, limitsRR=limits, xat=xat,
                 col=col, dots=TRUE, RRasPD=TRUE,
                 plottitles=plottitles,
                 xtitles=xtitles, ytitles=ytitles,
                 file="m3gUnitsFDRR")

# Plot the results by state with weighted mean
cigUnitLineplots(FD=simFD.3g.scen1, RR=simRR.3g.scen1,
                 periods=periods, units=stateList,
                 showUnits=TRUE, showMean=TRUE,
                 avg=avg, avgLab="US",
                 unitLabX=1995.125, labSet="tb", labK=1,
                 initialZero=TRUE, CI=FALSE,
                 limitsFD=limits, limitsRR=limits, xat=xat,
                 col=col, dots=TRUE, RRasPD=TRUE,
                 plottitles=plottitles,
                 xtitles=xtitles, ytitles=ytitles,
                 file="m3gUnitsMeanFDRR")

# Plot the results by state with weighted mean
cigUnitLineplots(FD=simFD.3g.scen1, RR=simRR.3g.scen1,
                 periods=periods, units=stateList,
                 showUnits=FALSE, showMean=TRUE,
                 avg=avg, avgLab="US",
                 unitLabX=1995.125, labSet="tb", labK=1,
                 initialZero=TRUE, CI=TRUE,
                 limitsFD=limits, limitsRR=limits, xat=xat,
                 col=col, dots=TRUE, RRasPD=TRUE,
                 plottitles=plottitles,
                 xtitles=xtitles, ytitles=ytitles,
                 file="m3gMeanFDRR")

# Plot the representative state results
cigUnitLineplots(FD=simFD.3g.scen1, RR=simRR.3g.scen1,
                 compFD=sim.fd3g, compRR=sim.rr3g,
                 periods=periods, units=stateList,
                 showUnits=FALSE, showMean=FALSE,
                 avg=avg, avgLab="Wgt", compLab="Rep",
                 unitLabX=1995.15, labSet="tb", labK=1,
                 initialZero=TRUE, CI=TRUE, compCI=TRUE,
                 limitsFD=limits, limitsRR=limits, xat=xat,
                 col=col, dots=TRUE, RRasPD=TRUE,
                 plottitles=plottitles,
                 xtitles=xtitles, ytitles=ytitles,
                 file="m3gCompFDRR")
             
# Plot the comparison of aggregate results and representative results
cigUnitLineplots(FD=simFD.3g.scen1, RR=simRR.3g.scen1,
                 compFD=sim.fd3g, compRR=sim.rr3g,
                 periods=periods, units=stateList,
                 showUnits=FALSE, showMean=TRUE,
                 avg=avg, avgLab="Wgt", compLab="Rep",
                 unitLabX=1995.15, labSet="tb", labK=1,
                 initialZero=TRUE, CI=TRUE, compCI=TRUE,
                 limitsFD=limits, limitsRR=limits, xat=xat,
                 col=col, dots=TRUE, RRasPD=TRUE,
                 plottitles=plottitles,
                 xtitles=xtitles, ytitles=ytitles,
                 file="m3gMeanCompFDRR")
             

# Plot dotplot of results for final period by state: all states
cigUnitDotplot(x=simFD.3g.scen1$units,
               units=stateFull, period=3,
               showRank=TRUE,
               vertmark=0,
               limits=c(-62, 5),
               labspace=0.35,
               at=seq(-60,0,10),
               col=col[1],
               xtitle="Change, Packs pc 3 years after +$.60 tax",
               toptitle="Change, Packs pc 3 years after +$.60 tax",
               file="m3gUnitFDdotplotALL")

# Plot dotplot of results for final period by state: top half
cigUnitDotplot(x=simFD.3g.scen1$units,
               units=stateFull, period=3,
               showRank=TRUE,
               showSelected=1:24,
               vertmark=0,
               limits=c(-62, 5),
               labspace=0.35,
               at=seq(-60,0,10),
               col=col[1],
               xtitle="Change, Packs pc 3 years after +$.60 tax",
               toptitle="Change, Packs pc 3 years after +$.60 tax",
               file="m3gUnitFDdotplotTOP")

# Plot dotplot of results for final period by state: bottomhalf
cigUnitDotplot(x=simFD.3g.scen1$units,
               units=stateFull, period=3,
               showRank=TRUE,
               showSelected=25:48,
               vertmark=0,
               limits=c(-62, 5),
               labspace=0.35,
               at=seq(-60,0,10),
               col=col[1],
               xtitle="Change, Packs pc 3 years after +$.60 tax",
               toptitle="Change, Packs pc 3 years after +$.60 tax",
               file="m3gUnitFDdotplotBOTTOM")
