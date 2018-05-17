###########################################################################
#CSSS 512
#Lab Session 4 - Modeling Nonstationary Time Series
#4/27/18

###########################################################################
rm(list=ls())

#Load libraries
library(tseries)           # For unit root tests
library(forecast)		   # For decompose()
library(lmtest)            # For Breusch-Godfrey LM test of serial correlation
library(urca)              # For estimating cointegration models
library(simcf)             # For counterfactual simulation via ldvsimev()
library(MASS)              # For mvrnorm()
library(RColorBrewer)      # For nice colors
library(Zelig)			   # For approval data
library(quantmod)		   # For creating lags
source("TSplotHelper.R")   # Helper function for counterfactual time series plots

###########################################################################
#Load data
#US Presidential approval data (Bush, Monthly, 2/2001--6/2006)
#Includes average oil price data ($/barrel?)
#
#Variable names:  month year approve disapprove unsure
#                 sept.oct.2001 iraq.war avg.price

data(approval)
attach(approval)

# Look at the time series
#pdf("tsapproval.pdf",width=6,height=3.25)
plot(approve,type="l",ylab="Percent Approving",xlab="Time",
     main = "US Presidential Approval")
lines(x=c(8,8),y=c(-1000,1000),col="red")
lines(x=c(26,26),y=c(-1000,1000),col="blue")
text("9/11",x = 10, y = 40, col="red",cex=0.7)
text("Iraq \n War",x = 28, y = 40, col="blue",cex=0.7) 
#dev.off()

#pdf("tsprice.pdf",width=6,height=3.25)
plot(avg.price,type="l",ylab="$ per Barrel",xlab="Time",
     main = "Average Price of Oil")
lines(x=c(8,8),y=c(-1000,1000),col="red")
lines(x=c(26,26),y=c(-1000,1000),col="blue")
text("9/11",x = 10, y = 175, col="red",cex=0.7)
text("Iraq \n War",x = 28, y = 175, col="blue",cex=0.7) 
#dev.off()

#Look at the ACF
#pdf("acfapproval.pdf",width=6,height=3.25)
acf(approve)
#dev.off()

#pdf("acfprice.pdf",width=6,height=3.25)
acf(avg.price)
#dev.off()

#Look at the PACF
#pdf("pacfapproval.pdf",width=6,height=3.25)
pacf(approve)
#dev.off()

#pdf("pacfprice.pdf",width=6,height=3.25)
pacf(avg.price)
#dev.off()

#Look at the decomposed time series
plot(decompose(ts(approve[12:59],freq=12)))
plot(decompose(ts(avg.price[12:59],freq=12)))


#Check for seasonality in approval
col <- brewer.pal(8, "RdYlGn")

#Gather the data (sort the number of deaths by month and year in a matrix)
appmat <- matrix(approve[12:59],nrow=12,ncol=length(approve[12:59])/12, byrow=FALSE)

#Repeat them as many times as needed
col <-  as.vector(t(matrix(col, nrow=length(col), ncol=ceiling(ncol(appmat)/length(col)))))

#Plot each year over the months
matplot(appmat, type="l", col=col, lty=1, xaxt="n", ylab="Percent Approving", xlab="Month",
        main=expression(paste("Monthly View of Presidential Approval in the US, 2002-2005")))
axis(1, at=1:12, labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
abline(a=0,b=0,lty="dashed")


#Check for seasonality in avg price
col <- brewer.pal(8, "RdYlBu")


pricemat <- matrix(avg.price[12:59],nrow=12,ncol=length(avg.price[12:59])/12, byrow=FALSE)
col <-  as.vector(t(matrix(col, nrow=length(col), ncol=ceiling(ncol(pricemat)/length(col)))))
matplot(pricemat, type="l", col=col, lty=1, xaxt="n", ylab="$ per Barrel", xlab="Month",
        main=expression(paste("Monthly View of Average Oil Price, 2002-2005")))
axis(1, at=1:12, labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
abline(a=0,b=0,lty="dashed")


#Check for a unit root in approval
PP.test(approve)
adf.test(approve)

#Check for a unit root in average price
PP.test(avg.price)
adf.test(avg.price)

###########################################################################
#Investigate some other (potentially) non-stationary time series data

#Simulated random walk
set.seed(1)
phony <- rnorm(length(approve))
for (i in 2:length(phony)){
    phony[i] <- phony[i-1] + rnorm(1) 
}

#Plot the data
plot(phony, type="l", col="red", ylab="y",xlab="Time",
     main = "Simulated Random Walk")

#Check the ACF
acf(phony)

#Check the PACF
pacf(phony)

#Check for seasonality 
col <- brewer.pal(8, "RdYlGn")
phonymat <- matrix(phony[12:59],nrow=12,ncol=length(phony[12:59])/12, byrow=FALSE)
matplot(phonymat, type="l", col=col, lty=1, xaxt="n", ylab="y", xlab="Time",
        main=expression(paste("Monthly View of Simulated Random Walk")))
axis(1, at=1:12, labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
abline(a=0,b=0,lty="dashed")

PP.test(phony)
adf.test(phony)


#US unemployment rate, Watson (2014)
unemployment <- read.csv("unemployment.csv", header=TRUE)

date <- unemployment$date
rate <- unemployment$unemployment_rate

#Plot the data
plot(date, rate, xlab="Month", ylab="Unempoyment Rate", main=expression(paste("Monthly Unemployment Rate in the United States, 1948-2013")))

#Check the ACF
acf(rate)

#Check the PACF
pacf(rate)

#Check for seasonality
col <- brewer.pal(8, "Blues")
unempmat <- matrix(rate[1:780],nrow=12,ncol=length(rate[1:780])/12, byrow=FALSE)
matplot(unempmat, type="l", col=col, lty=1, xaxt="n", ylab="y", xlab="Time",
        main=expression(paste("Monthly View of Unemployment in the US, 1948-2013")))
axis(1, at=1:12, labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
abline(a=0,b=0,lty="dashed")

#Perform a unit root test
PP.test(rate)
adf.test(rate)


#Crude oil prices, Hamilton (2008)
crude <- read.csv("crude_oil.csv", header=TRUE)

price <- crude$spot_price_fob

#Plot the data
plot(price, type="l", xlab="Month", ylab="Spot Price FOB ($ per Barrel)", main=expression(paste("Crude Oil Prices, 1986-2008")))
lines(x=c(189,190),y=c(-1000,1000),col="red")
lines(x=c(207,208),y=c(-1000,1000),col="blue")
text("9/11",x = 182, y = 40, col="red",cex=0.7)
text("Iraq \n War",x = 200, y = 40, col="blue",cex=0.7) 

#Check the ACF
acf(price)

#Check the PACF
pacf(price)

#Check for seasonality
col <- brewer.pal(8, "Reds")
oilmat <- matrix(price[1:264],nrow=12,ncol=length(rate[1:264])/12, byrow=FALSE)
#col <-  as.vector(t(matrix(col, nrow=length(col), ncol=ceiling(ncol(price)/length(col)))))
matplot(oilmat, type="l", col=col, lty=1, xaxt="n", ylab="y", xlab="Month",
        main=expression(paste("Monthly View of Crude Oil Prices, 1986-2008")))
axis(1, at=1:12, labels=c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
abline(a=0,b=0,lty="dashed")

#Perform a unit root test
PP.test(price)
adf.test(price)

###########################################################################
#Differencing time series data

#Back to the Presidential Approval dataset
#Consider the first difference of each variable
approveLag <- c(NA, approve[1:(length(approve)-1)])
approveLag2 <- as.vector(Lag(approve,k=1))

approveDiff <- approve - approveLag
approveLagDiff <- cbind(approve, approveLag, approveDiff)
approveLagDiff

avg.priceLag <- c(NA, avg.price[1:(length(avg.price)-1)])
avg.priceDiff <- avg.price - avg.priceLag

avg.priceLagDiff <- cbind(avg.price, avg.priceLag, avg.priceDiff)
avg.priceLagDiff


#Consider the second difference of each variable
approve2Lag <- c(NA, NA, approve[2:length(approve)-2])
approve2Lag2 <- as.vector(Lag(approve,k=2))

approve2Diff <- approve - approve2Lag
approveLagDiff <- cbind(approve, approveLag, approve2Lag, approveDiff, approve2Diff)
approveLagDiff

avg.price2Lag <- c(NA, NA, avg.price[2:length(avg.price)-2])
avg.price2Lag2 <- as.vector(Lag(avg.price,k=2))

avg.price2Diff <- avg.price - avg.price2Lag
avg.priceLagDiff <- cbind(avg.price, avg.priceLag, avg.price2Lag, avg.priceDiff, avg.price2Diff)
avg.priceLagDiff


#Look at the DIFFERENCED time series
#pdf("tsapprovalDiff.pdf",width=6,height=3.25)
plot(approveDiff,type="l",ylab="Change in Percent Approving",xlab="Time",
     main = "US Presidential Approval")
lines(x=c(8,8),y=c(-1000,1000),col="red")
lines(x=c(26,26),y=c(-1000,1000),col="blue")
text("9/11",x = 10, y = 15, col="red",cex=0.7)
text("Iraq \n War",x = 28, y = 15, col="blue",cex=0.7) 
#dev.off()

#pdf("tspriceDiff.pdf",width=6,height=3.25)
plot(avg.priceDiff,type="l",ylab="Change in $ per Barrel",xlab="Time",
     main = "Average Price of Oil")
lines(x=c(8,8),y=c(-1000,1000),col="red")
lines(x=c(26,26),y=c(-1000,1000),col="blue")
text("9/11",x = 10, y = -30, col="red",cex=0.7)
text("Iraq \n War",x = 28, y = -30, col="blue",cex=0.7) 
#dev.off()

#Look at the new ACF and PACF for approval
#pdf("acfapprovalDiff.pdf",width=6,height=3.25)
acf(approveDiff, na.action=na.pass)
#dev.off()

#pdf("pacfapprovalDiff.pdf",width=6,height=3.25)
pacf(approveDiff, na.action=na.pass)
#dev.off()


#Look at the new ACF and PACF for oil prices
#pdf("acfpriceDiff.pdf",width=6,height=3.25)
acf(avg.priceDiff, na.action=na.pass)
#dev.off()

#pdf("pacfpriceDiff.pdf",width=6,height=3.25)
pacf(avg.priceDiff, na.action=na.pass)
#dev.off()

# Check for a unit root in differenced time series

PP.test(as.vector(na.omit(approveDiff)))
adf.test(na.omit(approveDiff))

PP.test(as.vector(na.omit(avg.priceDiff)))
adf.test(na.omit(avg.priceDiff))


unempLag <- c(NA, rate[1:(length(rate)-1)])
unempLag2 <- as.vector(Lag(rate,k=1))

unempDiff <- rate - unempLag
unempLagDiff <- cbind(rate, unempLag, unempDiff)
unempLagDiff

plot(unempDiff,type="l",ylab="Change in Unemployment Rate",xlab="Time",
     main = "Unemployment in the US, 1948-2013")

#Look at the new ACF and PACF for unemployment
#pdf("acfapprovalDiff.pdf",width=6,height=3.25)
acf(unempDiff, na.action=na.pass)
#dev.off()

#pdf("pacfapprovalDiff.pdf",width=6,height=3.25)
pacf(unempDiff, na.action=na.pass)
#dev.off()

PP.test(as.vector(na.omit(unempDiff)))
adf.test(na.omit(unempDiff))

###########################################################################
#Estimation and model selection

## Model 1a:  ARIMA(0,1,0) model of approve
##

## Estimate an ARIMA(0,1,0) using arima
xcovariates <- cbind(sept.oct.2001, iraq.war, avg.price)
arima.res1a <- arima(approve, order = c(0,1,0),
                     xreg = xcovariates, include.mean = TRUE
                     )
print(arima.res1a)

# Extract estimation results from arima.res1a
pe.1a <- arima.res1a$coef                    # parameter estimates (betas)
se.1a <- sqrt(diag(arima.res1a$var.coef))    # standard errors
ll.1a <- arima.res1a$loglik                  # log likelihood at its maximum
sigma2hat.1a <- arima.res1a$sigma2           # standard error of the regression
aic.1a <- arima.res1a$aic                    # Akaike Information Criterion
resid.1a <- arima.res1a$resid                # residuals

###########################################################################
## Model 1b:  ARIMA(1,1,0) model of approve
##

## Estimate an ARIMA(1,1,0) using arima
xcovariates <- cbind(sept.oct.2001, iraq.war, avg.price)
arima.res1b <- arima(approve, order = c(1,1,0),
                     xreg = xcovariates, include.mean = TRUE
                     )
print(arima.res1b)

# Extract estimation results from arima.res1a
pe.1b <- arima.res1b$coef                    # parameter estimates (betas)
se.1b <- sqrt(diag(arima.res1b$var.coef))    # standard errors
ll.1b <- arima.res1b$loglik                  # log likelihood at its maximum
sigma2hat.1b <- arima.res1b$sigma2           # standard error of the regression
aic.1b <- arima.res1b$aic                    # Akaike Information Criterion
resid.1b <- arima.res1b$resid                # residuals



###########################################################################
## Model 1c:  ARIMA(2,1,2) model of approve
##

## Estimate an ARIMA(2,1,2) using arima
xcovariates <- cbind(sept.oct.2001, iraq.war, avg.price)
arima.res1c <- arima(approve, order = c(2,1,2),
                     xreg = xcovariates, include.mean = TRUE
                     )
print(arima.res1c)

# Extract estimation results from arima.res1a
pe.1c <- arima.res1c$coef                    # parameter estimates (betas)
se.1c <- sqrt(diag(arima.res1c$var.coef))    # standard errors
ll.1c <- arima.res1c$loglik                  # log likelihood at its maximum
sigma2hat.1c <- arima.res1c$sigma2           # standard error of the regression
aic.1c <- arima.res1c$aic                    # Akaike Information Criterion
resid.1c <- arima.res1c$resid                # residuals



###########################################################################
## Model 1d:  ARIMA(0,1,0) model of approve including a spurrious regressor
##

## Estimate an ARIMA(0,1,0) using arima
xcovariates <- cbind(sept.oct.2001, iraq.war, avg.price, phony)
arima.res1d <- arima(approve, order = c(0,1,0),
                     xreg = xcovariates, include.mean = TRUE
                     )
print(arima.res1d)

# Extract estimation results from arima.res1a
pe.1d <- arima.res1d$coef                    # parameter estimates (betas)
se.1d <- sqrt(diag(arima.res1d$var.coef))    # standard errors
ll.1d <- arima.res1d$loglik                  # log likelihood at its maximum
sigma2hat.1d <- arima.res1d$sigma2           # standard error of the regression
aic.1d <- arima.res1d$aic                    # Akaike Information Criterion
resid.1d <- arima.res1d$resid                # residuals


#Based on ACF, PACF, and AIC, let's select Model 1a to be Model 1
acf(approveDiff, na.action=na.pass)

pacf(approveDiff, na.action=na.pass)

arima.res1a$aic
arima.res1b$aic
arima.res1c$aic
arima.res1d$aic

###########################################################################
## What would happen if we used linear regression on a single lag of approval?
lm.res1e <- lm(approve ~ approveLag + sept.oct.2001 + iraq.war + avg.price)
print(summary(lm.res1e))

# linear regression with a spurious regressor?
lm.res1f <- lm(approve ~ approveLag + sept.oct.2001 + iraq.war + avg.price + phony)
print(summary(lm.res1f))

# Check LS result for serial correlation in the first or second order (null of no serial correlation)
bgtest(lm.res1e,1)
bgtest(lm.res1f,2)


###########################################################################
#Forecasting and simulation

#### Because Zelig's arima simulator isn't currently available,
#### Let's use ldvsimev().

# First we need the simulated parameters
sims <- 10000
pe <- arima.res1a$coef
vc <- arima.res1a$var.coef
simparams <- mvrnorm(sims, pe, vc)
simbetas <- simparams[, (1+sum(arima.res1a$arma[1:2])):ncol(simparams)]
#simphi <- simparams[, 1:arima.res1$arma[1]]  # if AR(1) or greater

# Choose counterfactual periods
periodsToSim <- seq(from=26, to=65, by=1) #March 2003 to June 2006
nPeriodsToSim <- length(periodsToSim)

## Create hypothetical covariates over these periods
# Start with factual values
# For ARIMA, need to enter these covariates in same order as model
model <- approve ~ sept.oct.2001 + iraq.war + avg.price
selectdata <- extractdata(model, approval, na.rm = TRUE)
perioddata <- selectdata[periodsToSim,]

# Treat factual data as the baseline counterfactual
xhyp0 <- subset(perioddata, select = -approve )

# Suppose no Iraq War occurred
xhyp <- subset(perioddata, select = -approve )
xhyp$iraq.war <- rep(0,nPeriodsToSim)

# Leave out phi because this specification has no AR component
# phi <- simphi

# Construction of prior lags depends on order of ARIMA
# Only need initialY for I(1)
initialY <- selectdata$approve[periodsToSim[1] - 1]		# Original level of the response (for differenced models)
            # if ARIMA(1,1,0), instead: selectdata$approve[periodsToSim[1] - 2]
# Only need lagY for AR(1) or higher
lagY <- selectdata$approve[periodsToSim[1] - 1]			# The prior levels of y (scalar or vector)
            # if ARIMA(1,1,0), append: - selectdata$approve[periodsToSim[1] - 2]

# Simulate expected values of Y out to periods
# given hypothetical values of X, and an initial level of Y
# No Iraq scenario
noIraq.ev1 <- ldvsimev(xhyp,               # The matrix of hypothetical x's
                       simbetas,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=NA,        # Column containing the constant;
                                           #  no constant because model is differenced
                       #phi=phi,           # estimated AR parameters; length must match lagY 
                       #lagY=lagY,         # lags of y, most recent last
                       transform="diff",   # "log" to undo log transformation,
                                           # "diff" to under first differencing
                                           # "difflog" to do both
                       initialY=initialY   # for differenced models, lag of level of y
                       )

# Iraq scenario (factual)
Iraq.ev1 <- ldvsimev(xhyp0,              # The matrix of hypothetical x's
                     simbetas,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=NA,        # Column containing the constant;
                                         #  no constant because model is differenced
                     #phi=phi,           # estimated AR parameters; length must match lagY 
                     #lagY=lagY,         # lags of y, most recent last
                     transform="diff",   # "log" to undo log transformation,
                                         # "diff" to under first differencing
                                         # "difflog" to do both
                     initialY=initialY   # for differenced models, lag of level of y
                     )


at.xaxis <- seq(from = -1, to = length(approve[1:25]) + nPeriodsToSim, by = 12)
lab.xaxis <- seq(from = 2001, by = 1, length.out = length(at.xaxis))

#pdf("noIraqARIMA.pdf",width=6,height=3.25)
ctrfactTS(observed = approve[1:25],
          predicted = noIraq.ev1$pe,
          lower = noIraq.ev1$lower,
          upper = noIraq.ev1$upper,
          #se = NULL,
          predicted0 = Iraq.ev1$pe,
          lower0 = Iraq.ev1$lower,
          upper0 = Iraq.ev1$upper,
          #se0 = NULL,
          factual = approve[26:65],
          at.xaxis = at.xaxis,
          lab.xaxis = lab.xaxis,
          #ylim = c(0, 100),
          main = "Ctrfactual Effect of No Iraq War from ARIMA(0,1,0)",
          xlab = "Year",
          ylab = "Presidential Approval (%)",
          col = "red",
          col0 = "blue")
#dev.off()         


###########################################################################
# What if we used an ARMA model in this case (ignored possible
# unit root)?  Would we get better long run estimates?

## Model 2a:  ARIMA(1,0,0) model of approve
##

## Estimate an ARIMA(1,0,0) using arima
xcovariates <- cbind(sept.oct.2001, iraq.war, avg.price)
arima.res2a <- arima(approve, order = c(1,0,0),
                     xreg = xcovariates, include.mean = TRUE
                     )
print(arima.res2a)

# Extract estimation results from arima.res1a
pe.2a <- arima.res2a$coef                    # parameter estimates (betas)
se.2a <- sqrt(diag(arima.res2a$var.coef))    # standard errors
ll.2a <- arima.res2a$loglik                  # log likelihood at its maximum
sigma2hat.2a <- arima.res2a$sigma2           # standard error of the regression
aic.2a <- arima.res2a$aic                    # Akaike Information Criterion
resid.2a <- arima.res2a$resid                # residuals

# First we need the simulated parameters
sims <- 10000
pe <- arima.res2a$coef
vc <- arima.res2a$var.coef
simparams <- mvrnorm(sims, pe, vc)
simbetas <- simparams[, (1+sum(arima.res2a$arma[1:2])):ncol(simparams)]
simphi <- simparams[, 1:arima.res2a$arma[1]]  # if AR(1) or greater

# Choose counterfactual periods
periodsToSim <- seq(from=26, to=65, by=1)
nPeriodsToSim <- length(periodsToSim)

## Create hypothetical covariates over these periods
# Start with factual values
# For ARIMA, need to enter these covariates in same order as model
model <- approve ~ sept.oct.2001 + iraq.war + avg.price
selectdata <- extractdata(model, approval, na.rm = TRUE)
perioddata <- selectdata[periodsToSim,]

# Treat factual data as the baseline counterfactual
xhyp0 <- subset(perioddata, select = -approve )

# Suppose no Iraq War occurred
xhyp <- subset(perioddata, select = -approve )
xhyp$iraq.war <- rep(0,nPeriodsToSim)

# Construction of prior lags depends on order of ARIMA
# Only need initialY for I(1)
# Only need lagY for AR(1) or higher
lagY <- selectdata$approve[periodsToSim[1] - 1]
            
# Simulate expected values of Y out to periods
# given hypothetical values of X, and an initial level of Y
# No Iraq scenario
noIraq.ev1 <- ldvsimev(xhyp,               # The matrix of hypothetical x's
                       simbetas,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=1,         # Column containing the constant;
                       phi=mean(simphi),   # estimated AR parameters; length must match lagY 
                       lagY=lagY           # lags of y, most recent last
                       )

# Iraq scenario (factual)
Iraq.ev1 <- ldvsimev(xhyp0,              # The matrix of hypothetical x's
                     simbetas,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=1,         # Column containing the constant;
                     phi=mean(simphi),   # estimated AR parameters; length must match lagY 
                     lagY=lagY           # lags of y, most recent last
                     )


at.xaxis <- seq(from = -1, to = length(approve[1:25]) + nPeriodsToSim, by = 12)
lab.xaxis <- seq(from = 2001, by = 1, length.out = length(at.xaxis))

#pdf("noIraqARMA.pdf",width=6,height=3.25)
ctrfactTS(observed = approve[1:25],
          predicted = noIraq.ev1$pe,
          lower = noIraq.ev1$lower,
          upper = noIraq.ev1$upper,
          #se = NULL,
          predicted0 = Iraq.ev1$pe,
          lower0 = Iraq.ev1$lower,
          upper0 = Iraq.ev1$upper,
          #se0 = NULL,
          factual = approve[26:65],
          at.xaxis = at.xaxis,
          lab.xaxis = lab.xaxis,
          #ylim = c(0, 100),
          main = "Ctrfactual Effect of No Iraq War (Red) from AR(1)",
          xlab = "Year",
          ylab = "Presidential Approval (%)",
          col = "red",
          col0 = "blue")
#dev.off()         


###########################################################################
#What if we used the linear regression model with a lagged DV?

lm.res3 <- lm(approve ~ approveLag + sept.oct.2001 + iraq.war + avg.price)
print(summary(lm.res3))

#Extract estimation results from arima.res1a
pe.3 <- coef(lm.res3)                # parameter estimates
vc.3 <- vcov(lm.res3)                # variance-covariance (of params)
se.3 <- sqrt(diag(vc.3))             # standard errors

#First we need the simulated parameters
sims <- 10000
simparams <- mvrnorm(sims, pe.3, vc.3)
# Pluck out the lag parameter
simphi <- simparams[,2]
# Get the rest of the betas
simbetas <- simparams[,c(1,3:ncol(simparams))]

#Choose counterfactual periods
periodsToSim <- seq(from=26, to=65, by=1)
nPeriodsToSim <- length(periodsToSim)

#Create hypothetical covariates over these periods
# Start with factual values (note I leave out the lag;
# ldvsimev() will construct it)
model <- approve ~ sept.oct.2001 + iraq.war + avg.price
selectdata <- extractdata(model, approval, na.rm = TRUE)
perioddata <- selectdata[periodsToSim,]

#Treat factual data as the baseline counterfactual
xhyp0 <- subset(perioddata, select = -approve )

#Suppose no Iraq War occurred
xhyp <- subset(perioddata, select = -approve )
xhyp$iraq.war <- rep(0,nPeriodsToSim)

#Construction of prior lags depends on order of ARIMA
#Only need initialY for I(1)
#Only need lagY for AR(1) or higher
lagY <- selectdata$approve[periodsToSim[1] - 1]
            
#Simulate expected values of Y out to periods
#given hypothetical values of X, and an initial level of Y
#No Iraq scenario
noIraq.ev1 <- ldvsimev(xhyp,               # The matrix of hypothetical x's
                       simbetas,           # The matrix of simulated betas
                       ci=0.95,            # Desired confidence interval
                       constant=1,         # Column containing the constant;
                       phi=mean(simphi),   # estimated AR parameters; length must match lagY 
                       lagY=lagY           # lags of y, most recent last
                       )

#Iraq scenario (factual)
Iraq.ev1 <- ldvsimev(xhyp0,              # The matrix of hypothetical x's
                     simbetas,           # The matrix of simulated betas
                     ci=0.95,            # Desired confidence interval
                     constant=1,         # Column containing the constant;
                     phi=mean(simphi),   # estimated AR parameters; length must match lagY 
                     lagY=lagY           # lags of y, most recent last
                     )


at.xaxis <- seq(from = -1, to = length(approve[1:25]) + nPeriodsToSim, by = 12)
lab.xaxis <- seq(from = 2001, by = 1, length.out = length(at.xaxis))

#pdf("noIraqLS.pdf",width=6,height=3.25)
ctrfactTS(observed = approve[1:25],
          predicted = noIraq.ev1$pe,
          lower = noIraq.ev1$lower,
          upper = noIraq.ev1$upper,
          #se = NULL,
          predicted0 = Iraq.ev1$pe,
          lower0 = Iraq.ev1$lower,
          upper0 = Iraq.ev1$upper,
          #se0 = NULL,
          factual = approve[26:65],
          at.xaxis = at.xaxis,
          lab.xaxis = lab.xaxis,
          #ylim = c(0, 100),
          main = "Ctrfactual Effect of No Iraq War (Red) from LS with Lagged DV",
          xlab = "Year",
          ylab = "Presidential Approval (%)",
          col = "red",
          col0 = "blue")
#dev.off()         


################################################
#Cointegration analysis

set.seed(123456)

# Generate cointegrated data
e1 <- rnorm(100)
e2 <- rnorm(100)
x <- cumsum(e1)
y <- 0.6*x + e2

#Run step 1 of the Engle-Granger two step
coint.reg <- lm(y ~ x -1)				#Estimate the cointegration vector by least squares with no constant
coint.err <- residuals(coint.reg)		#This gives us the cotingeration vector

#Check for stationarity of the cointegration vector
punitroot(adf.test(coint.err)$statistic, trend="nc")

#Make the lag of the cointegration error term
coint.err.lag <- coint.err[1:(length(coint.err)-2)]

#Make the difference of y and x
dy <- diff(y)
dx <- diff(x)

#And their lags
dy.lag <- dy[1:(length(dy)-1)]
dx.lag <- dx[1:(length(dx)-1)]

#Delete the first dy, because we are missing lags for this obs
dy <- dy[2:length(dy)]

#Estimate an Error Correction Model with LS
ecm1 <- lm(dy ~ coint.err.lag + dy.lag + dx.lag)
summary(ecm1)


#Alternatively, we can use the Johansen estimator
#Create a matrix of the cointegrated variables
cointvars <- cbind(y,x)
# Perform cointegration tests
coint.test1 <- ca.jo(cointvars,
                   ecdet = "const",
                   type="eigen",
                   K=2,
                   spec="longrun")

summary(coint.test1)


ecm.res1 <- cajorls(coint.test1,
                    r = 1,           # Cointegration rank
                    reg.number = 1)  # which variable(s) to put on LHS
#(column indexes of cointvars)

summary(ecm.res1)


########################################################
#For the presidential approval example, use an ECM equivalent to the
#ARIMA(1,0,1) model that we chose earlier

cointvars <- cbind(approve,avg.price)
ecm.test1 <- ca.jo(cointvars,
                   ecdet = "const",
                   type="eigen",
                   K=2,
                   spec="longrun",
                   dumvar=cbind(sept.oct.2001,iraq.war)
                   )

summary(ecm.test1)

ecm.res1 <- cajorls(ecm.test1,
                    r = 1,
                    reg.number = 1)

summary(ecm.res1$rlm)


########################################################
#Cointegration analysis with a spurious regressor
cointvars <- cbind(approve,avg.price,phony)

ecm.test1 <- ca.jo(cointvars,
                   ecdet = "const",
                   type="eigen",
                   K=2,
                   spec="longrun",
                   dumvar=cbind(sept.oct.2001,iraq.war)
                   )

summary(ecm.test1)

ecm.res1 <- cajorls(ecm.test1,
                    r = 1,
                    reg.number = 1)

summary(ecm.res1$rlm)

