# ARIMA estimation, fitiing, and interpretation
# Chris Adolph
rm(list=ls())

# Load libraries
library(forecast)     # For auto.arima and cross-validation
library(tseries)      # For unit root tests
library(lmtest)       # For Breusch-Godfrey LM test of serial correlation
library(RColorBrewer) # For nice colors
library(boot)
library(MASS)
library(simcf)


# ARIMA Cross-validation by rolling windows
# Adapted from Rob J Hyndman's code:
# http://robjhyndman.com/hyndsight/tscvexample/
#
# Could use further generalization, e.g. to seasonality
# Careful!  This can produce singularities using categorical covariates
arimaCV <- function(x, order, xreg, include.mean, forward=1, minper=50, seasonal=NULL) {
    require(forecast)
    if (!any(class(x)=="ts")) x <- ts(x)
    n <- length(x)
    mae <- matrix(NA, nrow=n-minper, ncol=forward)
    st <- tsp(x)[1]+(minper-2)   
    for(i in 1:(n-minper)) {
        xshort <- window(x, start=st+(i-minper+1), end=st+i)
        xnext <- window(x, start=st+(i+1), end=min(n, st+(i+forward)))
        xregshort <- window(xreg, start=st+(i-minper+1), end=st+i)
        xregnext <- window(xreg, start=st+(i+1), end=min(n, st+(i+forward)))
        fit <- Arima(xshort, order=order, seasonal=seasonal, xreg=xregshort, include.mean=include.mean)
        fcast <- forecast(fit, h=length(xnext), xreg=xregnext)
        mae[i,1:length(xnext)] <- abs(fcast[['mean']]-xnext)
    }
    colMeans(mae, na.rm=TRUE)
}



## Load data
##
## Democratic seat totals from the US House and Senate 
## Sources: Wikipedia, Bureau of Labor Statistics
##
## Variable names:  Congress StartYear DemSenateSeats DemHouseSeats Midterm DemPresident Unemployment

congress <- read.csv("congress.csv",header=TRUE)
attach(congress)

# Look at the time series
#pdf("tsDemSenateMaj.pdf",width=6,height=3.25)
plot(DemSenateMaj,type="l",ylab="seats",xlab="Time",
     main = "Democratic Senate Majority, 1963--2017")
#dev.off()

# Look at the ACF
#pdf("acfDemSenateMaj.pdf",width=6,height=3.25)
acf(DemSenateMaj)
#dev.off()

# Look at the PACF
#pdf("pacfDemSenateMaj.pdf",width=6,height=3.25)
pacf(DemSenateMaj)
#dev.off()

# Check for a unit root
PP.test(DemSenateMaj)
adf.test(DemSenateMaj)

### Remove period means to check structural break

periodMean <- (StartYear<1995)*mean(DemSenateMaj[StartYear<1995]) +
    (StartYear>1994)*mean(DemSenateMaj[StartYear>1994])
demeanDemSenateMaj <- DemSenateMaj - periodMean

# Look at the time series
#pdf("tsDemSenateMajDemean.pdf",width=6,height=3.25)
plot(demeanDemSenateMaj,type="l",ylab="seats",xlab="Time",
     main = "Demeaned Democratic Senate Majority, 1963--2017")
#dev.off()

# Look at the ACF
#pdf("acfDemSenateMajDemean.pdf",width=6,height=3.25)
acf(demeanDemSenateMaj)
#dev.off()

# Look at the PACF
#pdf("pacfDemSenateMajDemean.pdf",width=6,height=3.25)
pacf(demeanDemSenateMaj)
#dev.off()

# Check for a unit root
PP.test(demeanDemSenateMaj)
adf.test(demeanDemSenateMaj)

# Set rolling window length and look ahead period for cross-validation
minper <- 20
forward <- 3



#################################################################
## Model 1z:  AR(0) model of DemSenateMaj 
##

## Estimate an AR(0) using arima
xcovariates <- cbind(PartisanMidterm, PartisanUnem, Coattails, Pre1994)
order <- c(0,0,0)
arima.res1z <- arima(DemSenateMaj, order = order,
                     xreg = xcovariates, include.mean = TRUE
                     )
print(arima.res1z)

# Extract estimation results from arima.res1a
pe.1z <- arima.res1z$coef                    # parameter estimates (betas)
se.1z <- sqrt(diag(arima.res1z$var.coef))    # standard errors
ll.1z <- arima.res1z$loglik                  # log likelihood at its maximum
sigma2hat.1z <- arima.res1z$sigma2           # standard error of the regression
aic.1z <- arima.res1z$aic                    # Akaike Information Criterion
resid.1z <- arima.res1z$resid                # residuals

# Attempt at rolling window cross-validation (see caveats)
cv.1z <- arimaCV(DemSenateMaj, order=order, forward=forward,
                 xreg=xcovariates, include.mean=TRUE, minper=minper)





#################################################################
## Model 1a:  AR(1) model of DemSenateMaj 
##

## Estimate an AR(1) using arima
xcovariates <- cbind(PartisanMidterm, PartisanUnem, Coattails, Pre1994)
order <- c(1,0,0)
arima.res1a <- arima(DemSenateMaj, order = order,
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

# Attempt at rolling window cross-validation (see caveats)
cv.1a <- arimaCV(DemSenateMaj, order=order, forward=forward,
                 xreg=xcovariates, include.mean=TRUE, minper=minper)




#################################################################
## Model 1b:  AR(2) model of DemSenateMaj 
##

## Estimate an AR(1) using arima
xcovariates <- cbind(PartisanMidterm, PartisanUnem, Coattails, Pre1994)
order <- c(2,0,0)
arima.res1b <- arima(DemSenateMaj, order = order,
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

# Attempt at rolling window cross-validation (see caveats)
cv.1b <- arimaCV(DemSenateMaj, order=order, forward=forward,
                 xreg=xcovariates, include.mean=TRUE, minper=minper)




#################################################################
## Model 1c:  AR(3) model of DemSenateMaj 
##

## Estimate an AR(3) using arima
xcovariates <- cbind(PartisanMidterm, PartisanUnem, Coattails, Pre1994)
order <- c(3,0,0)
arima.res1c <- arima(DemSenateMaj, order = order,
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

# Attempt at rolling window cross-validation (see caveats)
cv.1c <- arimaCV(DemSenateMaj, order=order, forward=forward,
                 xreg=xcovariates, include.mean=TRUE, minper=minper)




#################################################################
## Model 1d:  MA(1) model of DemSenateMaj 
##

## Estimate an MA(1) using arima
xcovariates <- cbind(PartisanMidterm, PartisanUnem, Coattails, Pre1994)
order <- c(0,0,1)
arima.res1d <- arima(DemSenateMaj, order = order,
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

# Attempt at rolling window cross-validation (see caveats)
cv.1d <- arimaCV(DemSenateMaj, order=order, forward=forward,
                 xreg=xcovariates, include.mean=TRUE, minper=minper)



#################################################################
## Model 1e:  ARMA(1,1) model of DemSenateMaj 
##

## Estimate an ARMA(1,1) using arima
xcovariates <- cbind(PartisanMidterm, PartisanUnem, Coattails, Pre1994)
order <- c(1,0,1)
arima.res1e <- arima(DemSenateMaj, order = order,
                     xreg = xcovariates, include.mean = TRUE
                     )
print(arima.res1e)

# Extract estimation results from arima.res1a
pe.1e <- arima.res1e$coef                    # parameter estimates (betas)
se.1e <- sqrt(diag(arima.res1e$var.coef))    # standard errors
ll.1e <- arima.res1e$loglik                  # log likelihood at its maximum
sigma2hat.1e <- arima.res1e$sigma2           # standard error of the regression
aic.1e <- arima.res1e$aic                    # Akaike Information Criterion
resid.1e <- arima.res1e$resid                # residuals

# Attempt at rolling window cross-validation (see caveats)
cv.1e <- arimaCV(DemSenateMaj, order=order, forward=forward,
                 xreg=xcovariates, include.mean=TRUE, minper=minper)


#################################################################
## Model 1f:  ARIMA(0,1,0) model of DemSenateMaj 
##

## Estimate an ARIMA(0,1,0) using arima
xcovariates <- cbind(PartisanMidterm, PartisanUnem, Coattails, Pre1994)
order <- c(0,1,0)
arima.res1f <- arima(DemSenateMaj, order = order,
                     xreg = xcovariates, include.mean = TRUE
                     )
print(arima.res1f)

# Extract estimation results from arima.res1a
pe.1f <- arima.res1f$coef                    # parameter estimates (betas)
se.1f <- sqrt(diag(arima.res1f$var.coef))    # standard errors
ll.1f <- arima.res1f$loglik                  # log likelihood at its maximum
sigma2hat.1f <- arima.res1f$sigma2           # standard error of the regression
aic.1f <- arima.res1f$aic                    # Akaike Information Criterion
resid.1f <- arima.res1f$resid                # residuals

# Attempt at rolling window cross-validation (see caveats)
cv.1f <- arimaCV(DemSenateMaj, order=order, forward=forward,
                 xreg=xcovariates, include.mean=TRUE, minper=minper)


#################################################################
## Model 1g:  AR(1)AR(1)_3 model of DemSenateMaj
##

## Estimate an AR(1)AR(1)_3 using arima
xcovariates <- cbind(PartisanMidterm, PartisanUnem, Coattails, Pre1994)
order <- c(1,0,0)
arima.res1g <- arima(DemSenateMaj, order = order,
                     seasonal = list(order = c(1,0,0), period = 3),
                     xreg = xcovariates, include.mean = TRUE
                     )
print(arima.res1g)

# Extract estimation results from arima.res1a
pe.1g <- arima.res1g$coef                    # parameter estimates (betas)
se.1g <- sqrt(diag(arima.res1g$var.coef))    # standard errors
ll.1g <- arima.res1g$loglik                  # log likelihood at its maximum
sigma2hat.1g <- arima.res1g$sigma2           # standard error of the regression
aic.1g <- arima.res1g$aic                    # Akaike Information Criterion
resid.1g <- arima.res1g$resid                # residuals

# Attempt at rolling window cross-validation (see caveats)
cv.1g <- arimaCV(DemSenateMaj, order=order, seasonal = list(order = c(1,0,0), period = 3),
                 forward=forward, xreg=xcovariates, include.mean=TRUE, minper=minper)


mean(cv.1z)
mean(cv.1a)
mean(cv.1b)
mean(cv.1c)
mean(cv.1d)
mean(cv.1e)
mean(cv.1f)
mean(cv.1g)


# Interpret the model using custom code
sims <- 10000
simparam <- mvrnorm(sims, pe.1a, arima.res1a$var.coef)

xhyp <- rbind(c(-1, -(5.6-mean(Unemployment)), 0,0),
              c(0, -(5.6-mean(Unemployment)), 1,0),
              c(1, 5.6-mean(Unemployment), 0,0))

simphi <- simparam[,1] 
simbetas <- simparam[,2:ncol(simparam)]
lagY <- DemSenateMaj[length(DemSenateMaj)] # Hypothetical previous Y for simulation

# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev1a <- ldvsimev(xhyp,               # The matrix of hypothetical x's
                    simbetas,           # The matrix of simulated betas
                    ci=0.95,            # Desired confidence interval
                    constant=1,         # Column containing the constant
                                        # set to NA for no constant
                    phi=simphi,            # estimated AR parameters; length must match lagY 
                    lagY=lagY          # lags of y, most recent last
                    )


