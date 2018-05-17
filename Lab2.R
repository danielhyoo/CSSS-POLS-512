


###########################################################################
#CSSS 594
#Lab Session 2 - Temporal Concepts: Trends, Stochastic Processes, and Seasonality
#4/10/15

rm(list=ls())
library(tseries)
library(MASS)
set.seed(12321)

###########################################################################
#Deterministic Trends

###########################################################################

#Simulate a deterministic trend with noise, de-trend the data and plot the time series.
#Set the slope of the trend
b1 <- 3/2

#Set the intercept
b0 <- 2

#Set the number of periods
n <- 50
t <- seq(0,n)
y <- rep(0,n)

#Simulate the data
for (i in 1:length(t)){
  y[i] <- b1*t[i] + b0 + rnorm(1,0,15) #The rnorm gives us the noise with mean 0, variance 15
}

#Plot the data
par(mfrow=c(2,1))
plot(y,type="l", col="red",ylab="y",xlab="Time",main=expression(paste("Simulated Deterministic Trend, y=2+3/2t + Noise")))
abline(a=2,b=3/2,lty="dashed")

#Now de-trend the time series 
y.minus.tbeta <- rep(0,n)
for (i in 1:length(t)){
  y.minus.tbeta[i] <- y[i] - b1*t[i]
}
plot(y.minus.tbeta,type="l", col="red",ylab="y",xlab="Time",main=expression(paste("Detrended Time Series")))
abline(a=2,b=0,lty="dashed")

#Take a minute to inspect the residuals

###########################################################################
#Find the least squares estimate of the slope
slope1 <- lm(y~t)
slope1

#How does it compare to the true beta?

#Plot the data with the true beta and the estimated beta
plot(y,type="l", col="red",ylab="y",xlab="Time",main=expression(paste("Simulated Deterministic Trend y=2+3/2t + Noise")))
abline(a=2,b=3/2,lty="dashed")
abline(a=slope1$coefficients[1],b=slope1$coefficients[2],lty="dashed",col="green")

###########################################################################
#Simulate new data with a deterministic trend and serial correlation

#Set the slope
b1 <- 3/2

#Set the intercept
b0 <- 2

#Set phi
phi <- 0.33

#Set the number of periods
n <- 50
t <- seq(0,n)
y <- rep(0,n)

for (i in 2:length(t)){
  y[i] <- y[i-1]*phi + b1*t[i] + b0 + rnorm(1,0,15)
}

#Plot the data and also the de-trended time series
par(mfrow=c(2,1))
plot(y,type="l", col="red",ylab="y",xlab="Time",main=expression(paste("Simulated Deterministic Trend + Noise + Serial Correlation")))
abline(a=2,b=3/2,lty="dashed")

y.minus.tbeta <- rep(0,n)
for (i in 1:length(t)){
  y.minus.tbeta[i] <- y[i] - b1*t[i]
}

plot(y.minus.tbeta,type="l", col="red",ylab="y",xlab="Time",main=expression(paste("Detrended Time Series + Noise + Serial Correlation")))
abline(a=2,b=0,lty="dashed")

#Take a minute to inspect the residuals again

###########################################################################
#Compare the two sets of plots and discuss the differences between a deterministic trend and stochastic process.

#What are some issues that can arise when analyzing de-trended time series data using regression?

###########################################################################
#Seasonality 

###########################################################################
#Accidental Deaths in the United States from 1973-1978, (from P. J. Brockwell and R. A. Davis (1991))
setwd("~/desktop")
accidents <- read.csv("USAccDeaths.csv",header=TRUE)
attach(accidents)

plot(time, USAccDeaths,type="l",col="red",ylab="y",xlab="Year", main = expression(paste("Accidental Deaths in the United States from 1973-1978")))

#Simulate a time series with seasonal variation
#Assume the data is de-trended
b1 <- 0

#Set the intercept
b0 <- 2

#Set the number of periods
n <- 60			#Assume a one month period for 5 years
t <- seq(0,n)
y <- rep(0,n)

#Simulate the data
for (i in 1:n){
  y[i] <- b1*t[i] + b0 + rnorm(1,0,1)
}

#Introduce additive seasonality during the first three months of each year
a <- seq(1,60, by=12)
b <- seq(2,60, by=12)
c <- seq(3,60, by=12)
q <- sort(c(a,b,c))

for (i in q){
  y[i] <- y[i]+6 #Seasonality can be additive or multiplicative
}

#Plot the data
plot(y,type="l",col="red",ylab="y",xlab="Time", main = expression(paste("Simulated Time Series with Three Month Additive Seasonality")))

#R has a special class of objects that corresponds to time series data
#The ts function allows for you to create a time series object, use help(ts) for reference
ts.1 <- ts(y, start=c(2000,1), end=c(2005,12), frequency=12) #We are creating a time series of length 60 months that starts from Jan 2000 until Dec 2005
help(ts)
ts.1

quartz()
plot(ts.1,type="l",col="red",ylab="y",xlab="Time", main = expression(paste("Simulated Times Series with Three Month Additive Seasonality")))

#Now remove the seasonal variation with the decompose function, use help(decompose for reference)
rm.seas.1 <- decompose(ts.1,type="additive")
rm.seas.1
plot(rm.seas.1)

#Alternatively, find the mean for each month, then subtract the corresponding monthly mean from each observation
month.avg <- rep(NA, 12)
m <- seq(0,48,by=12)

for (i in 1:12){
  month.avg[i] <- mean(y[m+i]) #Find the monthly average
}     

month.avg <- rep(month.avg,5)
rm.seas.2 <- y-month.avg

quartz()
plot(rm.seas.2,type="l",col="red",ylab="y",xlab="Time", main = expression(paste("Simulated Times Series with Three Month Additive Seasonality Removed")))

###########################################################################
#AR Processes

###########################################################################

#Simulate an AR(1) process with phi of 0.5. Plot the data and examine the ACF and PACF

#Sample from an AR(1), phi_1 = 0.5, using arima.sim() 
y <- arima.sim(list(order = c(1,0,0), ar = 0.50, ma = NULL), n=1000)

#Plot the series against time
plot(y,type="l",col="red",ylab="y",xlab="Time",
     main = expression(paste("Simulated AR(1) process with ",phi[1]," = 0.50")))
abline(a=0,b=0,lty="dashed")

#Plot the ACF
acf(y, main = expression(paste("ACF of AR(1) process with ",phi[1]," = 0.50")))

#Plot the PACF
pacf(y, main = expression(paste("PACF of AR(1) process with ",phi[1]," = 0.50")))

###########################################################################
#Make some general observations about the AR(1) plot. What do we learn from the ACF and PACF?

###########################################################################
#Simulate several AR(1) processes with -1 < phi < 1. Plot the data and examine the ACF and PACF

#Sample from AR(1) with phi of 0.8
ar1.1 <- arima.sim(list(order = c(1,0,0), ar=0.8, ma=NULL),n=1000)
plot(ar1.1,type="l",col="red",ylab="y",xlab="Time",
     main = expression(paste("Simulated AR(1) process with ",phi[1]," = 0.8")))
abline(a=0,b=0,lty="dashed")

acf(ar1.1, main = expression(paste("ACF of AR(1) process with ",phi[1]," = 0.8")))
pacf(ar1.1, main = expression(paste("ACF of AR(1) process with ",phi[1]," = 0.8")))

#Sample from AR(1) with phi of 0.15
ar1.2 <- arima.sim(list(order = c(1,0,0), ar=0.15, ma=NULL),n=1000)
plot(ar1.2,type="l",col="red",ylab="y",xlab="Time",
     main = expression(paste("Simulated AR(1) process with ",phi[1]," = 0.15")))
abline(a=0,b=0,lty="dashed")

acf(ar1.2, main = expression(paste("ACF of AR(1) process with ",phi[1]," = 0.15")))
pacf(ar1.2, main = expression(paste("ACF of AR(1) process with ",phi[1]," = 0.15")))

#Sample from AR(1) with phi of 0.99
ar1.3 <- arima.sim(list(order = c(1,0,0), ar=0.99, ma=NULL),n=1000)
plot(ar1.3,type="l",col="red",ylab="y",xlab="Time",
     main = expression(paste("Simulated AR(1) process with ",phi[1]," = 0.99")))
abline(a=0,b=0,lty="dashed")

acf(ar1.3, main = expression(paste("ACF of AR(1) process with ",phi[1]," = 0.99")))
pacf(ar1.3, main = expression(paste("ACF of AR(1) process with ",phi[1]," = 0.99")))

###########################################################################
#Check for a unit root on one of the AR(1) processes

#Perform a Phillips-Perron test or Augmented Dickey-Fuller test
PP.test(ar1.1)
adf.test(ar1.1)

###########################################################################
#Simulate an AR(1) process with phi = 1. Plot the data and examine the ACF and PACF

#Set number of observations
n <- 1000

#Set phi
phi <- 1

#Set y 
ar1.4 <- rep(0,n)

#Simulate AR(1) process with unit root
for (i in 2:n){
  ar1.4[i] <- ar1.4[i-1] + rnorm(1)
}

#Plot the time series
plot(ar1.4,type="l",col="red",ylab="y",xlab="Time",
     main = expression(paste("Simulated AR(1) process with ",phi[1]," = 1.0"))
)
abline(a=0,b=0,lty="dashed")

#Plot the ACF
acf(ar1.4, main = expression(paste("ACF of AR(1) process with ",phi[1]," = 1.0")))

#Plot the PACF
pacf(ar1.4, main = expression(paste("PACF of AR(1) process with ",phi[1]," = 1.0")))


###########################################################################
#Make some general observations about the AR(1) plot. What do we learn from the ACF and PACF?

###########################################################################
#Perform a unit root test on the data

#Perform a Phillips-Perron test or Augmented Dickey-Fuller test
PP.test(ar1.4)
adf.test(ar1.4)

###########################################################################
#Simulate an AR(2) process with phi_1 = 0.5 and phi_2 = 0.2. Plot the data and inspect the ACF and PACF

ar2.1 <- arima.sim(list(order = c(2,0,0), ar = c(0.50,0.2), ma = NULL), n=1000)

#Plot the series against time
plot(ar2.1,type="l",col="red",ylab="y",xlab="Time",
     main = expression(paste("Simulated AR(2) process with ",phi[1]," = 0.5, ", phi[2]," =0.2")))
abline(a=0,b=0,lty="dashed")

#Plot the ACF
acf(ar2.1, main = expression(paste("ACF of AR(2) process with ",phi[1]," = 0.50, ", phi[2]," =0.2")))

#Plot the PACF
pacf(ar2.1, main = expression(paste("PACF of AR(2) process with ",phi[1]," = 0.50, ", phi[2]," =0.2")))


###########################################################################
#Try solving for the roots using the polyroot() function
polyroot(c(1,0.5,0.2))

#Is the time series stationary?

#Confirm results with a unit root test
PP.test(ar2.1)
adf.test(ar2.1)

###########################################################################
#Sample from an AR(2), phi_1 = 1.2, phi_2 = -0.2 and plot the ACF and PACF

#Set number of observations
n <- 1000

#Set phi
phi_1 <- 1.2
phi_2 <- -0.2

#Set y vector
ar2.2 <- rep(0,n)

#Simulate AR(2) process with unit root
for (i in 3:n){
  ar2.2[i] <- ar2.2[i-1]*phi_1 + ar2.2[i-2]*phi_2 + rnorm(1)
}

#Plot the time series
plot(ar2.2,type="l",col="red",ylab="y",xlab="Time",
     main = expression(paste("Simulated AR(2) process with ",phi[1]," = 1.2 ", phi[2]," =-0.2"))
)
abline(a=0,b=0,lty="dashed")

#Plot the ACF
acf(ar2.2, main = expression(paste("ACF of AR(2) process with ",phi[1]," = 1.2 ", phi[2]," =-0.2")))

#Plot the PACF
pacf(ar2.2, main = expression(paste("PACF of AR(2) process with ",phi[1]," = 1.2 ", phi[2]," =-0.2")))

#Again, what can we (and can we not) infer from the ACF and PACF?

###########################################################################
#Try to check whether process is stationary with a unit root test

#Try solving for the roots
polyroot(c(1,1.2,-0.2))

#Perform a Phillips-Perron or Augmented Dickey-Fuller test
adf.test(ar2.2)
PP.test(ar2.2)

###########################################################################
#MA Processes

###########################################################################

#Simulate several MA(q) processes. Plot the data and examine the ACFs and PACFs

#Sample from an MA(1), psi_1 = 0.5
ma1.1 <- arima.sim(list(order = c(0,0,1), ar = NULL, ma = 0.5), n=1000)

#Plot the data
plot(ma1.1,type="l",col="red",ylab="y",xlab="Time",
     main = expression(paste("Simulated MA(1) process with ",psi[1]," = 0.50"))
)
abline(a=0,b=0,lty="dashed")

#Plot the ACF
acf(ma1.1, main = expression(paste("ACF of MA(1) process with ",psi[1]," = 0.50")))

#Plot the PACF
pacf(ma1.1, main = expression(paste("PACF of MA(1) process with ",psi[1]," = 0.50")))

#Sample from MA(2) with psi_1 = 0.3 and psi_2=0.7
ma2.1 <- arima.sim(list(order=c(0,0,2), ar=NULL, ma=c(0.3,0.7)),n=1000)
plot(ma2.1,type="l",col="red",ylab="y",xlab="Time", main = expression(paste("Simulated MA(2) process with ",psi[1]," = 0.3 ",psi[2]," =0.7")))
abline(a=0,b=0,lty="dashed")
acf(ma2.1, main = expression(paste("ACF of MA(2) process with ",psi[1]," = 0.3 ",psi[2]," =0.7")))
pacf(ma2.1, main = expression(paste("ACF of MA(2) process with ",psi[1]," = 0.3 ",psi[2]," =0.7")))

#MA(5) with psi_1 = 0.3 and psi_2=0.7 and psi_3=0.5 and psi_4=0.7 and psi_5=1.2
ma5.1 <- arima.sim(list(order=c(0,0,5), ar=NULL, ma=c(0.3,0.7,0.5,0.7,1.2)),n=1000)
plot(ma5.1,type="l",col="red",ylab="y",xlab="Time", main = expression(paste("Simulated MA(5) process with ",psi[1]," = 0.3 ",psi[2]," =0.7 ",psi[3]," =0.5 ", psi[4]," =0.7 ", psi[5], " =1.2")))
abline(a=0,b=0,lty="dashed")
acf(ma5.1, main = expression(paste("ACF of MA(5) process with ",psi[1]," = 0.3 ",psi[2]," =0.7 ",psi[3]," =0.5 ", psi[4]," =0.7 ", psi[5], "=1.2")))
pacf(ma5.1, main = expression(paste("ACF of MA(5) process with ",psi[1]," = 0.3 ",psi[2]," =0.7 ",psi[3]," =0.5 ", psi[4]," =0.7 ", psi[5], "=1.2")))

###########################################################################
#What do we learn about the effect of past shocks in an MA(q) process from the ACFs and PACFs? How can we identify an AR versus an MA process from the ACF and PACF.

###########################################################################
#ARMA Processes

###########################################################################


arma1.1 <- arima.sim(list(order=c(1,0,1), ar=0.3, ma=0.5), n=1000)
acf(arma1.1)
pacf(arma1.1)






