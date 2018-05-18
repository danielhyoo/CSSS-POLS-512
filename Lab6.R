########################################################################################################
#CSSS 512
#Lab Session 5 - Panel Data Models with Variable Intercepts
#5/4/18

########################################################################################################
# Clear memory
rm(list=ls())

# Load libraries
library(nlme)      # Estimation of mixed effects models
library(lme4)      # Alternative package for mixed effects models
library(plm)       # Econometrics package for linear panel models
library(arm)       # Gelman & Hill code for mixed effects simulation
library(pcse)      # Calculate PCSEs for LS models (Beck & Katz)
library(tseries)   # For ADF unit root test
library(simcf)     # For panel functions and simulators

# Load Democracy data
# Variable names:
# COUNTRY	YEAR	BRITCOL 	CATH
# CIVLIB	EDT	ELF60	GDPW	MOSLEM
# NEWC    	OIL	POLLIB	REG	STRA
setwd("~/desktop")
data <- read.csv("democracy.csv",header=TRUE,na.strings=".")

# Create lags and differences now to correctly listwise delete
GDPWlag0 <- lagpanel(data$GDPW,data$COUNTRY,data$YEAR,1)			#Lag GDPW by using simcf and specifying the units and periods
GDPWdiff0 <- data$GDPW - GDPWlag0									#Compute the first difference
GDPWdifflag0 <- lagpanel(GDPWdiff0,data$COUNTRY,data$YEAR,1)		#Store the first lag of the first difference
GDPWdifflag20 <- lagpanel(GDPWdiff0,data$COUNTRY,data$YEAR,2)		#Store the second lag of the first difference

nrow(data)
ncol(data)
any(is.na(data))

# Listwise delete using only the data we need
selectdata <- na.omit(cbind(data$COUNTRY,
                            data$YEAR,
                            data$GDPW,
                            data$OIL,
                            data$REG,
                            data$EDT,
                            GDPWlag0,
                            GDPWdiff0,
                            GDPWdifflag0,
                            GDPWdifflag20))
selectdata <- as.data.frame(selectdata)
names(selectdata) <- c("COUNTRY","YEAR","GDPW","OIL",
                       "REG","EDT","GDPWlag","GDPWdiff", "GDPWdifflag", "GDPWdifflag2")

nrow(data)
nrow(selectdata)
nrow(data)

# Now attach the listwise deleted data
attach(selectdata)

selectdata[1:10,]

# Create a list of the unique country numbers
countrylist <- unique(COUNTRY)						#Generate a list of 113 country numbers
length(countrylist)									#Double check the number of countries

# Create a matrix of fixed effect dummies (for convenience)
fe <- makeFEdummies(COUNTRY)							#Create dummy variables for each of the 113 countries
ncol(fe)											#Double check the dimensions of the fe matrix
nrow(fe)

# Designate the number of units
n <- length(countrylist)

########################################################################################################
#Examine the time series plots, ACFs, and PACFs for panel data

setwd("~/desktop/plots")
# Diagnose time series
# Look at the time series for each country
for (i in 1:length(countrylist)) {									#Create a for loop from 1 to the number of countries (113)
    currcty <- countrylist[i]										#Make note of the country by number in the loop
    filename <- paste("tsGDPWcty",currcty,".pdf",sep="")			#Create the file name of the plot
    pdf(filename,width=6,height=3.25)								#Generate the PDF file
    plot(GDPW[COUNTRY==currcty],type="l",ylab="GDP",xlab="Time",	#Generate the plot of GDPW for the country by its number
         main = paste("Country",currcty) )
    dev.off()														#Turn off the PDF device
}

# Look at the ACF for each country
for (i in 1:length(countrylist)) {									#Create a for loop from 1 to the number of countries (113)
    currcty <- countrylist[i]										#Make note of the country by its number in the loop
    filename <- paste("acfGDPWcty",currcty,".pdf",sep="")			#Create the file name of the plot
    pdf(filename,width=6,height=3.25)								#Generate the PDF file
    acf(GDPW[COUNTRY==currcty])										#Generate the ACF plot of GDPW for the country by its number
    dev.off()														#Turn off the PDF device
}
   
# Look at the PACF
for (i in 1:length(countrylist)) {
    currcty <- countrylist[i]
    filename <- paste("pacfGDPWcty",currcty,".pdf",sep="")
    pdf(filename,width=6,height=3.25)
    pacf(GDPW[COUNTRY==currcty])									#Generate the ACF plot of GDPW for the country by its number
    dev.off()
}

# Check for a unit root in each country
PPtest.pvalues <- rep(0,113)											#Create empty vectors for PP test p-values
adftest.pvalues <- rep(0,113)											#Create empty vectors for adf test p-values

for (i in 1:length(countrylist)) {									#Create a for loop from 1 to the number of countries (113)
    currcty <- countrylist[i]										#Make note of the country by its number in the loop

    # Check PP unit root test, omitting errors due to short series
    curPP <- try(PP.test(GDPW[COUNTRY==currcty])$p.value)			#Find the p-value of the PP test for the country
    if (any(class(curPP)=="try-error")) curPP <- NA					#Make note if there is an error in the PP test, if so, fill with an NA
    PPtest.pvalues[i] <- curPP										#Store the p-value of the PP test in the PP test vector

    curadf <- try(adf.test(GDPW[COUNTRY==currcty])$p.value)			#Do the same with the adf test results
    if (any(class(curadf)=="try-error")) curadf <- NA
    adftest.pvalues[i] <- curadf
  }


#pdf("PPtest.pdf",width=6,height=3.25)
hist(PPtest.pvalues)          # Plot a histogram of the p-values
#dev.off()

#pdf("adftest.pdf",width=6,height=3.25)
hist(adftest.pvalues)         # Plot a histogram of the p-values
#dev.off()


# Repeat the ACF, PACF, and unit root tests on the differences data
setwd("~/desktop/diffplots")

# Look at the time series for each country in DIFFERENCES
for (i in 1:length(countrylist)) {
    currcty <- countrylist[i]
    filename <- paste("tsGDPWdiffcty",currcty,".pdf",sep="")
    pdf(filename,width=6,height=3.25)
    plot(GDPWdiff[COUNTRY==currcty],type="l",ylab="GDP",xlab="Time",
         main = paste("Country",currcty) )
    dev.off()
}

# Look at the ACF for each country in DIFFERENCES
for (i in 1:length(countrylist)) {
    currcty <- countrylist[i]
    filename <- paste("acfGDPWdiffcty",currcty,".pdf",sep="")
    pdf(filename,width=6,height=3.25)
    acf(GDPWdiff[COUNTRY==currcty])
    dev.off()
}
   
# Look at the PACF in DIFFERENCES
for (i in 1:length(countrylist)) {
    currcty <- countrylist[i]
    filename <- paste("pacfGDPWdiffcty",currcty,".pdf",sep="")
    pdf(filename,width=6,height=3.25)
    pacf(GDPWdiff[COUNTRY==currcty])
    dev.off()
}

# Check for a unit root in each country, differenced
PPtestdiff.pvalues <- rep(0,113)
adftestdiff.pvalues <- rep(0,113)
for (i in 1:length(countrylist)) {
    currcty <- countrylist[i]

    # Check PP unit root test, omitting errors due to short series
    curPPdiff <- try(PP.test(GDPWdiff[COUNTRY==currcty])$p.value)
    if (any(class(curPPdiff)=="try-error")) curPPdiff <- NA
    PPtestdiff.pvalues[i] <- curPPdiff

    curadfdiff <- try(adf.test(GDPWdiff[COUNTRY==currcty])$p.value)
    if (any(class(curadfdiff)=="try-error")) curadfdiff <- NA
    adftestdiff.pvalues[i] <- curadfdiff
  }

#pdf("PPtestdiff.pdf",width=6,height=3.25)
hist(PPtestdiff.pvalues)          # Plot a histogram of the p-values
#dev.off()

#pdf("adftestdiff.pdf",width=6,height=3.25)
hist(adftestdiff.pvalues)         # Plot a histogram of the p-values
#dev.off()


# Alternative panel based diagnostics available in the plm library
# (This package recently expanded to contain many many panel data tests
#  for serial correlation, fixed effects, and unit roots)

# First, create a plm data frame (special data frame that "knows" the
# unit variable and time variable
#pdata <- plm.data(selectdata, index=c("COUNTRY", "YEAR"))
pdata <- pdata.frame(selectdata, index=c("COUNTRY", "YEAR"))


# Do an panel unit root test on the undifferenced GDP data;
# there are many options; see ?purtest

# Note:  for some reason this isn't working
##purtest(GDPW~1, data=pdata, test="ips")

# Let's estimate some models, starting with random effects ARIMA


#########################################################################
#Random Effects Models

# Estimate a random effects AR(I)MA(p,q) model using lme (Restricted ML)
lme.res1 <- lme(# A formula object including the response,
                # the fixed covariates, and any grouping variables
                fixed = GDPWdiff ~ OIL + REG + EDT,			# i.e. response variable and explanatory variables 

                # The random effects component
                random = ~ 1 | COUNTRY,						# 1 indicates the intercept and COUNTRY indicates the grouping

                # The TS dynamics: specify the time & group variables,
                # and the order of the ARMA(p,q) process
                correlation = corARMA(form = ~ YEAR | COUNTRY,
                                      p = 1,  # AR(p) order
                                      q = 0   # MA(q) order
                                      ) 
                )

# Extract model results
pe.res1 <- fixed.effects(lme.res1)        # Point estimates of fixed effects
vc.res1 <- vcov(lme.res1)                 # Var-cov matrix of fixed effects estimates
se.res1 <- sqrt(diag(vc.res1))            # std erros of fixed effects estimates
re.res1 <- random.effects(lme.res1)       # "Estimated" random effects by group 
ll.res1 <- logLik(lme.res1)               # Log-likelihood at maximum
resid.res1 <- resid(lme.res1)             # Residuals
aic.res1 <- AIC(lme.res1)                 # Akaike Information Criterion

summary(lme.res1)

# Interpret the model using custom code
sims <- 1000										#Set the number of simulations
simbetas <- mvrnorm(sims,pe.res1,vc.res1)			#Sample the betas from a multivariate normal distribution

# Make matrix of hypothetical x's
formula <- GDPWdiff ~ OIL + REG + EDT				#Define the model
periods.out <- 50									#Set the number of prediction periods
xhyp <- cfMake(formula,selectdata,periods.out)		#Create a matrix of the mean values of the covariates
for (i in 1:periods.out) 							#Start a for loop from 1 to the number of periods to predict
  xhyp <- cfChange(xhyp, "EDT", x=mean(EDT,na.rm=TRUE)+sd(EDT,na.rm=TRUE),		#Use simcf to increase the EDT value by one std dev
                                scen=i)

xhyp

phi <- 0.25 # from model summary()
lagY <- mean(GDPWdiff) # Hypothetical previous change in Y for simulation
initialY <- mean(GDPW) # Hypothetical initial level of Y for simulation

# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev1 <- ldvsimev(xhyp,               # The matrix of hypothetical x's
                    simbetas,           # The matrix of simulated betas
                    ci=0.95,            # Desired confidence interval
                    constant=1,         # Column containing the constant
                                        # set to NA for no constant
                    phi=phi,            # estimated AR parameters; length must match lagY 
                    lagY=lagY,          # lags of y, most recent last
                    transform="diff",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                    initialY=initialY   # for differenced models, the lag of the level of y
                    )
    

# Simulate first differences in Y (on original level scale)
# out to periods.out given hypothetical future values of X, X0,
# and initial lags of the change in Y
sim.fd1 <- ldvsimfd(xhyp,               # The matrix of hypothetical x's
                    simbetas,           # The matrix of simulated betas
                    ci=0.95,            # Desired confidence interval
                    constant=1,         # Column containing the constant
                                        # set to NA for no constant
                    phi=phi,            # estimated AR parameters; length must match lagY 
                    lagY=lagY,          # lags of y, most recent last
                    transform="diff",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                    )


################################################################################
# Fixed Effects Models

# lme always wants to include random effects, so if we add fixed effects to the
# specification above, we will get a mixed effects model
#
# What if you just want fixed effects ARIMA, with no random effects?
# I'm not aware of a package that does this, but plm can come close:
# it can do fixed effects models with lagged dependent variables
#
# It also offers tests for serial correlation in these cases to check
# if you need the full ARMA toolkit: see pwartest(), pbgtest(), and pdwtest()

# Check for time invariant variables:
pvar(pdata)

# Note that we have to drop OIL in order to include FEs:  it is time invariant!

# "within" option tells plm to do fixed effects
# it does this by removing country means of GDPWdiff before estimation
plm.res2 <- plm(GDPWdiff ~ GDPWdifflag + REG + EDT, data = pdata, model="within")

# Some tests for serial correlation of errors (needed because we have a linear regression
# with lags of the dependent variable on the RHS

# the standard LM test (note we could specify order)
pbgtest(plm.res2)

# a robust test for small samples
# may not work at the moment?
# pwartest(plm.res2)

# Noting that the LM test finds serial correlation, we try adding a second lag
plm.res3 <- plm(GDPWdiff ~ GDPWdifflag + GDPWdifflag2 + REG + EDT, data = pdata, model="within")

# the standard LM test (note we could specify order)
pbgtest(plm.res3)

# Still borderline, but results are little changed...
# We will return to the single-lag specification to be consistent with the rest of the example,
# but in a more thorough analysis, we would want to get this right!


## Robust var-cov matrix alternatives for fixed effects models...
robust <- "None"   # Choose var-cov estimator here
if (robust=="None") vc <- vcov(plm.res2)
if (robust=="Arellano") vc <- vcovHC(plm.res2)   # Arellano (1987) heteroskedastic and serial correlation robust VC
if (robust=="BeckKatz") vc <- vcovBK(plm.res2)   # Beck and Katz (1995) panel corrected VC
if (robust=="DriscollKraay") vc <- vcovSCC(plm.res2)   # Driscoll and Kraay panel corrected VC

# Extract model results
pe.res2 <- coef(plm.res2)                                 # Point estimates of parameters
vc.res2 <- vc                                             # Var-cov matrix of point estimates
se.res2 <- sqrt(diag(vc.res2))                            # std erros of point estimates
tstat.res2 <- abs(pe.res2/se.res2)                        # t-statistics
df.res2 <- rep(plm.res2$df.residual, length(tstat.res2))  # residual degrees of freedom
pval.res2 <- 2*pt(tstat.res2, df.res2, lower.tail=FALSE)  # p-values
fe.res2 <- fixef(plm.res2)                                # the (removed) fixed effects by group 
resid.res2 <- resid(plm.res2)                             # Residuals

# Interpret the model using custom code
sims <- 1000												#Set the number of simulations
simparam <- mvrnorm(sims,pe.res2,vc.res2)					#Sample parameters from a multivariate normal dist
simphi <- simparam[,1]										#Pull off the simulated lag coefficient
# Put together the "constant" term (avg of the FEs, or a specific FE if you like)
# with the rest of the regressors
simbetas <- cbind(rep(mean(fe.res2), sims), simparam[,2:ncol(simparam)])

# Make matrix of hypothetical x's: drop OIL!
#                                  and now we need hypothetical countries!
#                                  and no intercept!

# Make matrix of hypothetical x's
# Note the formula below is a little different; we've removed the lag
# and will let the simulator handle lag structure
formula <- "GDPWdiff ~ REG + EDT"  
formula <- as.formula(formula)
  
periods.out <- 50										#Set the number of periods to predict
xhyp <- cfMake(formula,selectdata,periods.out)			#Create the initial matrix of mean of observed values of covariates
for (i in 1:periods.out) 
  xhyp <- cfChange(xhyp, "EDT", x=mean(EDT,na.rm=TRUE)+sd(EDT,na.rm=TRUE), scen=i)	#Use simcf to increase the EDT value by one std dev

phi <- mean(simphi) 
lagY <- mean(GDPWdiff) # Hypothetical previous change in Y for simulation
initialY <- mean(GDPW) # Hypothetical initial level of Y for simulation


# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev2 <- ldvsimev(xhyp,               # The matrix of hypothetical x's
                    simbetas,           # The matrix of simulated betas
                    ci=0.95,            # Desired confidence interval
                    constant=NA,        # NA indicates no constant!
                    phi=phi,            # estimated AR parameters; length must match lagY 
                    lagY=lagY,          # lags of y, most recent last
                    transform="diff",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                    initialY=initialY   # for differenced models, the lag of the level of y
                    )


# Simulate first differences in Y (on original level scale)
# out to periods.out given hypothetical future values of X, X0,
# and initial lags of the change in Y
sim.fd2 <- ldvsimfd(xhyp,               # The matrix of hypothetical x's
                    simbetas,           # The matrix of simulated betas
                    ci=0.95,            # Desired confidence interval
                    constant=1,         # Column containing the constant
                                        # set to NA for no constant
                    phi=phi,            # estimated AR parameters; length must match lagY 
                    lagY=lagY,          # lags of y, most recent last
                    transform="diff",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                    )




# Suppose we wanted Beck-Katz pcse's...
robust <- "BeckKatz"
if (robust=="BeckKatz") vc <- vcovBK(plm.res2)   # Beck and Katz (1995) panel corrected VC
vc.res2 <- vc                                             # Var-cov matrix of point estimates
se.res2 <- sqrt(diag(vc.res2))                            # std erros of point estimates
tstat.res2 <- abs(pe.res2/se.res2)                        # t-statistics
df.res2 <- rep(plm.res2$df.residual, length(tstat.res2))  # residual degrees of freedom
pval.res2 <- 2*pt(tstat.res2, df.res2, lower.tail=FALSE)  # p-values


# Interpret the model using custom code
sims <- 1000
simparam <- mvrnorm(sims,pe.res2,vc.res2)
# Pull off the simulated lag coefficient
simphi <- simparam[,1]
# Put together the "constant" term (avg of the FEs, or a specific FE if you like)
# with the rest of the regressors
simbetas <- cbind(rep(mean(fe.res2), sims), simparam[,2:ncol(simparam)])

# Make matrix of hypothetical x's: drop OIL!
#                                  and now we need hypothetical countries!
#                                  and no intercept!

# Make matrix of hypothetical x's
# Note the formula below is a little different; we've removed the lag
# and will let the simulator handle lag structure
formula <- "GDPWdiff ~ REG + EDT"  
formula <- as.formula(formula)
  
periods.out <- 50
xhyp <- cfMake(formula,selectdata,periods.out)
for (i in 1:periods.out) 
  xhyp <- cfChange(xhyp, "EDT", x=mean(EDT,na.rm=TRUE)+sd(EDT,na.rm=TRUE), scen=i)

phi <- mean(simphi) 
lagY <- mean(GDPWdiff) # Hypothetical previous change in Y for simulation
initialY <- mean(GDPW) # Hypothetical initial level of Y for simulation


# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev2bk <- ldvsimev(xhyp,               # The matrix of hypothetical x's
                    simbetas,           # The matrix of simulated betas
                    ci=0.95,            # Desired confidence interval
                    constant=NA,        # NA indicates no constant!
                    phi=phi,            # estimated AR parameters; length must match lagY 
                    lagY=lagY,          # lags of y, most recent last
                    transform="diff",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                    initialY=initialY   # for differenced models, the lag of the level of y
                    )


# Simulate first differences in Y (on original level scale)
# out to periods.out given hypothetical future values of X, X0,
# and initial lags of the change in Y
sim.fd2bk <- ldvsimfd(xhyp,               # The matrix of hypothetical x's
                    simbetas,           # The matrix of simulated betas
                    ci=0.95,            # Desired confidence interval
                    constant=1,         # Column containing the constant
                                        # set to NA for no constant
                    phi=phi,            # estimated AR parameters; length must match lagY 
                    lagY=lagY,          # lags of y, most recent last
                    transform="diff",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                    )



######
# Try LSDV:  less efficient route to fixed effects






################################################################################      
#Mixed Effects Models

# Estimate a mixed effects AR(I)MA(p,q) model using lme (Restricted ML)
lme.res3 <- lme(# A formula object including the response,
                # the fixed covariates, and the country fixed effects
                # (either as dummy variables or as a "factor" variable),
                # then remove the model intercept with - 1
                fixed = GDPWdiff ~ REG + EDT + fe - 1 ,
                             # NOTE:  I must drop OIL, which doesn't vary over
                             #        time for any country.
                             #        If I leave it in, I get a singularity error

                # The random effects component
                random = ~ 1 | COUNTRY,

                # The TS dynamics: specify the time & group variables,
                # and the order of the ARMA(p,q) process
                correlation = corARMA(form = ~ YEAR | COUNTRY,
                                      p = 1,  # AR(p) order
                                      q = 0   # MA(q) order
                                      ) 
                )

# Extract model results
pe.res3 <- fixed.effects(lme.res3)        # Point estimates of fixed effects
vc.res3 <- vcov(lme.res3)                 # Var-cov matrix of fixed effects estimates
se.res3 <- sqrt(diag(vc.res3))            # std erros of fixed effects estimates
re.res3 <- random.effects(lme.res3)       # "Estimated" random effects by group 
ll.res3 <- logLik(lme.res3)               # Log-likelihood at maximum
resid.res3 <- resid(lme.res3)             # Residuals
aic.res3 <- AIC(lme.res3)                 # Akaike Information Criterion


# Interpret the model using custom code
sims <- 1000									#Set the number of simulations
simbetas <- mvrnorm(sims,pe.res3,vc.res3)		#Sample from a multivariate normal dist.

# Make matrix of hypothetical x's: drop OIL!
#                                  and now we need hypothetical countries!
#                                  and no intercept!

# Make matrix of hypothetical x's
formula <- "GDPWdiff ~ REG + EDT - 1"			#Set the model less the intercept
selectdatafe <- cbind(selectdata,fe)			
fenames <- NULL
for (i in 1:ncol(fe)) {							#Create names for the fixed effects
  formula <- paste(formula,"+ fe",i," ",sep="")	#Unique fe name
  fenames <- c(fenames,paste("fe",i,sep=""))	#Combine with others
}
names(selectdatafe) <- c(names(selectdata),fenames)
formula <- as.formula(formula)
  
periods.out <- 50								#Set number of prediction periods
xhyp <- cfMake(formula,selectdatafe,periods.out)	#Create the initial matrix of mean values of covariates
for (i in 1:periods.out) 
  xhyp <- cfChange(xhyp, "EDT", x=mean(EDT,na.rm=TRUE)+sd(EDT,na.rm=TRUE),	#Increase EDT by one standard deviation
                                scen=i)

phi <- 0.25 # from model summary()
lagY <- mean(GDPWdiff) # Hypothetical previous change in Y for simulation
initialY <- mean(GDPW) # Hypothetical initial level of Y for simulation



# Simulate expected values of Y (on original level scale)
# out to periods.out given hypothetical future values of X,
# initial lags of the change in Y, and an initial level of Y
sim.ev3 <- ldvsimev(xhyp,               # The matrix of hypothetical x's
                    simbetas,           # The matrix of simulated betas
                    ci=0.95,            # Desired confidence interval
                    constant=NA,        # NA indicates no constant!
                    phi=phi,            # estimated AR parameters; length must match lagY 
                    lagY=lagY,          # lags of y, most recent last
                    transform="diff",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                    initialY=initialY   # for differenced models, the lag of the level of y
                    )
    

# Simulate first differences in Y (on original level scale)
# out to periods.out given hypothetical future values of X, X0,
# and initial lags of the change in Y
sim.fd3 <- ldvsimfd(xhyp,               # The matrix of hypothetical x's
                    simbetas,           # The matrix of simulated betas
                    ci=0.95,            # Desired confidence interval
                    constant=NA,        # Column containing the constant
                                        # set to NA for no constant
                    phi=phi,            # estimated AR parameters; length must match lagY 
                    lagY=lagY,          # lags of y, most recent last
                    transform="diff",   # "log" to undo log transformation,
                                        # "diff" to under first differencing
                                        # "difflog" to do both
                    )
    


#####################################################################
# Make some graphs

# Random effects model of change in GDP given increase in EDT
#pdf("simfdREupEDT.pdf",width=5,height=4.5)
plot.new()
par(usr=c(1,50,-15000,15000))			#Set extremes of the plotting region
axis(1,at=seq(1,50,10))					#Range of x axis on plot
axis(2,at=seq(-15000,15000,5000))		#Range of y axis on plot
title(xlab="Time",ylab="Expected cumulative change in constant GDP ($ pc)",main="Random effects ARIMA(1,1,0)") 

# Make the x-coord of a confidence envelope polygon
xpoly <- c(1:periods.out,				#lower envelope
           rev(1:periods.out),			#upper envelope
           1)

# Make the y-coord of a confidence envelope polygon
ypoly <- c(sim.fd1$lower.cumulative,		#lower envelope
           rev(sim.fd1$upper.cumulative),	#upper envelope
           sim.fd1$lower.cumulative[1])

# Choose the color of the polygon 
col <- "lightblue"

# Plot the polygon first, before the points & lines
polygon(x=xpoly,
        y=ypoly,
        col=col,
        border=FALSE
        )

# Plot the fitted line
lines(x=1:periods.out,y=sim.fd1$pe.cumulative,col="darkblue")
lines(x=c(0,50),y=c(0,0),lty="solid")

#dev.off()




# Fixed effects model of change in GDP given increase in EDT
#pdf("simfdFEupEDT.pdf",width=5,height=4.5)
plot.new()
par(usr=c(1,50,-15000,15000))
axis(1,at=seq(1,50,10))
axis(2,at=seq(-15000,15000,5000))
title(xlab="Time",ylab="Expected cumulative change in constant GDP ($ pc)",main="Fixed effects ARIMA(1,1,0)") 

# Make the x-coord of a confidence envelope polygon
xpoly <- c(1:periods.out,
           rev(1:periods.out),
           1)

# Make the y-coord of a confidence envelope polygon
ypoly <- c(sim.fd2$lower.cumulative,
           rev(sim.fd2$upper.cumulative),
           sim.fd2$lower.cumulative[1])

# Choose the color of the polygon 
col <- "lightgreen"

# Plot the polygon first, before the points & lines
polygon(x=xpoly,
        y=ypoly,
        col=col,
        border=FALSE
        )

# Plot the fitted line
lines(x=1:periods.out,y=sim.fd2$pe.cumulative,col="darkgreen")
lines(x=c(0,50),y=c(0,0),lty="solid")

#dev.off()




# Fixed effects model of change in GDP given increase in EDT, with PCSEs
#pdf("simfdFEPCSEupEDT.pdf",width=5,height=4.5)
plot.new()
par(usr=c(1,50,-15000,15000))
axis(1,at=seq(1,50,10))
axis(2,at=seq(-15000,15000,5000))
title(xlab="Time",ylab="Expected cumulative change in constant GDP ($ pc)",main="Fixed effects ARIMA(1,1,0)") 

# Make the x-coord of a confidence envelope polygon
xpoly <- c(1:periods.out,
           rev(1:periods.out),
           1)

# Make the y-coord of a confidence envelope polygon
ypoly <- c(sim.fd2bk$lower.cumulative,
           rev(sim.fd2bk$upper.cumulative),
           sim.fd2bk$lower.cumulative[1])

# Choose the color of the polygon 
col <- "lightgreen"

# Plot the polygon first, before the points & lines
polygon(x=xpoly,
        y=ypoly,
        col=col,
        border=FALSE
        )

# Plot the fitted line
lines(x=1:periods.out,y=sim.fd2bk$pe.cumulative,col="darkgreen")
lines(x=c(0,50),y=c(0,0),lty="solid")

#dev.off()



##########################
# Mixed effects model of change in GDP given increase in EDT
#pdf("simfdMEupEDT.pdf",width=5,height=4.5)
plot.new()
par(usr=c(1,50,-15000,15000))
axis(1,at=seq(1,50,10))
axis(2,at=seq(-15000,15000,5000))
title(xlab="Time",ylab="Expected cumulative change in constant GDP ($ pc)",main="Mixed effects ARIMA(1,1,0)") 

# Make the x-coord of a confidence envelope polygon
xpoly <- c(1:periods.out,
           rev(1:periods.out),
           1)

# Make the y-coord of a confidence envelope polygon
ypoly <- c(sim.fd3$lower.cumulative,
           rev(sim.fd3$upper.cumulative),
           sim.fd3$lower.cumulative[1])

# Choose the color of the polygon 
col <- "pink"

# Plot the polygon first, before the points & lines
polygon(x=xpoly,
        y=ypoly,
        col=col,
        border=FALSE
        )

# Plot the fitted line
lines(x=1:periods.out,y=sim.fd3$pe.cumulative,col="red")
lines(x=c(-10,1100),y=c(0,0),lty="solid")

#dev.off()
