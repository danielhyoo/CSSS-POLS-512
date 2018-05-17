


# logistics: lab sessions, office hours, homeworks
# logistics: when the course is over, you should be able to do the following
# working with time series and panel data objects

# r refresher and time series objects


data(AirPassengers)
AP <- AirPassengers
AP

class(AP)

start(AP); end(AP); frequency(AP)

plot(AP, ylab="Passengers (1000's)")

layout(1:2)
plot(aggregate(AP))
boxplot(AP ~ cycle(AP))


www <- "https://raw.githubusercontent.com/svkerr/R_Files/master/TimeSeries/Maine.dat"
Maine.month <- read.table(www, header = TRUE)

attach(Maine.month)
class(Maine.month)
Maine.month

Maine.month.ts <- ts(unemploy, start = c(1996, 1), freq = 12)
Maine.month.ts

Maine.annual.ts <- aggregate(Maine.month.ts)/12
Maine.annual.ts

layout(1:2)
plot(Maine.month.ts, ylab="unemployed (%)")
plot(Maine.annual.ts, ylab="unemployed (%)")

Maine.Feb <- window(Maine.month.ts, start = c(1996,2), freq = TRUE)
Maine.Aug <- window(Maine.month.ts, start = c(1996,8), freq = TRUE)
Feb.ratio <- mean(Maine.Feb) / mean(Maine.month.ts)
Aug.ratio <- mean(Maine.Aug) / mean(Maine.month.ts)

Maine.Feb
Maine.month.ts


www <- "https://raw.githubusercontent.com/dallascard/Introductory_Time_Series_with_R_datasets/master/global.dat"
Global <- scan(www)


Global.ts <- ts(Global, st = c(1856, 1), end = c(2005, 12), fr = 12)

Global.annual <- aggregate(Global.ts, FUN = mean)
plot(Global.ts)
plot(Global.annual)

New.series <- window(Global.ts, start=c(1970, 1), end=c(2005, 12))
New.time <- time(New.series)
plot(New.series); abline(reg=lm(New.series ~ New.time))


www <- "https://raw.githubusercontent.com/kaybenleroll/dublin_r_workshops/master/ws_timeseries_201309/cbe.dat"
CBE <- read.table(www, header=T)

CBE[1:4,]

class(CBE)

Elec.ts <- ts(CBE[, 3], start = 1958, freq = 12)
Beer.ts <- ts(CBE[, 2], start = 1958, freq = 12)
Choc.ts <- ts(CBE[, 1], start = 1958, freq = 12)
plot(cbind(Elec.ts, Beer.ts, Choc.ts))



library(foreign)
library(tidyverse)

setwd("/Users/danielyoo/CSSS-POLS-512/Labs")
data<-read.csv("Lab1data.csv", header=T)  



# Data Frames

is.data.frame(data) #Yes!
is.matrix(data) #No
is.character(data$Year)
data$Year<-as.character(data$Year)

# Data Frames
names(data)

# Data Frames
summary(data)

# Data Frames
head(unique(data$country)) # observations on 174 countries
head(tapply(data$country, data$Year, length))
head(tapply(data$Year, data$country, length))

# Data Frames
data<-na.omit(data) # listwise deletion!!

dim(data)

attach(data)



# Data Frames
plot(polity2, GDP.per.capita.PPP.current.international, ylab="Polity2", xlab="GDP per capita")


# Data Frames
data$democracy[data$polity2>0]<-1
data$democracy[data$polity2<0|data$polity2==0]<-0
summary(data$democracy)


# Data Frames
data$democracy.2<-rep(NA, length(data$polity2)) # 1305

for (i in 1:length(data$polity2)) {
  if (data$polity2[i]>0) data$democracy.2[i]<-1
  else data$democracy.2[i]<-0
}

head(cbind(data$democracy, data$democracy.2))

# Data Frames

#19. Subset the data frame to show only country name and GDP per capita

#20. Rearrange the columns of the data frame ascending by polity score

#21. Show only values of GDP per capita for South Africa from 2002 to 2008

#22. Create a new variable that takes the first letter of the country and attaches it to the year of observation

#23. Find the mean of GDP per capita for each year of observation

# Data Frames
library(tidyverse)
head(select(data, country, GDP.per.capita.PPP.current.international))
head(data[, c(1,3)])
head(data.frame(data$country, data$GDP.per.capita.PPP.current.international))


# Data Frames
head(arrange(data, polity2))
head(data[order(data$polity2),])

# Data Frames
head(filter(data, country==c("South Africa"), Year>=2002 & Year<=2008))
head(subset(data, data$country==c("South Africa") & data$Year>=2002 & Year<=2008))

# Data Frames
head(mutate(data, paste(substring(data$country, 1, 1), data$Year, sep="")))

# Data Frames
data%>%
  group_by(Year)%>%
  summarize(mean(GDP.per.capita.PPP.current.international, na.rm=T)
  )

