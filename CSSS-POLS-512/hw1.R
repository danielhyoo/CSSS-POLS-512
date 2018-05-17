

# Homework 1

data <- read.csv("mysterytsUW.csv", header=T)

colnames(data)

# How to do this question?

# Assume stationary process

# Could have some combination of
# Deterministic trend
# Seasonality (additive with montly data)
# Autoregressive process
# Moving average process

# Box-Jenkins Method
# Be familiar with generic forms and processes of temporality
# Study these realizations in the data, using diagnostics 
# Make guess, assess and reiterate

# Main diagnostics will be ACF, PACF, decomposition, visual inspection
# Should know how to construct ts objects and use the decompose function

# If there is a deterministic trend, indicate evidence, describe it, remove it
# If there is seasonality, describe the seasonal cycle and remove the seasonal means from the data
# If there is autoregressive process, describe the order and likely sign and magnitude
# If there is moving average process, describe the order and likely sign and magnitude

attach(data)
# Inspect first time series

series <- e

plot(series, type="l", col="purple")

# Any thoughts?

# Construct a ts object
ts.series <- ts(series, frequency=12)

# Is there seasonality?
decompose(ts.series)
plot(decompose(ts.series))

acf(series)
pacf(series)









