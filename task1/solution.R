# Time Series Analysis II

rm(list = ls())
dev.off()

library(BVAR)
library(dplyr)
library(forecast)
library(ggplot2)
library(lmtest)
library(readxl)
library(sandwich)
library(tseries)

setwd("D:\code_github\r")

#######################################
##### Ex. 3: Multivariate Process #####
#######################################

# 2.b) --------------------------------
# Create matrix A
A1 <- matrix(c(0.9, -0.6, 0.4, 0.8), 2, 2) 
A2 <- matrix(c(0.2, 0.5, 0.2, -0.4), 2, 2)
I2 <- diag(1, 2, 2)
N  <- matrix(0, 2 ,2)
A <- rbind(cbind(A1, A2), cbind(I2, N))

# Compute eigenvalues:
eigenvalues <- abs(eigen(A)$values)
eigenvalues

eigenvalues <= 1 # all eigenvalues are smaller than 1 except for the first one

# Or compute roots:
roots <- 1/abs(eigen(A)$values)
roots

roots >= 1 # all roots are larger than one except for the first one

# The process is not stationary - there is a unit root.
# A root near 1 in the AR polynomial suggests that the data should be differenced before fitting an ARMA model.

#######################################
##### Ex. 3: Consumer Price Index #####
#######################################

# 3.a) --------------------------------
data <- read_excel("consumer_price_index.xlsx")[4:255, 2]

cpi <- ts(data, start = c(2000, 1), frequency = 12) # monthly data

plot(cpi, main = "Swiss Consumer Price Index",
     ylab = "Value (Base: December 2020 = 100)",
     xlab = "Time")

# 3.b) --------------------------------
# The CPI shows a clear upward trend between 2000 and 2008.
# Moreover, there appears to be seasonality (cycle that repeats itself).
# Therefore, the CPI is not stationary.

plot(decompose(cpi)) # Decomposition of the time series

# 3.c) --------------------------------

## Alternative 1
cpi_mom_growth_rate <- diff(log(cpi)) # eliminates trend

plot(cpi_mom_growth_rate, main = "Swiss Consumer Price Index - Month-on-Month Growth Rate",
     ylab = "Month-on-Month Growth Rate",
     xlab = "Time")
# It is stationary.

## Alternative 2
cpi_yoy_growth_rate <- diff(log(cpi), lag = 12) # eliminates seasonality

plot(cpi_yoy_growth_rate, main = "Swiss Consumer Price Index - Year-on-Year Growth Rate",
     ylab = "Year-on-Year Growth Rate",
     xlab = "Time")
# There is not stationary.

## Alternative 3
# To make this time series stationary, both trend and seasonality need to be eliminated:
cpi_stationary <- log(cpi) %>%
  diff(lag = 1) %>%  # eliminates trend
  diff(lag = 12) # eliminates seasonality

plot(cpi_stationary, main = "Swiss Consumer Price Index - Growth Rate (Log Difference)",
     ylab = "Log Difference (Growth Rate)",
     xlab = "Time")

# The time series looks stationary.

# 3.d) --------------------------------
acf(cpi_stationary, main = "Autocorrelation Function of CPI Growth Rate") # lag 1 = after one year (12 months)
# The ACF plot shows a significant negative spike at lag 1.

# 3.e) --------------------------------
pacf(cpi_stationary)
# The PACF plot shows a significant (negative) spike at lag 1, which could suggest using an AR(1) model.

# 3.f) --------------------------------
arma21 <- arima(cpi_stationary, order = c(2, 0, 1), method = "CSS", include.mean = F)
summary(arma21)

#######################################
########### Ex. 4: Oil Price ##########
#######################################

# 4.a) --------------------------------
#data(package = "BVAR") # to check the available data within a package
data <- fred_qd
data$TIME <- as.Date(rownames(data)) # create a date variable from columns
attach(data) # data to search path, i.e. colnames can be addressed directly

#XData <- data.frame(TIME,OILPRICEx)
plot(OILPRICEx ~ TIME, type = "l", main = "Oil Price Index")

# 4.b) --------------------------------
DIFF <- diff(OILPRICEx, 1)
plot(DIFF~TIME[-1], type = "l")

# 4.c) --------------------------------
# estimate AR(4)
AR4 <- lm(DIFF ~ -1 + lag(DIFF, 1) + lag(DIFF, 2) + lag(DIFF, 3) + lag(DIFF, 4))
summary(AR4)

arima <- arima(DIFF, order = c(4, 0, 0), method = "CSS", include.mean = F)
summary(arima)

# classic standard errors
coeftest(AR4, vcov = vcovHC(AR4, type = "const"))

# (heteroskedasticity) robust standard errors
coeftest(AR4, vcov = vcovHC(AR4, type = "HC1"))

# (heteroskedasticity) robust standard errors
coeftest(AR4, vcov = vcovHC(AR4, type = "HC3"))

# NeweyWest standard errors (Heteroskedasticity and Autocorrelation consistent s.e.)
coeftest(AR4, vcov = NeweyWest(AR4))

# 4.d) --------------------------------
newtime <- c(TIME[2:length(TIME)],
             seq(as.Date(TIME[length(TIME)]), by = "3 months", length.out = 20))
newdiff <- c(DIFF, rep(NA, 20))
new.T <- length(newtime)
n.periods <- 20
forecast <- predict(arima, n.ahead = n.periods)

# 4.e) --------------------------------
# Forecasted values and confidence interval
forecast.val  <- forecast$pred
forecast.loci <- forecast$pred - 1.96*forecast$se
forecast.upci <- forecast$pred + 1.96*forecast$se

plot(newtime, newdiff, type="l")
points(newtime[(new.T - n.periods + 1):new.T], forecast.val, type = "l", col = 2)
points(newtime[(new.T - n.periods + 1):new.T], forecast.loci, type = "l", col = 4, lty = 2)
points(newtime[(new.T - n.periods + 1):new.T], forecast.upci , type = "l", col = 4, lty = 2)

detach(data)
