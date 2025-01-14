# Time Series Analysis II

rm(list = ls())
dev.off()

library(readxl)
library(vars)
library(aTSA)
library(urca)

setwd("D:\code_github\r")

#######################################
######## Ex. 2: US Macroeconomy #######
#######################################

data <- read_excel("gdp_cons_inv.xlsx", range = "B2:D237", col_names = FALSE)
time <- seq(as.Date("1947-01-01"), as.Date("2005-10-01"), by = "quarters")

X <- log(data)
colnames(X) <- c("Y", "C", "I")

# 2.a) --------------------------------
y_range <- range(c(X$Y, X$C, X$I), na.rm = TRUE) # overall y-axis range

plot(X$Y ~ time, type = "l", col = "blue", lwd = 2, 
     main = "US Macro Variables (1950-2005)", 
     xlab = "Time", ylab = "Value", ylim = y_range)
lines(X$C ~ time, col = "red", lwd = 2)
lines(X$I ~ time, col = "green", lwd = 2)

legend("topleft", legend = c("Real GDP", "Real Private Consumption", "Real Gross Investment"), 
       col = c("blue", "red", "green"), lty = 1, lwd = 2)

# All three variables exhibit a trend, so they are not stationary.

# 2.b) --------------------------------
aTSA::adf.test(X$Y)
aTSA::adf.test(X$C)
aTSA::adf.test(X$I)

# The adf.test from aTSA tests for the null hypothesis of a unit root and incorporates:
# 1. A linear model without drift or trend
# 2. A linear model with drift but no trend
# 3. A linear model with both drift and trend

# There is a visible trend from the plot in (a),
# so to find out whether the time series is trend-stationary or difference-stationary, type 3 model is used.

# Under type 3:
# - H0 implies that the time series is difference-stationary
# - H1 implies that the time series is trend-stationary

# The p-values from the test output of the third type are all higher than 5%,
# The null at this significance level cannot be rejected. This implies the three series are difference-stationary.

# Therefore, all three series are not stationary.

# 2.c) --------------------------------
# GDP (Y_t) and Consumption (C_t)
# Before the estimation, the cointegrating vector is normalized: beta = (1, -beta_2)'.
# Then beta_2 can be estimated with an OLS regression:
summary(lm(X$Y ~ X$C)) # estimate of X$C is 0.94644 (= beta_2)

# Therefore: beta = (1, -0.946)

# 2.d) --------------------------------
# The non-stationary time series are cointegrated if there is a linear combination of them that is stationary.

# Thus, the interpretation of the cointegrating vector beta is that
#     1 * X$Y - 0.946 * X$C
# is the linear combination of these two series that is stationary.

# This plot confirms that this linear combination is stationary (at least graphically):
plot(time, 1 * X$Y - 0.946 * X$C, type = "l",
     xlab = "Time", main = "Linear combination of the two series: X$Y - 0.946 * X$C")

# Residual-based tests for cointegration (Granger and Engle)
# --> can only be used if there is at most one cointegrating vector

# Engle and Granger propose testing the no-cointegration hypothesis with a unit root test (without a constant and trend),
# using the estimated cointegrating residual û_t:

error <- residuals(lm(X$Y ~ X$C)) # estimated cointegrating residual û_t
aTSA::adf.test(error) # H0: "unit root" can be rejected, thus the residuals are stationary.

# Therefore, the residuals are stationary, meaning the two variables are cointegrated.

# 2.e) --------------------------------
# Lag-order selection

VARselect(y = cbind(X$Y, X$C), type = "const") # VARselect(cbind(X$Y, X$C))
# AIC: 4 --> use AIC
# BIC: 2

# The AIC suggests p = 4 lags for the VAR. Thus, include 4-1 = 3 lags in the VECM.

aic_p <- 4

# 2.f) --------------------------------

# Prepare the variables (slide 6.60)
n <- dim(cbind(X$Y, X$C))[1] # number of observations
gdp <- X$Y
cons <- X$C

# Difference terms
dgdp <- diff(gdp)
dcons <- diff(cons)

# Generate the lagged variables (need 4-1 = 3 lags) 
dgdp_lagged <- sapply(seq(aic_p-1), dplyr::lag, x = dgdp)
colnames(dgdp_lagged) <- paste("dgdp_lag", seq(aic_p-1), sep="")

dcons_lagged <- sapply(seq(aic_p-1), dplyr::lag, x = dcons)
colnames(dcons_lagged) <- paste("dcons_lag", seq(aic_p-1), sep="")

# Error lagged
error_lag1 <- error[1:(n-1)] # error was calculated in (d): residuals(lm(X$Y ~ X$C))

# Create the dataset with all variables and exclude first 3 observations so there are no NAs
diff_dat <- as.data.frame(cbind(dgdp, dcons, dgdp_lagged, dcons_lagged, error_lag1)[aic_p:(n-1), ])

# Estimation of the VECM
vecm.reg1 <- summary(lm(dgdp ~ error_lag1 + dgdp_lag1 + dgdp_lag2 + dgdp_lag3 + dcons_lag1 + dcons_lag2 + dcons_lag3, data = diff_dat))
vecm.reg2 <- summary(lm(dcons ~ error_lag1 + dgdp_lag1 + dgdp_lag2 + dgdp_lag3 + dcons_lag1 + dcons_lag2 + dcons_lag3, data = diff_dat))

# Speed of adjustment: alpha
vecm.reg1 # alpha1 (= estimate for error_lag1) is -0.102053 (significant: **)
vecm.reg2 # alpha2 (= estimate for error_lag1) is 0.011346 (not significantly different from zero)

# The speed of adjustment coefficients capture the reactions of GDP and Consumption to disequilibrium
# (a deviation from the long-term equilibrium).
# Here the condition for convergence is satisfied (one alpha is positive, one negative).
# However, alpha2 is not significantly different from zero.

# In disequilibrium, the series of the first regression (for GDP) is too high and
# needs to decrease (alpha1 is negative) to return to the equilibrium.

# The series of the second regression (for Consumption) is too low (in disequlibrium) and
# needs to increase (alpha2 is positive) to return to the equilibrium.

# The speed of adjustment is higher for alpha1 (corresponding to GDP) because
# its absolute value is higher (0.102 > 0.011). Therefore, GDP reacts faster than Consumption.

# 2.g) --------------------------------
# Lag-order selection

VARselect(X) # now use all three series
# AIC: 4

# No, the optimal lag order for the VAR is still 4 and for the VECM still 4-1 = 3.
aic_p <- 4

# 2.h) --------------------------------

# Johansen's Trace Test
vecm_trace_test <- ca.jo(X, type = "trace", K = aic_p) # K is the lag order of the series (levels) in the VAR
summary(vecm_trace_test)

# The trace test suggests that the rank is 2 (using 5% significance level). (Or r = 1 (using 1% significance level).)

# Johansen's Eigenvalue Test
vecm_eigen_test <- ca.jo(X, type = "eigen", K = aic_p)
summary(vecm_eigen_test)

# The eigenvalue test also suggests that the rank is r = 2 (using 5% significance level).

# 2.i) --------------------------------

# Results with r = 2:
vecm.r2 <- cajorls(vecm_eigen_test, r = 2)

# Estimate for beta
vecm.r2$beta

# Note normalization set: 
# Y=1 & C=0 for the first vector
# Y=0 & C=1 for the second vector

# Beta contains the cointegrating vectors. These columns represents the coefficients of
# the variables that should be used in a linear combination so that it is stationary:

# 1 * GDP - 0.777 * Investment
# 1 * Consumption - 0.830 * Investment

# The existence of this cointegrating vector implies that deviations from the long-term equilibrium will be corrected over time.
# If the variables deviate from their equilibrium relationship, they will adjust in the long run to restore that relationship.

plot(time, 1 * X$Y - 0.777 * X$I, type = "l",
     xlab = "Time", main = "Linear combination: 1 * GDP - 0.777 * Investment")

plot(time, 1 * X$C - 0.83 * X$I, type = "l",
     xlab = "Time", main = "Linear combination: 1 * Consumption - 0.830 * Investment")

# -------------------------------------
# Using 1% significance level
# -------------------------------------

# Results with r = 1:
vecm.r1 <- cajorls(vecm_eigen_test, r = 1)

# Estimate for beta
vecm.r1$beta

# Beta contains the cointegrating vector. This column represents the coefficients of
# the variables that should be used in a linear combination so that it is stationary:

# 1 * GDP - 1.399 * Consumption + 0.385 * Investment

# The existence of this cointegrating vector implies that deviations from the long-term equilibrium will be corrected over time.
# If the variables deviate from their equilibrium relationship, they will adjust in the long run to restore that relationship.

plot(time, 1 * X$Y - 1.399 * X$C + 0.385 * X$I, type = "l",
     xlab = "Time", main = "Linear combination of the three series: 1 * X$Y - 1.399 * X$C + 0.385 * X$I")

# 2.j) --------------------------------

# Pi captures the long-run equilibrium relationships.
vecm_eigen_test@PI

# Impose r=2, Pi can be computed by multiplying alpha with beta':
vecm.r2$rlm$coefficients[1:2,] # alpha'
vecm.r2$beta # beta

# Multiply alpha with beta' to get Pi
t(vecm.r2$rlm$coefficients[1:2,]) %*% t(vecm.r2$beta) # Pi = alpha * beta'

#######################################
######### Ex. 3: Factor models ########
#######################################

rm(list = ls()) 

load("sectors.Rdata")

# First de-mean each variable (two alternative ways)
sectors2 <- sectors - matrix(colMeans(sectors), nrow = dim(sectors)[1], ncol = dim(sectors)[2], byrow = T)
sectors <- apply(sectors, 2, function(x) x - mean(x))

sectors - sectors2

# 3.b) --------------------------------

h <- eigen(var(sectors))$vectors[,1]
h

o_sectors <- order(h)
sorted_sectors <- h[o_sectors]
names(sorted_sectors) <- colnames(sectors)[o_sectors]
sorted_sectors

# The coefficients are positive mostly for the services sector and are negative for the primary sector and for energy.

# 3.c) --------------------------------

# OLS regression of the 13th sectoral shares on the estimated factor loadings without a constant.
# For instance, for the first state in the matrix:
lm(sectors[1,] ~ h - 1)
# Estimate all the values with a single regression:
lm(t(sectors) ~ h - 1)
# The value is positive and large for MA and negative for WY. This means that this states have a very different structure. 
# MA has a strong service sector and Wyoming a strong primary and energy sector.

# 3.e) --------------------------------

plot(eigen(var(sectors))$values/sum(eigen(var(sectors))$values))
# The decline is more or less continuous.
# The share explained by the second factor is large (about 25%).
# The 6th factor still explains a non-negligible part of the variation.
