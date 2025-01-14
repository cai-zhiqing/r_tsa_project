# Time Series Analysis II

rm(list = ls())
dev.off()

library(readxl)
library(vars)
library(aTSA)
library(expm)
library(ggplot2)

setwd("D:\code_github\r")

#######################################
#### Ex. 1: Advertisement and Sales ###
#######################################

data <- read_excel("advertisement_sales.xlsx", range = "B2:C55", col_names = FALSE)
time <- seq(as.Date("1907-01-01"), as.Date("1960-01-01"), by = "years")

X <- log(data)
colnames(X) <- c("AE", "S")

# 1.a) --------------------------------
plot(X$AE ~ time, type = "l", main = "Advertisement Expenditures")
plot(X$S ~ time, type = "l", main = "Sales")

# Both series display a high dependence on previous values, but they are not exploding processes.
# Suspicion of a unit root.

# Events from that time might have influenced the two series:
# - Economic downturn during the Great Depression of the 1930s
# - World War II (1939 - 1945)

# 1.b) --------------------------------
# Lag order selection
VARselect(X)
# AIC: 4 --> never selects too short a model
# BIC: 1 --> more parsimonious

# Use a VAR(2):
VARselect(X, lag.max = 6)
# HQ: 2
# - HQ is strongly consistent
# - HQ is not asymptotically efficient

# 1.c) --------------------------------
# Estimate reduced-form VAR
VAR <- VAR(X, p = 2, type = "const") # Include a constant
summary(VAR)

# 1.e) --------------------------------
Sigma <- cov(residuals(VAR)) # or var(residuals(VAR))
Sigma

# Alternative: summary(VAR)$covres
# However, there is a small difference in the numbers. See https://stackoverflow.com/questions/23554920/different-covariance-matrices-in-var-model
# It is probably due to a different correction for the degrees of freedom.

# 1.g) --------------------------------
B <- t(chol(Sigma))
B


# 1.i) --------------------------------
# Orthogonalized impulse response function

# Compute companion matrix
Phi <- t(cbind(VAR$varresult$AE$coefficients, 
              VAR$varresult$S$coefficients)[-length(VAR$varresult$AE$coefficients),])

companion <- rbind(Phi, cbind(diag(2), 0, 0))
companion

# Check whether the series are stationary 
# abs(eigen(companion)$values) # all eigenvalues below 1
# 1/abs(eigen(companion)$values) # all roots above 1

# Compute the OIRF
IRF <- matrix(NA, 11, 4)
for (i in 1:11){
  # postmultiply the IRF by B
  IRF[i, 1] <- ((companion %^% (i-1))[(1:2), (1:2)]%*% B )[1, 1] # Variable AE, Shock AE
  IRF[i, 2] <- ((companion %^% (i-1))[(1:2), (1:2)]%*% B )[2, 1] # Variable S, Shock AE
  IRF[i, 3] <- ((companion %^% (i-1))[(1:2), (1:2)]%*% B )[1, 2] # Variable AE, Shock S
  IRF[i, 4] <- ((companion %^% (i-1))[(1:2), (1:2)]%*% B )[2, 2] # Variable S, Shock S
}

plot(IRF[, 1] ~ seq(0, 10), type = "l", xlab = "period", main = "Effect of an advertisement shock on advertisement expenditures")
plot(IRF[, 2] ~ seq(0, 10), type = "l", xlab = "period", main = "Effect of an advertisement shock on sales")
plot(IRF[, 3] ~ seq(0, 10), type = "l", xlab = "period", main = "Effect of a sales shock on advertisement expenditures") # here we can see our zero restriction (b12 = 0)
plot(IRF[, 4] ~ seq(0, 10), type = "l", xlab = "period", main = "Effect of a sales shock on sales")

# compare to
plot(irf(VAR, ortho = TRUE))

# 1.j) --------------------------------
# Forecast error variance decomposition
fevd <- fevd(VAR, n.ahead = 10)
fevd

plot(fevd$AE[, 1], type = "l", main = "Effect of an advertisement shock on advertisement expenditures")
plot(fevd$AE[, 2], type = "l", main = "Effect of a sales shock on advertisement expenditures") # here we can see our zero restriction (b12 = 0)

plot(fevd$S[, 1], type = "l", main = "Effect of an advertisement shock on sales")
plot(fevd$S[, 2], type = "l", main = "Effect of a sales shock on sales")

# Alternative
plot(fevd, col = 1:2)

#######################################
########### Ex. 2: Bootstrap ##########
#######################################

rm(list = ls())
load("housing.Rdata")

# 2.a) --------------------------------
# Assess stationarity
plot(permit$houst, type = "l")
plot(permit$permit, type = "l")
plot(permit$realln, type = "l")
# The variables houst and permit are probably stationary, but realln is not.

# Take the first difference of the log for realln (third column in permit)
Y <- data.frame(permit = permit[-1, 2], houst = permit[-1, 1], realln = diff(log(permit[, 3]))) 
# Note that we have to drop the first observation for the other time series.

# Determine the optimal number of lags
VARselect(Y, lag.max = 8)
LL <- VARselect(Y, lag.max = 8)$selection["SC(n)"] # store the value from BIC

# Estimate the reduced-form VAR with the lag order chosen by BIC

FrequentistVAR <- VAR(Y, p = LL, type = "const")
coef(FrequentistVAR)

# B is the Choleski decomposition of the variance-covariance matrix of the residuals:
B <- t(chol(var(residuals(FrequentistVAR))))
B

# Alternative: OLS by hand
T <- dim(Y)[1]

YY <- as.matrix(Y[(LL+1):T, ])
X <- as.matrix(cbind(1, Y[LL:(T-1), ], Y[(LL-1):(T-2), ], Y[(LL-2):(T-3), ],Y[1:(T-LL), ]))
A <- solve(t(X) %*% X) %*% t(X) %*% YY
t(A)

B2 <- t(chol(var(YY - X %*% A)))

# The two approaches arrive at the same matrix B
round(B - B2, 10)

# Find the companion matrix so that we can compute the IRF

# Using the VAR() results
AA <- t(cbind(FrequentistVAR$varresult$permit$coefficients, 
              FrequentistVAR$varresult$houst$coefficients, 
              FrequentistVAR$varresult$realln$coefficients)[-13,])
# Note the position of the intercept

# Using OLS by hand
AA2  <- t(A)[, -1]                        
# Note the position of the intercept

# The two approaches result in the same coefficient matrix for the companion matrix
AA/ AA2

companion <- rbind(AA, cbind(diag(9), 0, 0, 0)) # the companion matrix will be 12x12
companion

# To be sure that the code is correct, you can compare the IRF(3) to the one computed using the function irf()
irf(FrequentistVAR, ortho = F, n.ahead = 3, response = "permit")$irf # look at the last row only
(companion %^% (3))[(1:3), (1:3)][1,] # They are the same

# Compute the orthogonalized IRF of houst to a permit shock for h = 0, 1, ..., 30.
# Therefore, we only need entry (2, 1).

IRF <- matrix(NA, 31, 1)
for (i in 1:31){
  # we need to postmultiply the IRF by B
  IRF[i] <- ((companion %^% (i-1))[(1:3), (1:3)]%*% B )[2, 1] # second variable (houst), first shock (permit)
}

# 2.b) --------------------------------
# Residual Bootstrap
reps <- 100
boot_IRF <- matrix(NA, 31, reps)

e <- residuals(FrequentistVAR) # We want to sample rows from this 

# Separate the matrix A into one matrix for each lag (these are the same for all bootstrap iterations)
A0 <- A[1,] 
A1 <- A[2:4,]
A2 <- A[5:7,]
A3 <- A[8:10,]
A4 <- A[11:13,]

set.seed(123)
for (r in 1:reps) {
  y_star <- Y*NA # this is the longer vector Y
  y_star[1:LL, ] <- Y[1:LL, ] # set initial conditions 
  y_star <- as.matrix(y_star)
  sample_row <- sample(1:dim(e)[1], size = dim(e)[1], replace = TRUE) # sample of rows
  e_star <- e[sample_row, ] # sample the residuals
  
  for (i in (LL+1):dim(y_star)[1]) { # recursively compute y
    y_star[i,] <- A0 + (y_star[(i-1),]) %*% A1 + (y_star[(i-2),]) %*%  A2 + (y_star[(i-3),]) %*% A3 +  t(y_star[(i-4),]) %*% A4 + e_star[i-LL,]
  }
  
  T2 <- dim(y_star)[1]
  
  # Compute the IRFs as above:
  # Create matrix of xs
  x_star <- as.matrix(cbind(1, y_star[LL:(T-1), ], y_star[(LL-1):(T-2),], y_star[(LL-2):(T-3),],y_star[(1):(T-LL),]))
  # Estimate the VAR(4) and save coefficients
  Aboot <- solve(t(x_star)%*%x_star)%*%t(x_star)%*%y_star[-(1:LL),]
  # Apply cholesky to residuals to derive B
  Bboot <-t(chol(var(y_star[-(1:LL),] - x_star %*% Aboot)))
  # Create companion matrix
  AA <- t(Aboot)[, 2:dim(Aboot)[1]]     # remove the constant 
  companion <- rbind(AA, cbind(diag(9), 0, 0,0))
  # Construct the IRFs
  for (i in 1:31){
    boot_IRF[i,r] <- ((companion %^% (i-1))[(1:3), (1:3)]%*% Bboot)[2,1]
  }
}

# bootstrapped standard errors 
se <- apply(boot_IRF, 1, sd)

# compute CI
lb <- IRF - 1.96 * se 
ub <- IRF + 1.96 * se

# Alternative (percentile bootstrap)
apply(boot_IRF, 1, quantile, c(0.025, 0.975))

# 2.c) --------------------------------
# Pairwise Bootstrap
reps <- 100
boot_IRF <- matrix(NA, 31, reps)
dat <- cbind(YY, X) # matrix containing the entire data, since we want to draw pairs

for (r in 1:reps) {
  sample_row <- sample(nrow(dat), size = dim(dat)[1], replace = TRUE)
  dat_star <- dat[sample_row, ] # sample with replacement
  y_star <- dat_star[, 1:3] 
  x_star <- dat_star[, 4:dim(dat_star)[2]]
  
  
  # Same steps as above to compute the IRF:
  A <- solve(t(x_star)%*%x_star)%*%t(x_star)%*%y_star
  
  B <-t(chol(var(y_star - x_star %*% A)))
  
  AA <- t(A)[, 2:dim(A)[1]]     # remove the constant 
  companion <- rbind(AA, cbind(diag(9), 0, 0,0))
  
  for (i in 1:31){
    boot_IRF[i,r] <- ((companion %^% (i-1))[(1:3), (1:3)]%*% B)[2,1]
  }
}

# bootstrapped standard errors 
se <- apply(boot_IRF, 1, sd)

# compute CI
lb <- IRF - 1.96 * se 
ub <- IRF + 1.96 * se

# 2.d) --------------------------------
dat <- as.data.frame(cbind("Horizon" = 0:30, IRF, lb, ub))

# Plot IRF with CI
ggplot(data = dat, aes(Horizon)) +
  geom_line(aes(y = IRF), linewidth = 0.5, color = "black") +
  geom_line(aes(y = lb), linewidth = 0.5, color = "red", linetype = "dashed") +
  geom_line(aes(y = ub), linewidth = 0.5, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "IRF") +
  theme(axis.title = element_text(face = "bold"))

#######################################
#### Ex. 3: Long-run restrictions #####
#######################################

rm(list = ls())
load("labor.Rdata")

# 3.a) --------------------------------
ggplot(data = dat, aes(Date)) +
  geom_line(aes(y = PROD), size = 0.5, color = "black") +
  geom_line(aes(y = TECHI), size = 0.5, color = "blue") +  
  geom_line(aes(y = (HOURS+7.4)*5), size = 0.5, color = "red") +  
  scale_y_continuous(name = "Labor Productivity and Tech. Index", sec.axis = sec_axis(trans=~./5-7.4, name="Hours Worked")) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  labs(title = "Labor Productivity (Black), Technology Index (Blue) and  Hours Worked (red)") +
  theme(axis.title = element_text(face = "bold"))

# 3.b) --------------------------------
Prod <- dat$PROD
Tech <- dat$TECHI
Hours <- dat$HOURS

# For Prod and Tech, Type 2 of the test below seems relevant. 
# For Hours, Type 1 of the test below seems relevant.

aTSA::adf.test(dat$PROD)
aTSA::adf.test(dat$TECHI)
aTSA::adf.test(dat$HOURS)
# For all three series, we cannot reject that the series has a unit root. 

# 3.c) --------------------------------
# Note that all variables are already in log!
dP <- diff(Prod)
dT <- diff(Tech)
dH <- diff(Hours) 

aTSA::adf.test(dP)
aTSA::adf.test(dT)
aTSA::adf.test(dH)

# 3.d) --------------------------------
maxL <- 10
Y <- cbind(dT, dP, dH)

InfoCrits <- matrix(0, maxL, 2)

for (ii in 1:maxL){
  YY <- Y[(maxL+1-ii):dim(Y)[1], ]
  
  # Estimate VAR
  estVAR <- VAR(YY, p = ii, type = "const")
  print(estVAR$obs)
  # AIC
  InfoCrits[ii, 1] <- AIC(estVAR)
  
  # BIC
  InfoCrits[ii, 2] <- BIC(estVAR)
}

IC <- cbind(seq(1, maxL), InfoCrits)
IC <- data.frame(Lags = IC[, 1], AIC = IC[, 2], BIC = IC[, 3])

ggplot(data = IC, aes(Lags)) +
  geom_line(aes(y = AIC), size = 0.5, color = "red") +
  geom_line(aes(y = BIC), size = 0.5, color = "blue") +
  theme_minimal() +
  ylab("AIC / BIC") +
  labs(title = "Values of Information Criteria for Different Lags", subtitle = "Red = AIC, Blue = BIC") +
  theme(axis.title = element_text(face = "bold")) +
  scale_x_continuous(breaks = seq(1,maxL))

LL <- 2

# 3.e) --------------------------------
# The variables being ordered as Tech, Prod, and Hours.
# In the long run, shocks to labor productivity and hours worked do not affect technology.
# Further, shocks to hours do not affect labor productivity in the long run. 

# In the 3x3 matrix of accumulated long run effects (C in the slides) is
C <- matrix(c("*","*","*",0,"*","*",0,0,"*"),3,3)
C

# 3.f) --------------------------------
Z <- cbind(dT, dP, dH)
T <- nrow(Z)

# Estimate preferred VAR:
Y <- as.matrix(Z[(LL+1):T,])
X <- as.matrix(cbind(1, Z[LL:(T-1), ], Z[(LL-1):(T-2),]))

allA <- solve(t(X)%*%X) %*% t(X)%*%Y
A <- t(allA)[, -1] # remove the constant

Res <- Y - X%*% allA

# Covariance of Residuals
BB <- t(Res)%*%(Res)
BB <- BB/(T-(2*LL)-1) # This is our estimate for Sigma (adjusted for df)
# Alternatively
BB2 <- var(Res)

BB / BB2 # There is a small difference (different DF correction)

# A(1) = (I_m - A_1 - A_2 - ... - A_p)
A1 <- diag(1,3)
for (ii in 1:LL){
  A1 <- A1 - A[1:3,(3*ii-2):(3*ii)]
}

# A(1)^(-1)
A1inv <- solve(A1)

# A(1)^(-1) BB' A(1)^(-1)'
CC <- A1inv %*% BB %*% t(A1inv)

# Identification: restrict three entries in CC to be 0
C <- t(chol(CC))

# Structural effects
Bhat <- A1 %*% C

# construct companion:
bigA <- rbind(A, cbind(diag(3), 0, 0, 0 ))

IRFT <- matrix(NA, LL+24, 3)
IRFP <- IRFT
IRFH <- IRFT
for (i in 1:dim(IRFT)[1]){
  IRFT[i, ] <- ((bigA %^% (i-1))[(1:3), (1:3)]%*% Bhat )[,1]    # column for technology shock
  IRFP[i, ] <- ((bigA %^% (i-1))[(1:3), (1:3)]%*% Bhat )[,2]    # column for productivity 
  IRFH[i, ] <- ((bigA %^% (i-1))[(1:3), (1:3)]%*% Bhat )[,3]    # column for hours shock
}

CIRFT <- rbind(cumsum(IRFT[, 1]), cumsum(IRFT[, 2]), cumsum(IRFT[, 3]))
CIRFP <- rbind(cumsum(IRFP[, 1]), cumsum(IRFP[, 2]), cumsum(IRFP[, 3]))
CIRFH <- rbind(cumsum(IRFH[, 1]), cumsum(IRFH[, 2]), cumsum(IRFH[, 3]))

# IRF Technology ----------------------
Horizon <- seq(0,25)
PlotIRFT <- data.frame(Horizon, t(CIRFT))

ggplot(data = PlotIRFT, aes(Horizon)) +
  geom_line(aes(y = X1), size = 0.5, color = "black") +
  geom_line(aes(y = X2), size = 0.5, color = "blue") +
  geom_line(aes(y = X3), size = 0.5, color = "red") +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  labs(title = "Technology (black), Labor Productivity (blue), Hours (red)") +
  theme(axis.title = element_text(face = "bold"))

# IRF Labor Productivity  -------------
PlotIRFP <- data.frame(Horizon, t(CIRFP))

ggplot(data = PlotIRFP, aes(Horizon)) +
  geom_line(aes(y = X1), size = 0.5, color = "black") +
  geom_line(aes(y = X2), size = 0.5, color = "blue") +
  geom_line(aes(y = X3), size = 0.5, color = "red") +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  labs(title = "Technology (black), Labor Productivity (blue), Hours (red)") +
  theme(axis.title = element_text(face = "bold"))

# IRF Hours  --------------------------
PlotIRFH <- data.frame(Horizon, t(CIRFH))

ggplot(data = PlotIRFH, aes(Horizon)) +
  geom_line(aes(y = X1), size = 0.5, color = "black") +
  geom_line(aes(y = X2), size = 0.5, color = "blue") +
  geom_line(aes(y = X3), size = 0.5, color = "red") +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  labs(title = "Technology (black), Labor Productivity (blue), Hours (red)") +
  theme(axis.title = element_text(face = "bold"))

# Interpretation ----------------------
# In response to a positive investment-specific technology shock, the investment-
# specific technology as well as the labor productivity increase on impact 
# and converge then to their long-run response.

# The investment-specific technology reacts strongly between period 1 to 5. 
# Hours worked shows the typical hump-shaped response which is positive on impact, a
# further strong increase until period six followed by a convergence to the long-run effect.

# A positive labor-productivity shock only increases labor-productivity on impact which is already close to the
# long-run response. Investment-specific technology slightly decreases on impact but the response
# converges fast to zero, which is the condition set. Hours worked decreases slightly on impact, the
# response becomes positive after two quarters and then converges to the long-run effect.

#######################################
# Ex. 4: Stationarity & Cointegration #
#######################################

rm(list = ls())
dev.off()

# 4.a) --------------------------------
A1 <- matrix(c(0.5, 0, 0.2, 1), 2, 2)
abs(eigen(A1)$values)
# The vector is not stationary (slide 6.49)
pi.mat <- -(diag(2) - A1)
# What is the rank of pi?
Matrix::rankMatrix(pi.mat) # check the rank
# The rank is 1, therefore y1 and y2 are cointegrated.
pi.mat
# The speed of adjustment coefficients are -0.5 and 0,
# the cointegrating vector is c(1, -0.4).

# 4.b) --------------------------------
A1 <- matrix(c(0.5,-0.2,0.2,1),2,2)
abs(eigen(A1)$values)
# The vector is stationary.

# 4.c) --------------------------------
A1 <- matrix(c(0.5,0,0.2,0.5),2,2)
A2 <- matrix(c(0.5,0,0.2,0.1),2,2)
bigA <- cbind(rbind(A1,diag(2)),rbind(A2,matrix(0,2,2)))
abs(eigen(bigA)$values)
# The vector is not stationary. There is one unit root.
# To be cointegrated, both y1 and y2 must have a unit root.
# y2 is not affected by y1 and it has no unit root. Thus, y2 is I(0) and y1 is I(1).
# By definition, they are not cointegrated.

# 4.d) --------------------------------
A1 <- matrix(c(0.5,0,0.2,0.5),2,2)
A2 <- matrix(c(0.5,0,0.2,0.5),2,2)
abs(eigen(cbind(rbind(A1,diag(2)),rbind(A2,matrix(0,2,2))))$values)
# Not stationary. There are two unit roots for two variables. Both are I(1).
pi.mat <- -(diag(2)-A1-A2)
Matrix::rankMatrix(pi.mat)
# The rank of pi is 1. They are cointegrated.
pi.mat
# The two variables are integrated with one cointegrating relationship.
# The cointegrating vector beta is (0, 1) and the speed of adjustment is (0.4, 0).
