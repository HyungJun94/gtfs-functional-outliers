# 00_smoke_test.R
# This script verifies that package functions work correctly on your machine.

# Load required packages
library(gtfsOutliers) 
library(plotly)

## Setting ---------------------------------------------------------------------
n <- 200      # number of observation 
c <- 0.025    # contamination ratio
p <- 50       # number of obs. pt.
t <- seq(0, 1, len = p)

## Execution and Test Plotting -------------------------------------------------
par(mfrow = c(1, 3))

# 1. Generate peak outlier synthetic data (Model 1)
X <- simul_1(n, p, t, c)
plot_mod(X$raw, main = 'Jump Outliers', n = n, c = c)

# 2. Add random intercept
X.RI <- add_intercept(X, method = 'uniform')
plot_mod(X.RI$raw, main = 'With Random Intercept', n = n, c = c)


x <- X$raw

plot(t,x[1,], type='l', ylim= c(min(x), max(x)))
for(i in 1:n){
  lines(t, x[i,], col='grey')
}


temp <- gtfs_preprocess(x,t , 'fourier', 15)

cent <- temp$x.ic

plot(t,cent[1,], type='l', ylim= c(min(cent), max(cent)))
for(i in 1:n){
  lines(t, cent[i,], col='grey')
}

