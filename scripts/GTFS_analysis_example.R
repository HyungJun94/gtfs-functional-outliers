# scripts/run_analysis_simulation.R
# Test and simulation run for checking the GTFS package performance

# 1. Load your official package
library(gtfsOutliers)
library(mvtnorm)
library(fields)

## Settings ---------------------------------------------------------------------
n <- 200            # number of observations
p <- 50             # number of observation points
t <- seq(0, 1, len = p)
c <- 0.025          # contamination ratio
oind <- 196:200     # actual outlier indices


## GTFS test run ----------------------------------------------------------------
x <- simul_6(n,p,t,c)
x <- add_intercept(x, method='normal_adapt')  # add intercept
x <- x$raw

system.time(
  trial.C <- GTFS(x,t,alpha=0.05, alpha.p= 0.2,n.init=10,iter.max=10,
                  basis='bspline', nbasis=31,
                  v.prop=0.9, nharm=15,
                  refined=T, gumbel.cutoff=0.01,
                  coeffs.type = 'algorithm')   # algorithm -> Algorithm 1 & practical -> practical ver.
)

trial.C$outlier
length(trial.C$outlier)
trial.C$step.iter

cat('Original df: ', trial.C$d, ' & Reduced df: ',trial.C$dd, '\n')


## plot
par(mfrow=c(1,2))
plot(trial.C$delta.H, ylab='delta_i(H)', main='delta_i(H) plot')
points(oind, trial.C$delta.H[oind], col='red')
abline(h = trial.C$cutoff, col='blue', lwd=3)

# Check by histogram of delta
hist(trial.C$delta.H,probability = T,nclass=n/5, ylim=c(0, quantile(trial.C$density,0.99)+0.1),
     main='Hist. of delta & actual dist')
lines(seq(0,50, len=300), trial.C$density, lwd=2, col='orange')
abline(v = trial.C$cutoff, lwd=3, col='blue')



## Scheme1 test run -------------------------------------------------------------
x <- simul_1(n,p,t,c)
x <- add_intercept(x, method='normal_adapt')  # add intercept
x <- x$raw

system.time(
  trial <- GTFS_scheme1(x,t,alpha=0.05,n.init=10,iter.max=10,
                        basis='bspline', nbasis=31,
                        v.prop=0.9, nharm=15,
                        refined=T)
)


trial$outlier
length(trial$outlier)


## plot
par(mfrow=c(1,2))
plot(trial$delta.H, ylab='delta_i(H)', main='delta_i(H) plot')
points(oind, trial$delta.H[oind], col='red')
abline(h = trial$cutoff, col='blue', lwd=3)


# Check by histogram of delta
hist(trial$delta.H,probability = T,nclass=n/5, ylim=c(0,max(trial$density)+2),
     main='Hist. of delta & actual dist')
lines(seq(0,2, len=300), trial$density, lwd=3, col='orange')
abline(v = trial$cutoff, lwd=3, col='blue')




## Scheme2 test run -------------------------------------------------------------
x <- simul_31(n,p,t,c)
x <- add_intercept(x, method='normal_adapt')  # add intercept
x <- x$raw

system.time(
  trial.D <- GTFS_scheme2(x,t,alpha=0.05,n.init=10,iter.max=10,
                          basis='bspline', nbasis=31,
                          v.prop=0.9, nharm=15,
                          refined=T)
)

trial.D$outlier
length(trial.D$outlier)

## plot
par(mfrow=c(1,2))
plot(trial.D$delta.H, ylab='delta_i(H)', main='delta_i(H) plot')
points(oind, trial.D$delta.H[oind], col='red')
abline(h = trial.D$cutoff, col='blue', lwd=3)


# Check by histogram of delta
hist(trial.D$delta.H,probability = T,nclass=n/5, ylim=c(0, max(trial.D$density)+0.1),
     main='Hist. of delta & actual dist')
lines(seq(0,50, len=300), trial.D$density, lwd=2, col='orange')
abline(v = trial.D$cutoff, lwd=3, col='blue')








