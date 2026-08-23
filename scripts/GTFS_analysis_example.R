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
  trial.C <- GTFS_outlier(scheme='scheme1',   # 'GTFS', 'GTFS(P)', 'scheme1', 'scheme2' 
                          x,t, 
                          basis='bspline', nbasis=31, nharm=20, 
                          v.prop=0.9,alpha=0.05, q.star='default',
                          n.init=5,iter.max=10,
                          refined = T)   
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




