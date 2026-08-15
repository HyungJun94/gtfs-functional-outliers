#' Generate Synthetic Peak Outliers (Model 1)
#'
#' This function generates functional synthetic data with peak outliers,
#' replicating the simulation design proposed by Arribas-Gil and Romo (2014).
#' It returns both noise-free raw data and smoothed data using a Matern covariance structure.
#'
#' @param n Integer. Total number of samples.
#' @param p Integer. Number of time points (dimensions).
#' @param t Numeric vector. Time evaluation points.
#' @param c Numeric. Contamination rate (between 0 and 1).
#'
#' @return A list containing two matrices: \code{raw} (noise-free) and \code{smooth} (Matern-smoothed).
#'
#' @references Arribas-Gil, A., & Romo, J. (2014). Shape outlier detection and visualization 
#' for functional data: the outliergram. \emph{Biostatistics}, 15(4), 603-619.
#'
#' @export
simul_1 <- function(n, p, t, c) {
  X <- matrix(0, nrow = n, ncol = p)
  X.s <- matrix(0, nrow = n, ncol = p)
  
  C0 <- abs(outer(t, t, FUN = '-')) # covariance function
  C <- exp(-C0) 
  
  for(i in 1:n) { # normal data
    X[i, ] <- 4 * t + mvtnorm::rmvnorm(1, rep(0, p), sigma = C)
  }
  for(i in (n * (1 - c) + 1):n) { # outlier 
    u <- rbinom(1, 1, 1 / 2)
    mu <- runif(1, 0.25, 0.75)
    X[i, ] <- 4 * t + (-1)^u * 1.8 + (1 / sqrt(2 * pi * 0.01)) * exp(-(t - mu)^2 / 0.02) + mvtnorm::rmvnorm(1, rep(0, p), sigma = C)
  }
  
  # smoothing
  C2 <- apply(C0, c(1, 2), fields::Matern, range = 1, nu = 3 / 2)
  for(i in 1:n) { # normal data
    X.s[i, ] <- 4 * t + mvtnorm::rmvnorm(1, rep(0, p), sigma = C2)
  }
  for(i in (n * (1 - c) + 1):n) { # outlier 
    u <- rbinom(1, 1, 1 / 2)
    mu <- runif(1, 0.25, 0.75)
    X.s[i, ] <- 4 * t + (-1)^u * 1.8 + (1 / sqrt(2 * pi * 0.01)) * exp(-(t - mu)^2 / 0.02) + mvtnorm::rmvnorm(1, rep(0, p), sigma = C2)
  }
  
  out <- list(raw = X, smooth = X.s)
  return(out)
}


#' Generate Synthetic Jump Outliers (Model 2)
#'
#' This function generates functional synthetic data with jump outliers,
#' replicating the simulation design proposed by Dai, Wenlin , et al. (2020).
#' It returns both noise-free raw data and smoothed data using a Matern covariance structure.
#'
#' @param n Integer. Total number of samples.
#' @param p Integer. Number of time points (dimensions).
#' @param t Numeric vector. Time evaluation points.
#' @param c Numeric. Contamination rate (between 0 and 1).
#'
#' @return A list containing two matrices: \code{raw} (noise-free) and \code{smooth} (Matern-smoothed).
#'
#' @references Dai, Wenlin, et al. (2020) "Functional outlier detection and taxonomy by sequential transformations." 
#' \emph{Computational Statistics & Data Analysis} 149, 106960.
#'
#' @export
simul_2 <- function(n,p,t,c){
  X <- matrix(0, nrow=n,ncol=p) # noise-free function
  X.s <- matrix(0, nrow=n,ncol=p) # noise-contaminated function
  C0 = abs(outer(t, t, FUN = '-')) # covariance function
  #C0 = abs(outer(1:p, 1:p, FUN = '-')) # covariance function
  C = exp(-C0) 
  
  for(i in 1:n){
    X[i,] <- 4*t + rmvnorm(1, rep(0,p), sigma = C)
  }
  for(i in (n*(1-c)+1):n){
    U <- runif(1)
    X[i,] <- 4*t + 3*(t>U) + rmvnorm(1, rep(0,p), sigma = C) 
  }
  # smoothing
  C2 <-apply(C0,c(1,2), Matern, range=1, nu=3/2)
  for(i in 1:n){
    X.s[i,] <- 4*t + rmvnorm(1, rep(0,p), sigma = C2)
  }
  for(i in (n*(1-c)+1):n){
    U <- runif(1)
    X.s[i,] <- 4*t + 3*(t>U) + rmvnorm(1, rep(0,p), sigma = C2) 
  }
  
  
  out <- list(raw = X, smooth = X.s)
  return(out)
}


#' Generate Synthetic Slope Outliers (Model 3)
#'
#' This function generates functional synthetic data with slope outliers,
#' replicating the simulation design proposed by Dai, Wenlin , et al. (2020).
#' It returns both noise-free raw data and smoothed data using a Matern covariance structure.
#'
#' @param n Integer. Total number of samples.
#' @param p Integer. Number of time points (dimensions).
#' @param t Numeric vector. Time evaluation points.
#' @param c Numeric. Contamination rate (between 0 and 1).
#'
#' @return A list containing two matrices: \code{raw} (noise-free) and \code{smooth} (Matern-smoothed).
#'
#' @references Dai, Wenlin, et al. (2020) "Functional outlier detection and taxonomy by sequential transformations." 
#' \emph{Computational Statistics & Data Analysis} 149, 106960.
#'
#' @export
simul_3 <- function(n,p,t,c){
  X <- matrix(0, nrow=n,ncol=p) # noise-free function
  X.s <- matrix(0, nrow=n,ncol=p) # noise-contaminated function
  C0 = abs(outer(t, t, FUN = '-'))/0.3 # covariance function
  C = 0.1*exp(-C0) 
  
  for(i in 1:n){ # normal data
    A <- rnorm(1, 0, sd=2)
    B <- rexp(1,1)
    X[i,] <- A + B*atan(t) + rmvnorm(1, rep(0,p), sigma = C) 
  }
  for(i in (n*(1-c)+1):n){ # outlier 
    X[i,] <- 1-2*atan(t) + rmvnorm(1, rep(0,p), sigma = C) 
  }
  # smoothing
  C2 <-apply(C0,c(1,2), Matern, range=1, nu=3/2)
  C2 <- 0.1*C2
  for(i in 1:n){ # normal data
    A <- rnorm(1, 0, sd=2)
    B <- rexp(1,1)
    X.s[i,] <- A +B*atan(t) + rmvnorm(1, rep(0,p), sigma = C2) 
  }
  for(i in (n*(1-c)+1):n){ # outlier 
    X.s[i,] <- 1-2*atan(t) + rmvnorm(1, rep(0,p), sigma = C2) 
  }
  out <- list(raw = X, smooth = X.s)
  return(out)
}


#' Generate Synthetic altered Slope Outliers (Model 3-1)
#'
#' This function generates altered version of the slope outliers proposed by Dai, Wenlin , et al. (2020).
#' Non-outliers are redesigned to share a common mean function.
#' It returns both noise-free raw data and smoothed data using a Matern covariance structure.
#'
#' @param n Integer. Total number of samples.
#' @param p Integer. Number of time points (dimensions).
#' @param t Numeric vector. Time evaluation points.
#' @param c Numeric. Contamination rate (between 0 and 1).
#'
#' @return A list containing two matrices: \code{raw} (noise-free) and \code{smooth} (Matern-smoothed).
#'
#' @references Dai, Wenlin, et al. (2020) "Functional outlier detection and taxonomy by sequential transformations." 
#' \emph{Computational Statistics & Data Analysis} 149, 106960.
#'
#' @export
simul_31 <- function(n,p,t,c){
  X <- matrix(0, nrow=n,ncol=p) # noise-free function
  X.s <- matrix(0, nrow=n,ncol=p) # noise-contaminated function
  C0 = abs(outer(t, t, FUN = '-'))/0.3 # covariance function
  C = 0.1*exp(-C0) 
  
  for(i in 1:n){ # normal data
    A <- rnorm(1, 0, sd=2)
    #B <- rexp(1,1)
    X[i,] <- A + 2*atan(t) + rmvnorm(1, rep(0,p), sigma = C) 
  }
  for(i in (n*(1-c)+1):n){ # outlier 
    X[i,] <- 1-2*atan(t) + rmvnorm(1, rep(0,p), sigma = C) 
  }
  # smoothing
  C2 <-apply(C0,c(1,2), Matern, range=1, nu=3/2)
  C2 <- 0.1*C2
  for(i in 1:n){ # normal data
    A <- rnorm(1, 0, sd=2)
    B <- rexp(1,1)
    X.s[i,] <- A +2*atan(t) + rmvnorm(1, rep(0,p), sigma = C2) 
  }
  for(i in (n*(1-c)+1):n){ # outlier 
    X.s[i,] <- 1-2*atan(t) + rmvnorm(1, rep(0,p), sigma = C2) 
  }
  out <- list(raw = X, smooth = X.s)
  return(out)
}


#' Generate Synthetic Phase Outliers (Model 4)
#'
#' This function generates functional synthetic data with phase variation-induced outliers,
#' replicating the simulation design proposed by Arribas-Gil and Romo (2014).
#' It returns both noise-free raw data and smoothed data using a Matern covariance structure.
#'
#' @param n Integer. Total number of samples.
#' @param p Integer. Number of time points (dimensions).
#' @param t Numeric vector. Time evaluation points.
#' @param c Numeric. Contamination rate (between 0 and 1).
#'
#' @return A list containing two matrices: \code{raw} (noise-free) and \code{smooth} (Matern-smoothed).
#'
#' @references Arribas-Gil, A., & Romo, J. (2014). Shape outlier detection and visualization 
#' for functional data: the outliergram. \emph{Biostatistics}, 15(4), 603-619.
#'
#' @export
simul_4 <- function(n,p,t,c){
  X <- matrix(0, nrow=n,ncol=p) # noise-free function
  X.s <- matrix(0, nrow=n,ncol=p) # noise-free function
  
  C0 = abs(outer(t, t, FUN = '-'))/0.3 # covariance function
  C = 0.3*exp(-C0) 
  
  for(i in 1:n){ # normal functions (first n(1-c) functions)
    X[i,] <- 30*t*(1-t)^(3/2) + rmvnorm(1, rep(0,p), sigma = C) # original curve
  }
  for(i in (n*(1-c)+1):n){ # outliers (last n*c functions)
    X[i,] <- 30*t^(3/2)*(1-t) + rmvnorm(1, rep(0,p), sigma = C) 
  }
  # smoothing
  C2 <-apply(C0,c(1,2), Matern, range=1, nu=3/2)
  C2 <- 0.3*C2
  for(i in 1:n){ # normal functions (first n(1-c) functions)
    X.s[i,] <- 30*t*(1-t)^(3/2) + rmvnorm(1, rep(0,p), sigma = C2) # original curve
  }
  for(i in (n*(1-c)+1):n){ # outliers (last n*c functions)
    X.s[i,] <- 30*t^(3/2)*(1-t) + rmvnorm(1, rep(0,p), sigma = C2) 
  }
  
  out <- list(raw = X, smooth = X.s)
  return(out)
}

#' Generate Synthetic Frequency Outliers (Model 6)
#'
#' This function generates functional synthetic data with frequency outliers,
#' inspired by the simulation design proposed by Harris, Trevor, et al.
#' It returns both noise-free raw data and smoothed data using a Matern covariance structure.
#'
#' @param n Integer. Total number of samples.
#' @param p Integer. Number of time points (dimensions).
#' @param t Numeric vector. Time evaluation points.
#' @param c Numeric. Contamination rate (between 0 and 1).
#'
#' @return A list containing two matrices: \code{raw} (noise-free) and \code{smooth} (Matern-smoothed).
#'
#' @references Harris, Trevor, et al. (2021) "Elastic depths for detecting shape anomalies in functional data." 
#' \emph{Technometrics} 63.4, 466-476.
#'
#' @export
simul_6 <- function(n,p,t,c){
  X <- matrix(0, nrow=n,ncol=p) # noise-free function
  X.s <- matrix(0, nrow=n,ncol=p) # noise-contaminated function
  
  #C0 = abs(outer(1:p, 1:p, FUN = '-')) # covariance function
  C0 = abs(outer(t, t, FUN = '-')) # covariance function
  C1 <<- exp(-C0/0.5)
  
  for(i in 1:n){ # normal data
    X[i,] <- sin(2*pi*t) + 4*t + rmvnorm(1, rep(0,p), sigma = C1)  + rnorm(1)
  }
  for(i in (n*(1-c)+1):n){ # outlier
    X[i,] <- sin(12*pi*t) + 4*t + rmvnorm(1, rep(0,p), sigma = C1)  + rnorm(1)
  }
  # smoothing
  C2 <-apply(C0,c(1,2), Matern, range=1, nu=3/2)
  for(i in 1:n){ # normal data
    X.s[i,] <- sin(2*pi*t) + 4*t + rmvnorm(1, rep(0,p), sigma = C2)  + rnorm(1)
  }
  for(i in (n*(1-c)+1):n){ # outlier
    X.s[i,] <- sin(12*pi*t) + 4*t + rmvnorm(1, rep(0,p), sigma = C2)  + rnorm(1)
  }
  
  out <- list(raw = X, smooth = X.s)
  return(out)
}


