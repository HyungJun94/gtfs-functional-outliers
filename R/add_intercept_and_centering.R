#' Add Random Intercepts to Synthetic Data
#'
#' @param tmp List. Output from synthetic generator containing raw and smooth matrices.
#' @param method Character. Distribution method for random intercepts ('uniform', 'normal', 'uniform01', 'normal_adapt').
#' @return A list containing modified raw and smooth matrices.
#' @export
add_intercept <- function(tmp, method='uniform') {
  n <- nrow(tmp$raw)
  X <- tmp$raw
  X.s <- tmp$smooth
  
  if(method == 'uniform') {
    X.range <- abs(max(X) - min(X))
    U <- runif(n, -X.range / 2, X.range / 2)
    X <- X + U
    
    X.s.range <- abs(max(X.s) - min(X.s))
    U.s <- runif(n, -X.s.range / 2, X.s.range / 2)
    X.s <- X.s + U.s
  } else if(method == 'normal') {
    RN <- rnorm(n, sd = 2)
    X <- X + RN
    RN.s <- rnorm(n, sd = 2)
    X.s <- X.s + RN.s
  } else if(method == 'uniform01') {
    U <- runif(n, -1, 1)
    X <- X + U
    U.s <- runif(n, -1, 1)
    X.s <- X.s + U.s
  } else if(method == 'normal_adapt') {
    X.range <- abs(max(X) - min(X))
    RN <- rnorm(n, 0, sd = X.range / 4)
    X <- X + RN 
    
    X.s.range <- abs(max(X.s) - min(X.s))
    RN.s <- rnorm(n, 0, sd = X.s.range / 4)
    X.s <- X.s + RN.s
  }
  
  out <- list(raw = X, smooth = X.s)
  return(out)
}

#' Perform Integral Centering on Functional Data
#'
#' @param X Matrix. N x p matrix of functional data.
#' @return Centered functional matrix.
#' @export
integral_centering <- function(X) {
  t(apply(X, 1, function(x) { x - mean(x) }))
}
