#' Cutoff value calculation function
#'
#'
#' @author Hyungjun Lim
#' @import coga
#' @name scheme1_cutoff
NULL

#' Scheme 1 cutoff calculation 
#'
#' @param d number of eigenfunction
#' @param e.c.pc fpca result
#' @param alpha expected false positive rate
#' @param score.H GTFS score 
#' @param adjust eigenvalue adjustment described in the manuscript 
#' 
#' @return A list containing:
#' \item{cutoff}{Numeric. The theoretical outlier detection threshold determined by the Chi-Square distribution.}
#' \item{density}{Numeric vector. The density estimate for the outlyingness statistics used for empirical distribution tracking.}
#' @export
scheme1_cutoff <- function(d, e.c.pc , alpha, score.H, adjust=T){
  
  # d : number of eigenfunctions
  # ec.pc : fpca result
  # alpha : expected false positive rate
  # score.H : outlyingness score estimates
  # adjust : whether to adjust eigenvalues 
  
  eval <- e.c.pc$values
  
  # adjusting constant theta 
  if(adjust==T){
    # finding median
    median_finding <- function(x){pcoga(x, rep(1/2,d), 1/2/eval[1:d]) -0.5}
    medi <- uniroot(median_finding, interval = c(0,100))$root    # theoretical median  
    theta <- median(score.H)/medi                                # adjusting constant
    eval <- eval * theta                                         # adjusted eigenvalues
  }
  
  # find cutoff
  cutoff_finding <- function(x){pcoga(x, rep(1/2,d), 1/2/eval[1:d]) - (1-alpha)}
  cutoff <- uniroot(cutoff_finding, interval = c(0,100))$root
  
  # density
  density <- dcoga(seq(0,2, len=300), rep(1/2, d), 1/2/eval[1:d])      # density values for visualization 
  
  # output
  result <- list(cutoff = cutoff, 
                 density = density)
  return(result)
}

