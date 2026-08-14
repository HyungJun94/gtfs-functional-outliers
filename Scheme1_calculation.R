#' GTFS variant (Scheme1)
#'
#' This file contains the implementation of the GTFS variant 'Scheme1' 
#' for functional shape outlier detection.
#'
#' @author Hyungjun Lim
#' @import fda
#' @import coga 
#' @name scheme1_core
NULL


#' Score calculation
#'
#' @param xfd input matrix converted into fd object
#' @param H subset index
#' @param N sample size
#' @param nbasis number of basis functions for smoothing
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' 
#' @return A list containing the statistic (delta.H), number of incorporated eigenfunctions (r), FPCA raw result (fpca), and integral-centered function (e.c).
#' @export
delta.calculation <- function(xfd, H, N, nbasis, nharm, v.prop){
  # Final delta statistic calculation 
  xfd.sub <- xfd$fd[H]                      # final subset fd objects
  mean.sub <- fda::mean.fd(xfd.sub)              # mean function estimate 
  xfd.c <- xfd                              # e_i = X_i - mu0_hat
  xfd.c$fd$coefs <- xfd.c$fd$coefs - matrix(rep(mean.sub$coefs,N), ncol=N)  # update e_i coefficients
  
  # centering 
  ei.bar <- fda::inprod(xfd.c$fd)                # vector containing <e_i,1> for all i=1,...,N
  e.c <- xfd.c                              # centered e_i = X_i - mu0_hat   
  e.c$fd$coefs <- xfd.c$fd$coefs - matrix(rep(ei.bar,rep(nbasis,N)),ncol=N)  # update centered e_i coefficients
  
  # Evaluate eigenvalues 
  e.c.pc <- fda::pca.fd(e.c$fd[H], nharm=nharm)              # conduct FPCA
  eval <- e.c.pc$values                                 # eigenvalues
  r <- min(which(cumsum(eval)/sum(eval) >= v.prop))     # number of eigenvalues explaining pre-specified proportion of total variation
  
  # Find cutoff D-version
  e.c.scores <- fda::inprod(e.c$fd, e.c.pc$harmonics[1:r])
  delta.H <- apply(e.c.scores,1, function(x){sum(x^2)})  
  
  out <- list(delta.H =delta.H, r=r, fpca=e.c.pc, e.c = e.c)
  return(out)
}


#' Finding a threshold incorporating the chi-square convolution 
#'
#' @param r number of incorporated eigenfunctions
#' @param ec.pc FPCA result obtained from centered data
#' @param alpha expected false positive rate for threshold 
#' @param delta.H calculated scores 
#' @param H subset indice
#' @param adjust eigenvalue adjustment described in the manuscript
#' @param re.adjust additional adjustment option
#' 
#' @return list containing the cutoff value and the density estimate for the scores
#' @export
find_cutoff <- function(r, ec.pc , alpha, delta.H, H, adjust=T, re.adjust=F){
  
  eval <- ec.pc$values
  
  seq0 <- seq(0,1,len=7000)                                      # sequence for a grid search
  
  # adjusting constant theta 
  if(adjust==T){
    # finding median
    i <- 1                                                       # sequence index 
    prob <- 0
    while(prob <= 0.5){
      i = i + 1
      prob = coga::pcoga(seq0[i], rep(1/2, r), 1/2/eval[1:r])  
    }
    medi <- seq0[i]                                              # theoretical median  
    if(re.adjust==T){
      theta <- median(delta.H[H])/medi                           # adjusting constant
    }else{
      theta <- median(delta.H)/medi                              # adjusting constant
    }
    eval <- eval * theta                                         # adjusted eigenvalues
  }
  
  # find cutoff   (coga : R package for convolution of gamma(chi-square) distribution)
  i <- 1 
  prob <- 0                                                
  while(prob <= 1-alpha){
    i = i +1  
    prob = coga::pcoga(seq0[i], rep(1/2, r), 1/2/eval[1:r])                  # find a cutoff point that suits the significance level 
  }
  cutoff <- seq0[i]                                                    # cutoff point
  density <- coga::dcoga(seq(0,2, len=300), rep(1/2, r), 1/2/eval[1:r])      # density values for visualization 
  result <- list(cutoff = cutoff, 
                 density = density)
  return(result)
}


#' Iteration step function 
#'
#' @param xfd input matrix converted into fd object
#' @param N sample size
#' @param H subset index
#' @param h subset size 
#' @param nbasis number of basis functions for smoothing
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param init.pca initial FPCA estimates 
#'  
#' @return A list containing the clean subset indices (H1), sum of calculated scores (sum.delta.H), calculated scores (delta.H).
#' @export
Iter.step.delta <- function(xfd,N,H,h,nbasis,v.prop=0.9, nharm=10, init.pca){

  xfd.sub <- xfd$fd[H]              # subset fd objects
  mean.sub <- fda::mean.fd(xfd.sub)      # mean estimate from the subset
  xfd.c <- xfd                      # centered x subset fd objects   
  xfd.c$fd$coefs <- xfd.c$fd$coefs - matrix(rep(mean.sub$coefs,N), ncol=N)  # update centered x subset coefficients
  
  ei.bar <- fda::inprod(xfd.c$fd)        # vector containing <e_i,1> for all i=1,...,N
  e.c <- xfd.c                      # centered e_i = X_i - mu0_hat 
  e.c$fd$coefs <- xfd.c$fd$coefs - matrix(rep(ei.bar,rep(nbasis,N)),ncol=N)  # update centered e_i coefficients
  
  
  eval <- init.pca$values
  efunc <- init.pca$harmonics
  
  r <- min(which(cumsum(eval)/sum(eval) >= v.prop))  # number of eigenvalues explaining pre-specified proportion of total variation
  e.c.pc.scores <- fda::inprod(e.c$fd, efunc[,1:r])
  if(r == 1){
    delta.H <- (e.c.pc.scores[,1])^2
  }else{
    delta.H <- apply(e.c.pc.scores[,1:r],1, function(x){sum(x^2)})  
  }
  
  H1 <- order(delta.H)[1:h]               # updated subset (with the smallest sum of delta)
  sum.delta.H <- sum(delta.H[H1])         # sum of delta of updated subset  
  list(H1=H1, sum.delta.H=sum.delta.H, delta.H=delta.H)
}





#' GTFS scheme1 calculation function 
#'
#' @param x N x p functional data
#' @param t sampling time grid points
#' @param alpha expected false positive rate for threshold 
#' @param n.init number of random multi-starts
#' @param iter.max maximum number of iteration
#' @param basis smoothing basis choice: 'bspline' or 'fourier'
#' @param nbasis number of basis functions for smoothing
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param refined whether to conduct one-step refinement described in manuscript 
#'  
#' @return A list containing the clean subset (H), calculated statistic (delta.H), indices labeled outliers (outlier), cutoff value (cutoff), density estimate for scores (density).
#' @export
GTFS_scheme1 <- function(x,t,alpha=0.05,n.init=5,iter.max=10,
                         basis='bspline', nbasis=31,
                         v.prop=0.9, nharm=10,
                         refined = T){
  
  # scaling (for convenience)
  x <- x/sd(x)
  
  # Dataset size & settings
  N <- nrow(x)           # number of obs in the data 
  p <- length(t)         # number of observations per functional object
  h <- floor(N/2)+1      # size of subset for mean function estimation
  basis.c <- fda::create.constant.basis(c(0, 1), names="const") # constant basis for mean function approximation
  
  # convert to fd object (basis function selection)
  if(basis=='bspline'){
    basis0 <- fda::create.bspline.basis(c(0,1), nbasis=nbasis)
  }else if(basis=='fourier'){
    basis0 <- fda::create.fourier.basis(c(0,1), nbasis=nbasis)
  }
  lambdas <- 10^(-(1:10))  # regularization parameter 
  gcvs <- c()              # Generalized Cross Validation Scores
  for(i in 1:10){
    tmp.fd <- fda::smooth.basis(t,t(x),fda::fdPar(basis0,2, lambdas[i]))
    gcvs <- c(gcvs, mean(tmp.fd$gcv))
  }
  mypar <- fda::fdPar(basis0,2, lambda=lambdas[which.min(gcvs)]) 
  xfd <- fda::smooth.basis(t,t(x),mypar)  # smoothed functional objects 
  
  # initial eigenfunction 
  MDP.result <- MDP(2,h,x,N,p)                                               # MDP result
  delta.cal <- delta.calculation(xfd, MDP.result$Hmdp, N, nbasis, nharm, 
                                 v.prop=v.prop)                              # calculate delta statistic with MDP subset
  init.pca <- delta.cal$fpca                                                 # FPCA result
  
  # clean subset extraction step 
  delta.sum.vec <- rep(0,n.init)           # vector of outlyingness measure
  H.out.mat <- matrix(0,nrow=n.init,h)     # vector of final subset indices 
  for(k in 1:n.init){                      # multiple initial starting (n.init) 
    rand.init <- sample(1:N, 2, replace = F)                             # random initial subset indices (cardinality=2)
    H.init <- Iter.step.delta(xfd,N, rand.init,h,
                              nbasis, v.prop=v.prop, nharm=nharm,
                              init.pca)                                  # undated subset indices
    H <- H.init$H1 
    n.iter <- 1                                                          # update iteration count
    while(n.iter <= iter.max ){                                          # repeat until the maximum iteration
      H.updated <- Iter.step.delta(xfd,N,rand.init,h,
                                   nbasis, v.prop=v.prop, nharm=nharm,
                                   init.pca)                             # undated subset indices
      if(sum(H == H.updated$H1) == h){                                   # break when subset members don't change by update
        H.out <- H.updated  
        break
      }else{
        H = H.updated$H1                                  
        if(n.iter == iter.max){
          H.out <- H.updated
        }
        n.iter <- n.iter + 1  
      }
    }
    H.out.mat[k,] <- H.out$H1               # updated subset indices 
    delta.sum.vec[k] <- H.out$sum.delta.H   # updated sum of deltas
  }
  ind.H2 <- which.min(delta.sum.vec)        # choose the subset with the smallest sum of deltas among the multiple starting results
  H2 <- H.out.mat[ind.H2,]                  # final subset indices
  delta.sum <- delta.sum.vec[ind.H2]        # final sum of deltas 
  
  # statistic calculation 
  delta.result <- delta.calculation(xfd, H2, N, nbasis, nharm, v.prop=v.prop)
  delta.H <- delta.result$delta.H
  r <- delta.result$r
  e.c <- delta.result$e.c
  ec.pc <- delta.result$fpca
  
  # Find cutoff according to pre-specified significance level alpha
  cutoff.result <- find_cutoff(r, ec.pc, alpha, delta.H, H2,
                               adjust=T)
  cutoff <- cutoff.result$cutoff
  outlier <- which(delta.H > cutoff) 
  
  result <- list(H=H2,delta.H=delta.H, delta.sum=delta.sum,
                 outlier=outlier, cutoff=cutoff,
                 density = cutoff.result$density)
  
  
  # Refined procedure
  if(refined==T){
    cutoff.result <- find_cutoff(r, ec.pc, alpha/2, delta.H, H2,
                                 adjust=T)
    cutoff <- cutoff.result$cutoff
    HR <- which(delta.H < cutoff)           # indices that are classified as non-outlying observations
    
    # calculate statistic
    delta.result.HR <- delta.calculation(xfd, HR, N, nbasis, nharm, v.prop=v.prop)
    delta.HR <- delta.result.HR$delta.H
    r.HR <- delta.result.HR$r
    e.c.HR <- delta.result.HR$e.c
    ec.pc.HR <- delta.result.HR$fpca
    
    # Find cutoff according to pre-specified significance level alpha
    cutoff.result.HR <- find_cutoff(r.HR, ec.pc.HR, alpha, delta.HR, HR,
                                    adjust=T)
    cutoff.HR <- cutoff.result.HR$cutoff
    density.HR <- cutoff.result.HR$density
    outlier.HR <- which(delta.HR > cutoff.HR) 
    
    result <- list(H=HR,delta.H=delta.HR,
                   outlier=outlier.HR, cutoff=cutoff.HR,
                   density = density.HR)
  }
  
  return(result)
}