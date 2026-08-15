#' GTFS variant (Scheme2)
#'
#' This file contains the implementation of the GTFS variant 'Scheme2' 
#' for functional shape outlier detection.
#'
#' @author Hyungjun Lim
#' @import fda
#' @name scheme1_core
NULL


#' Score calculation
#'
#' @param xfd input matrix converted into fd object
#' @param H subset index
#' @param N sample size
#' @param nbasis number of basis functions for smoothing
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param adjust eigenvalue adjustment described in the manuscript
#' @param re.adjust additional adjustment option
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' 
#' @return A list containing:
#' \item{delta.H}{Numeric vector of the calculated robust outlyingness scores for all N observations.}
#' \item{r}{Integer. The number of principal components explaining the pre-specified proportion of total variation.}
#' \item{fpca}{The raw Functional Principal Component Analysis (FPCA) object returned from \code{pca.fd}.}
#' @export
D.calculation <- function(xfd, H, N, nbasis, nharm, adjust=T, re.adjust=F , v.prop){
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
  eval <- e.c.pc$values                                 # eigen values
  r <- min(which(cumsum(eval)/sum(eval) >= v.prop))     # number of eigenvalues explaining pre-specified proportion of total variation
  
  # Find cutoff D-version
  e.c.scores <- fda::inprod(e.c$fd, e.c.pc$harmonics[1:r])
  delta.H <- apply(e.c.scores,1, function(x){sum(x^2/eval[1:r])})  
  
  if(adjust == T){
    # adjusting constant theta 
    medi = qchisq(0.5, d=r)                           # theoretical median  
    if(re.adjust == T){
      theta <- median(delta.H[H])/medi                         # adjusting constant
    }else{
      theta <- median(delta.H)/medi                         # adjusting constant
    }
    eval <- eval * theta  
  }
  
  # Re-evaluate delta with adjusted eigenvalues
  delta.H <- apply(e.c.scores,1, function(x){sum(x^2/eval[1:r])}) 
  
  out <- list(delta.H =delta.H, r=r, fpca=e.c.pc)
  return(out)
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
#' @return A list containing:
#' \item{H1}{Integer vector containing the updated subset indices of size h with the smallest sum of delta.}
#' \item{sum.delta.H}{Numeric. The sum of outlyingness scores within the updated subset H1.}
#' \item{delta.H}{Numeric vector of the updated outlyingness scores for all N observations.}
#' @export
Iter.step.D <- function(xfd,N,H,h,nbasis,v.prop=0.9, nharm=10, init.pca){
  # xfd : x converted to fd object
  # N : sample size 
  # H : subset index
  # h : subset size 
  # basis.c : constant basis function
  
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
  e.c.pc.scores <- fda::inprod(e.c$fd, efunc)
  if(r == 1){
    delta.H <- (e.c.pc.scores[,1])^2/eval[1]
  }else{
    delta.H <- apply(e.c.pc.scores[,1:r],1, function(x){sum(x^2/eval[1:r])})  # divided by eigenvalues
  }
  
  H1 <- order(delta.H)[1:h]               # updated subset (with the smallest sum of delta)
  sum.delta.H <- sum(delta.H[H1])         # sum of delta of updated subset  
  list(H1=H1, sum.delta.H=sum.delta.H, delta.H=delta.H)
}







#' GTFS scheme2 calculation function 
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
#' @return A list containing:
#' \item{H}{Integer vector. The final optimized clean subset indices of size h used for robust estimation.}
#' \item{delta.H}{Numeric vector. The calculated robust functional outlyingness scores for all N observations.}
#' \item{outliers}{Integer vector. The final detected outlier indices whose scores exceed the cutoff threshold.}
#' \item{cutoff}{Numeric. The theoretical outlier detection threshold determined by the Chi-Square distribution.}
#' \item{density}{Numeric vector. The density estimate for the outlyingness statistics used for empirical distribution tracking.}
#' @export
GTFS_scheme2 <- function(x,t,alpha=0.05,n.init=5,iter.max=10,
                         basis='bspline', nbasis=31,
                         v.prop=0.9, nharm=10,
                         refined = T){
  # x : N x P data
  # t : sampling grid points
  # h : subsample size
  # n.init : number of random initial start
  # iter.max : maximum number of iteration
  # basis: orthonormal basis for functional approximation ('bspline'/ 'fourier')
  # v.prop: proportion of variation for d 
  # nharm: number of eigenfunctions to be estimated in FPCA
  # refined: whether to conduct one-step refinement step 
  
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
  MDP.result <- MDP(2,h,x,N,p)                                           # MDP result
  init.pca <- D.calculation(xfd, MDP.result$Hmdp, N, nbasis, nharm, F,
                            v.prop=v.prop)   # FPCA
  init.pca <- init.pca$fpca                                              # eigenfunction
  
  # clean subset extraction step 
  delta.sum.vec <- rep(0,n.init)           # vector of outlyingness measure
  H.out.mat <- matrix(0,nrow=n.init,h)     # vector of final subset indices 
  for(k in 1:n.init){                      # multiple initial starting (n.init) 
    rand.init <- sample(1:N, 2)                                        # random initial subset indices (cardinality=2)
    #t0 <- system.time(
    H.init <- Iter.step.D(xfd,N,rand.init,h,
                          nbasis, v.prop=v.prop, nharm=nharm,
                          init.pca)                                    # undated subset indices
    #)
    #print(t0)
    H <- H.init$H1  
    n.iter <- 1                                                        # update iteration count
    while(n.iter <= iter.max ){                                        # until iteration reaches pre-specified maximum
      #print(n.iter)
      H.updated <- Iter.step.D(xfd,N,H,h,
                               nbasis, v.prop=v.prop, nharm=nharm,
                               init.pca)                               # update subset indices
      if(sum(H == H.updated$H1) == h){                                 # break when subset members don't change by update
        H.out <- H.updated                                   
        n.iter <- n.iter + 1
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
  D.result <- D.calculation(xfd, H2, N, nbasis, nharm, T,v.prop=v.prop)
  delta.H <- D.result$delta.H
  r <- D.result$r
  
  # Refined procedure
  if(refined==T){
    cutoff <- qchisq(1-alpha/2, df=r) 
    HR <- which(delta.H < cutoff)           # indices that are classified as non-outlying observations
    
    D.result <- D.calculation(xfd, HR, N, nbasis, nharm, adjust=T, re.adjust =F,v.prop=v.prop)
    delta.HR <- D.result$delta.H
    r.HR <- D.result$r
    
    # Find cutoff according to pre-specified significance level alpha
    cutoff.HR <- qchisq(1-alpha, df=r.HR)  
    density.HR <- dchisq(seq(0,50,len=300), df= r.HR)
    outlier.HR <- which(delta.HR > cutoff.HR) 
    
    result <- list(H=HR,delta.H=delta.HR,
                   outlier=outlier.HR, cutoff=cutoff.HR,
                   density = density.HR)
  }else{
    # Find cutoff according to pre-specified significance level alpha
    cutoff <- qchisq(1-alpha, df=r)  
    density <- dchisq(seq(0,50,len=300), df= r)
    outlier <- which(delta.H > cutoff) 
    
    result <- list(H=H2,delta.H=delta.H, delta.sum=delta.sum,
                   outlier=outlier, cutoff=cutoff,
                   density = density)
  }
  
  return(result)
}

