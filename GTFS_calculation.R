#' Generalized Trimmed Functional Score (GTFS)
#'
#' This file contains the implementation of the GTFS algorithm 
#' for robust functional outlier detection.
#'
#' @author Hyungjun Lim
#' @import fda
#' @import ordinal
#' @name GTFS_core
NULL


#' C-statistic Calculation for GTFS
#'
#' @param xfd input matrix converted into fd object
#' @param H subset index
#' @param N sample size
#' @param nbasis number of basis functions for smoothing
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param adjust eigenvalue adjustment described in the manuscript
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' @param alpha.p Gumbel quantile for selection cutoff
#' @param gumbel.cutoff sample size dependent gumbel cutoff value described in the manuscript
#' 
#' @return A list containing the statistic (delta.H), number of incorporated eigenfunctions (d), FPCA raw result (e.c.pc), and number of selected eigenfunctions (dd).
#' @export
C.calculation <- function(xfd, H, N, nbasis, nharm, adjust=T, re.adjust=F , v.prop, alpha.p=0.2, gumbel.cutoff=0.01){
  # Final delta statistic calculation 
  xfd.sub <- xfd$fd[H]                      # final subset fd objects
  mean.sub <- fda::mean.fd(xfd.sub)         # mean function estimate 
  xfd.c <- xfd                              # e_i = X_i - mu0_hat
  xfd.c$fd$coefs <- xfd.c$fd$coefs - matrix(rep(mean.sub$coefs,N), ncol=N)  # update e_i coefficients
  
  # centering 
  ei.bar <- fda::inprod(xfd.c$fd)           # vector containing <e_i,1> for all i=1,...,N
  e.c <- xfd.c                              # centered e_i = X_i - mu0_hat   
  e.c$fd$coefs <- xfd.c$fd$coefs - matrix(rep(ei.bar,rep(nbasis,N)),ncol=N)  # update centered e_i coefficients
  
  # Evaluate eigenvalues 
  e.c.pc <- fda::pca.fd(e.c$fd[H], nharm=nharm)         # conduct FPCA
  eval <- e.c.pc$values                                 # eigen values
  d <- min(which(cumsum(eval)/sum(eval) >= v.prop))     # number of eigenvalues explaining pre-specified proportion of total variation
  
  # coefficients
  e.c.scores <- fda::inprod(e.c$fd, e.c.pc$harmonics[1:d])        # N x r matrix of FPC scores. 
  max.scores <- apply(e.c.scores, 2, function(x){max(x^2)})       # max squared score per column 
  max.scores <- max.scores/eval[1:d]                              # max squared score divided by eval
  cN <- 2*ordinal::qgumbel(1-alpha.p) + 2*log(N) - log(log(N)) - 2* log(gamma(1-gumbel.cutoff))    # cutoff value
  
  if( sum(max.scores> cN)  == 0){                     
    coeffs <- c(max.scores == max(max.scores))     
    coeffs.ind <- which(coeffs!=0)
    dd <- 1                                                   # dd: number of filtered eigenfunctions
  }else{
    coeffs <- c(max.scores > cN)/ sum(max.scores > cN)
    coeffs.ind <- which(coeffs!=0)
    dd <- sum(max.scores > cN)
  }
  # Find cutoff 
  delta.H <- apply(e.c.scores,1, function(x){temp = x^2/eval[1:d];  return(sum(coeffs* temp))})  
  
  if(adjust == T){
    # adjusting constant theta    
    medi = qgamma(0.5, shape=dd/2, scale=2/dd)                 # theoretical median  
    if(re.adjust == T){ 
      theta <- median(delta.H[H])/medi                         # adjusting constant
    }else{
      theta <- median(delta.H)/medi                            # adjusting constant
    }
    eval <- eval * theta  
  }
  
  # Re-evaluate delta with adjusted eigenvalues
  delta.H <- apply(e.c.scores,1, function(x){sum(coeffs*x^2/eval[1:d])})  
  
  out <- list(delta.H =delta.H, d=d, fpca=e.c.pc, dd=dd)
  return(out)
}


#' Iteration step function 
#'
#' @param xfd input matrix converted into fd object
#' @param H subset index
#' @param N sample size
#' @param h subset size 
#' @param nbasis number of basis functions for smoothing
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param init.pca initial FPCA estimates 
#' @param alpha.p Gumbel quantile for selection cutoff
#' @param gumbel.cutoff sample size dependent gumbel cutoff value described in the manuscript
#' @param init.coeffs initial selection coefficients
#' @param coeffs.type 'practical' for GTFS(P) and 'algorithm' for GTFS
#'  
#' @return A list containing the clean subset indices (H1), sum of calculated statistics (sum.delta.H), calculated statistics (delta.H), and selection coefficients (coeffs)
#' @export
Iter.step.C <- function(xfd,N,H,h,nbasis,v.prop=0.9, nharm=10, init.pca, alpha.p=0.2, gumbel.cutoff=0.01, 
                        init.coeffs=NULL, coeffs.type='practical'){
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
  e.c <- xfd.c                           # centered e_i = X_i - mu0_hat 
  e.c$fd$coefs <- xfd.c$fd$coefs - matrix(rep(ei.bar,rep(nbasis,N)),ncol=N)  # update centered e_i coefficients
  
  eval <- init.pca$values
  efunc <- init.pca$harmonics
  
  d <- min(which(cumsum(eval)/sum(eval) >= v.prop))       # number of eigenvalues explaining pre-specified proportion of total variation
  e.c.pc.scores <- fda::inprod(e.c$fd, efunc[1:d])        # N x r FPC scores. 
  
  if(coeffs.type=='practical'){                           # practical scheme : update coefficients (weights) on every H update
    max.scores <- apply(e.c.pc.scores, 2, function(x){max(x^2)})
    max.scores <- max.scores/eval[1:d]
    cN <- 2*qgumbel(1-alpha.p) + 2*log(N) - log(log(N)) - 2* log(gamma(1-gumbel.cutoff))
    if( sum(max.scores> cN)  == 0){
      coeffs <- c(max.scores == max(max.scores))
    }else{
      coeffs <- c(max.scores > cN)/ sum(max.scores > cN)
    }
    delta.H <- apply(e.c.pc.scores,1, function(x){sum(coeffs*x^2/eval[1:d])})  # divided by eigenvalues
  }else if(coeffs.type == 'algorithm'){                    # algorithm scheme : fix coefficients (weights) regardless of H update
    coeffs <- init.coeffs
    delta.H <- apply(e.c.pc.scores,1, function(x){sum(coeffs*x^2/eval[1:d])})  # divided by eigenvalues
  }
  
  H1 <- order(delta.H)[1:h]               # updated subset (with the smallest sum of delta)
  sum.delta.H <- sum(delta.H[H1])         # sum of delta of updated subset  
  list(H1=H1, sum.delta.H=sum.delta.H, delta.H=delta.H, coeffs= coeffs)
}


#' GTFS calculation function 
#'
#' @param x N x p functional data
#' @param t sampling time grid points
#' @param alpha expected false positive rate for threshold 
#' @param alpha.p Gumbel quantile for selection cutoff
#' @param n.init number of random multi-starts
#' @param iter.max maximum number of iteration
#' @param basis smoothing basis choice: 'bspline' or 'fourier'
#' @param nbasis number of basis functions for smoothing
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param refined whether to conduct one-step refinement described in manuscript 
#' @param gumbel.cutoff sample size dependent gumbel cutoff value described in the manuscript
#' @param coeffs.type 'practical' for GTFS(P) and 'algorithm' for GTFS
#'  
#' @return A list containing the clean subset (H2), calculated statistic (delta.H), sum of calculated statistics (delta.sum), indices labeled outliers (outlier), cutoff value (cutoff), density estimate for statistic (density), number of used and selected eigenfunctions (d & dd), and number of iteration (step.iter).
#' @export
GTFS <- function(x,t,alpha=0.05, alpha.p=0.2,n.init=5,iter.max=10,
                 basis='bspline', nbasis=31,
                 v.prop=0.9, nharm=20,
                 refined = T, gumbel.cutoff=0.01, 
                 coeffs.type='practical'){
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
    tmp.fd <- fda::smooth.basis(t,t(x),fdPar(basis0,2, lambdas[i]))
    gcvs <- c(gcvs, mean(tmp.fd$gcv))
  }
  mypar <- fda::fdPar(basis0,2, lambda=lambdas[which.min(gcvs)]) 
  xfd <- fda::smooth.basis(t,t(x),mypar)  # smoothed functional objects 
  
  # initial eigenfunction 
  MDP.result <- MDP(2,h,x,N,p)                                           # MDP result
  init.pca <- C.calculation(xfd, MDP.result$Hmdp, N, nbasis, nharm, F,
                            v.prop=v.prop,alpha.p, gumbel.cutoff=gumbel.cutoff)   # FPCA
  init.pca <- init.pca$fpca                                              # eigenfunction
  
  # clean subset extraction step 
  delta.sum.vec <- rep(0,n.init)           # vector of outlyingness measure
  H.out.mat <- matrix(0,nrow=n.init,h)     # vector of final subset indices 
  n.iter.vec <- rep(0,n.init)              # vector of interation numbers
  for(k in 1:n.init){                      # multiple initial starting (n.init) 
    rand.init <- sample(1:N, 2)                                        # random initial subset indices (cardinality=2)
    H.init <- Iter.step.C(xfd,N,rand.init,h,
                          nbasis, v.prop=v.prop, nharm=nharm,
                          init.pca,alpha.p, gumbel.cutoff=gumbel.cutoff)     # undated subset indices
    H <- H.init$H1  
    init.coeffs <- H.init$coeffs
    n.iter <- 1                                                        # update iteration count
    while(n.iter <= iter.max ){                                        # until iteration reaches pre-specified maximum
      H.updated <- Iter.step.C(xfd,N,H,h,
                               nbasis, v.prop=v.prop, nharm=nharm,
                               init.pca,alpha.p, gumbel.cutoff=gumbel.cutoff,
                               init.coeffs = init.coeffs, coeffs.type=coeffs.type)             # update subset indices
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
    if(coeffs.type=='algorithm'){                                             # last weight update in case of fixed weight scheme
      H.out <- Iter.step.C(xfd,N,rand.init,h,
                           nbasis, v.prop=v.prop, nharm=nharm,
                           init.pca,alpha.p, gumbel.cutoff=gumbel.cutoff)     # undated subset indices
    }
    n.iter.vec[k] <- n.iter
    H.out.mat[k,] <- H.out$H1               # updated subset indices 
    delta.sum.vec[k] <- H.out$sum.delta.H   # updated sum of deltas
  }
  ind.H2 <- which.min(delta.sum.vec)        # choose the subset with the smallest sum of deltas among the multiple starting results
  H2 <- H.out.mat[ind.H2,]                  # final subset indices
  delta.sum <- delta.sum.vec[ind.H2]        # final sum of deltas 
  step.iter <- n.iter.vec[ind.H2]           # number of iterations for H2
  
  
  # statistic calculation 
  C.result <- C.calculation(xfd, H2, N, nbasis, nharm, T,v.prop=v.prop,alpha.p, gumbel.cutoff=gumbel.cutoff)
  delta.H <- C.result$delta.H
  d <- C.result$d
  dd <- C.result$dd
  
  # Refined procedure
  if(refined==T){
    cutoff <- qgamma(1-alpha, shape = dd/2, scale= 2/dd)
    HR <- which(delta.H < cutoff)           # indices that are classified as non-outlying observations
    
    C.result <- C.calculation(xfd, HR, N, nbasis, nharm, adjust=T, re.adjust =F,v.prop=v.prop,alpha.p, gumbel.cutoff=gumbel.cutoff)
    delta.HR <- C.result$delta.H
    d.HR <- C.result$d
    dd.HR <- C.result$dd
    
    # Find cutoff according to pre-specified significance level alpha
    cutoff.HR <- qgamma(1-alpha, shape = dd.HR/2, scale= 2/dd.HR)
    density.HR <- dgamma(seq(0,50,len=300),shape = dd.HR/2, scale = 2/dd.HR)
    outlier.HR <- which(delta.HR > cutoff.HR) 
    
    result <- list(H=HR,delta.H=delta.HR,
                   outlier=outlier.HR, cutoff=cutoff.HR,
                   density = density.HR, d=d.HR, dd= dd.HR, step.iter = step.iter)
  }else{
    # Find cutoff according to pre-specified significance level alpha
    cutoff <- qgamma(1-alpha, shape = dd/2, scale= 2/dd)
    density <- dgamma(seq(0,50,len=300), shape = dd/2, scale= 2/dd)
    outlier <- which(delta.H > cutoff) 
    
    result <- list(H=H2,delta.H=delta.H, delta.sum=delta.sum,
                   outlier=outlier, cutoff=cutoff,
                   density = density, d=d, dd=dd, step.iter = step.iter)
  }
  return(result)
}



