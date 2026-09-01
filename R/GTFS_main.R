#' Generalized Trimmed Functional Score (GTFS)
#'
#' This file contains the implementation of the GTFS algorithm 
#' for functional shape outlier detection.
#'
#' @author Hyungjun Lim
#' @import fda
#' @import ordinal
#' @name GTFS_core
NULL


#' Preprocessing for GTFS including integral centering
#'
#' @param x N x p functional data
#' @param t sampling time grid points
#' @param basis  smoothing basis ('bspline' or 'fourier')
#' @param nbasis number of basis functions for smoothing
#' 
#' @return A list containing:
#' \item{xfd}{x converted into fd object}
#' \item{x.ic}{Matrix. xfd curves evaluated on measurement points t}
#' @export
gtfs_preprocess <- function(x, t, basis='bspline', nbasis=20){
 
  # scaling (for convenience) 
  x <- x/sd(x)
  
  # Dataset size & settings
  N <- nrow(x)           # number of obs in the data 
  p <- length(t)         # number of observations per functional object
  
  # convert to fd object (basis function selection)
  if(basis=='bspline'){
    basis0 <- create.bspline.basis(c(0,1), nbasis=nbasis)
  }else if(basis=='fourier'){
    basis0 <- create.fourier.basis(c(0,1), nbasis=nbasis)
  }
  
  # smoothness penalty grid search
  lambdas <- 10^(-(1:10))  # regularization parameter
  gcvs <- c()              # Generalized Cross Validation Scores
  for(i in 1:10){
    tmp.fd <- smooth.basis(t,t(x),fdPar(basis0,2, lambdas[i]))
    gcvs <- c(gcvs, mean(tmp.fd$gcv))
  }
  mypar <- fdPar(basis0, 2, lambda=lambdas[which.min(gcvs)]) 
  xfd <- smooth.basis(t,t(x),mypar)  # smoothed functional objects 
  
  
  ###############################################################
  
  # integral centering 
  x.int <- fda::inprod(xfd$fd)                           # vector containing <e_i,1> for all i=1,...,N
  xfd$fd$coefs <- xfd$fd$coefs - matrix(rep(x.int,rep(nbasis,N)),ncol=N)  # update centered e_i coefficients
  
  
  # integral centering
  #rangeval <- xfd$fd$basis$rangeval  # 0 1
  #interval.length <- diff(rangeval)  # 1
  
  #basis.c <- create.constant.basis(rangeval) # constant basis for mean function
  #fd.const <- fd(matrix(1,1,1), basis.c)     # 1 everywhere
  
  #x.int <- fda::inprod(xfd$fd)           # <the ith curve, 1>
  #curve.int <- x.int/interval.length     # divided by range
  #curve.int.fd <- fd(matrix(curve.int, 1,N), basis.c)  # N constant functions
  
  #xfd$fd <- xfd$fd - curve.int.fd   # integral-centering
  
  ###############################################################
  
  # back to matrix form 
  x.ic <- t(eval.fd(t, xfd$fd))
  
  # output 
  out <- list(xfd = xfd, x.ic = x.ic)
  return(out)
}



#' C-statistic Calculation for GTFS
#'
#' @param scheme GTFS variants selection ('GTFS', 'GTFS(P)', 'scheme1', 'scheme2')
#' @param xfd input matrix converted into fd object
#' @param N sample size
#' @param H clean subset 
#' @param nbasis number of basis functions for smoothing
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' @param q.star Gumbel cutoff value described in the manuscript ('default' gives 2/sqrt(N) as described in the manuscript)
#' 
#' @return A list containing:
#' \item{score.H}{Numeric vector of the calculated robust outlyingness scores for all N observations.}
#' \item{fpca}{The raw Functional Principal Component Analysis (FPCA) object returned from \code{pca.fd}.}
#' \item{d}{Integer. The number of principal components explaining the pre-specified proportion of total variation.}
#' \item{d.select}{Integer. The number of filtered/selected eigenfunctions after thresholding with the Gumbel cutoff.}
#' @export
C.calculation <- function(scheme = 'GTFS', 
                          xfd, N, H, 
                          nbasis, nharm, v.prop,
                          q.star=0.1){
  
  # mean centering
  xfd.sub <- xfd$fd[H]                      # final subset fd objects
  mean.sub <- fda::mean.fd(xfd.sub)              # mean function estimate 
  e.c <- xfd                              # e_i = X_i - mu0_hat
  e.c$fd$coefs <- e.c$fd$coefs - matrix(rep(mean.sub$coefs,N), ncol=N)  # update e_i coefficients
  
  # FPCA 
  e.c.pc <- fda::pca.fd(e.c$fd[H], nharm=nharm)              # conduct FPCA
  eval <- e.c.pc$values                                 # eigenvalues
  d <- min(which(cumsum(eval)/sum(eval) >= v.prop))     # d eigenvalues explaining v.prop of total variation
  
  # coefficients
  e.c.scores <- fda::inprod(e.c$fd, e.c.pc$harmonics[1:d])        # N x r matrix of FPC scores. 
  
  # GTFS score calculation for each weighting scheme
  if(scheme %in% c('GTFS', 'GTFS(P)')){
    
    max.scores <- apply(e.c.scores, 2, function(x){max(x^2)})       # max squared score per column 
    max.scores <- max.scores/eval[1:d]                              # max squared score divided by eval
    cN <- 2*ordinal::qgumbel(1-q.star) + 2*log(N) - log(log(N)) - 2*log(gamma(0.5))    # Gumbel cutoff
    
    if( sum(max.scores > cN)  == 0 ){                         # equal-weight fallback
      coeffs <- rep(1/d, d)     
      d.select <- d                                           # d.select: number of filtered eigenfunctions
    }else{
      select.index <- c(max.scores > cN)                      # direction selection
      d.select <- sum(select.index)                           # number of selected directions
      coeffs <- select.index / sum(select.index)              # coefficients
    }
    # evaluate GTFS score 
    score.H <- apply(e.c.scores,1, function(x){temp = x^2/eval[1:d];  return(sum(coeffs* temp))})  
    
    # eigenvalue adjustment described in the manuscript
    medi = qgamma(0.5, shape=d.select/2, scale=2/d.select)        # theoretical median  
    theta <- median(score.H)/medi                                 # adjusting constant
    eval <- eval * theta  
    
    # Re-evaluate score with adjusted eigenvalues
    score.H <- apply(e.c.scores,1, function(x){sum(coeffs*x^2/eval[1:d])})  
    
  }else if(scheme == 'scheme1'){
    
    # Scheme 1 weighting  
    d.select <- d 
    coeffs <- rep(1/d , d)
    #score.H <- apply(e.c.scores,1, function(x){sum(coeffs*x^2)})  
    score.H <- apply(e.c.scores,1, function(x){sum(x^2)})  # xd scale for convenience 
    # no adjustment step for Scheme 1 (eigenvalue not used)
    
  }else if(scheme == 'scheme2'){
    
    # Scheme 2 weighting 
    d.select <- d 
    coeffs <- rep(1/d, d) 
    score.H <- apply(e.c.scores,1, function(x){sum(coeffs*x^2/eval[1:d])})  
    
    # eigenvalue adjustment described in the manuscript
    medi = qgamma(0.5, shape=d/2, scale=2/d)                      # theoretical median  
    theta <- median(score.H)/medi                                 # adjusting constant
    eval <- eval * theta  
    
    # Re-evaluate delta with adjusted eigenvalues
    score.H <- apply(e.c.scores,1, function(x){sum(coeffs*x^2/eval[1:d])}) 
    
  }
  
  # Output calculation 
  out <- list(score.H = score.H, 
              fpca=e.c.pc, 
              d=d, 
              d.select=d.select)
  return(out)
}

#' Iteration step function 
#'
#' @param scheme GTFS variants selection ('GTFS', 'GTFS(P)', 'scheme1', 'scheme2')
#' @param xfd input matrix converted into fd object
#' @param N sample size
#' @param H clean subset 
#' @param h prespecified clean subset size 
#' @param nbasis number of basis functions for smoothing
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' @param q.star Gumbel cutoff value described in the manuscript ('default' gives 2/sqrt(N) as described in the manuscript)
#' @param init.pca initial FPCA estimates 
#' @param init.coeffs initial selection coefficients
#'  
#' @return A list containing:
#' \item{H1}{Integer vector containing the updated subset indices of size h with the smallest sum of delta.}
#' \item{score.H}{Numeric vector of the updated outlyingness scores for all N observations.}
#' \item{coeffs}{Numeric vector of the updated weights/coefficients applied to each eigenfunction score.}
#' @export
Iter.step.C <- function(scheme ='GTFS',
                        xfd, N, H, h,
                        nbasis, nharm=10,v.prop=0.9, q.star=0.1,
                        init.pca, init.coeffs=NULL){
  
  # mean centering with H
  xfd.sub <- xfd$fd[H]              # subset fd objects
  mean.sub <- fda::mean.fd(xfd.sub)      # mean estimate from the subset
  e.c <- xfd                        # mean-centered x subset fd objects   
  e.c$fd$coefs <- e.c$fd$coefs - matrix(rep(mean.sub$coefs,N), ncol=N)  # update centered x subset coefficients
  
  # FPCA components and scores evaluation with H
  eval <- init.pca$values
  efunc <- init.pca$harmonics
  d <- min(which(cumsum(eval)/sum(eval) >= v.prop))    # d eigenvalues explaining v.prop of total variation
  e.c.scores <- fda::inprod(e.c$fd, efunc[1:d])        # N x d FPC scores. 
  
  # C-step
  if(scheme == 'GTFS(P)'  || scheme == 'GTFS' && is.null(init.coeffs)){          
    
    # GTFS(P) : update coefficients (weights) on every H update
    max.scores <- apply(e.c.scores, 2, function(x){max(x^2)})       # max squared score per column 
    max.scores <- max.scores/eval[1:d]                              # max squared score divided by eval
    cN <- 2*ordinal::qgumbel(1-q.star) + 2*log(N) - log(log(N)) - 2*log(gamma(0.5))    # Gumbel cutoff
    if( sum(max.scores > cN)  == 0 ){                         # equal-weight fallback
      coeffs <- rep(1/d, d)     
      d.select <- d                                           # d.select: number of filtered eigenfunctions
    }else{
      select.index <- c(max.scores > cN)                      # direction selection
      d.select <- sum(select.index)                           # number of selected directions
      coeffs <- select.index / sum(select.index)              # coefficients
    }
    # evaluate GTFS score 
    score.H <- apply(e.c.scores,1, function(x){temp = x^2/eval[1:d];  return(sum(coeffs* temp))})
    
  }else if(scheme == 'GTFS'){                   
    
    # GTFS : fixed coefficients (weights) 
    coeffs <- init.coeffs
    # evaluate GTFS score 
    score.H <- apply(e.c.scores,1, function(x){sum(coeffs*x^2/eval[1:d])})  
    
  }else if(scheme == 'scheme1'){
    
    coeffs <- rep(1/d, d)
    # evaluate GTFS score 
    #score.H <- apply(e.c.scores,1, function(x){sum(coeffs*x^2)})
    score.H <- apply(e.c.scores,1, function(x){sum(x^2)}) # xd scale 
    
  }else if(scheme == 'scheme2'){
    
    coeffs <- rep(1/d, d)
    # evaluate GTFS score 
    score.H <- apply(e.c.scores,1, function(x){sum(coeffs*x^2/eval[1:d])})  
    
  }
  
  # Output 
  H1 <- order(score.H)[1:h]                # updated subset (with the smallest sum of delta)
  
  out <-list(H1=H1,  
             score.H=score.H, 
             coeffs= coeffs)
  return(out)
}



#' GTFS calculation function 
#'
#' @param scheme GTFS variants selection ('GTFS', 'GTFS(P)', 'scheme1', 'scheme2')
#' @param x N x p functional data
#' @param t sampling time grid points
#' @param basis  smoothing basis ('bspline' or 'fourier')
#' @param nbasis number of basis functions for smoothing
#' @param nharm number of eigenfunctions to be estimated in FPCA
#' @param v.prop proportion of total variation expected to be explained by d eigenfunctions
#' @param alpha expected false positive rate
#' @param q.star Gumbel cutoff value described in the manuscript ('default' gives 2/sqrt(N) as described in the manuscript)
#' @param n.init number of random multi-starts
#' @param iter.max maximum number of iteration
#' @param refined whether to conduct one-step refinement described in manuscript 
#'  
#' @return A list containing:
#' \item{H}{Integer vector. The final optimized clean subset indices of size h used for robust estimation.}
#' \item{score.H}{Numeric vector. The calculated robust functional outlyingness scores for all N observations.}
#' \item{cutoff}{Numeric. The theoretical outlier detection threshold determined by the Chi-Square distribution.}
#' \item{outlier}{Integer vector. The final detected outlier indices whose scores exceed the cutoff threshold.}
#' \item{density}{Numeric vector. The density estimate for the outlyingness statistics used for empirical distribution tracking.}
#' \item{d}{Integer. The original number of estimated eigenfunctions explaining the pre-specified variation proportion.}
#' \item{d.select}{Integer. The reduced number of filtered eigenfunctions actively used after Gumbel thresholding.}
#' \item{step.iter}{Integer. The number of concentration iterations taken until the subset convergence criteria was met.}
#' @export

GTFS_outlier <- function(scheme='GTFS',
                         x,t, 
                         basis='bspline', nbasis=31, nharm=20, 
                         v.prop=0.9,alpha=0.05, q.star='default',
                         n.init=5,iter.max=10,
                         refined = T){
  
  # scheme : GTFS variants selection 
  # x : N x p matrix of N functional data observed at p grid points
  # t : measurement points 
  # basis : smoothing basis ('bspline' or 'fourier')
  # nbasis : number of basis to be used
  # nharm : number of eigenfunctions to be used
  # v.prop : proportion of total variation explained by eigenfunctions
  # alpha : expected false positive rate
  # q.star : quantile for Gumbel cutoff ('default' : 2/sqrt(N))
  # n.init : number of random initial start
  # iter.max : maximum number of iteration for C-step 
  # refined : refinement described in manuscript
  
  # setting  
  N <- nrow(x)           # sample size
  p <- ncol(x)           # length of measurement points 
  h <- floor(N/2)+1      # size of subset for mean function estimation
  
  # Gumbel cutoff 
  if(q.star == 'default'){
    q.star <- 2/sqrt(N)
  }else(
    q.star <- q.star 
  )
  
  # preprocessing (with Integral centering)
  x.p <- gtfs_preprocess(x, t, basis=basis, nbasis = nbasis)
  xfd <- x.p$xfd    # integral centered fd object
  x.ic <- x.p$x.ic  # integral centered x matrix
  
  # initial FPCA (by MDP)
  MDP.result <- MDP(2,h,x.ic,N,p)                    # initial clean subset via MDP                
  init.pca <- C.calculation(scheme, 
                            xfd, N, MDP.result$Hmdp, 
                            nbasis, nharm, v.prop,
                            q.star)             
  init.pca <- init.pca$fpca                          # initial FPCA via MDP
  
  # clean subset extraction step 
  score.sum.vec <- rep(0,n.init)           # vector for scores
  H.out.mat <- matrix(0,n.init,h)          # vector for candidate clean subset indices 
  n.iter.vec <- rep(0,n.init)              # vector of iteration numbers
  for(k in 1:n.init){                      # multiple initial starting (n.init) 
    rand.init <- sample(1:N, 2)            # random initial subset indices (cardinality=2)
    H.init <- Iter.step.C(scheme,
                          xfd,N,rand.init,h,
                          nbasis, nharm=nharm, v.prop, q.star,
                          init.pca, init.coeffs=NULL)             
    
    H.cand <- H.init$H1                         
    init.coeffs <- H.init$coeffs
    n.iter <- 1                                     # update iteration count
    while(n.iter <= iter.max ){                     # until iteration reaches pre-specified maximum
      H.updated <- Iter.step.C(scheme,
                               xfd, N, H.cand,h,
                               nbasis, nharm, v.prop, q.star,
                               init.pca, init.coeffs)          # update subset indices
      
      if(setequal(H.cand, H.updated$H1)){                      # stop when subset converges
        H.out <- H.updated                                   
        n.iter <- n.iter + 1
        break
      }else{
        H.cand = H.updated$H1                                  
        if(n.iter == iter.max){
          H.out <- H.updated
        }
        n.iter <- n.iter + 1  
      }
    }
    n.iter.vec[k] <- n.iter                 # iteration number save
    H.out.mat[k,] <- H.out$H1               # updated subset indices 
    
    # scores update with H output 
    candidate.eval <- Iter.step.C(scheme,
                                  xfd, N, H.out$H1, h, 
                                  nbasis, nharm, v.prop, q.star,
                                  init.pca, NULL) 
    
    score.sum.vec[k] <- sum(candidate.eval$score.H[H.out$H1])  # updated sum of scores
  }
  
  # last clean subset H 
  ind.H2 <- which.min(score.sum.vec)        # choose the subset with the smallest sum of deltas among the multiple starting results
  H2 <- H.out.mat[ind.H2,]                  # final subset indices
  score.sum <- score.sum.vec[ind.H2]        # final sum of deltas 
  step.iter <- n.iter.vec[ind.H2]           # number of iterations for H2
  
  # statistic calculation 
  C.result <- C.calculation(scheme,
                            xfd, N, H2,
                            nbasis, nharm, v.prop , 
                            q.star)
  
  score.H <- C.result$score.H
  fpca.H <- C.result$fpca
  d.H <- C.result$d
  d.select.H <- C.result$d.select
  
  # Refined procedure
  if(refined==T){
    if(scheme %in% c('GTFS', 'GTFS(P)')){
      cutoff <- qgamma(1-alpha/2, shape = d.select.H/2, scale= 2/d.select.H)
    }else if(scheme == 'scheme1'){
      cutoff.search <- scheme1_cutoff(d.H, fpca.H, alpha/2, score.H, T)
      cutoff <- cutoff.search$cutoff 
    }else if(scheme == 'scheme2'){
      cutoff <- qgamma(1-alpha/2, shape = d.H/2, scale= 2/d.H)
    }
    
    # nonoutliers
    HR <- which(score.H < cutoff)          
    
    C.result <- C.calculation(scheme, 
                              xfd, N, HR,
                              nbasis, nharm, v.prop,
                              q.star)
    
    score.HR <- C.result$score.H
    fpca.HR <- C.result$fpca
    d.HR <- C.result$d
    d.select.HR <- C.result$d.select
    
    # Find cutoff according to pre-specified significance level alpha
    if(scheme %in% c('GTFS', 'GTFS(P)')){
      cutoff.HR <- qgamma(1-alpha, shape = d.select.HR/2, scale= 2/d.select.HR)
      density.HR <- dgamma(seq(0,50,len=300),shape = d.select.HR/2, scale = 2/d.select.HR)
      outlier.HR <- which(score.HR > cutoff.HR) 
    }else if(scheme == 'scheme1'){
      cutoff.search.HR <- scheme1_cutoff(d.HR, fpca.HR, alpha, score.HR, T)
      cutoff.HR <- cutoff.search.HR$cutoff
      density.HR <- cutoff.search.HR$density
      outlier.HR <- which(score.HR > cutoff.HR)
    }else if(scheme == 'scheme2'){
      cutoff.HR <- qgamma(1-alpha, shape = d.HR/2, scale= 2/d.HR)
      density.HR <- dgamma(seq(0,50,len=300), shape = d.HR/2, scale = 2/d.HR)
      outlier.HR <- which(score.HR > cutoff.HR) 
    }
    
    result <- list(H = HR,
                   score.H=score.HR,
                   cutoff=cutoff.HR,
                   outlier=outlier.HR, 
                   density = density.HR, 
                   d=d.HR, 
                   d.select= d.select.HR, 
                   step.iter = step.iter)
  }else{
    
    if(scheme %in% c('GTFS', 'GTFS(P)')){
      cutoff.H <- qgamma(1-alpha, shape = d.select.H/2, scale= 2/d.select.H)
      density.H <- dgamma(seq(0,50,len=300),shape = d.select.H/2, scale = 2/d.select.H)
      outlier.H <- which(score.H > cutoff.H) 
    }else if(scheme == 'scheme1'){
      cutoff.search.H <- scheme1_cutoff(d.H, fpca.H, alpha, score.H, T)
      cutoff.H <- cutoff.search.H$cutoff
      density.H <- cutoff.search.H$density
      outlier.H <- which(score.H > cutoff.H)
    }else if(scheme == 'scheme2'){
      cutoff.H <- qgamma(1-alpha, shape = d.H/2, scale= 2/d.H)
      density.H <- dgamma(seq(0,50,len=300), shape = d.H/2, scale = 2/d.H)
      outlier.H <- which(score.H > cutoff.H) 
    }
    
    result <- list(H = H2,
                   score.H=score.H, 
                   cutoff=cutoff.H,
                   outlier=outlier.H, 
                   density = density.H, 
                   d=d.H, 
                   d.select=d.select.H, 
                   step.iter = step.iter)
  }
  return(result)
}



