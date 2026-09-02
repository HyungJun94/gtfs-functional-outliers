#' Least Trimmed Functional Score (LTFS) Benchmark Toolkit
#'
#' This file contains the official benchmark implementation of the LTFS algorithm
#' for functional outlier detection, as proposed by Ren, et. al. (2017) Biometrika.
#'
#' @author Ren, Haojie, Nan Chen, and Changliang Zou
#' @import QRM
#' @import mvoutlier
#' @import MASS
#' @import Matrix
#' @import splines
#' @import fda
#' @name LTFS_benchmark
NULL



#' Concentration Step (C-step) for LTFS Algorithm
#'
#' Implements Step 1 to Step 3 in Algorithm 1 of Ren, Chen, and Zou (2017)
#' to iteratively update the clean subset based on scaled PC scores.
#'
#' @param H0_LTS Integer vector. Current subset indices.
#' @param h Integer. Subset size.
#' @param X Matrix. The \code{N x p} functional/high-dimensional data matrix.
#' @param N Integer. Total sample size.
#' @param ev_cpp List. Contains eigen elements: \code{vectors} (eigenvectors matrix) and \code{values} (eigenvalues vector).
#' @param d Integer. Number of principal component vectors to retain.
#'
#' @return A list containing:
#' \text{H1} {The updated subset indices of size h with lowest distances.}
#' \text{sumDH} {The sum of distances within the updated subset.}
#' \text{DH} {The accumulated distance vector for all N observations.}
#' \text{d} {The number of principal components used.}
#'
#' @references Ren, H., Chen, N., & Zou, C. (2017). Projection-based outlier detection 
#' in functional data. \emph{Biometrika}, 104(2), 411-423.
#'
#' @export
Cstep<-function(H0_LTS,h,X,N,ev_cpp,d){
  C=X[H0_LTS,] # initial subset
  Cbar=apply(C,2,mean) # sample mean of initial subset 
  DH=matrix(0,N,1) 
  for (i in 1:N){
    temp1=0
    for (j in 1:d){
      temp2=(X[i,]-Cbar)%*%ev_cpp$vectors[,j] # PC score
      temp2=temp2**2/ev_cpp$values[j] # scaled PC score
      temp1=temp1+temp2 
    }
    DH[i]=temp1 # D_1(H_0) of ith obs
  }
  H1=order(DH)[1:h]   
  sumDH=sum(DH[H1])
  list(H1=H1,sumDH=sumDH,DH=DH,d=d)
}





#' Least Trimmed Functional Score (LTFS) Algorithm
#'
#' @param mdp Integer vector. Robust initial subset indices (e.g., from MDP algorithm).
#' @param m1 Integer. Number of random initial subsets of size 2.
#' @param m2 Integer. Number of best candidate subsets to keep for the final C-step.
#' @param C Matrix. The \code{N x p} functional/high-dimensional data matrix.
#' @param N Integer. Total sample size.
#' @param h Integer. The final target clean subset size (typically \code{[N/2] + 1}).
#'
#' @return A list containing:
#' \text{Hltfs} {The final optimized LTFS subset indices.}
#' \text{DH} {The robust functional outlyingness score for all N observations.}
#' \text{dltfs} {The number of eigenfunctions explaining up to 90\% of variation.}
#'
#' @references Ren, H., Chen, N., & Zou, C. (2017). Projection-based outlier detection 
#' in functional data. \emph{Biometrika}, 104(2), 411-423.
#'
#' @export
LTFS<-function(mdp,m1,m2,C,N,h){
  # mdp: MDP index set
  # m1: number of initial subsets of size 2
  # m2: 
  # C:
  # N: sample size
  # h: subset size [N/2]+1
  
  cpp=cov(C[mdp,])   
  while(det(cpp)==0){
    a=setdiff(1:N,mdp)[1]
    mdp=sort(c(mdp,a))
    cpp=cov(C[mdp,])
  }
  ev_cpp=eigen(cpp) # eigen decomposition result of MDP sets
  
  d=0
  per=0
  while(per<0.90){
    d=d+1
    per=per+ev_cpp$values[d]/sum(ev_cpp$values)
  }
  
  H0_LTS=matrix(0,2,m1)
  H0_LTS=apply(H0_LTS,2,function(X){sort(sample(N,2))})    #initial subset for MDP
  H0_LTS=t(H0_LTS) #dimension m1*k
  
  H2_LTS=matrix(0,m1,h)
  sumDH2=matrix(0,m1,1)
  for(i in 1:m1){ # repeat Step1~Step2 m1 times
    H1_LTS=Cstep(H0_LTS[i,],h,C,N,ev_cpp,d)$H1
    temp=Cstep(H1_LTS,h,C,N,ev_cpp,d)
    H2_LTS[i,]=temp$H1
    sumDH2[i]=temp$sumDH 
  }
  
  ord_H2=order(sumDH2)[1:m2]
  H3=H2_LTS[ord_H2,]
  sumDH3=sumDH2[ord_H2]
  
  for(i in 1:m2){
    sum_dh0=1000000
    sum_dh1=sumDH3[i]
    h1=H3[i,]
    
    dk=1
    while((sum_dh0-sum_dh1)>0 & dk<=10){
      sum_dh1=sum_dh0
      temp=Cstep(h1,h,C,N,ev_cpp,d)
      h1=temp$H1
      sum_dh0=temp$sumDH
      dk=dk+1
    }
    H3[i,]=h1
    sumDH3[i]=sum_dh1
  }
  tt=which.min(sumDH3)
  Hltfs=H3[tt,]
  
  cpp=cov(C[Hltfs,])
  ev_cpp=eigen(cpp)
  d=0
  per=0
  while(per<0.90){
    d=d+1
    per=per+ev_cpp$values[d]/sum(ev_cpp$values)
  }
  Cbar=apply(C[Hltfs,],2,mean)
  DH=matrix(0,N,1)
  #scores.mat <- matrix(0, N, d)
  
  
  for (i in 1:N){
    temp1=0
    for (j in 1:d){
      temp2=(C[i,]-Cbar)%*%ev_cpp$vectors[,j] 
      temp2=temp2**2/ev_cpp$values[j] 
      
      #scores.mat[i,j] <- temp2**2/ev_cpp$values[j]
      temp1=temp1+temp2
    }
    DH[i]=temp1
  }
  
  list(Hltfs=Hltfs, # LTFS subset index
       DH=DH, # outlyingness measure
       dltfs=d # number of eigen vectors explain up to 90% of total variation
  ) 
  
}

#' One-Step Refined Estimation for LTFS Algorithm
#'
#' @param N Integer. Total sample size.
#' @param od Integer vector. Indices of preliminary outliers detected by LTFS.
#' @param allC Matrix. The complete \code{N x p} functional data matrix.
#' @return A list containing the refined distance statistics and updated dimension d.
#' @export
Refined<-function(N,od,allC){
  nod=length(od)
  C=allC[-od,]
  if(nod==0){C=allC}
  newN=N-nod
  
  cbar=apply(C,2,mean)
  cpp=cov(C)
  ev_cpp=eigen(cpp)
  
  d=0
  per=0
  while(per<0.90){
    d=d+1
    per=per+ev_cpp$values[d]/sum(ev_cpp$values)
  }
  
  stat=matrix(0,N,1)
  for (i in 1:N){
    temp1=0
    for (j in 1:d){
      temp2=(allC[i,]-cbar)%*%ev_cpp$vectors[,j]
      temp2=temp2**2/ev_cpp$values[j]
      temp1=temp1+temp2
    }
    stat[i]=temp1
  }
  list(stat=stat,d=d)
}



#' Full LTFS Pipeline Wrapper for Outlier Detection
#'
#' This master function runs the entire LTFS pipeline including functional basis 
#' expansion transformation (Fourier or B-spline), robust MDP initialization, 
#' core LTFS optimization, and one-step refinement thresholding.
#'
#' @param x Matrix. An \code{N x p} raw functional data matrix.
#' @param N Integer. Total sample size.
#' @param p Integer. Number of time grid evaluation points.
#' @param nbasis Integer. Number of basis functions to be used for smoothing.
#' @param basis Character. Type of basis expansion, either \code{'fourier'} or \code{'bspline'}.
#' @param alpha Numeric. Significance level for outlier thresholding (default is 0.05).
#' @param plot Logical. If \code{TRUE}, generates before-and-after diagnostic diagnostic distance plots.
#'
#' @return A list containing:
#' \item{outlier}{Indices of preliminary outliers detected before refinement.}
#' \item{refined}{Indices of final outliers confirmed after the refinement step.}
#' \item{LTFS_result}{The raw output list from the underlying \code{LTFS} call.}
#'
#' @references Ren, H., Chen, N., & Zou, C. (2017). Projection-based outlier detection 
#' in functional data. \emph{Biometrika}, 104(2), 411-423.
#'
#' @export
LTFS_outlier <- function(x,N,p,nbasis,basis='fourier',alpha=0.05, plot=F){
  
  # x: NxP dataset
  # n: sample size 
  # p: number of time grid points
  # nbasis: number of basis functions to be used
  # alpha: significance level
  
  h=floor(N/2)+1 # subset siez
  
  if(basis=='fourier'){
    fb=fourier(seq(0,1,1/p)[-1],nbasis=nbasis)  #fourier basis functions K=15
    fb_lse=solve(t(fb)%*%fb)%*%t(fb)    #lse: (fb'fb)^{-1}fb'
    C=fb_lse%*%t(x)
    C=t(C)
  }else if(basis=='bspline'){
    bs = create.bspline.basis(c(0,1), nbasis=nbasis) # bspline basis
    bs = eval.basis(seq(0,1,len=p),bs)
    bs_lse=solve(t(bs)%*%bs)%*%t(bs)   
    C=bs_lse%*%t(x)
    C=t(C)
  }
  
  #MDP_result
  m1=100
  MDP_result=MDP(m1,h,x,N,p)
  #LTFS result
  m11=100
  m2=10
  LTFS_result=LTFS(MDP_result$Hmdp,m11,m2,C,N,h)
  
  ##Refined LTFS,delta=alpha/2    
  DHltfs=LTFS_result$DH
  dltfs=LTFS_result$dltfs
  c=median(DHltfs)/qchisq(0.5,dltfs)  #scale parameter of Hltfs
  
  de_und=qchisq(1-alpha/2,dltfs)
  de_od=which(DHltfs>(de_und*c)) #outlier
  
  re_result=Refined(N,de_od,C)
  re_und=qchisq(1-alpha,re_result$d)
  re_c=median(re_result$stat[-de_od])/qchisq(0.5,re_result$d)  
  re_od=which(re_result$stat>(re_und*re_c))
  
  # plot
  if(plot==T){
    par(mfrow=c(1,2))
    aa <- LTFS_result$DH
    plot(1:200, aa, xlab='index', ylab='D(H)' , main='outlyingness distance')
    points(196:200,aa[196:200], col='red')
    
    # refined
    aa <- re_result$stat
    plot(1:200, aa, xlab='index', ylab='D(H)' , main='outlyingness distance')
    points(196:200,aa[196:200], col='red')
  }
  out <- list(outlier=de_od, refined=re_od, LTFS_result=LTFS_result)
  return(out)
}



