#' Minimum Diagonal Product (MDP) Algorithm
#'
#' Implements the MDP algorithm for robust initial subset selection in high-dimensional data,
#' as proposed by Ro, Zou, Wang, and Yin (2015).
#'
#' @param m1 Integer. Number of random initial subsets of size 2.
#' @param h Integer. The size of the clean subset (typically \code{floor(N/2) + 1}).
#' @param X Matrix. The \code{N x p} high-dimensional data matrix.
#' @param N Integer. Number of samples (rows of X).
#' @param p Integer. Number of variables/dimensions (columns of X).
#'
#' @return A list containing:
#' \text{Hmdp} {The final robust clean subset index vector of size h.}
#' \text{dis} {The modified Mahalanobis distances calculated using the robust subset.}
#' \text{trR2} {The estimated trace variation metric.}
#'
#' @references Ro, K., Zou, C., Wang, Z., & Yin, G. (2015). Outlier detection for 
#' high-dimensional data. \emph{Biometrika}, 102(3), 589–599.
#'
#' @export
MDP<-function(m1,h,X,N,p){
  
  H0=matrix(0,2,m1) # initial index subsets of size 2 
  H0=apply(H0,2,function(X){sort(sample(N,2))})    #initial subset for MDP
  H0=t(H0) #dimension m1*2
  
  H0_LTS=matrix(0,m1,h)  #save the best MDP subset 
  detD_mdp=matrix(0,m1,1) #save the best MDP detD
  trR2_mdp=matrix(0,m1,1) #save the best  MDP trR2
  dis_mdp=matrix(0,m1,N)
  
  for(i in 1:m1){
    Y=X[H0[i,],]
    Ybar=apply(Y,2,mean) # sample mean with two obs
    S1=cov(Y) # sample covariance of two obs
    D=diag(S1) 
    detD=prod(D) # product of diag terms of S1
    dis=matrix(0,N,1) # modified Mahalanobis dist 
    for (j in 1:N){
      temp2=as.matrix(X[j,]-Ybar)
      dis[j]=t(temp2/D)%*%temp2
    }
    nn=sort(order(dis)[1:h]) # indices of h lowest dist
    crit=100 # 
    
    k=1 # max iteration
    while(crit!=0 & k<10){
      Y=X[nn,]
      Ybar=apply(Y,2,mean) # sample mean of samples with index in [nn]
      S2=cov(Y) # sample covariance function
      D1=diag(S2) # diagonal terms
      
      detD=prod(D1) # product of diagonals
      dis=matrix(0,N,1) # re-calculate the distance
      for (j in 1:N){
        temp2=as.matrix(X[j,]-Ybar)
        dis[j]=t(temp2/D1)%*%temp2
      }
      nn2=sort(order(dis)[1:h]) # new indices with h lowest dist
      crit=sum(abs(nn2-nn)) # number of different indices in [nn] and [nn2]
      nn=nn2 # set [nn2] as a new [nn]
      k=k+1 # add iteration
    }
    ER=cor(X[nn,]) 
    trR2=sum(diag(ER%*%ER))-p^2/h
    
    H0_LTS[i,]=nn 
    detD_mdp[i]=detD 
    trR2_mdp[i]=trR2
    dis_mdp[i,]=dis
    
  }
  
  loc_mdp=which.min(detD_mdp)
  list(Hmdp=H0_LTS[loc_mdp,], # MDP index set (what we need)
       dis=dis_mdp[loc_mdp,],
       trR2=trR2_mdp[loc_mdp])
  
}




#' Refined Minimum Diagonal Product (ReMDP) Algorithm
#'
#' Implements the one-step refinement step for the MDP algorithm (ReMDP)
#' as proposed by Ro, Zou, Wang, and Yin (2015).
#'
#' @param MDP_result List. The output object from the \code{MDP} function.
#' @param alpha Numeric. Significance level for the cutoff threshold.
#' @param X Matrix. The \code{N x p} high-dimensional data matrix.
#' @param N Integer. Number of samples (rows of X).
#' @param p Integer. Number of variables/dimensions (columns of X).
#' @param h Integer. The robust subset size used in the MDP step.
#' @param An Integer. The baseline index separating true normal data and synthetic outliers for performance tracking.
#'
#' @return A list containing:
#' \text{od} {Indices of final detected outliers.}
#' \text{NNod} {Number of detected true outliers (indices less than An).}
#' \text{Uod} {Number of missed true outliers (or false negatives).}
#'
#' @references Ro, K., Zou, C., Wang, Z., & Yin, G. (2015). Outlier detection for 
#' high-dimensional data. \emph{Biometrika}, 102(3), 589–599.
#'
#' @export
ReMDP<-function(MDP_result,alpha,X,N,p,h,An){
  
  delta=alpha/2
  cpn=1+(MDP_result$trR2+p^2/h)/p^1.5	 
  und=p+qnorm(1-delta)*sqrt(2*cpn*MDP_result$trR2)
  scale=median(MDP_result$dis)/p
  od=which(MDP_result$dis>(und*scale))
  nod=length(od)
  
  Y=X[-od,]
  if(nod==0){Y=X}
  newN=N-nod
  Ybar=apply(Y,2,mean)
  
  S1=cov(Y)#
  D=diag(S1)
  detD=prod(D)
  #D=diag(diag(S1))
  #detD=det(D*newN)
  #sD=solve(D)
  
  dis=matrix(0,N,1)
  for (i in 1:N){
    temp2=as.matrix(X[i,]-Ybar)
    dis[i]=t(temp2/D)%*%temp2
  }
  
  ER=cor(Y)
  trRAW=sum(diag(ER%*%ER))
  cpn=1+trRAW/p^1.5
  trR2=trRAW-p^2/newN
  und=p+qnorm(1-alpha)*sqrt(2*cpn*trR2)
  scale=1+sqrt(2*trR2)/(1-delta)/p/sqrt(2*pi)*exp(-qnorm(1-delta)^2/2)
  reod=which(dis>(und*scale))
  
  reNNod=length(reod[reod<An])
  reUod=length(setdiff(An:N,reod))
  
  list(od=reod,NNod=reNNod,Uod=reUod)
  
}


