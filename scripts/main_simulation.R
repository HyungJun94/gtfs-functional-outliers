# Main simulation code

######################################################
# packages
required_packages <- c("fda", "mvtnorm", "rainbow", "chemometrics", "fda.usc", "roahd", 
                       "pcaPP", "TeachingDemos", "plyr", "robustbase", "depthTools", 
                       "bootstrap", "cluster", "ks", "mrfDepth", "fdasrvf", "fdaoutlier", 
                       "FUNTA", "coga")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

library(fda)
library(mvtnorm)
library(rainbow)
library(chemometrics) # robust mahalanobis distance
library(fda.usc) # outlier depth trim
library(roahd) # outliergram

library(pcaPP)
library(TeachingDemos)
library(plyr)
library(robustbase)
library(depthTools)
library(bootstrap)
library(cluster)
library(ks)

library(mrfDepth) # DO & FastMUOD
library(fdasrvf) # elastic depth
library(fdaoutlier) # extremal depth & MUOD
library(FUNTA)

library(coga)  # for cutoff of out method

#' Simulation Run Function for Shape Outliers
#'
#' @param N Integer. Sample size.
#' @param p Integer. Number of time points.
#' @param cont Numeric. Contamination ratio.
#' @param N_sim Integer. Number of simulation iterations.
#' @param outlier_type Character. Outlier model scenario ('Peak', 'Jump', etc.)
#' @param intercept Logical. Whether to add random intercept.
#' @export
Sim.shape.outlier <- function(N=200,p=50,cont=0.025,
                              N_sim=100,N_method=15,
                              outlier_type=c('Peak','Jump','Slope1',
                                             'Slope2','Phase','Frequency'),
                              smooth=F,seed=1,
                              intercept=F){
  
  t <- seq(0,1,len=p)  # generate time sequence
  N_type <- length(outlier_type) # number of outliers
  oind = (N*(1-cont)+1):N
  
  for(j in 1:N_type){
    # data storage
    TPR.mat <- matrix(0,nrow=N_sim,ncol=N_method)
    FPR.mat <- matrix(0,nrow=N_sim,ncol=N_method)
    ACC.mat <- matrix(0,nrow=N_sim,ncol=N_method)
    iteration.mat <- matrix(0, nrow=N_sim, ncol=4)
    
    for(i in 1:N_sim){
      # data generation 
      if(outlier_type[j] == 'Peak'){
        X <- simul_1(N,p,t,cont) # Peak
      }else if(outlier_type[j] == 'Jump'){
        X <- simul_2(N,p,t,cont) # Jump
      }else if(outlier_type[j] == 'Slope1'){
        X <- simul_3(N,p,t,cont) # Slope 1
      }else if(outlier_type[j] == 'Slope2'){
        X <- simul_31(N,p,t,cont) # Slope 2
      }else if(outlier_type[j] == 'Phase'){
        X <- simul_4(N,p,t,cont) # Phase
      }else if(outlier_type[j] == 'Frequency'){
        X <- simul_6(N,p,t,cont) # Frequency
      }
      
      # adding random intercept
      if(intercept==T){
        X <- add_intercept(X,method='normal_adapt')
      }
      
      # whether to use smooth synthetic data
      if(smooth==T){
        x <- X$smooth
      }else{
        x <- X$raw
      }
      
      
      # temporary storage
      TPR.temp <- c()
      FPR.temp <- c()
      ACC.temp <- c()
      method.names <- c()
      
      # functional boxplot [Sun and Genton, 2011] (MBD) 
      fbp <- fda::fbplot(t(x), method = 'MBD', plot=F)
      out.temp <- fbp$outpoint # detected outliers
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'FOBox')
      #print('Fbox')
      
      # Robust Mahalanobis distance [Hyndman & Shang, 2010] (not depth)
      rdist <- Moutlier(x,quantile=0.993, plot=F)
      out.temp <- which(rdist$rd> rdist$cutoff)
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'RMD')
      #print('rmd')
      
      # Integrated squared error  [Hyndman and Ullah 2007] (not depth)
      Yfds <- fds(t, t(x), xname = "time", yname = "Simulated value")
      FO.HU <- foutliers(Yfds, method='HUoutliers')
      out.temp <- FO.HU$outliers
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'ISE')
      #print('ise')
      
      # Outliergram [Arribas-Gil & Romo, 2014] (MBD & MEI)
      Yfdata <- fData(t,x)
      og1 <- outliergram(Yfdata, display = F)
      out.temp <- og1$ID_outliers
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'Outgram')
      #print('ogram')
      
      # Total variation depth [Huang & Sun ,2019]
      tvd1 <- fdaoutlier::total_variation_depth(x)
      c <- quantile(tvd1$mss,1/4) - 3*IQR(tvd1$mss)
      out.temp <- which(tvd1$mss <= c)
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'TVD(MSS)')
      #print('tvd')
      
      # MUOD index [Azcorra et al, 2018/Ojo et al., 2021] (shape)
      muod1 <- fdaoutlier::muod(x, cut_method = c('boxplot'))
      if(is.vector(muod1)==F){
        out.temp <- NULL
      }else if(is.atomic(muod1$outliers) ==T){
        out.temp <- NULL
      }else{
        out.temp <- muod1$outliers$shape
      }
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'MUOD(shape)')
      #print('muod s')
      
      # MUOD index [Ojo et al., 2021] (amplitude)
      muod1 <- fdaoutlier::muod(x, cut_method = c('boxplot'))
      if(is.vector(muod1)==F){
        out.temp <- NULL
      }else if(is.atomic(muod1$outliers) ==T){
        out.temp <- NULL
      }else{
        out.temp <- muod1$outliers$shape
      }
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'MUOD(amp)')
      #print('muod a')
      
      ########################################################
      
      # Least Trimmed Functional Set 
      LTFS1 <- LTFS_outlier(x,N,p,15,'bspline',0.05,F)
      out.temp <- LTFS1$outlier
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'LTFS')
      #print('ltfs s')
      
      # Least Trimmed Functional Set [Refined]
      LTFS1 <- LTFS_outlier(x,N,p,15,'bspline',0.05,F)
      out.temp <- LTFS1$refined
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'LTFS.re')
      #print('ltfs r')
      
      ########################################################
      
      gumbel.cutoff = 1/sqrt(N)*2
      
      # GTFS scheme 1 [0.05]
      trial0 <- GTFS_scheme1(x,t,alpha=0.05,n.init=10,iter.max=10,
                             basis='bspline', nbasis=31,
                             v.prop=0.9, nharm=15,
                             refined=T)
      out.temp <- trial0$outlier
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'scheme1 0.05')
      #print('new 0.05')
      
      # GTFS scheme 1 [0.01]
      trial0 <- GTFS_scheme1(x,t,alpha=0.01,n.init=10,iter.max=10,
                             basis='bspline', nbasis=31,
                             v.prop=0.9, nharm=15,
                             refined=T)
      out.temp <- trial0$outlier
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'scheme1 0.01')
      #print('new 0.01')
      
      ########################################################
      
      # GTFS scheme 2 [0.05]
      trial0 <- GTFS_scheme2(x,t,alpha=0.05,n.init=10,iter.max=10,
                             basis='bspline', nbasis=31,
                             v.prop=0.9, nharm=15,
                             refined=T)
      out.temp <- trial0$outlier
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'scheme2 0.05')
      #print('new 0.05')
      
      # GTFS scheme 2 [0.01]
      trial0 <- GTFS_scheme2(x,t,alpha=0.01,n.init=10,iter.max=10,
                             basis='bspline', nbasis=31,
                             v.prop=0.9, nharm=15,
                             refined=T)
      out.temp <- trial0$outlier
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'scheme2 0.01')
      #print('new 0.01')
      
      
      ########################################################
      
      # GTFS (Algorithm 1) [0.05]
      trial0 <- TRIAL.fd.C(x,t,alpha=0.05, alpha.p= gumbel.cutoff,n.init=10,iter.max=10,
                           basis='bspline', nbasis=31,
                           v.prop=0.9, nharm=15,
                           refined=T, 
                           coeffs.type = 'algorithm')
      out.temp <- trial0$outlier
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'GTFS 0.05')
      #print('new 0.01')
      
      iteration.mat[i,1] <- trial0$step.iter
      
      # GTFS (Algorithm 1) [0.01]
      trial0 <- TRIAL.fd.C(x,t,alpha=0.01, alpha.p= gumbel.cutoff,n.init=10,iter.max=10,
                           basis='bspline', nbasis=31,
                           v.prop=0.9, nharm=15,
                           refined=T,
                           coeffs.type = 'algorithm')
      out.temp <- trial0$outlier
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'GTFS 0.01')
      #print('new 0.01')
      
      iteration.mat[i,2] <- trial0$step.iter
      
      ########################################################
      
      # GTFS (practical ver.) [0.05]
      trial0 <- TRIAL.fd.C(x,t,alpha=0.05, alpha.p= gumbel.cutoff,n.init=10,iter.max=10,
                           basis='bspline', nbasis=31,
                           v.prop=0.9, nharm=15,
                           refined=T,
                           coeffs.type = 'practical')
      out.temp <- trial0$outlier
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'GTFS.P 0.05')
      #print('new 0.01')
      
      iteration.mat[i,3] <- trial0$step.iter
      
      
      # GTFS (practical ver.) [0.01]
      trial0 <- TRIAL.fd.C(x,t,alpha=0.01, alpha.p= gumbel.cutoff,n.init=10,iter.max=10,
                           basis='bspline', nbasis=31,
                           v.prop=0.9, nharm=15,
                           refined=T,
                           coeffs.type = 'practical')
      out.temp <- trial0$outlier
      TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
      FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
      ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
      method.names <- c(method.names, 'GTFS.P 0.01')
      #print('new 0.01')
      
      iteration.mat[i,4] <- trial0$step.iter
      
      
      
      
      if(i == 1){models <<- method.names}
      TPR.mat[i,] <- TPR.temp
      FPR.mat[i,] <- FPR.temp
      ACC.mat[i,] <- ACC.temp
      
      print(i) # print progress
    }
    cat('Simulation of [',outlier_type[j],'/N=',N,'/rho=',cont, '] is finished.','\n')
    
    colnames(TPR.mat) <- models
    colnames(FPR.mat) <- models
    colnames(ACC.mat) <- models
    
    # save the result of sim for each outlier type
    write.csv(TPR.mat, paste('sim_raw_trp_',outlier_type[j],'_N',N,'_rho',cont,'_260000.csv', sep=''))
    write.csv(FPR.mat, paste('sim_raw_frp_',outlier_type[j],'_N',N,'_rho',cont,'_260000.csv', sep=''))
    write.csv(ACC.mat, paste('sim_raw_acc_',outlier_type[j],'_N',N,'_rho',cont,'_260000.csv', sep=''))
    write.csv(iteration.mat, paste('sim_raw_iter_',outlier_type[j],'_N',N,'_rho',cont,'_260000.csv', sep=''))
  }
}


###########################################################################
# full simulation

# simulation for rho = 0.025
set.seed(1)
system.time(
  for(N in c(200,600,1000)){
    for(r in c(0.025)){
      Sim.shape.outlier(N=N,p=50,cont=r,
                        N_sim=100,N_method=17,
                        outlier_type=c('Peak','Jump',
                                       'Slope1','Slope2',
                                       'Phase','Frequency'),
                        smooth=F,seed=10, intercept = T)
    }
  }
)



models <- c('FOBox', 'RMD','ISE','Outgram','TVS(MSS)',
            'MUOD(shape)','MUOD(amp)','LTFS','LTFS(RE)',
            'scheme1 0.05', 'scheme1 0.01', 'scheme2 0.05', 'scheme2 0.01',
            'GTFS 0.05', 'GTFS 0.01', 'GTFS.P 0.05', 'GTFS.P 0.01')


###########################################################################
# To table

for(N in c(200,600,1000)){
  for(r in c(0.025,0.05,0.1)){
    outlier_type=c('Peak','Jump','Slope1',
                   'Slope2','Phase','Frequency' )#, 'Local')
    N_type = length(outlier_type)
    models #<- method.names
    N_method = length(models)
    
    
    # summary statistics
    TPR.result <- matrix(0, nrow= N_type,ncol=N_method)
    FPR.result <- matrix(0, nrow= N_type,ncol=N_method)
    rownames(TPR.result)<- outlier_type
    rownames(FPR.result)<- outlier_type
    colnames(TPR.result) <- models
    colnames(FPR.result) <- models
    
    for(j in 1:N_type){
      TPR.temp <- read.csv(paste('sim_raw_trp_',outlier_type[j],'_N',N,'_rho',r,'_260000.csv', sep=''))[,-1]
      FPR.temp <- read.csv(paste('sim_raw_frp_',outlier_type[j],'_N',N,'_rho',r,'_260000.csv', sep=''))[,-1]
      
      TPR.result[j,] <-   round(apply(TPR.temp, 2, mean),5)
      FPR.result[j,] <-   round(apply(FPR.temp, 2, mean),5)
    }
    
    print("##############################################")
    print(paste('Simul: N: ',N,' rho: ',r,  sep=''))
    print(TPR.result)
    print(FPR.result)
  }
}

