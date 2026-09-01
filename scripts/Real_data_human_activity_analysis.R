# Section 4. Real data analysis 

library(fda)

################################################################################
## Human activity recognition ################################

## The whole dataset can be downloaded from the following URL
## (https://archive.ics.uci.edu/dataset/240/human+activity+recognition+using+smartphones)


################################################################################


# human activity data
tab0 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/train/Inertial Signals/body_acc_x_train.txt')
tab1 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/test/Inertial Signals/body_acc_x_test.txt')
tab0 <- rbind(tab0,tab1)
dim(tab0)  # 10299 x 128 



# activity label
lab0 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/train/y_train.txt')
lab1 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/test/y_test.txt')
lab0 <- rbind(lab0,lab1)
dim(lab0)   # 10299 x 1 
table(lab0)




####################################
### Main analysis
####### Outliergram, ISE, Fuctionalboxplot, MUOD, TVD failed
#######  iteration 100 / subsample size 400 / contamination ratio 0.05


set.seed(1)
for(type in 1:3){
  N = 100
  cont = 0.05
  p = 128
  t <- seq(0,1,len=p)
  
  
  N_sim = 10
  N_method = 4
  
  TPR.mat <- matrix(0,nrow=N_sim,ncol=N_method)
  FPR.mat <- matrix(0,nrow=N_sim,ncol=N_method)
  ACC.mat <- matrix(0,nrow=N_sim,ncol=N_method)
  
  
  # Read data
  lab0 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/train/y_train.txt')
  lab1 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/test/y_test.txt')
  lab0 <- rbind(lab0,lab1)
  #dim(dat0)
  
  if(type==1){
    dat0 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/train/Inertial Signals/total_acc_y_train.txt')
    dat1 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/test/Inertial Signals/total_acc_y_test.txt')
    dat0 <- rbind(dat0,dat1)
    
    activity1 = 2
    activity2 = 4
  }else if(type==2){
    dat0 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/train/Inertial Signals/total_acc_x_train.txt')
    dat1 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/test/Inertial Signals/total_acc_x_test.txt')
    dat0 <- rbind(dat0,dat1)
    
    activity1 = 3
    activity2 = 6
  }else if(type==3){
    dat0 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/train/Inertial Signals/total_acc_z_train.txt')
    dat1 <- read.table('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/test/Inertial Signals/total_acc_z_test.txt')
    dat0 <- rbind(dat0,dat1)
    
    activity1 = 1
    activity2 = 6
  }
  
  # iteration subsampling
  for(i in 1:N_sim){
    
    # temporary storage
    TPR.temp <- c()
    FPR.temp <- c()
    ACC.temp <- c()
    method.names <- c()
    
    # sample data
    x <- dat0[sample(which(lab0 ==activity1), N*(1-cont)), ]
    outind <- sample( which(lab0 == activity2) ,N*(cont))
    x <- rbind(x, dat0[outind,])
    x <- as.matrix(x)
    rownames(x) <- 1:N
    
    oind <- (N- (N*cont) +1):N
    
    ########################################################
    
    # GTFS [0.05]
    trial0 <-  GTFS_outlier(scheme='GTFS',
                            x,t, 
                            basis='fourier', nbasis=15, nharm=15, 
                            v.prop=0.9,alpha=0.05, q.star='default',
                            n.init=10,iter.max=10,
                            refined = T)   
    out.temp <- trial0$outlier
    out.temp
    TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
    FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
    ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
    method.names <- c(method.names, 'GTFS 0.05')
    #print('new 0.01')
    
    
    ########################################################
    
    # GTFS - practical [0.05]
    trial1 <-  GTFS_outlier(scheme='GTFS(P)',
                            x,t, 
                            basis='fourier', nbasis=15, nharm=15, 
                            v.prop=0.9,alpha=0.05, q.star='default',
                            n.init=10,iter.max=100,
                            refined = T)   
    out.temp <- trial1$outlier
    out.temp
    TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
    FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
    ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
    method.names <- c(method.names, 'GTFS.P 0.05')
    #print('new 0.01')
    
    
    
    ########################################################
    
    # Least Trimmed Functional Set 
    LTFS1 <- LTFS_outlier(x,N,p,15,'fourier',0.05,F)
    out.temp <- LTFS1$outlier
    TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
    FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
    ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
    method.names <- c(method.names, 'LTFS')
    
    
    
    # MUOD index [Ojo et al., 2021] (magnitude)
    muod1 <- fdaoutlier::muod(x, cut_method = c('boxplot'))
    if(is.vector(muod1)==F){
      out.temp <- NULL
    }else if(is.atomic(muod1$outliers) ==T){
      out.temp <- NULL
    }else{
      out.temp <- muod1$outliers$magnitude
    }
    TPR.temp <- c(TPR.temp, sum(out.temp %in% oind)/(N*cont))
    FPR.temp <- c(FPR.temp, (length(out.temp)-sum(out.temp %in% oind))/(N*(1-cont))  )
    ACC.temp <- c(ACC.temp, (sum(out.temp %in% oind) + (N*(1-cont))- sum(!(out.temp %in% oind)) )/N)
    method.names <- c(method.names, 'MUOD(mag)')
    #print('muod a')
    
    
    
    if(i == 1){models <<- method.names}
    TPR.mat[i,] <- TPR.temp
    FPR.mat[i,] <- FPR.temp
    ACC.mat[i,] <- ACC.temp
    
    print(i) # print progress
    
  }
  
  if(type==1){
    write.csv(TPR.mat, '~/human+activity+recognition+using+smartphones/UCI HAR Dataset/TPR_y24.csv')
    write.csv(FPR.mat, '~/human+activity+recognition+using+smartphones/UCI HAR Dataset/FPR_y24.csv')
  }else if(type==2){
    write.csv(TPR.mat, '~/human+activity+recognition+using+smartphones/UCI HAR Dataset/TPR_x36.csv')
    write.csv(FPR.mat, '~/human+activity+recognition+using+smartphones/UCI HAR Dataset/FPR_x36.csv')
  }else if(type==3){
    write.csv(TPR.mat, '~/human+activity+recognition+using+smartphones/UCI HAR Dataset/TPR_z16.csv')
    write.csv(FPR.mat, '~/human+activity+recognition+using+smartphones/UCI HAR Dataset/FPR_z16.csv')
  }
  
  
}


for(type in 1:3){
  if(type==1){
    TPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/TPR_y24.csv')
    FPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/FPR_y24.csv')
    TPR.mat <- TPR.mat[,-1]
    FPR.mat <- FPR.mat[,-1]
    colnames(TPR.mat)  <- method.names
    colnames(FPR.mat)  <- method.names
    tab.r <- rbind(apply(TPR.mat,2,mean),
                   apply(FPR.mat,2,mean))
    #colnames(tab.r)  <- method.names
    rownames(tab.r)  <- c('TPR', 'FPR')
    cat('\n','Y coord 2 vs 4 : ','\n')
    print(round(tab.r,3))
    
    par(mfrow=c(1,2))
    boxplot(TPR.mat, main='[TPR] Y 2 vs 4')
    boxplot(FPR.mat, main='[FPR] Y 2 vs 4')
    abline(h=0.05, col='blue')
    
    
    
    plot(1:4,apply(TPR.mat,2,mean), ylim= c(0,1), pch=15, main='[TPR] Y 2 vs 4')
    points(1:4, apply(TPR.mat,2,mean) + 2* apply(TPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=25)
    points(1:4, apply(TPR.mat,2,mean) - 2* apply(TPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=24)
    
    
    plot(1:4,apply(FPR.mat,2,mean), ylim= c(0,1), pch=15, main='[FPR] Z 2 vs 4')
    points(1:4, apply(FPR.mat,2,mean) + 2* apply(FPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=25)
    points(1:4, apply(FPR.mat,2,mean) - 2* apply(FPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=24)
    abline(h=0.05, col='blue')
  }else if(type==2){
    if(type==2){
      TPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/TPR_x36.csv')
      FPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/FPR_x36.csv')
      TPR.mat <- TPR.mat[,-1]
      FPR.mat <- FPR.mat[,-1]
      colnames(TPR.mat)  <- method.names
      colnames(FPR.mat)  <- method.names
      tab.r <- rbind(apply(TPR.mat,2,mean),
                     apply(FPR.mat,2,mean))
      tab.r
      #colnames(tab.r)  <- method.names
      rownames(tab.r)  <- c('TPR', 'FPR')
      cat('\n','X 3 vs 6 : ','\n')
      print(round(tab.r,3))
      
      par(mfrow=c(1,2))
      boxplot(TPR.mat, main='[TPR] X 3 vs 6')
      boxplot(FPR.mat, main='[FPR] X 3 vs 6')
      abline(h=0.05, col='blue')
      
      
      
      plot(1:4,apply(TPR.mat,2,mean), ylim= c(0,1), pch=15, main='[TPR] X 3 vs 6')
      points(1:4, apply(TPR.mat,2,mean) + 2* apply(TPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=25)
      points(1:4, apply(TPR.mat,2,mean) - 2* apply(TPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=24)
      
      
      plot(1:4,apply(FPR.mat,2,mean), ylim= c(0,1), pch=15, main='[FPR] X 3 vs 6')
      points(1:4, apply(FPR.mat,2,mean) + 2* apply(FPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=25)
      points(1:4, apply(FPR.mat,2,mean) - 2* apply(FPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=24)
      abline(h=0.05, col='blue')
    }
  }else if(type==3){
    TPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/TPR_z16.csv')
    FPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/FPR_z16.csv')
    TPR.mat <- TPR.mat[,-1]
    FPR.mat <- FPR.mat[,-1]
    colnames(TPR.mat)  <- method.names
    colnames(FPR.mat)  <- method.names
    tab.r <- rbind(apply(TPR.mat,2,mean),
                   apply(FPR.mat,2,mean))
    #colnames(tab.r)  <- method.names
    rownames(tab.r)  <- c('TPR', 'FPR')
    cat('\n','Z 1 vs 6 : ','\n')
    print(round(tab.r,3))
    
    par(mfrow=c(1,2))
    boxplot(TPR.mat, main='[TPR] Z 1 vs 6')
    boxplot(FPR.mat, main='[FPR] Z 1 vs 6')
    abline(h=0.05, col='blue')
    
    
    plot(1:4,apply(TPR.mat,2,mean), ylim= c(0,1), pch=15, main='[TPR] Z 1 vs 6')
    points(1:4, apply(TPR.mat,2,mean) + 2* apply(TPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=25)
    points(1:4, apply(TPR.mat,2,mean) - 2* apply(TPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=24)
    
    
    plot(1:4,apply(FPR.mat,2,mean), ylim= c(0,1), pch=15, main='[FPR] Z 1 vs 6')
    points(1:4, apply(FPR.mat,2,mean) + 2* apply(FPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=25)
    points(1:4, apply(FPR.mat,2,mean) - 2* apply(FPR.mat,2,function(x)(sd(x)/sqrt(N_sim))) , col='red', pch=24)
    abline(h=0.05, col='blue')
    
  }
}







##########################################
# Summarize to table 


TPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/TPR_x36.csv')
FPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/FPR_x36.csv')
TPR.mat <- TPR.mat[,-1]
FPR.mat <- FPR.mat[,-1]
colnames(TPR.mat)  <- c('GTFS', 'GTFS(P)', 'LTFS', 'MUOD')
colnames(FPR.mat)  <- c('GTFS', 'GTFS(P)', 'LTFS', 'MUOD')
tab.r <- cbind(apply(TPR.mat,2,mean),
               apply(FPR.mat,2,mean))

T1 <- TPR.mat %>% 
  gather( key = 'method', value = 'TPR',) %>% 
  mutate( activity = 'X36')
F1 <- FPR.mat %>% 
  gather( key = 'method', value = 'FPR',) %>% 
  mutate( activity = 'X36')


TPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/TPR_y24.csv')
FPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/FPR_y24.csv')
TPR.mat <- TPR.mat[,-1]
FPR.mat <- FPR.mat[,-1]
colnames(TPR.mat)  <- c('GTFS', 'GTFS(P)', 'LTFS', 'MUOD')
colnames(FPR.mat)  <- c('GTFS', 'GTFS(P)', 'LTFS', 'MUOD')
tab.r <- cbind(tab.r,
               apply(TPR.mat,2,mean),
               apply(FPR.mat,2,mean))

T2 <- TPR.mat %>% 
  gather( key = 'method', value = 'TPR',) %>% 
  mutate( activity = 'Y24')
F2 <- FPR.mat %>% 
  gather( key = 'method', value = 'FPR',) %>% 
  mutate( activity = 'Y24')

TPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/TPR_z16.csv')
FPR.mat <- read.csv('~/human+activity+recognition+using+smartphones/UCI HAR Dataset/FPR_z16.csv')
TPR.mat <- TPR.mat[,-1]
FPR.mat <- FPR.mat[,-1]
colnames(TPR.mat)  <- c('GTFS', 'GTFS(P)', 'LTFS', 'MUOD')
colnames(FPR.mat)  <- c('GTFS', 'GTFS(P)', 'LTFS', 'MUOD')
tab.r <- cbind(tab.r,
               apply(TPR.mat,2,mean),
               apply(FPR.mat,2,mean))

T3 <- TPR.mat %>% 
  gather( key = 'method', value = 'TPR',) %>% 
  mutate( activity = 'Z16')
F3 <- FPR.mat %>% 
  gather( key = 'method', value = 'FPR',) %>% 
  mutate( activity = 'Z16')

colnames(tab.r) <- rep(c('TPR', 'FPR'), 3)
tab.r <- round(tab.r,4)
tab.r


library(xtable)
xtable(tab.r, digits=3)


# Visualization

tab.tpr <- rbind(T1,T2,T3)
tab.tpr <- data.frame(tab.tpr)

tab.fpr <- rbind(F1,F2,F3)
tab.fpr <- data.frame(tab.fpr)

# for interval calculation
mean_2se <- function(x){
  temp <- matrix(c(mean(x), 
                   max(0,mean(x)- 2*sd(x)/sqrt(length(x))), 
                   min(1,mean(x)+ 2*sd(x)/sqrt(length(x)))),1,3)
  colnames(temp) <- c('y', 'ymin','ymax')
  return(data.frame(temp))
}




##########################################
# True-positive rate plot
P.tpr.1 <- tab.tpr %>%
  mutate(acc = as.factor(activity) ) %>%
  filter(acc == 'X36') %>%
  ggplot(aes(method, TPR, col=method)) + 
  geom_errorbar(stat = "summary", fun.data = mean_2se, width = 0.8, show.legend=F) +
  geom_point(size = 3, stat = "summary", fun = mean, aes(shape=method), show.legend=F) +
  scale_shape_manual(values = c(15:18)) +
  theme_classic() +
  theme(text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(0,1))+ 
  xlab(NULL) + ylab(NULL)

P.tpr.2 <- tab.tpr %>%
  mutate(acc = as.factor(activity) ) %>%
  filter(acc == 'Y24') %>%
  ggplot(aes(method, TPR, col=method)) + 
  geom_errorbar(stat = "summary", fun.data = mean_2se, width = 0.8, show.legend=F) +
  geom_point(size = 3, stat = "summary", fun = mean, aes(shape=method), show.legend=F) +
  scale_shape_manual(values = c(15:18)) +
  theme_classic() +
  theme(text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(0,1))+ 
  xlab(NULL) + ylab(NULL)

P.tpr.3 <- tab.tpr %>%
  mutate(acc = as.factor(activity) ) %>%
  filter(acc == 'Z16') %>%
  ggplot(aes(method, TPR, col=method)) + 
  geom_errorbar(stat = "summary", fun.data = mean_2se, width = 0.8, show.legend=F) +
  geom_point(size = 3, stat = "summary", fun = mean, aes(shape=method), show.legend=F) +
  scale_shape_manual(values = c(15:18)) +
  theme_classic() +
  theme(text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(0,1))+ 
  xlab(NULL) + ylab(NULL)

P.tpr.1
P.tpr.2
P.tpr.3



# False-positive rate plot

P.fpr.1 <- tab.fpr %>%
  mutate(acc = as.factor(activity) ) %>%
  filter(acc == 'X36') %>%
  ggplot(aes(method, FPR, col=method)) + 
  geom_errorbar(stat = "summary", fun.data = mean_2se, width = 0.8, show.legend=F) +
  geom_point(size = 3, stat = "summary", fun = mean, aes(shape=method), show.legend=F) +
  geom_abline(intercept=0.05,slope=0, col='gray', size=1, alpha=0.8) + 
  scale_shape_manual(values = c(15:18)) +
  theme_classic() +
  theme(text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(0,0.15)) + 
  xlab(NULL) + ylab(NULL)

P.fpr.2 <- tab.fpr %>%
  mutate(acc = as.factor(activity) ) %>%
  filter(acc == 'Y24') %>%
  ggplot(aes(method, FPR, col=method)) + 
  geom_errorbar(stat = "summary", fun.data = mean_2se, width = 0.8, show.legend=F) +
  geom_point(size = 3, stat = "summary", fun = mean, aes(shape=method), show.legend=F) +
  geom_abline(intercept=0.05,slope=0, col='gray', size=1, alpha=0.8) + 
  scale_shape_manual(values = c(15:18)) +
  theme_classic() +
  theme(text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(0,0.15)) + 
  xlab(NULL) + ylab(NULL)

P.fpr.3 <- tab.fpr %>%
  mutate(acc = as.factor(activity) ) %>%
  filter(acc == 'Z16') %>%
  ggplot(aes(method, FPR, col=method)) + 
  geom_errorbar(stat = "summary", fun.data = mean_2se, width = 0.8, show.legend=F) +
  geom_point(size = 3, stat = "summary", fun = mean, aes(shape=method), show.legend=F) +
  geom_abline(intercept=0.05,slope=0, col='gray', size=1, alpha=0.8) + 
  scale_shape_manual(values = c(15:18)) +
  theme_classic() +
  theme(text = element_text(size = 20)) + 
  scale_y_continuous(limits = c(0,0.15)) + 
  xlab(NULL) + ylab(NULL)


P.fpr.1
P.fpr.2
P.fpr.3


##########################################
# Multi-grid plots
plot_grid(P.tpr.1,P.tpr.2,P.tpr.3,
          P.fpr.1,P.fpr.2,P.fpr.3, nrow=2)
