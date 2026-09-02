### GTFS Tecator data analysis

######################################################
# packages
required_packages <- c("fda", "mvtnorm", "rainbow", "chemometrics", "fda.usc", "roahd", 
                       "pcaPP", "TeachingDemos", "plyr", "robustbase", "depthTools", 
                       "bootstrap", "cluster", "ks", "mrfDepth", "fdasrvf", "fdaoutlier", 
                       "FUNTA", "coga", "fsemipar")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

###########################################################################
# packages
library(fda)
library(mvtnorm)
library(rainbow)
library(chemometrics) # robust mahalanobis distance
#library(fda.usc) # outlier depth trim
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


###########################################################################
# Load data (R package "fsemipar")

#install.packages('fsemipar')
library(fsemipar)

# load data
data(Tecator)
names(Tecator)

###########################################################################
# Explorative Data Aanlysis

spectra0 <- Tecator$absor.spectra
dim(spectra0)  # 215 x 100 
# near-infrared absorbance spectra observed at 100 time points
# wavelengths in the range of 850-1050nm.

wavelength <- seq(850,1050,len=100) # wavelength sequence



###########################################################################
# Setting 

x <- spectra0
N <- nrow(x)
p <- ncol(x)
t <- seq(0,1,len=p)


###########################################################################
# plot 215 spectra
par(mfrow=c(1,1))
plot(wavelength, spectra0[1,], type='l', 
     ylim= c(min(spectra0)-0.1, max(spectra0)+0.1),
     main='Tecator', ylab="absorbance", xlab='wavelength')
for(i in 1:nrow(spectra0)){
  lines(wavelength, spectra0[i,], col=rgb(0.2,0.2,0.2,0.2))
}




###########################################################################
# Run analysis


set.seed(100)

method_out <- list()  # List containing the outlier-classified by competing methods


# GTFS - practical
trial1 <- GTFS_outlier(scheme='GTFS(P)',
                       x,t, 
                       basis='bspline', nbasis=15, nharm=15, 
                       v.prop=0.9,alpha=0.05, q.star='default',
                       n.init=20,iter.max=20,
                       refined = T)   
method_out[['GTFS(P)']] <- trial1$outlier



# GTFS - Algorithm 1
trial0 <- GTFS_outlier(scheme='GTFS',
                       x,t, 
                       basis='bspline', nbasis=15, nharm=15, 
                       v.prop=0.9,alpha=0.05, q.star='default',
                       n.init=20,iter.max=20,
                       refined = T) 

method_out[['GTFS']] <- trial0$outlier




# functional boxplot [Sun and Genton, 2011] (MBD) 
fbp <- fda::fbplot(t(x), method = 'MBD', plot=F)
method_out[['FBox']] <- fbp$outpoint 

# Integrated squared error  [Hyndman and Ullah 2007] (not depth)
Yfds <- fds(t, t(x), xname = "time", yname = "Simulated value")
FO.HU <- foutliers(Yfds, method='HUoutliers')

method_out[['ISE']] <- FO.HU$outliers


# Outliergram [Arribas-Gil & Romo, 2014] (MBD & MEI)
Yfdata <- fData(t,x)
og1 <- outliergram(Yfdata, display = F)

method_out[['Outliergram']] <- og1$ID_outliers


# Total variation depth, MSS [Huang & Sun ,2019]
tvd1 <- fdaoutlier::total_variation_depth(x)
c <- quantile(tvd1$mss,1/4) - 3*IQR(tvd1$mss)
method_out[['TVD (MSS)']] <- which(tvd1$mss <= c)




# Robust Mahalanobis distance [Hyndman & Shang, 2010]
rdist <- Moutlier(x,quantile=0.993, plot=F)
method_out[['RMD']] <- which(rdist$rd> rdist$cutoff)



# MUOD index [Ojo et al., 2021] (shape)
muod1 <- fdaoutlier::muod(x, cut_method = c('boxplot'))
if(is.vector(muod1)==F){
  out.temp <- NULL
}else if(is.atomic(muod1$outliers) ==T){
  out.temp <- NULL
}else{
  out.temp1 <- muod1$outliers$shape
  out.temp2 <- muod1$outliers$magnitude
}

method_out[['MUOD (shape)']] <- out.temp1
method_out[['MUOD (mag)']] <- out.temp2




# Least Trimmed Functional Set 
LTFS1 <- LTFS_outlier(x,N,p,15,'bspline',0.05,F)
method_out[['LTFS']] <- LTFS1$outlier

# Least Trimmed Functional Set [Refined]
method_out[['LTFS.re']] <- LTFS1$refined



###########################################################################
# Result summary 

# classified outlier index
method_out

# number of outliers from each model
sapply(method_out, length)

# Benchmark models
#benchmark_methods <- names(method_out)[-c(1,2)]
benchmark_methods <- c("FBox","ISE","Outliergram", "RMD",
                       "TVD (MSS)","MUOD (mag)","MUOD (shape)","LTFS" ,"LTFS.re"    )

# outliers labeled by benchmark models
benchmark_outliers <- sort(as.numeric(unique(unlist(method_out[benchmark_methods]))))
benchmark_outliers

# outlier detected only by GTFS
GTFS_only_outliers <-   method_out[['GTFS']] [  -which( method_out[['GTFS']] %in% benchmark_outliers  ) ]
GTFS_only_outliers


temp.num <- length(GTFS_only_outliers)

num_gtfs_only <- c(num_gtfs_only, temp.num)


###########################################################################
# Plot

library(tidyverse)

curve_wide <- as.data.frame(spectra0)
curve_wide$id <- seq_len(N)

curve_long <- curve_wide %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "grid",
    values_to = "value"
  ) %>%
  mutate(
    grid_index = as.integer(gsub("V", "", grid)),
    wavelength = wavelength[grid_index],
    GTFS_selected = id %in% method_out[['GTFS']],
    GTFS_only = id %in% GTFS_only_outliers
  )

curve_long


##########################
# Figure 5-1.

fig_51 <- ggplot() +
  geom_line(
    data = curve_long %>% dplyr::filter(!GTFS_selected),
    aes(x = wavelength, y = value, group = id),
    color = "grey75",
    alpha = 0.35,
    linewidth = 0.3,
    linetype=5
  ) +
  geom_line(
    data = curve_long %>% dplyr::filter(GTFS_selected),
    aes(x = wavelength, y = value, group = id),
    color = "skyblue",
    alpha = 0.85,
    linewidth = 0.7,
    linetype = 1
  ) + geom_line(
    data = mean_gtfs,
    aes(x = wavelength, y = mean_value),
    color = "darkblue",
    linewidth = 2.0,
    linetype=1
  )  +
  geom_line(
    data = mean_rest,
    aes(x = wavelength, y = mean_value),
    color = "black",
    linewidth = 2,
    linetype=5
  ) +
  theme_bw() + xlab('') + ylab('') 


print(fig_51)


##########################
# Figure 5-2.


fig_52 <- ggplot() +
  geom_line(
    data = curve_long %>% dplyr::filter(!GTFS_only),
    aes(x = wavelength, y = value, group = id),
    color = "grey80",
    alpha = 0.4,
    linewidth = 0.3
  ) +
  geom_line(
    data = curve_long %>% dplyr::filter(GTFS_only),
    aes(x = wavelength, y = value, group = id,
        col=as.factor(id)),
    #color = "red",
    alpha = 0.80,
    linewidth = 1.4,
    linetype=1
  )  +
  scale_colour_manual(values=c('red','blue','darkgreen'))+
  theme_bw() +
  xlab('') + ylab('') +
  theme(legend.position = 'none')

print(fig_52)



##########################
# Figure 5-3.


dat2 <- data.frame(id = 1:nrow(x), 
                   protein= Tecator$protein, 
                   fat = Tecator$fat, 
                   moisture = Tecator$moisture)
dat2 <- as.tibble(dat2)
dat2


score_vec <- primary_gtfs$score.H
if (is.null(score_vec)) score_vec <- rep(NA_real_, N)

cutoff_val <- primary_gtfs$cutoff
if (is.null(cutoff_val)) cutoff_val <- NA_real_

chem_tbl <- as.tibble(chem_tbl)

selected_tbl <- dat2 %>%
  mutate(
    GTFS_score= trial0$score.H,
    GTFS_selected= id %in% method_out[['GTFS']],
    detected_by_nonGTFS = id %in% benchmark_outliers,
    GTFS_only = GTFS_selected & !detected_by_nonGTFS
  ) %>%
  arrange(desc(GTFS_score)) %>%
  as.tibble()

head(selected_tbl)

score_df <- selected_tbl %>%
  mutate(
    label_group = case_when(
      GTFS_only ~ "GTFS-only",
      GTFS_selected ~ "GTFS-selected",
      TRUE ~ "Others"
    ),
    GTFS_only_id = ifelse(GTFS_only, paste0("ID ", id), "Others")
  )

head(score_df)

gtfs_only_ids <- paste0("ID ", GTFS_only_outliers)


col_vals <- c("Others" = "grey60",  setNames(c('red','blue','darkgreen'), gtfs_only_ids) )
shape_vals <- c("Others" = 16,      setNames(c(15,17,19),  gtfs_only_ids))
size_vals <- c("Others" = 1.2,      setNames(c(3.6,3.6,3.6),  gtfs_only_ids) )


fig_53 <- ggplot(score_df, aes(x = id, y = GTFS_score)) +
  geom_point(aes(color = GTFS_only_id, size=GTFS_only_id,
                 shape= GTFS_only_id),  alpha=0.55) +
  geom_hline(yintercept = cutoff_val, color = "black", 
             linewidth = 0.8, linetype=4) +
  scale_color_manual(values=col_vals) +
  scale_size_manual(values=size_vals) +
  scale_shape_manual(values=shape_vals) +
  theme_bw() +
  xlab('') + ylab('') +
  theme(
    legend.position = "none"
  )


print(fig_53)





##########################
# Figure 6. Violin plot 

chem_long2 <- selected_tbl %>%
  mutate(
    GTFS_group = ifelse(GTFS_selected, "GTFS-selected", "Not selected"),
    GTFS_only_id = ifelse(GTFS_only, paste0("ID ", id), "Others")
  ) %>%
  dplyr::select(id, GTFS_group, protein, fat, moisture, GTFS_only_id) %>%
  pivot_longer(
    cols = c(protein, fat, moisture),
    names_to = "Variable",
    values_to = "Value"
  )




plot_one_biochem <- function(var_name){
  
  ggplot(chem_long2 %>% dplyr::filter(Variable == var_name),
         aes(x = GTFS_group, y = Value)) +
    geom_violin(
      aes(fill = GTFS_group),
      alpha = 0.25,
      trim = FALSE,
      color = "grey40"
    ) +
    geom_boxplot(
      width = 0.13,
      outlier.shape = NA,
      alpha = 0.55
    ) +
    geom_jitter(
      aes(color = GTFS_only_id,
          shape = GTFS_only_id,
          size = GTFS_only_id,
          alpha = GTFS_only_id),
      width = 0.2
    ) +
    scale_color_manual(values = col_vals) +
    scale_shape_manual(values = shape_vals) +
    scale_size_manual(values = size_vals) +
    scale_fill_manual(values = c(
      "Not selected" = "grey85",
      "GTFS-selected" = "tomato"
    )) +
    scale_alpha_manual(values = c(rep(0.7, length(unique(chem_long2$GTFS_only_id)) - 1), 0.35)) +
    theme_bw() +
    labs(
      x = "",
      y = '', # varname
      fill = "",
      color = "",
      shape = "",
      size = "",
      alpha = ""
    ) +
    theme(
      legend.position = "none"
    )
}

fig_61  <- plot_one_biochem("protein")
fig_62  <- plot_one_biochem("fat")
fig_63  <- plot_one_biochem("moisture")

par(mfrow=c(1,3))
print(fig_61)
print(fig_62)
print(fig_63)
















