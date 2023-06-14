library(fda)
library(tidyverse)
library(mvtnorm)
library(fdaoutlier)

result_df <- data.frame(FAST_Detected = NA, tvd_detected = NA, 
                        bagplot_detected = NA)


#False alarms already set in offline partial period study

################
#Model 3 - GP Residuals, Magnitude Anomaly
#################

set.seed(100)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- sample(c(-1,1), 3, TRUE) * runif(3, 0.01, 0.1)
  y[1:500, i] <- y[1:500,i] + a[1] + a[2]*((time[1:500])/lt)
  a[3]*exp(-(time[1:500]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- 1.9

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  #use anomaly and null together for offline detection
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  
  fast_run <- fast_end_of_day(time, xfd, 2, combined_fd_smooth[(n_test+1)],
                              threshold = threshold)  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1

  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  print(j)
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}


new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test)
result_df <- rbind(result_df, new_result_df_row)

################
#Model 4 - GP Residuals, Polynomial Anomaly, Random Start
#################
#not tested as fixed anomaly length
new_result_df_row <- c(NA, NA, NA) 
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 5- Sinusoidal Anomaly
######################

set.seed(200)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.5)
  y[1:500, i] <- y[1:500,i] + a[1]*(sin(2*pi*time[1:500]/50)) +
    a[2]*(cos(2*pi*time[1:500]/50))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- 1.9

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  #use anomaly and null together for offline detection
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  
  fast_run <- fast_end_of_day(time, xfd, 2, combined_fd_smooth[(n_test+1)],
                              threshold = threshold)  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  print(j)
}

#compute ADD as in Tartakovsky (2014)

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test)
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 6 - Abrupt Magnitude Anomaly
######################

set.seed(200)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  y[1:500, i] <- y[1:500,i] - x_underlying[1:500]
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- 1.9

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  #use anomaly and null together for offline detection
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  
  fast_run <- fast_end_of_day(time, xfd, 2, combined_fd_smooth[(n_test+1)],
                              threshold = threshold)  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  print(j)
}

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test)
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 7 - Change In Shape New Covariance
######################

set.seed(2000)

covar2 <- squared_exp_covar(1:500, 15, 0.3)  #for shape anomaly

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
noise_anomaly <- rmvnorm(100, mean = rep(0, 500), sigma = covar2) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  y[1:500, i] <- y[1:500,i] - noise[1:500, i] + noise_anomaly[, i]
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- 1.9

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  #use anomaly and null together for offline detection
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  
  fast_run <- fast_end_of_day(time, xfd, 2, combined_fd_smooth[(n_test+1)],
                              threshold = threshold)  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  print(j)
}

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test)
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 8 - Translation in time
######################

set.seed(2000)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  translation <- sample( 1:10, 1)
  y[1:500, i] <- y[1:500,i] - x_underlying[1:500] + x_underlying[(101 - translation):(200-translation)]
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- 1.9

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  #use anomaly and null together for offline detection
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  
  fast_run <- fast_end_of_day(time, xfd, 2, combined_fd_smooth[(n_test+1)],
                              threshold = threshold)  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
  }

  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  print(j)
}

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test)
result_df <- rbind(result_df, new_result_df_row)

############
#Model 9 - T Residuals, Sinusoidal Anomaly
############

set.seed(200)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.5)
  y[1:500, i] <- y[1:500,i] + a[1]*(sin(2*pi*time[1:500]/50)) +
    a[2]*(cos(2*pi*time[1:500]/50))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- 1.9

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvt(n_train, df=2, sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  #use anomaly and null together for offline detection
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  
  fast_run <- fast_end_of_day(time, xfd, 2, combined_fd_smooth[(n_test+1)],
                              threshold = threshold)
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
  }
  
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  print(j)
}

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test)
result_df <- rbind(result_df, new_result_df_row)

#####################

result_df_offline_full <- result_df[-1,]
#add false alarms from other offline study
result_df_offline_full <- rbind(result_df_offline_partial[1:2,], result_df_offline_full)
rownames(result_df_offline_full) <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6",
                                         "Model 7", "Model 8", "Model 9")
