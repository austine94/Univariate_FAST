library(fda)
library(tidyverse)
library(mvtnorm)
library(fdaoutlier)

result_df <- data.frame(FAST_Detected = NA, tvd_detected = NA, 
                        bagplot_detected = NA,
                        FAST_Delay = NA, FAST_Delay_sd = NA)

##################
#Model 1 - GP Residuals, No Anomalies
################

set.seed(3000)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 15, 0.01)
lt <- length(time)

mss_scale_factor <- 1.05 #as in Huang (2019) but reduced to allow pfa up to 0.05
tvd_scale_factor <- 0.5
bagplot_inflation <- 2.15 #control pfa at 0.05

n_train <- 100
n_test <- 500

x_underlying <- exp(-time/lt) + exp(-time^2/lt^2) -  sinh(-0.5+time/500) - cosh(-0.5+time/500)

fast_detected <- 0
fast_detected_after_change <- 500
fast_delay <- rep(NA, 500)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  #training data for each rep of the simulation
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  fast_run <- fast_end_of_day(time, xfd, 4, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies as in Huang (2019)
  if(length(tvd_run$outliers) > 0){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + length(tvd_run$outliers)
  }
  
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(!is.na(bag_run[1])){
    bagplot_detected <- bagplot_detected + length(bag_run)
  }
  
  
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected / n_test, tvd_detected / (n_test*(n_train+1)), 
                       bagplot_detected / (n_test*(n_train+1)),
                       fast_ADD, sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

################
#Model 2 - T Residuals, No Anomaly
#################

set.seed(10)

x_underlying <- exp(-time/lt) + exp(-time^2/lt^2) -  sinh(-0.5+time/500) - cosh(-0.5+time/500)

fast_detected <- 0
fast_detected_after_change <- 500
fast_delay <- rep(NA, 500)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
y <- y + noise
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  #training data for each rep of the simulation
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvt(n_train, df=2, sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  fast_run <- fast_end_of_day(time, xfd, 4, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  if(length(tvd_run$outliers) > 0){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(!is.na(bag_run[1])){
    bagplot_detected <- bagplot_detected + 1
  }
  
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected / n_test, tvd_detected / (n_test*(n_train+1)), 
                       bagplot_detected / (n_test*(n_train+1)),
                       fast_ADD, sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

################
#Model 3 - GP Residuals, Magnitude Anomaly
#################

set.seed(200)

x_underlying <- exp(-time/lt) + exp(-time^2/lt^2) -  sinh(-0.5+time/500) - cosh(-0.5+time/500)

fast_detected <- 0
fast_detected_after_change <- 500
fast_delay <- rep(NA, 500)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- sample(c(-1,1), 3, TRUE) * runif(3, 0.01, 0.1)
  y[100:200, i] <- y[100:200,i] + a[1] + a[2]*((time[100:200]-time[100])/lt)
  a[3]*exp(-(time[100:200] - time[100]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  fast_run <- fast_end_of_day(time, xfd, 4, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test,
                       fast_ADD, sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

################
#Model 4 - GP Residuals, Polynomial Anomaly, Random Start
#################

set.seed(30)

x_underlying <- exp(-time/lt) + exp(-time^2/lt^2) -  sinh(-0.5+time/500) - cosh(-0.5+time/500)

fast_detected <- 0
fast_detected_after_change <- 500
fast_delay <- rep(NA, 500)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
start_times <- sample(101:491, n_test, TRUE)
end_times <- sapply(start_times, function(x) min((x+99), lt))
for(i in 1:ncol(y)){
  a <- sample(c(-1,1), 3, TRUE) * runif(3, 0.01, 0.1)
  y[(start_times[i]:end_times[i]), i] <- y[(start_times[i]:end_times[i]),i] + a[3] +
    a[1]*((start_times[i]:end_times[i]) / (lt)) +
    a[2]*exp(-(start_times[i]:end_times[i]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  fast_run <- fast_end_of_day(time, xfd, 4, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > start_times[i]){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - (start_times[i] - 1))
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test,
                       fast_ADD, sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 5- Sinusoidal Anomaly
######################

set.seed(1213)

x_underlying <- exp(-time/lt) + exp(-time^2/lt^2) -  sinh(-0.5+time/500) - cosh(-0.5+time/500)

fast_detected <- 0
fast_detected_after_change <- 500
fast_delay <- rep(NA, 500)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.3)
  y[101:200, i] <- y[101:200,i] + a[1]*(sin(2*pi*time[101:200]/50 - sin(2*pi*time[100]))) +
    a[2]*(cos(2*pi*time[101:200]/50) - cos(2*pi*time[100]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  fast_run <- fast_end_of_day(time, xfd, 4, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test,
                       fast_ADD, sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 6 - Abrupt Magnitude Anomaly
######################

set.seed(200)

x_underlying <- exp(-time/lt) + exp(-time^2/lt^2) -  sinh(-0.5+time/500) - cosh(-0.5+time/500)

fast_detected <- 0
fast_detected_after_change <- 500
fast_delay <- rep(NA, 500)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  y[175:200, i] <- y[175:200,i] - x_underlying[175:200]
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  fast_run <- fast_end_of_day(time, xfd, 4, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test,
                       fast_ADD, sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 7 - Change In Shape New Covariance
######################

set.seed(2000)

covar2 <- squared_exp_covar(101:200, 35, 0.3)  #for shape anomaly

x_underlying <- exp(-time/lt) + exp(-time^2/lt^2) -  sinh(-0.5+time/500) - cosh(-0.5+time/500)

fast_detected <- 0
fast_detected_after_change <- 500
fast_delay <- rep(NA, 500)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
noise_anomaly <- rmvnorm(n_test, mean = rep(0, 100), sigma = covar2) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  y[101:200, i] <- y[101:200,i] - noise[101:200, i] + noise_anomaly[, i]
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  fast_run <- fast_end_of_day(time, xfd, 4, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test,
                       fast_ADD, sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 8 - Translation in time
######################

set.seed(8000)

x_underlying <- exp(-time/lt) + exp(-time^2/lt^2) -  sinh(-0.5+time/500) - cosh(-0.5+time/500)

fast_detected <- 0
fast_detected_after_change <- 500
fast_delay <- rep(NA, 500)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  translation <- sample( 1:10, 1)
  y[101:200, i] <- y[101:200,i] - x_underlying[101:200] + x_underlying[(101 - translation):(200-translation)]
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  fast_run <- fast_end_of_day(time, xfd, 4, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{  
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test,
                       fast_ADD, sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

############
#Model 9 - T Residuals, Sinusoidal Anomaly
############

set.seed(333)

x_underlying <- exp(-time/lt) + exp(-time^2/lt^2) -  sinh(-0.5+time/500) - cosh(-0.5+time/500)

fast_detected <- 0
fast_detected_after_change <- 500
fast_delay <- rep(NA, 500)

tvd_detected <- 0
bagplot_detected <- 0

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.3)
  y[101:200, i] <- y[101:200,i] + a[1]*(sin(2*pi*time[101:200]/50 - sin(2*pi*time[100]))) +
    a[2]*(cos(2*pi*time[101:200]/50) - cos(2*pi*time[100]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)
y_eval <- eval.fd(time, yfd$fd)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvt(n_train, df=2, sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_eval <- eval.fd(time, xfd$fd)
  
  fast_run <- fast_end_of_day(time, xfd, 4, yfd$fd[j], threshold = threshold)
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1
    }
  }
  
  combined_fd <- cbind(x_eval, y_eval[,j])  #add anomaly to null obs
  tvd_run <- tvdmss(t(combined_fd), mss_scale_factor, tvd_scale_factor, 0.5) 
  #we use both metrics to detect anomalies
  if(101 %in% tvd_run$outliers){ #if at least one anomaly detected
    tvd_detected <- tvd_detected + 1
  }
  
  combined_fd_smooth <- fd(cbind(xfd$fd$coefs, yfd$fd$coefs[,j]), bas)
  bag_run <- fbagplot(combined_fd_smooth, bagplot_inflation)
  if(101 %in% bag_run){
    bagplot_detected <- bagplot_detected + 1
  }
  
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}
#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected / n_test, tvd_detected / n_test, 
                       bagplot_detected / n_test,
                       fast_ADD, sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

#####################

result_df_order4_comparison <- result_df[-1,c(-4,-5)]
rownames(result_df_order4_comparison) <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6",
                                           "Model 7", "Model 8", "Model 9")
print("Table 15: Done")