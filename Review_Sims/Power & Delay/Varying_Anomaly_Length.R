library(fda)
library(tidyverse)
library(mvtnorm)

s_len_vec <- seq(5, 100, by = 5)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

################
#Model 3 - GP Residuals, Magnitude Anomaly
#################

set.seed(10)

fast_detected_by_s <- rep(NA, length(s_len_vec))
fast_delay_by_s <- fast_delay_sd_by_s <- rep(NA, length(s_len_vec))

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

for(k in 1:length(s_len_vec)){ #iterate over each anomaly region length
  fast_detected <- 0
  fast_detected_after_change <- 0
  fast_delay <- rep(NA, 100)
  
  end_period <- 200
  start_period <- end_period - s_len_vec[k] + 1 #place anomalies at end of time
  
  #generate the 100  test observations
  y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
  noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
  y <- y + noise
  for(i in 1:ncol(y)){
    a <- runif(3, 0.01, 0.1)
    y[(start_period:end_period), i] <- y[(start_period:end_period),i] + a[3] + a[1]*((time[(start_period:end_period)] - start_period) / (lt)) +
      a[2]*exp(-time[(start_period:end_period)])
  }
  bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
  yfd <- smooth.basis(time, y, bas)
  
  threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )
  
  for(j in 1:ncol(y)){
    
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      
      fast_detected <- fast_detected + 1
      if(fast_run$detection_time > start_period){  #ADD as in Tartakovsky
        fast_delay[j] <- (fast_run$detection_time - start_period)
      }else{
        fast_delay[j] <- 0 
        fast_detected_after_change <- fast_detected_after_change + 1    }
    }
    print(c(k,j))
  }
  if(all(is.na(fast_delay))){
    fast_delay <- rep(Inf, n_test)
  }
  
  #compute ADD as in Tartakovsky (2014)
  fast_detected_by_s[k] <- fast_detected/n_test
  fast_delay_by_s[k] <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 
  fast_delay_sd_by_s[k] <- sd(fast_delay, na.rm = TRUE)
  
}
model3_detected <- fast_detected_by_s
model3_ADD <- fast_delay_by_s 
model3_ADD_sd <- fast_delay_sd_by_s


################
#Model 5 - GP Residuals, Sinusoidal Anomaly
#################

set.seed(400)

fast_detected_by_s <- rep(NA, length(s_len_vec))
fast_delay_by_s <- fast_delay_sd_by_s <- rep(NA, length(s_len_vec))

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

for(k in 1:length(s_len_vec)){ #iterate over each anomaly region length
  fast_detected <- 0
  fast_detected_after_change <- 0
  fast_delay <- rep(NA, 100)
  
  end_period <- 200
  start_period <- end_period - s_len_vec[k] + 1 #place anomalies at end of time
  
  #generate the 100  test observations
  y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
  noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
  y <- y + noise
  for(i in 1:ncol(y)){
    a <- runif(2, 0.1, 0.5)
    y[(start_period:end_period), i] <- y[(start_period:end_period),i] + a[1]*(sin(2*pi*time[(start_period:end_period)]/50) - sin(2*pi*start_period)) +
      a[2]*(cos(2*pi*time[(start_period:end_period)]/50) - cos(2*pi*start_period))
  }
  bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
  yfd <- smooth.basis(time, y, bas)
  
  threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )
  
  for(j in 1:ncol(y)){
    
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      
      fast_detected <- fast_detected + 1
      if(fast_run$detection_time > start_period){  #ADD as in Tartakovsky
        fast_delay[j] <- (fast_run$detection_time - start_period)
      }else{
        fast_delay[j] <- 0 
        fast_detected_after_change <- fast_detected_after_change + 1    }
    }
    print(c(k,j))
  }
  if(all(is.na(fast_delay))){
    fast_delay <- rep(Inf, n_test)
  }
  
  #compute ADD as in Tartakovsky (2014)
  fast_detected_by_s[k] <- fast_detected/n_test
  fast_delay_by_s[k] <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 
  fast_delay_sd_by_s[k] <- sd(fast_delay, na.rm = TRUE)
  
}
model5_detected <- fast_detected_by_s
model5_ADD <- fast_delay_by_s 
model5_ADD_sd <- fast_delay_sd_by_s

################
#Model 7 - GP Residuals, Covariance Change in Shape
#################

set.seed(400)

fast_detected_by_s <- rep(NA, length(s_len_vec))
fast_delay_by_s <- fast_delay_sd_by_s <- rep(NA, length(s_len_vec))

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

for(k in 1:length(s_len_vec)){ #iterate over each anomaly region length
  fast_detected <- 0
  fast_detected_after_change <- 0
  fast_delay <- rep(NA, 100)
  
  end_period <- 200
  start_period <- end_period - s_len_vec[k] + 1 #place anomalies at end of time
  covar2 <- squared_exp_covar((start_period:end_period), 15, 0.3)  #for shape anomaly
  
  #generate the 100  test observations
  y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
  noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
  y <- y + noise
  noise_anomaly <- rmvnorm(100, mean = rep(0, s_len_vec[k]), sigma = covar2) %>% t()
  for(i in 1:ncol(y)){
    y[(start_period:end_period), i] <- y[(start_period:end_period),i] - noise[(start_period:end_period), i] + noise_anomaly[, i]
  }
  bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
  yfd <- smooth.basis(time, y, bas)
  
  threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )
  
  for(j in 1:ncol(y)){
    
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      
      fast_detected <- fast_detected + 1
      if(fast_run$detection_time > start_period){  #ADD as in Tartakovsky
        fast_delay[j] <- (fast_run$detection_time - start_period)
      }else{
        fast_delay[j] <- 0 
        fast_detected_after_change <- fast_detected_after_change + 1    }
    }
    print(c(k,j))
  }
  if(all(is.na(fast_delay))){
    fast_delay <- rep(Inf, n_test)
  }
  
  #compute ADD as in Tartakovsky (2014)
  fast_detected_by_s[k] <- fast_detected/n_test
  fast_delay_by_s[k] <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 
  fast_delay_sd_by_s[k] <- sd(fast_delay, na.rm = TRUE)
  
}
model7_detected <- fast_detected_by_s
model7_ADD <- fast_delay_by_s 
model7_ADD_sd <- fast_delay_sd_by_s

################
#Model 9 - T Residuals, Sinusoidal Anomaly
#################

set.seed(10000)

fast_detected_by_s <- rep(NA, length(s_len_vec))
fast_delay_by_s <- fast_delay_sd_by_s <- rep(NA, length(s_len_vec))

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

for(k in 1:length(s_len_vec)){ #iterate over each anomaly region length
  fast_detected <- 0
  fast_detected_after_change <- 0
  fast_delay <- rep(NA, 100)
  
  end_period <- 200
  start_period <- end_period - s_len_vec[k] + 1 #place anomalies at end of time
  
  #generate the 100  test observations
  y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
  noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
  y <- y + noise
  for(i in 1:ncol(y)){
    a <- runif(2, 0.1, 0.5)
    y[(start_period:end_period), i] <- y[(start_period:end_period),i] + a[1]*(sin(2*pi*time[(start_period:end_period)]/50) - sin(2*pi*start_period)) +
      a[2]*(cos(2*pi*time[(start_period:end_period)]/50) - cos(2*pi*start_period))
  }
  bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
  yfd <- smooth.basis(time, y, bas)
  
  threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )
  
  for(j in 1:ncol(y)){
    
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    noise <- rmvt(n_train, df=2, sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      
      fast_detected <- fast_detected + 1
      if(fast_run$detection_time > start_period){  #ADD as in Tartakovsky
        fast_delay[j] <- (fast_run$detection_time - start_period)
      }else{
        fast_delay[j] <- 0 
        fast_detected_after_change <- fast_detected_after_change + 1    }
    }
    print(c(k,j))
  }
  if(all(is.na(fast_delay))){
    fast_delay <- rep(Inf, n_test)
  }
  
  #compute ADD as in Tartakovsky (2014)
  fast_detected_by_s[k] <- fast_detected/n_test
  fast_delay_by_s[k] <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 
  fast_delay_sd_by_s[k] <- sd(fast_delay, na.rm = TRUE)
  
}
model9_detected <- fast_detected_by_s
model9_ADD <- fast_delay_by_s 
model9_ADD_sd <- fast_delay_sd_by_s

######################

varying_s_df <- data.frame(model3 = model3_detected, model5 = model5_detected, 
           model7 = model7_detected, model9 = model9_detected)
row.names(varying_s_df) <- as.character(s_len_vec)

