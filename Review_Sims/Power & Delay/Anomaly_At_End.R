library(fda)
library(tidyverse)
library(mvtnorm)

result_df <- data.frame(FAST_Detected = NA, FAST_Delay = NA, FAST_Delay_sd = NA)

##################
#Model 1 - GP Residuals, No Anomalies
################

set.seed(10)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  #training data for each rep of the simulation
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 0){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 0)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  print(j)
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected/n_test, fast_ADD, 
                       sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

################
#Model 2 - T Residuals, No Anomaly
#################

set.seed(10)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
y <- y + noise
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  #training data for each rep of the simulation
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvt(n_train, df=2, sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(fast_run$detection_time[1])){    #if detected
    
    fast_detected <- fast_detected + 1
    if(fast_run$detection_time > 0){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 0)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  print(j)
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected/n_test, fast_ADD, 
                       sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

################
#Model 3 - GP Residuals, Magnitude Anomaly
#################

set.seed(500)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- sample(c(-1,1), 3, TRUE) * runif(3, 0.01, 0.1)
  y[401:500, i] <- y[401:500,i] + a[3] + a[1]*((time[401:500]) / (lt)) +
    a[2]*exp(-time[401:500])
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
    if(fast_run$detection_time > 400){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 400)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change -1     }
  }
  print(j)
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected/n_test, fast_ADD, 
                       sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

################
#Model 4 - GP Residuals, Magnitude Anomaly, Random Start
#################

#Not tested as only focussed on anomalies at the end.

new_result_df_row <- rep(NA, 3)
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 5- Sinusoidal Anomaly
######################

set.seed(20)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.5)
  y[401:500, i] <- y[401:500,i] + a[1]*(sin(2*pi*time[401:500]/50) - sin(2*pi*time[400])) +
    a[2]*(cos(2*pi*time[401:500]/50) - cos(2*pi*time[400]))
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
    if(fast_run$detection_time > 400){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 400)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  print(c(j, fast_run$detection_time))
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected/n_test, fast_ADD, 
                       sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 6 - Abrupt Magnitude Anomaly
######################

set.seed(200)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  y[475:500, i] <- y[475:500,i] - x_underlying[475:500]
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
    if(fast_run$detection_time > 400){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 400)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  print(j)
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected/n_test, fast_ADD, 
                       sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 7 - Change In Shape New Covariance
######################

set.seed(2000)

covar2 <- squared_exp_covar(401:500, 15, 0.3)  #for shape anomaly

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
noise_anomaly <- rmvnorm(100, mean = rep(0, 100), sigma = covar2) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  y[401:500, i] <- y[401:500,i] - noise[401:500, i] + noise_anomaly[, i]
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
    if(fast_run$detection_time > 400){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 400)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change -1    }
  }
  print(j)
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected/n_test, fast_ADD, 
                       sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

#####################
#Model 8 - Translation in time
######################

set.seed(2000)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  translation <- sample(1:6, 1)
  y[401:500, i] <- y[401:500,i] - x_underlying[401:500] + x_underlying[(401 - translation):(500-translation)]
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
    if(fast_run$detection_time > 400){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 400)
    }else{  
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  print(j)
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected/n_test, fast_ADD, 
                       sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

############
#Model 9 - T Residuals, Sinusoidal Anomaly
############

set.seed(200)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- 100
fast_delay <- rep(NA, 100)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.5)
  y[401:500, i] <- y[401:500,i] + a[1]*(sin(2*pi*time[401:500]/50) - sin(2*pi*time[400])) +
    a[2]*(cos(2*pi*time[401:500]/50) - cos(2*pi*time[400]))
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
    if(fast_run$detection_time > 400){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 400)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1
    }
  }
  print(j)
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}
#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected/n_test, fast_ADD, 
                       sd(fast_delay, na.rm = TRUE))
result_df <- rbind(result_df, new_result_df_row)

#####################

result_df_anomaly_at_end <- result_df[-1,]
rownames(result_df_anomaly_at_end) <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6",
                         "Model 7", "Model 8", "Model 9")