library(fda)
library(tidyverse)
library(mvtnorm)

result_df <- data.frame("1%" = NA, "5%"= NA, "10%" = NA)

##################
#Model 3 - GP Residuals, Magnitude Anomalies
################

set.seed(10)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 500
n_contaminated <- c(1, 5, 10) #1, 5, and 10% percent contamination

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 3)

for(k in 1:length(n_contaminated)){
  #generate the 100  test observations
  #also store the contaminated observations to add to the training data at each iteration
  y <- matrix(rep(x_underlying, (n_test + n_contaminated[k]*n_test)), nrow = length(time), 
              ncol = (n_test + n_contaminated[k]*n_test), byrow = FALSE)
  noise <- rmvnorm((n_test + n_contaminated[k]*n_test), mean = rep(0, length(time)), sigma = covar) %>% t()
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
  
  for(j in 1:n_test){
    #training data for each rep of the simulation
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    x[,((n_train - n_contaminated[k] + 1):n_train)] <- {
        y_eval[,((n_test + ((j-1)*n_contaminated[k] + 1)):(n_test + j*n_contaminated[k]))]
      }
    #above line adds contaminated functions to the training data
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      fast_detected[k] <- fast_detected[k] + 1
    }
     
  }
  
}
result_df <- rbind(result_df, fast_detected)

##################
#Model 4 - GP Residuals, Varying Start Magnitude Anomalies
################

set.seed(150)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 500
n_contaminated <- c(1, 5, 10) #1, 5, and 10% percent contamination

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 3)

for(k in 1:length(n_contaminated)){
  #generate the 100  test observations
  #also store the contaminated observations to add to the training data at each iteration
  y <- matrix(rep(x_underlying, (n_test + n_contaminated[k]*n_test)), nrow = length(time), 
              ncol = (n_test + n_contaminated[k]*n_test), byrow = FALSE)
  noise <- rmvnorm((n_test + n_contaminated[k]*n_test), mean = rep(0, length(time)), sigma = covar) %>% t()
  y <- y + noise
  start_times <- sample(101:491, (n_test + n_contaminated[k]*n_test), TRUE)
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
  
  for(j in 1:n_test){
    #training data for each rep of the simulation
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    x[,((n_train - n_contaminated[k] + 1):n_train)] <- {
      y_eval[,((n_test + ((j-1)*n_contaminated[k] + 1)):(n_test + j*n_contaminated[k]))]
    }
    #above line adds contaminated functions to the training data
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      fast_detected[k] <- fast_detected[k] + 1
    }
     
  }
  
}
result_df <- rbind(result_df, fast_detected)

##################
#Model 5 - GP Residuals, Sinusoidal Anomalies
################

set.seed(190)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 500
n_contaminated <- c(1, 5, 10) #1, 5, and 10% percent contamination

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 3)

for(k in 1:length(n_contaminated)){
  #generate the 100  test observations
  #also store the contaminated observations to add to the training data at each iteration
  y <- matrix(rep(x_underlying, (n_test + n_contaminated[k]*n_test)), nrow = length(time), 
              ncol = (n_test + n_contaminated[k]*n_test), byrow = FALSE)
  noise <- rmvnorm((n_test + n_contaminated[k]*n_test), mean = rep(0, length(time)), sigma = covar) %>% t()
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
  
  for(j in 1:n_test){
    #training data for each rep of the simulation
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    x[,((n_train - n_contaminated[k] + 1):n_train)] <- {
      y_eval[,((n_test + ((j-1)*n_contaminated[k] + 1)):(n_test + j*n_contaminated[k]))]
    }
    #above line adds contaminated functions to the training data
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      fast_detected[k] <- fast_detected[k] + 1
    }
     
  }
  
}
result_df <- rbind(result_df, fast_detected)

##################
#Model 6 - Abrupt Magnitude Anomalies
################

set.seed(30)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 500
n_contaminated <- c(1, 5, 10) #1, 5, and 10% percent contamination

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 3)

for(k in 1:length(n_contaminated)){
  #generate the 100  test observations
  #also store the contaminated observations to add to the training data at each iteration
  y <- matrix(rep(x_underlying, (n_test + n_contaminated[k]*n_test)), nrow = length(time), 
              ncol = (n_test + n_contaminated[k]*n_test), byrow = FALSE)
  noise <- rmvnorm((n_test + n_contaminated[k]*n_test), mean = rep(0, length(time)), sigma = covar) %>% t()
  y <- y + noise
  for(i in 1:ncol(y)){
    y[175:200, i] <- y[175:200,i] - x_underlying[175:200]
  }
  bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
  yfd <- smooth.basis(time, y, bas)
  y_eval <- eval.fd(time, yfd$fd)
  
  threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )
  
  for(j in 1:n_test){
    #training data for each rep of the simulation
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    x[,((n_train - n_contaminated[k] + 1):n_train)] <- {
      y_eval[,((n_test + ((j-1)*n_contaminated[k] + 1)):(n_test + j*n_contaminated[k]))]
    }
    #above line adds contaminated functions to the training data
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      fast_detected[k] <- fast_detected[k] + 1
    }
     
  }
  
}
result_df <- rbind(result_df, fast_detected)

##################
#Model 7 - Change in shape covariance
################

set.seed(40)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 500
n_contaminated <- c(1, 5, 10) #1, 5, and 10% percent contamination

covar2 <- squared_exp_covar(101:200, 35, 0.3)  #for shape anomaly

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 3)

for(k in 1:length(n_contaminated)){
  #generate the 100  test observations
  #also store the contaminated observations to add to the training data at each iteration
  y <- matrix(rep(x_underlying, (n_test + n_contaminated[k]*n_test)), nrow = length(time), 
              ncol = (n_test + n_contaminated[k]*n_test), byrow = FALSE)
  noise <- rmvnorm((n_test + n_contaminated[k]*n_test), mean = rep(0, length(time)), sigma = covar) %>% t()
  noise_anomaly <- rmvnorm((n_test + n_contaminated[k]*n_test), mean = rep(0, 100), sigma = covar2) %>% t()
  y <- y + noise
  for(i in 1:ncol(y)){
    y[101:200, i] <- y[101:200,i] - noise[101:200, i] + noise_anomaly[, i]
  }
  bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
  yfd <- smooth.basis(time, y, bas)
  y_eval <- eval.fd(time, yfd$fd)
  
  threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )
  
  for(j in 1:n_test){
    #training data for each rep of the simulation
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    x[,((n_train - n_contaminated[k] + 1):n_train)] <- {
      y_eval[,((n_test + ((j-1)*n_contaminated[k] + 1)):(n_test + j*n_contaminated[k]))]
    }
    #above line adds contaminated functions to the training data
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      fast_detected[k] <- fast_detected[k] + 1
    }
     
  }
  
}
result_df <- rbind(result_df, fast_detected)

##################
#Model 8 - Translation
################

set.seed(20000)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 500
n_contaminated <- c(1, 5, 10) #1, 5, and 10% percent contamination

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 3)

for(k in 1:length(n_contaminated)){
  #generate the 100  test observations
  #also store the contaminated observations to add to the training data at each iteration
  y <- matrix(rep(x_underlying, (n_test + n_contaminated[k]*n_test)), nrow = length(time), 
              ncol = (n_test + n_contaminated[k]*n_test), byrow = FALSE)
  noise <- rmvnorm((n_test + n_contaminated[k]*n_test), mean = rep(0, length(time)), sigma = covar) %>% t()
  y <- y + noise
  for(i in 1:ncol(y)){
    translation <- sample( 1:10, 1)
    y[101:200, i] <- y[101:200,i] - x_underlying[101:200] + x_underlying[(101 - translation):(200-translation)]
  }
  bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
  yfd <- smooth.basis(time, y, bas)
  y_eval <- eval.fd(time, yfd$fd)
  
  threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )
  
  for(j in 1:n_test){
    #training data for each rep of the simulation
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    x[,((n_train - n_contaminated[k] + 1):n_train)] <- {
      y_eval[,((n_test + ((j-1)*n_contaminated[k] + 1)):(n_test + j*n_contaminated[k]))]
    }
    #above line adds contaminated functions to the training data
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      fast_detected[k] <- fast_detected[k] + 1
    }
     
  }
  
}
result_df <- rbind(result_df, fast_detected)

##################
#Model 9 - Heavy Tailed Residuals, Sinusoidal
################

set.seed(123)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 500
n_contaminated <- c(1, 5, 10) #1, 5, and 10% percent contamination

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 3)

for(k in 1:length(n_contaminated)){
  #generate the 100  test observations
  #also store the contaminated observations to add to the training data at each iteration
  y <- matrix(rep(x_underlying, (n_test + n_contaminated[k]*n_test)), nrow = length(time), 
              ncol = (n_test + n_contaminated[k]*n_test), byrow = FALSE)
  noise <- rmvt((n_test + n_contaminated[k]*n_test), df = 2,sigma = covar) %>% t()
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
  
  for(j in 1:n_test){
    #training data for each rep of the simulation
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    x[,((n_train - n_contaminated[k] + 1):n_train)] <- {
      y_eval[,((n_test + ((j-1)*n_contaminated[k] + 1)):(n_test + j*n_contaminated[k]))]
    }
    #above line adds contaminated functions to the training data
    noise <- rmvt(n_train, df=2, sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      fast_detected[k] <- fast_detected[k] + 1
    }
     
  }
  
}
result_df <- rbind(result_df, fast_detected)

#########

result_df_contaminated <- result_df[-1,] / n_test
rownames(result_df_contaminated) <- c("Model 3", "Model 4", "Model 5", "Model 6",
                                      "Model 7", "Model 8", "Model 9")
colnames(result_df_contaminated) <- c("1%", "5%", "10%")
print("Contaminated Sim: Done")
