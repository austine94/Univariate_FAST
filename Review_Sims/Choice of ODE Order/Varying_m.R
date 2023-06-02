library(fda)
library(tidyverse)
library(mvtnorm)

result_df <- data.frame(Order1 = NA, Order2 = NA, Order3 = NA, Order4 = NA, Order5 = NA)

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

fast_detected <- rep(0, 5)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 8)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  #training data for each rep of the simulation
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  for(i in 1:5){
    fast_run <- fast(time, xfd, i, yfd$fd[j], threshold)
    if( !is.na(fast_run$detection_time[1])){
      fast_detected[i] <- fast_detected[i] + 1
    }
  }
  print(j)
}

new_result_df_row <- fast_detected/n_test
result_df <- rbind(result_df, new_result_df_row)

##################
#Model 3 - GP Residuals, Polynomial Anomalies
################

set.seed(40)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 5)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(4, 0.5, 1.5)
  b <- runif(1, 0, 0.5)
  y[101:200, i] <- y[101:200,i] + b[1] + a[1]*((time[101:200] - time[100]) / (lt)) +
    a[2]*((time[101:200] - time[100])^2/(lt)^2) + a[3]*((time[101:200] - time[100])^3/(lt)^3)
  + a[4]*((time[101:200] - time[100])^4/(lt)^4)
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 8)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  #training data for each rep of the simulation
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  for(i in 1:5){
    fast_run <- fast(time, xfd, i, yfd$fd[j], threshold)
    if( !is.na(fast_run$detection_time[1])){
      fast_detected[i] <- fast_detected[i] + 1
    }
  }
  print(j)
}

new_result_df_row <- fast_detected/n_test

result_df <- rbind(result_df, new_result_df_row)

##################
#Model 5 - GP Residuals, Sinusoidal Anomalies
################

set.seed(10)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 5)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.5)
  y[101:200, i] <- y[101:200,i] + a[1]*(sin(2*pi*time[101:200]/50) - sin(2*pi*time[100])) +
    a[2]*(cos(2*pi*time[101:200]/50) - cos(2*pi*time[100]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 8)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  #training data for each rep of the simulation
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  for(i in 1:5){
    fast_run <- fast(time, xfd, i, yfd$fd[j], threshold)
    if( !is.na(fast_run$detection_time[1])){
      fast_detected[i] <- fast_detected[i] + 1
    }
  }
  print(j)
}

new_result_df_row <- fast_detected/n_test

result_df <- rbind(result_df, new_result_df_row)

##################
#Model 7 - Change In Shape New Covariance
######################

set.seed(2000)

covar2 <- squared_exp_covar(101:200, 10, 0.3)  #for shape anomaly

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 5)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  y[101:200, i] <- y[101:200,i] - noise[101:200, i] + noise_anomaly[, i]
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 8)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  #training data for each rep of the simulation
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  for(i in 1:5){
    fast_run <- fast(time, xfd, i, yfd$fd[j], threshold)
    if( !is.na(fast_run$detection_time[1])){
      fast_detected[i] <- fast_detected[i] + 1
    }
  }
  print(j)
}

new_result_df_row <- fast_detected/n_test
result_df <- rbind(result_df, new_result_df_row)

####################
#Model 9 - T Residuals, Polynomial Anomalies
################

set.seed(40)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 5)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(4, 0.5, 1.5)
  b <- runif(1, 0, 0.5)
  y[101:200, i] <- y[101:200,i] + b[1] + a[1]*((time[101:200] - time[100]) / (lt)) +
    a[2]*((time[101:200] - time[100])^2/(lt)^2) + a[3]*((time[101:200] - time[100])^3/(lt)^3)
  + a[4]*((time[101:200] - time[100])^4/(lt)^4)
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 8)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

for(j in 1:ncol(y)){
  #training data for each rep of the simulation
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvt(n_train, df=2,sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  for(i in 1:5){
    fast_run <- fast(time, xfd, i, yfd$fd[j], threshold)
    if( !is.na(fast_run$detection_time[1])){
      fast_detected[i] <- fast_detected[i] + 1
    }
  }
  print(j)
}

new_result_df_row <- fast_detected/n_test

result_df <- rbind(result_df, new_result_df_row)


#####

result_df <- result_df[-1,]

