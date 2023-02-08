library(fda)
library(tidyverse)
library(mvtnorm)

######################
#GP Residuals
######################

#######################

#Polynomial Anomalies

########################

set.seed(10)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

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
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

detected_one <- detected_two <- 0
delay_one <- delay_two <- rep(NA, ncol(y))

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  test_one <- fast(time, xfd, 1, yfd$fd[j], threshold = threshold)
  test_two <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(test_one$detection_time[1])){    #if detected
    
    detected_one <- detected_one + 1
    if(test_one$detection_time > 100){  #CADD
      delay_one[j] <- (test_one$detection_time - 100)
    }else{
      delay_one[j] <- 1
    }
  }
  if( !is.na(test_two$detection_time[1])){
    
    detected_two <- detected_two + 1
    if(test_two$detection_time > 100){ 
      delay_two[j] <- (test_two$detection_time - 100)
    }else{
      delay_two[j] <- 1
    }
  }
  print(j)
}

power_one <- detected_one / ncol(y)
power_two <- detected_two / ncol(y)

cadd_one <- mean(delay_one, na.rm = TRUE)
cadd_two <- mean(delay_two, na.rm = TRUE)

sd_cadd_one <- sd(delay_one, na.rm = TRUE)
sd_cadd_two <- sd(delay_two, na.rm = TRUE)

gp_polynomial_df <- data.frame(power_one = power_one, power_two = power_two,
                            cadd_one = cadd_one, cadd_two = cadd_two,
                            sd_cadd_one = sd_cadd_one, sd_cadd_two = sd_cadd_two)

##########################################

#sinusoidal anomaly

##########################################

set.seed(100)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.5)
  y[101:200, i] <- y[101:200,i] + + a[1]*(sin(2*pi*time[101:200]/50) - sin(2*pi*time[100])) +
    a[2]*(cos(2*pi*time[101:200]/50) - cos(2*pi*time[100]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

detected_one <- detected_two <- 0
delay_one <- delay_two <- rep(NA, ncol(y))

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  test_one <- fast(time, xfd, 1, yfd$fd[j], threshold = threshold)
  test_two <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(test_one$detection_time[1])){    #if detected
    
    detected_one <- detected_one + 1
    if(test_one$detection_time > 100){  #CADD
      delay_one[j] <- (test_one$detection_time - 100)
    }else{
      delay_one[j] <- 1
    }
  }
  if( !is.na(test_two$detection_time[1])){
    
    detected_two <- detected_two + 1
    if(test_two$detection_time > 100){ 
      delay_two[j] <- (test_two$detection_time - 100)
    }else{
      delay_two[j] <- 1
    }
  }
  print(j)
}

power_one <- detected_one / ncol(y)
power_two <- detected_two / ncol(y)

cadd_one <- mean(delay_one, na.rm = TRUE)
cadd_two <- mean(delay_two, na.rm = TRUE)

sd_cadd_one <- sd(delay_one, na.rm = TRUE)
sd_cadd_two <- sd(delay_two, na.rm = TRUE)

gp_sinusoidal_df <- data.frame(power_one = power_one, power_two = power_two,
                            cadd_one = cadd_one, cadd_two = cadd_two,
                            sd_cadd_one = sd_cadd_one, sd_cadd_two = sd_cadd_two)

###########################################

#Noise Anomaly

###########################################


set.seed(500)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  #start at 175 so the process smoothly drops to zero first
  y[175:200, i] <- y[175:200,i] - x_underlying[175:200] 
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

detected_one <- detected_two <- 0
delay_one <- delay_two <- rep(NA, ncol(y))

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  test_one <- fast(time, xfd, 1, yfd$fd[j], threshold = threshold)
  test_two <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(test_one$detection_time[1])){    #if detected
    
    detected_one <- detected_one + 1
    if(test_one$detection_time > 100){  #CADD
      delay_one[j] <- (test_one$detection_time - 100)
    }else{
      delay_one[j] <- 1
    }
  }
  if( !is.na(test_two$detection_time[1])){
    
    detected_two <- detected_two + 1
    if(test_two$detection_time > 100){ 
      delay_two[j] <- (test_two$detection_time - 100)
    }else{
      delay_two[j] <- 1
    }
  }
  print(j)
}

power_one <- detected_one / ncol(y)
power_two <- detected_two / ncol(y)

cadd_one <- mean(delay_one, na.rm = TRUE)
cadd_two <- mean(delay_two, na.rm = TRUE)

sd_cadd_one <- sd(delay_one, na.rm = TRUE)
sd_cadd_two <- sd(delay_two, na.rm = TRUE)

gp_noise_df <- data.frame(power_one = power_one, power_two = power_two,
                       cadd_one = cadd_one, cadd_two = cadd_two,
                       sd_cadd_one = sd_cadd_one, sd_cadd_two = sd_cadd_two)
###################################

###################################

###################
#T Process Residuals
###################


#######################

#Polynomial Anomalies

########################

set.seed(10)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

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
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

detected_one <- detected_two <- 0
delay_one <- delay_two <- rep(NA, ncol(y))

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvt(n_train, df=  2, sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  test_one <- fast(time, xfd, 1, yfd$fd[j], threshold = threshold)
  test_two <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(test_one$detection_time[1])){    #if detected
    
    detected_one <- detected_one + 1
    if(test_one$detection_time > 100){  #CADD
      delay_one[j] <- (test_one$detection_time - 100)
    }else{
      delay_one[j] <- 1
    }
  }
  if( !is.na(test_two$detection_time[1])){
    
    detected_two <- detected_two + 1
    if(test_two$detection_time > 100){ 
      delay_two[j] <- (test_two$detection_time - 100)
    }else{
      delay_two[j] <- 1
    }
  }
  print(j)
}

power_one <- detected_one / ncol(y)
power_two <- detected_two / ncol(y)

cadd_one <- mean(delay_one, na.rm = TRUE)
cadd_two <- mean(delay_two, na.rm = TRUE)

sd_cadd_one <- sd(delay_one, na.rm = TRUE)
sd_cadd_two <- sd(delay_two, na.rm = TRUE)

t_polynomial_df <- data.frame(power_one = power_one, power_two = power_two,
                              cadd_one = cadd_one, cadd_two = cadd_two,
                              sd_cadd_one = sd_cadd_one, sd_cadd_two = sd_cadd_two)

##########################################

#sinusoidal anomaly

##########################################

set.seed(100)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df = 2, sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.5)
  y[101:200, i] <- y[101:200,i] + + a[1]*(sin(2*pi*time[101:200]/50) - sin(2*pi*time[100])) +
    a[2]*(cos(2*pi*time[101:200]/50) - cos(2*pi*time[100]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

detected_one <- detected_two <- 0
delay_one <- delay_two <- rep(NA, ncol(y))

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvt(n_train, df = 2, sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  test_one <- fast(time, xfd, 1, yfd$fd[j], threshold = threshold)
  test_two <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(test_one$detection_time[1])){    #if detected
    
    detected_one <- detected_one + 1
    if(test_one$detection_time > 100){  #CADD
      delay_one[j] <- (test_one$detection_time - 100)
    }else{
      delay_one[j] <- 1
    }
  }
  if( !is.na(test_two$detection_time[1])){
    
    detected_two <- detected_two + 1
    if(test_two$detection_time > 100){ 
      delay_two[j] <- (test_two$detection_time - 100)
    }else{
      delay_two[j] <- 1
    }
  }
  print(j)
}

power_one <- detected_one / ncol(y)
power_two <- detected_two / ncol(y)

cadd_one <- mean(delay_one, na.rm = TRUE)
cadd_two <- mean(delay_two, na.rm = TRUE)

sd_cadd_one <- sd(delay_one, na.rm = TRUE)
sd_cadd_two <- sd(delay_two, na.rm = TRUE)

t_sinusoidal_df <- data.frame(power_one = power_one, power_two = power_two,
                              cadd_one = cadd_one, cadd_two = cadd_two,
                              sd_cadd_one = sd_cadd_one, sd_cadd_two = sd_cadd_two)

###########################################

#Noise Anomaly

###########################################


set.seed(50)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df = 2, sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.5, 1)
  y[175:200, i] <- y[175:200,i] - x_underlying[175:200] 
  
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

detected_one <- detected_two <- 0
delay_one <- delay_two <- rep(NA, ncol(y))

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvt(n_train, df = 2, sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  test_one <- fast(time, xfd, 1, yfd$fd[j], threshold = threshold)
  test_two <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(test_one$detection_time[1])){    #if detected
    
    detected_one <- detected_one + 1
    if(test_one$detection_time > 100){  #CADD
      delay_one[j] <- (test_one$detection_time - 100)
    }else{
      delay_one[j] <- 1
    }
  }
  if( !is.na(test_two$detection_time[1])){
    
    detected_two <- detected_two + 1
    if(test_two$detection_time > 100){ 
      delay_two[j] <- (test_two$detection_time - 100)
    }else{
      delay_two[j] <- 1
    }
  }
  print(j)
}

power_one <- detected_one / ncol(y)
power_two <- detected_two / ncol(y)

cadd_one <- mean(delay_one, na.rm = TRUE)
cadd_two <- mean(delay_two, na.rm = TRUE)

sd_cadd_one <- sd(delay_one, na.rm = TRUE)
sd_cadd_two <- sd(delay_two, na.rm = TRUE)

t_noise_df <- data.frame(power_one = power_one, power_two = power_two,
                         cadd_one = cadd_one, cadd_two = cadd_two,
                         sd_cadd_one = sd_cadd_one, sd_cadd_two = sd_cadd_two)

##################################

##################################

######################
#No underlying shape
######################

#######################

#Polynomial Anomalies

########################

set.seed(10)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- rep(3, length(time))


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
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

detected_one <- detected_two <- 0
delay_one <- delay_two <- rep(NA, ncol(y))

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  test_one <- fast(time, xfd, 1, yfd$fd[j], threshold = threshold)
  test_two <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(test_one$detection_time[1])){    #if detected
    
    detected_one <- detected_one + 1
    if(test_one$detection_time > 100){  #CADD
      delay_one[j] <- (test_one$detection_time - 100)
    }else{
      delay_one[j] <- 1
    }
  }
  if( !is.na(test_two$detection_time[1])){
    
    detected_two <- detected_two + 1
    if(test_two$detection_time > 100){ 
      delay_two[j] <- (test_two$detection_time - 100)
    }else{
      delay_one[j] <- 1
    }
  }
  print(j)
}

power_one <- detected_one / ncol(y)
power_two <- detected_two / ncol(y)

cadd_one <- mean(delay_one, na.rm = TRUE)
cadd_two <- mean(delay_two, na.rm = TRUE)

sd_cadd_one <- sd(delay_one, na.rm = TRUE)
sd_cadd_two <- sd(delay_two, na.rm = TRUE)

constant_polynomial_df <- data.frame(power_one = power_one, power_two = power_two,
                                     cadd_one = cadd_one, cadd_two = cadd_two,
                                     sd_cadd_one = sd_cadd_one, sd_cadd_two = sd_cadd_two)

##########################################

#sinusoidal anomaly

##########################################

set.seed(100)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- rep(3, length(time))

y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.5)
  y[101:200, i] <- y[101:200,i] + + a[1]*(sin(2*pi*time[101:200]/50) - sin(2*pi*time[100])) +
    a[2]*(cos(2*pi*time[101:200]/50) - cos(2*pi*time[100]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

detected_one <- detected_two <- 0
delay_one <- delay_two <- rep(NA, ncol(y))

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  test_one <- fast(time, xfd, 1, yfd$fd[j], threshold = threshold)
  test_two <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(test_one$detection_time[1])){    #if detected
    
    detected_one <- detected_one + 1
    if(test_one$detection_time > 100){  #CADD
      delay_one[j] <- (test_one$detection_time - 100)
    }else{
      delay_two[j] <- 1
    }
  }
  if( !is.na(test_two$detection_time[1])){
    
    detected_two <- detected_two + 1
    if(test_two$detection_time > 100){ 
      delay_two[j] <- (test_two$detection_time - 100)
    }else{
      delay_two[j] <- 1
    }
  }
  print(j)
}

power_one <- detected_one / ncol(y)
power_two <- detected_two / ncol(y)

cadd_one <- mean(delay_one, na.rm = TRUE)
cadd_two <- mean(delay_two, na.rm = TRUE)

sd_cadd_one <- sd(delay_one, na.rm = TRUE)
sd_cadd_two <- sd(delay_two, na.rm = TRUE)

constant_sinusoidal_df <- data.frame(power_one = power_one, power_two = power_two,
                                     cadd_one = cadd_one, cadd_two = cadd_two,
                                     sd_cadd_one = sd_cadd_one, sd_cadd_two = sd_cadd_two)

###########################################

#Noise Anomaly

###########################################


set.seed(500)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 100

x_underlying <- rep(0, length(time))

y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
#no anomaly to add as process already constant
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )

detected_one <- detected_two <- 0
delay_one <- delay_two <- rep(NA, ncol(y))

for(j in 1:ncol(y)){
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  test_one <- fast(time, xfd, 1, yfd$fd[j], threshold = threshold)
  test_two <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
  
  if( !is.na(test_one$detection_time[1])){    #if detected
    
    detected_one <- detected_one + 1
    if(test_one$detection_time > 100){  #CADD
      delay_one[j] <- (test_one$detection_time - 100)
    }else{
      delay_one[j] <- 1
    }
  }
  if( !is.na(test_two$detection_time[1])){
    
    detected_two <- detected_two + 1
    if(test_two$detection_time > 100){ 
      delay_two[j] <- (test_two$detection_time - 100)
    }else{
      delay_two[j] <- 1
    }
  }
  print(j)
}

power_one <- detected_one / ncol(y)
power_two <- detected_two / ncol(y)

cadd_one <- mean(delay_one, na.rm = TRUE)
cadd_two <- mean(delay_two, na.rm = TRUE)

sd_cadd_one <- sd(delay_one, na.rm = TRUE)
sd_cadd_two <- sd(delay_two, na.rm = TRUE)

constant_noise_df <- data.frame(power_one = power_one, power_two = power_two,
                                cadd_one = cadd_one, cadd_two = cadd_two,
                                sd_cadd_one = sd_cadd_one, sd_cadd_two = sd_cadd_two)

##################################

#Table 2

#####################################

table_2a_df <- data.frame(gp_power = c(gp_polynomial_df$power_one, gp_sinusoidal_df$power_one, gp_noise_df$power_one),
                          gp_cadd = c(gp_polynomial_df$cadd_one, gp_sinusoidal_df$cadd_one, gp_noise_df$cadd_one),
                          gp_sd_cadd = c(gp_polynomial_df$sd_cadd_one, gp_sinusoidal_df$sd_cadd_one, gp_noise_df$sd_cadd_one),
                          t_power = c(t_polynomial_df$power_one, t_sinusoidal_df$power_one, t_noise_df$power_one),
                          t_cadd = c(t_polynomial_df$cadd_one, t_sinusoidal_df$cadd_one, t_noise_df$cadd_one),
                          t_sd_cadd = c(t_polynomial_df$sd_cadd_one, t_sinusoidal_df$sd_cadd_one, t_noise_df$sd_cadd_one),
                          constant_power = c(constant_polynomial_df$power_one, constant_sinusoidal_df$power_one, constant_noise_df$power_one),
                          constant_cadd = c(constant_polynomial_df$cadd_one, constant_sinusoidal_df$cadd_one, constant_noise_df$cadd_one),
                          constant_sd_cadd = c(constant_polynomial_df$sd_cadd_one, constant_sinusoidal_df$sd_cadd_one, constant_noise_df$sd_cadd_one),
                          row.names = c("Polynomial", "Sinusoidal", "Loss of Shape"))

table_2b_df <- data.frame(gp_power = c(gp_polynomial_df$power_two, gp_sinusoidal_df$power_two, gp_noise_df$power_two),
                          gp_cadd = c(gp_polynomial_df$cadd_two, gp_sinusoidal_df$cadd_two, gp_noise_df$cadd_two),
                          gp_sd_cadd = c(gp_polynomial_df$sd_cadd_two, gp_sinusoidal_df$sd_cadd_two, gp_noise_df$sd_cadd_two),
                          t_power = c(t_polynomial_df$power_two, t_sinusoidal_df$power_two, t_noise_df$power_two),
                          t_cadd = c(t_polynomial_df$cadd_two, t_sinusoidal_df$cadd_two, t_noise_df$cadd_two),
                          t_sd_cadd = c(t_polynomial_df$sd_cadd_two, t_sinusoidal_df$sd_cadd_two, t_noise_df$sd_cadd_two),
                          constant_power = c(constant_polynomial_df$power_two, constant_sinusoidal_df$power_two, constant_noise_df$power_two),
                          constant_cadd = c(constant_polynomial_df$cadd_two, constant_sinusoidal_df$cadd_two, constant_noise_df$cadd_two),
                          constant_sd_cadd = c(constant_polynomial_df$sd_cadd_two, constant_sinusoidal_df$sd_cadd_two, constant_noise_df$sd_cadd_two),
                          row.names = c("Polynomial", "Sinusoidal", "Loss of Shape"))

table2a <- kbl(table_2a_df, digits = 2) %>% kable_classic_2(full_width = F) %>% kable_styling(font_size = 30)
save_kable(table2a, file = "./Results/table2a.png")
table2b <- kbl(table_2b_df, digits = 2) %>% kable_classic_2(full_width = F) %>% kable_styling(font_size = 30)
save_kable(table2b, file = "./Results/table2b.png")

