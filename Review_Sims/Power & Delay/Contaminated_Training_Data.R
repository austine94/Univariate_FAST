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
n_test <- 100
n_contaminated <- c(1, 5, 10)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- rep(0, 3)

for(k in 1:length(n_contaminated)){
  #generate the 100  test observations
  #also store the contaminated observations to add to the training data at each iteration
  y <- matrix(rep(x_underlying, (n_test + n_contaminated[k]*n_test)), nrow = length(time), 
              ncol = n_test, byrow = FALSE)
  noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
  y <- y + noise
  bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
  yfd <- smooth.basis(time, y, bas)
  y_eval <- eval.fd(time, yfd$fd)
  
  threshold <- qnorm( 1 - (0.05 / (2*(length(time)-1)) ) )
  
  for(j in 1:ncol(y)){
    #training data for each rep of the simulation
    x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
    x[,((n_test - n_contaminated[k] + 1):n_test)] <- {
                             y_eval[,((n_test + ((j-1)*n_contaminated + 1)):(n_test + j*n_contaminated))]}
    #above line adds contaminated functions to the training data
    noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
    x <- x + noise
    xfd <- smooth.basis(time, x, bas)
    
    fast_run <- fast(time, xfd, 2, yfd$fd[j], threshold = threshold)
    
    if( !is.na(fast_run$detection_time[1])){    #if detected
      fast_detected[k] <- fast_detected[k] + 1
    }
    print(c(k,j))
  }
  
}
result_df <- rbind(result_df, fast_detected)