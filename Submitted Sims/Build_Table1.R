##################
#This Script Will Build Table 1 - Just run it in its entirety and it will save
#the table in the results folder.
#################

#############
#Wrapper for FAST
#############

fast <- function(time, training_fd, m, new_observation, threshold, basis = NULL, penalty_order = 1,
                 penalty_param = 0){
  
  if(is.na(threshold)){
    stop("threshold value must be a non-negative numeric when choosing the manual method")
  }
  
  fast_run <- pdafast_manual_threshold(time, training_fd, m, threshold,
                                       basis, penalty_order, penalty_param, new_observation)
  
  output_list <- list(detection_time = fast_run$detection_time, threshold = fast_run$threshold,
                      training_fd = training_fd, new_observation = new_observation)
  class(output_list) <- "FAST"             #S3 output for plot function
  
  return(output_list)
  
  
}

##########################
#Inner function for FAST
###########################


pdafast_manual_threshold <- function(time, fd, m, threshold, basis = NULL, penalty_order = 1,
                                     penalty_param = 0, new_observation){
  
  #time is the vector of points the fd object is observed over
  #fd is the functional data object of class fdSmooth
  #threshold is the threshold for the CUSUM test
  #basis is the basis object used for the pda weights
  #penalty param and penalty order are the smoothing parameters for fitting the basis functions
  #new_observation is the new fd obs to test
  
  #Now begin the CUSUM function:
  
  if(!is.numeric(time)){
    stop("time must be a vector of numerics")
  }
  
  if(class(fd) != "fdSmooth"){
    stop("Class of fd must be fdSmooth")
  }
  
  if(is.fd(new_observation) != TRUE){
    stop("Must be functional data object")
  }
  
  if(round(m) != m | m < 1){
    stop("m needs to be a positive integer")
  }
  
  n <- length(time)
  
  if( !is.numeric(threshold) | threshold <= 0){
    stop("threshold must be a positive numeric")
  }
  
  if( !is.null(basis) ){
    
    if(class(basis) != "basisfd"){
      stop("basis must either be null, or of class basisfd")
    }
    
    if( all(basis$rangeval != fd$fd$basis$rangeval)){
      stop("Basis range must match the range of the basis for the fd object")
    }
    
  }else{
    
    basis = fd$fd$basis   #use the original basis
  }
  
  if(round(penalty_order) != penalty_order | penalty_order < 0){
    stop("penalty_order needs to be a positive integer")
  }
  
  if( !is.numeric(penalty_param) | penalty_param < 0){
    stop("penalty_param must be a positive numeric")
  }
  
  #first, fit the ODE
  
  ODE <- ode_fit(time, fd, m, basis = basis, penalty_order = penalty_order, penalty_param = penalty_param)
  
  #Now the DE coeffs, and the residuals, are obtained move to the fitting procedure
  
  eval_fd_mat <- eval.fd(time, fd$fd, ODE$ODE)  #This evaluates the residual functions
  
  #Next scale the evaluated residuals to obtain null hypothesis observations
  
  eval_mat <- eval_mat_scale(eval_fd_mat)
  
  #We can now test this against a new observation:
  #Calculate the absolute one step change, and standardise
  
  obs_step_change <- eval.fd(time, new_observation, Lfdobj = ODE$ODE)
  obs_step_change <- obs_step_change %>% diff() %>% abs()
  obs_step_change <- obs_step_change^2
  
  obs_step_change <- (obs_step_change - eval_mat$mean_value) #standardise to mean 0
  
  cusum <- cumsum(obs_step_change) %>% abs()
  len_cusum <- length(cusum)
  cusum <- cusum / (sqrt(1:len_cusum) * eval_mat$sd_value_bm)   #standardise to unit variance
  #Now perform the detection
  if( sum(cusum > threshold) == 0 ){
    detection = NA
  }else{
    detection <- min(which(cusum > threshold))  #first time detection occurs
  }
  
  detection <- detection + 1 #add one as we start testing at time 2
  
  return(list(detection_time = detection, threshold = threshold))
  
}


##################
#ODE Fitting Function
##################

ode_fit <- function(time, fd, m, basis = NULL, penalty_order = 1, penalty_param = 0){
  
  #time is the range of time the function is observed over
  #fd is the functional data object to have pda fitted to it
  #m is the order of the fd object
  #basis is the basis system used to smooth the coeff functions
  #penalty_order is the order of smoothing penalty to use
  #penalty_param is the value of the penalty used to control smoothing
  
  #function returns the ODE, of type Lfd, and the residuals of type "fd"
  
  if(!is.numeric(time)){
    stop("time must be a vector of numerics")
  }
  
  if(class(fd) != "fdSmooth"){
    stop("Class of fd must be fdSmooth")
  }
  
  if(round(m) != m | m < 1){
    stop("m needs to be a positive integer")
  }
  
  if( !is.null(basis) ){
    
    if(class(basis) != "basisfd"){
      stop("basis must either be null, or of class basisfd")
    }
    
    if( all(basis$rangeval != fd$fd$basis$rangeval)){
      stop("Basis range must match the range of the basis for the fd object")
    }
    
  }else{
    
    basis = fd$fd$basis   #use the original basis
  }
  
  if(round(penalty_order) != penalty_order | penalty_order < 0){
    stop("penalty_order needs to be a positive integer")
  }
  
  if( !is.numeric(penalty_param) | penalty_param < 0){
    stop("penalty_param must be a positive numeric")
  }
  
  basis <- fd(NULL, basis)    #this ensures the basis system is of class "fdsmooth"
  
  fdpar <- fdPar(basis, penalty_order, penalty_param)  #Add smoothing penalty
  
  #Create list for basis of weight functions
  
  bwt_list <- vector("list", m)
  bwt_list <- lapply(bwt_list, function(x) fdpar)  #m basis systems for m coeffs
  
  fd_list <- list(fd$fd)   #Functional Data object to have DE fitted to it
  
  pda <- pda.fd(fd_list, bwt_list)  #apply pda
  
  fitted_equation_coefs <- vector("list", m)   #For the fitted equation coeffs to be stored
  
  #now fill the coeffs list to provide the list of m ODE coeffs fitted using pda
  
  for(i in 1:m){
    fitted_equation_coefs[[i]] <- pda$bwtlist[[i]]$fd  #Each entry is of type fd
  }
  
  #now create the linear ODE functional data object:
  
  Ldobj <- Lfd(m, fitted_equation_coefs)  #Linear differential object
  
  return(list(ODE = Ldobj, residuals = pda$resfdlist[[1]]) )
  
}

#####################
#Test Statistic Scaling
######################

#this function scales the matrix of evaluated residual functions to have zero mean and unit variance
#at each point in time

eval_mat_scale <- function(eval_mat){
  
  step_change <- diff(eval_mat)  #size of one step change
  
  step_change <- (step_change)^2   #squared score function
  
  mean_value <- apply(step_change, 1, mean)  #takes the mean at each point in time
  
  null_mat <- sweep(step_change, 1, mean_value, "-")  #centre each timepoint
  
  test_stat <- apply(null_mat, 2, cumsum)   #cumsum
  
  n <- nrow(null_mat)
  
  scaled_test_stat <- sweep(test_stat, 1, sqrt(1:n), "/" )  #scale to N(0, sigma^2)
  
  sd_value_bm <- apply(scaled_test_stat, 1, sd)   #sd for the N(, sigma^2)
  
  rescaled_test_stat <- sweep(scaled_test_stat, 1, sd_value_bm, "/")
  
  return(list(null_mat = null_mat, mean_value = mean_value,
              sd_value_bm = sd_value_bm, rescaled_test_stat = rescaled_test_stat))
  
}

########
#Function for BIC
########

pda_bic <- function(time, fd, m, basis = NULL, penalty_order = 1, penalty_param = 0){
  #time is the obs region
  #fd is the vanilla fd object
  #m is the ODE order
  
  if(!is.numeric(time)){
    stop("time must be a vector of numerics")
  }
  
  if( !inherits(fd, "fdSmooth") ){
    stop("fd must either be of class fdSmooth")
  }
  
  if(round(m) != m | m < 1){
    stop("Error m: m needs to be a positive integer")
  }
  
  if( !is.null(basis) ){
    
    if(!inherits(basis, "basisfd")){
      stop("basis must either be null, or of class basisfd")
    }
    
    if( all(basis$rangeval != fd$fd$basis$rangeval)){
      stop("Basis range must match the range of the basis for the fd object")
    }
    
  }else{
    
    basis = fd$fd$basis   #use the original basis
  }
  
  if(round(penalty_order) != penalty_order | penalty_order < 0){
    stop("penalty_order needs to be a positive integer")
  }
  
  if( !is.numeric(penalty_param) | penalty_param < 0){
    stop("penalty_param must be a positive numeric")
  }
  
  #first fit the ODE
  
  len_t <- length(time)
  n_obs <- ncol(fd$y)
  
  ODE <- ode_fit(time, fd, m, basis = basis, penalty_order = penalty_order,
                 penalty_param = penalty_param)
  residual_functions <- eval.fd(time, ODE$residuals)  #This evaluates the residual functions
  scaled_residual_functions <- residual_functions
  for(i in 1:nrow(residual_functions)){  #scale to unit variance for llhd expression to hold
    scaled_residual_functions[i,] <- (residual_functions[i,]) / sd(residual_functions[i,])
  }
  #mean zero unit variance
  SSE <- scaled_residual_functions^2 %>% sum()   #sum of square error
  
  likelihood <- (-n_obs/2)*(log(SSE/n_obs))
  BIC <- m * log(n_obs) - 2*likelihood
  return(c(BIC, SSE))
}    


############
#Function for squared exponential covariance
############

squared_exp_covar <- function(time, lengthscale, amplitude){
  
  if(!is.numeric(time)){
    stop("time must be a vector of numerics")
  }
  
  if( !(is.numeric(lengthscale)) | lengthscale < 0){
    stop("lengthscale must be a non-negative numeric")
  }
  
  if( !(is.numeric(amplitude)) | amplitude <= 0){
    stop("amplitude must be a positive numeric")
  }
  
  len_t <- length(time)
  
  covar <- matrix(0, nrow = len_t, ncol = len_t)
  
  if(lengthscale == 0){
    diag(covar) <- rep(amplitude, len_t)
  }else{
    for(i in 1:len_t){
      for(j in 1:i){
        covar[i,j] <- covar[j,i] <- amplitude * exp(-0.5 * (( (time[i] - time[j]) / lengthscale)^2) )
      }
    }
  }
  
  
  return(covar)
  
}

##########
#Function for matern covariance
##########

matern_covar <- function(time, rho, nu, amplitude){
  
  #matern covariance matrix observed over interval of time
  #note that nu is the parameter controlling differentiability, so it is nu - 1 times diff
  
  if(!is.numeric(time)){
    stop("time must be a vector of numerics")
  }
  
  if( !(is.numeric(rho)) | rho <= 0){
    stop("rho must be a positive numeric")
  }
  
  if( !(is.numeric(nu)) | nu <= 0){
    stop("nu must be a positive numeric")
  }
  
  
  if( !(is.numeric(amplitude)) | amplitude <= 0){
    stop("amplitude must be a positive numeric")
  }
  
  len_t <- length(time)
  
  covar <- matrix(0, nrow = len_t, ncol = len_t)
  
  for(i in 1:len_t){
    for(j in 1:i){
      
      distance <- abs(i-j)
      covar[i,j] <- covar[j,i] <- amplitude * (2^(1-nu) / gamma(nu)) *
        ( sqrt(2*nu) * distance / rho)^nu *
        besselK(sqrt(2*nu) * distance / rho, nu)
    }
  }
  
  diag(covar) <- 1
  
  return(covar)
  
}

###########
#FAST Plot
##########

fast_plot <- function(fast_run){
  
  if(class(fast_run) != "FAST"){
    stop("Must input an object of class FAST")
  }
  
  #plot training data and overlay the new observation
  
  if( is.na(fast_run$detection_time[1]) ){  #if no anomaly
    
    plot(fast_run$new_observation)
    lines(fast_run$training_fd, col = "grey")
    lines(fast_run$new_observation, col = "green", lty = 1, lwd = 3)
    
  }else{   #if anomaly detected add detection time
    
    plot(fast_run$new_observation)
    abline(v = fast_run$detection_time, col = "red", lty = 2)
    lines(fast_run$training_fd, col = "grey")
    lines(fast_run$new_observation, col = "red", lty = 1, lwd = 3)
    
    
  }
  
}

###############
#FAST End of Day
###############

fast_end_of_day <- function(time, training_fd, m, new_observation, threshold, basis = NULL, penalty_order = 1,
                            penalty_param = 0){
  
  if(is.na(threshold)){
    stop("threshold value must be a non-negative numeric when choosing the manual method")
  }
  
  fast_run <- pdafast_manual_threshold_end_of_day(time, training_fd, m, threshold,
                                                  basis, penalty_order, penalty_param, new_observation)
  
  output_list <- list(detection_time = fast_run$detection_time, threshold = fast_run$threshold,
                      training_fd = training_fd, new_observation = new_observation)
  class(output_list) <- "FAST"             #S3 output for plot function
  
  return(output_list)
  
  
}


pdafast_manual_threshold_end_of_day <- function(time, fd, m, threshold, basis = NULL, penalty_order = 1,
                                                penalty_param = 0, new_observation){
  
  #time is the vector of points the fd object is observed over
  #fd is the functional data object of class fdSmooth
  #threshold is the threshold for the CUSUM test
  #basis is the basis object used for the pda weights
  #penalty param and penalty order are the smoothing parameters for fitting the basis functions
  #new_observation is the new fd obs to test
  
  #Now begin the CUSUM function:
  
  if(!is.numeric(time)){
    stop("time must be a vector of numerics")
  }
  
  if(class(fd) != "fdSmooth"){
    stop("Class of fd must be fdSmooth")
  }
  
  if(is.fd(new_observation) != TRUE){
    stop("Must be functional data object")
  }
  
  if(round(m) != m | m < 1){
    stop("m needs to be a positive integer")
  }
  
  n <- length(time)
  
  if( !is.numeric(threshold) | threshold <= 0){
    stop("threshold must be a positive numeric")
  }
  
  if( !is.null(basis) ){
    
    if(class(basis) != "basisfd"){
      stop("basis must either be null, or of class basisfd")
    }
    
    if( all(basis$rangeval != fd$fd$basis$rangeval)){
      stop("Basis range must match the range of the basis for the fd object")
    }
    
  }else{
    
    basis = fd$fd$basis   #use the original basis
  }
  
  if(round(penalty_order) != penalty_order | penalty_order < 0){
    stop("penalty_order needs to be a positive integer")
  }
  
  if( !is.numeric(penalty_param) | penalty_param < 0){
    stop("penalty_param must be a positive numeric")
  }
  
  #first, fit the ODE
  
  ODE <- ode_fit(time, fd, m, basis = basis, penalty_order = penalty_order, penalty_param = penalty_param)
  
  #Now the DE coeffs, and the residuals, are obtained move to the fitting procedure
  
  eval_fd_mat <- eval.fd(time, fd$fd, ODE$ODE)  #This evaluates the residual functions
  
  #Next scale the evaluated residuals to obtain null hypothesis observations
  
  eval_mat <- eval_mat_scale(eval_fd_mat)
  
  #We can now test this against a new observation:
  #Calculate the absolute one step change, and standardise
  
  obs_step_change <- eval.fd(time, new_observation, Lfdobj = ODE$ODE)
  obs_step_change <- obs_step_change %>% diff() %>% abs()
  obs_step_change <- obs_step_change^2
  
  obs_step_change <- (obs_step_change - eval_mat$mean_value) #standardise to mean 0
  
  cusum <- cumsum(obs_step_change) %>% abs()
  len_cusum <- length(cusum)
  cusum <- cusum / (sqrt(1:len_cusum) * eval_mat$sd_value_bm)   #standardise to unit variance
  #Now perform the detection
  if( cusum[len_cusum] > threshold){
    detection <- 500
  }else{
    detection <- NA
  }
  
  return(list(detection_time = detection, threshold = threshold))
  
}
library(fda)
library(tidyverse)
library(mvtnorm)

result_df <- data.frame(FAST_Detected = NA, FAST_Delay = NA, FAST_Delay_sd = NA)

##################
#Model 1 - GP Residuals, No Anomalies
################

set.seed(3000)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 500

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- n_test
fast_delay <- rep(NA, n_test)

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
  
}
if(all(is.na(fast_delay))){
  fast_delay <- rep(Inf, n_test)
}

#compute ADD as in Tartakovsky (2014)
fast_ADD <- mean(fast_delay, na.rm = TRUE) / (fast_detected_after_change / n_test) 

new_result_df_row <- c(fast_detected/n_test, fast_ADD, 
                       sd(fast_delay, na.rm = TRUE))
print(fast_detected)
result_df <- rbind(result_df, new_result_df_row)
################
#Model 2 - T Residuals, No Anomaly
#################

set.seed(10)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- n_test
fast_delay <- rep(NA, n_test)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
y <- y + noise
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

set.seed(200)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- n_test
fast_delay <- rep(NA, n_test)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- sample(c(-1,1), 3, TRUE) * runif(3, 0.01, 0.1)
  y[100:200, i] <- y[100:200,i] + a[1] + a[2]*((time[100:200]-time[100])/lt)
  a[3]*exp(-(time[100:200] - time[100]))
}
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
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
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
#Model 4 - GP Residuals, Polynomial Anomaly, Random Start
#################

set.seed(30)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- n_test
fast_delay <- rep(NA, n_test)

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
    if(fast_run$detection_time > start_times[i]){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - (start_times[i] - 1))
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
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
#Model 5- Sinusoidal Anomaly
######################

set.seed(1213)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- n_test
fast_delay <- rep(NA, n_test)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.3)
  y[101:200, i] <- y[101:200,i] + a[1]*(sin(2*pi*time[101:200]/50) - sin(2*pi*time[100])) +
    a[2]*(cos(2*pi*time[101:200]/50) - cos(2*pi*time[100]))
}
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
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
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
fast_detected_after_change <- n_test
fast_delay <- rep(NA, n_test)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  y[175:200, i] <- y[175:200,i] - x_underlying[175:200]
}
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
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  
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

set.seed(20000)

covar2 <- squared_exp_covar(101:200, 35, 0.3)  #for shape anomaly

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- n_test
fast_delay <- rep(NA, n_test)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
noise_anomaly <- rmvnorm(n_test, mean = rep(0, 100), sigma = covar2) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  y[101:200, i] <- y[101:200,i] - noise[101:200, i] + noise_anomaly[, i]
}
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
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
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

set.seed(8000)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- n_test
fast_delay <- rep(NA, n_test)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
yfd <- smooth.basis(time, y, bas)
y <- eval.fd(time, yfd$fd)
for(i in 1:ncol(y)){
  translation <- sample(1:10, 1)
  y[101:200, i] <- y[101:200,i] - x_underlying[101:200] + x_underlying[(101 - translation):(200-translation)]
}
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
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{  
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1    }
  }
  print(fast_run$detection_time)
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

set.seed(333)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

fast_detected <- 0
fast_detected_after_change <- n_test
fast_delay <- rep(NA, n_test)

#generate the 100  test observations
y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
y <- y + noise
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.3)
  y[101:200, i] <- y[101:200,i] + a[1]*(sin(2*pi*time[101:200]/50) - sin(2*pi*time[100])) +
    a[2]*(cos(2*pi*time[101:200]/50) - cos(2*pi*time[100]))
}
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
    if(fast_run$detection_time > 100){  #ADD as in Tartakovsky
      fast_delay[j] <- (fast_run$detection_time - 100)
    }else{
      fast_delay[j] <- 0 
      fast_detected_after_change <- fast_detected_after_change - 1
    }
  }
  
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
result_df_order2 <- result_df[-1,]
rownames(result_df_order2) <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6",
                                "Model 7", "Model 8", "Model 9")
result_df_order2 <- result_df_order2[,-3]
colnames(result_df_order2) <- c("Power", "ADD")
result_df_order2 <- result_df_order2 %>%  mutate_if(is.numeric, round, digits = 2)
result_df_order2[1:2, 2] <- "-"
#####################
#return Table 4.1
library(gridExtra)
S4.1_df <- t(result_df_order2) %>% as.data.frame()
pdf("./Results/Table1.pdf", height=11, width=20)
grid.table(S4.1_df)
dev.off()
