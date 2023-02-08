####Script for examples for "Sequential Detection of Emergent Anomalous Structures in
####Functional Data

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
