library(fda)
library(tidyverse)
library(mvtnorm)

n_reps <- 100
n_obs <- 100

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)
bas <- create.bspline.basis(range(time), norder = 12, nbasis = 60)

m_vec <- 1:5  #orders to consider

##############
#Order 1 Function
###############
set.seed(1)

BIC_vals <- matrix(NA, nrow = n_reps, ncol = length(m_vec))

for(j in 1:n_reps){ #add noise to underlying shape and create fd object
  a <- runif(1, 1, 3)
  x_underlying <- a*exp(-a*time/500) 
  x <- matrix(rep(x_underlying, n_obs), nrow = length(time), ncol = n_obs, byrow = FALSE)
  noise <- rmvnorm(n_obs, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  #compute BIC and select order
  BIC_vals[j,] <- sapply(m_vec, function(x) pda_bic(time, xfd, x)[1])
  print(j)
}
order_selected <- apply(BIC_vals, 1, which.min)
correct_choice_model1 <- sum(order_selected == 1) / n_reps

##############
#Order 2 - Sinusoidal
##############

set.seed(2)

BIC_vals <- matrix(NA, nrow = n_reps, ncol = length(m_vec))

for(j in 1:n_reps){
  
  a <- runif(2, 1, 3)
  
  x_underlying <- a[1]*sin(2 * pi * time / 200) + a[2]*cos(2 * pi * time / 200)
  
  x <- matrix(rep(x_underlying, n_obs), nrow = length(time), ncol = n_obs, byrow = FALSE)
  noise <- rmvnorm(n_obs, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  BIC_vals[j,] <- sapply(m_vec, function(x) pda_bic(time, xfd, x)[1])
  print(j)
}
order_selected <- apply(BIC_vals, 1, which.min)

correct_choice_model2 <- sum(order_selected == 2) / n_reps

##############
#Order 3 - Sinusoidal and Polynomial
##############

set.seed(3)


BIC_vals <- matrix(NA, nrow = n_reps, ncol = length(m_vec))

for(j in 1:n_reps){
  a <- runif(3, 1, 3)
  x_underlying <- a[1]*(time/250)^3 + a[2]*((time-100)/500)^5 + a[3]*sin(2*pi*time/40)  #solution to a third order Cauchy-Euler ODE
  x <- matrix(rep(x_underlying, n_obs), nrow = length(time), ncol = n_obs, byrow = FALSE)
  noise <- rmvnorm(n_obs, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  BIC_vals[j,] <- sapply(m_vec, function(x) pda_bic(time, xfd, x)[1])
  print(j)
}
order_selected <- apply(BIC_vals, 1, which.min)

correct_choice_model3 <- sum(order_selected == 3) / n_reps

##############
#Order 4 - Mixed
##############

set.seed(4)
covar <- squared_exp_covar(1:length(time), 40, 0.3)

BIC_vals <- matrix(NA, nrow = n_reps, ncol = length(m_vec))

for(j in 1:n_reps){
  a <- runif(4, 1, 3)
  x_underlying <- a[1]*sin(2*pi*time/50) + a[2]*cos(2*pi*time/20)  + a[3]*sinh(-0.5+time/500) + a[4]*cosh(-0.5+time/500)
  #x_underlying <- a[1]*sin(2 * pi * time / 400) + a[2]*cos(2 * pi * time / 400) + a[3] * exp(-2*pi*time/400) * sin(2 * pi * time / 400) +  a[4] * exp(-2*pi*time/400) * cos(2 * pi * time / 400)
  x <- matrix(rep(x_underlying, n_obs), nrow = length(time), ncol = n_obs, byrow = FALSE)
  noise <- rmvnorm(n_obs, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  BIC_vals[j,] <- sapply(m_vec, function(x) pda_bic(time, xfd, x)[1])
  print(j)
}
order_selected <- apply(BIC_vals, 1, which.min)

correct_choice_model4 <- sum(order_selected == 4) / n_reps

##############
#Order 2 - Sinusoidal Heavy Tails
##############

set.seed(5)

BIC_vals <- matrix(NA, nrow = n_reps, ncol = length(m_vec))

for(j in 1:n_reps){
  
  a <- runif(2, 1, 3)
  
  x_underlying <- a[1]*sin(2 * pi * time / 200) + a[2]*cos(2 * pi * time / 200)
  
  x <- matrix(rep(x_underlying, n_obs), nrow = length(time), ncol = n_obs, byrow = FALSE)
  noise <- rmvt(n_obs, df=2, sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  BIC_vals[j,] <- sapply(m_vec, function(x) pda_bic(time, xfd, x)[1])
  print(j)
}
order_selected <- apply(BIC_vals, 1, which.min)

correct_choice_model5 <- sum(order_selected == 2) / n_reps


########################
#Order 4 but combination of 2 sinusoidal order 2's, showing a possible reason for error
#####################
set.seed(6)
covar <- squared_exp_covar(1:length(time), 40, 0.3)

BIC_vals <- matrix(NA, nrow = n_reps, ncol = length(m_vec))

for(j in 1:n_reps){
  a <- runif(4, 1, 3)
  x_underlying <- a[1]*sin(2 * pi * time / 400) + a[2]*cos(2 * pi * time / 400) + a[3] * exp(-2*pi*time/400) * sin(2 * pi * time / 200) +  a[4] * exp(-2*pi*time/400)* cos(2 * pi * time / 200)
  x <- matrix(rep(x_underlying, n_obs), nrow = length(time), ncol = n_obs, byrow = FALSE)
  noise <- rmvnorm(n_obs, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  BIC_vals[j,] <- sapply(m_vec, function(x) pda_bic(time, xfd, x)[1])
  print(j)
}
order_selected <- apply(BIC_vals, 1, which.min)

correct_choice_model6 <- sum(order_selected == 4) / n_reps
chosen_as_order2 <- sum(order_selected == 2) / n_reps

########################
result_df_bic <- data.frame(correct_order = c(correct_choice_model1,
                                              correct_choice_model2,
                                              correct_choice_model3,
                                              correct_choice_model4,
                                              correct_choice_model5))
rownames(result_df_bic) <- c("Order 1", "Order 2", "Order 3", "Order 4", 
                             "Order 2 HT")
#correct_choice_model1
#correct_choice_model2
#correct_choice_model3
#correct_choice_model4
#correct_choice_model5
#correct_choice_model6
#chosen_as_order2
