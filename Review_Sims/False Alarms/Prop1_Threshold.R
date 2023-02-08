library(mvtnorm)
library(fda)
library(tidyverse)
library(kableExtra)


#######

#False Alarm Simulation

#######################

#Gaussian Process Residuals

##############
set.seed(601)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

n_train <- 100
n_test <- 100

x_list <- vector(mode = "list", n_test)

bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)

for(i in 1:length(x_list)){
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_list[[i]] <- xfd
  print(i)
}

y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
yfd <- smooth.basis(time, y, bas)

quantiles <- c(0.01, 0.05, 0.1, 0.2)
thresholds <- qnorm( 1 - (quantiles / (2*(length(time)-1)) ) )

false_one_results <- false_two_results <- rep(0, length(thresholds))


for(i in 1:length(thresholds)){
  false_one <- false_two <- 0
  
  for(j in 1:ncol(y)){
    
    test_one <- fast(time, x_list[[j]], 1, yfd$fd[j], 
                      threshold = thresholds[i])
    test_two <- fast(time, x_list[[j]], 2, yfd$fd[j], 
                     threshold = thresholds[i])
    
    if( !is.na(test_one$detection_time[1])){    #if detected
        false_one <- false_one + 1
    }
    if( !is.na(test_two$detection_time[1])){
        false_two <- false_two + 1
    }
    
    print(c(i,j))
  }
  
  false_one_results[i] <- false_one / ncol(y)
  false_two_results[i] <- false_two / ncol(y)
  
}

gp_df <- data.frame(false_one = false_one_results, false_two = false_two_results)

#######################

#T-Process residuals

##############

set.seed(10)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

n_train <- 100
n_test <- 100

x_list <- vector(mode = "list", n_test)

for(i in 1:length(x_list)){
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvt(n_train, df=2, sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_list[[i]] <- xfd
  print(i)
}

y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvt(n_test, df=2, sigma = covar) %>% t()
y <- y + noise
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

quantiles <- c(0.01, 0.05, 0.1, 0.2)
thresholds <- qnorm( 1 - (quantiles / (2*(length(time)-1)) ) )

false_one_results <- false_two_results <- rep(0, length(thresholds))


for(i in 1:length(thresholds)){
  false_one <- false_two <- 0
  
  for(j in 1:ncol(y)){
    
    test_one <- fast(time, x_list[[j]], 1, yfd$fd[j], 
                     threshold = thresholds[i])
    test_two <- fast(time, x_list[[j]], 2, yfd$fd[j], 
                     threshold = thresholds[i])
    
    if( !is.na(test_one$detection_time[1])){    #if detected
      false_one <- false_one + 1
    }
    if( !is.na(test_two$detection_time[1])){
      false_two <- false_two + 1
    }
    
    print(c(i,j))
  }
  
  false_one_results[i] <- false_one / ncol(y)
  false_two_results[i] <- false_two / ncol(y)
  
}

t_df <- data.frame(false_one = false_one_results, false_two = false_two_results)


#######################

#No underlying shape

##############

set.seed(100)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

x_underlying <- rep(3, length(time))

n_train <- 100
n_test <- 100

x_list <- vector(mode = "list", n_test)

for(i in 1:length(x_list)){
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  x_list[[i]] <- xfd
  print(i)
}

y <- matrix(rep(x_underlying, n_test), nrow = length(time), ncol = n_test, byrow = FALSE)
noise <- rmvnorm(n_test, mean = rep(0, length(time)), sigma = covar) %>% t()
y <- y + noise
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

quantiles <- c(0.01, 0.05, 0.1, 0.2)
thresholds <- qnorm( 1 - (quantiles / (2*(length(time)-1)) ) )

false_one_results <- false_two_results <- rep(0, length(thresholds))


for(i in 1:length(thresholds)){
  false_one <- false_two <- 0
  
  x <- matrix(rep(x_underlying, n_train), nrow = length(time), ncol = n_train, byrow = FALSE)
  noise <- rmvnorm(n_train, mean = rep(0, length(time)), sigma = covar) %>% t()
  x <- x + noise
  xfd <- smooth.basis(time, x, bas)
  
  for(j in 1:ncol(y)){
    
    test_one <- fast(time, x_list[[j]], 1, yfd$fd[j], 
                      threshold = thresholds[i])
    test_two <- fast(time, x_list[[j]], 2, yfd$fd[j], 
                     threshold = thresholds[i])
    
    if( !is.na(test_one$detection_time[1])){    #if detected
      false_one <- false_one + 1
    }
    if( !is.na(test_two$detection_time[1])){
      false_two <- false_two + 1
    }
    
    print(c(i,j))
  }
  
  false_one_results[i] <- false_one / ncol(y)
  false_two_results[i] <- false_two / ncol(y)
  
}

constant_df <- data.frame(false_one = false_one_results, false_two = false_two_results)

#######
#Table 1
#######

table_df <- data.frame(Scenario1 = gp_df$false_two, Scenario2 = t_df$false_two, 
                       Scenario3 = constant_df$false_two,
                       row.names = c("0.01", "0.05", "0.1", "0.2")) 
table_1 <- kbl(table_df) %>% kable_classic_2(full_width = F) %>% kable_styling(font_size = 30)
save_kable(table_1, file = "./Results/table1.png")
