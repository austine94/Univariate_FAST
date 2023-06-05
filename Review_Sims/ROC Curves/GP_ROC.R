library(fda)
library(tidyverse)
library(mvtnorm)

theme_idris <- function() {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey20"),
    panel.border =  element_rect(fill = NA,
                                 colour = "grey20"),
  )
}

############

#Polynomial Outlier

####################

set.seed(10)

time <- 1:1000
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)

n_train <- 100
n_test <- 100

x_list <- vector("list", n_test)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

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
for(i in 1:ncol(y)){
  a <- runif(4, 0.5, 1.5)
  b <- runif(1, 0, 0.5)
  y[101:200, i] <- y[101:200,i] + b[1] + a[1]*((time[101:200] - time[100]) / (lt)) +
    a[2]*((time[101:200] - time[100])^2/(lt)^2) + a[3]*((time[101:200] - time[100])^3/(lt)^3)
  + a[4]*((time[101:200] - time[100])^4/(lt)^4)
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

quantiles <- c(1e-12, 1e-3, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05)

thresholds <- qnorm(1 - quantiles / (2 * (length(time) - 1)))

detected_one_results <- detected_two_results <- rep(0, length(thresholds))
delay_one_results <- delay_two_results <- rep(0, length(thresholds))

for(i in 1:length(thresholds)){
  detected_one <- detected_two <- 0

  for(j in 1:ncol(y)){

    test_one <- fast(time, x_list[[j]], 1, yfd$fd[j],
                     threshold = thresholds[i])
    test_two <- fast(time, x_list[[j]], 2, yfd$fd[j],
                     threshold = thresholds[i])

    if( !is.na(test_one$detection_time[1])){    #if detected
        detected_one <- detected_one + 1
    }
    if( !is.na(test_two$detection_time[1])){
        detected_two <- detected_two + 1
    }
    print(c(i,j))
  }
  detected_one_results[i] <- detected_one / ncol(y)
  detected_two_results[i] <- detected_two / ncol(y)
}

detected_one_results <- c(detected_one_results, 0)   #for start of ROC Curve
detected_two_results <- c(detected_two_results, 0)

df_quantiles <- c(quantiles, 0)

df_polynomial <- data.frame(false_one = df_quantiles, false_two = df_quantiles,
                            detected_one = detected_one_results, detected_two = detected_two_results)

ggplot(df_polynomial) + geom_line(aes(x = false_one, y = detected_one, col = "Order One")) +
  geom_line(aes(x = false_two, y = detected_two, col = "Order Two")) +
  labs(x = "False Alarms", y = "Anomalies Detected") +
  ylim(0,1) +
  theme_idris()

##################################

#Sinusoidal

##################################


set.seed(100)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 1000

x_list <- vector("list", n_test)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

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
for(i in 1:ncol(y)){
  a <- runif(2, 0.1, 0.5)
  y[101:200, i] <- y[101:200,i] + + a[1]*(sin(2*pi*time[101:200]/50) - sin(2*pi*time[100])) +
    a[2]*(cos(2*pi*time[101:200]/50) - cos(2*pi*time[100]))
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

quantiles <- c(1e-12, 1e-3, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05)
thresholds <- qnorm(1 - quantiles / (2 * (length(time) - 1)))


detected_one_results <- detected_two_results <- rep(0, length(thresholds))
delay_one_results <- delay_two_results <- rep(0, length(thresholds))

for(i in 1:length(thresholds)){
  detected_one <- detected_two <- 0

  for(j in 1:ncol(y)){

    test_one <- fast(time, x_list[[j]], 1, yfd$fd[j],
                     threshold = thresholds[i])
    test_two <- fast(time, x_list[[j]], 2, yfd$fd[j],
                     threshold = thresholds[i])

    if( !is.na(test_one$detection_time[1])){    #if detected
      detected_one <- detected_one + 1
    }
    if( !is.na(test_two$detection_time[1])){
      detected_two <- detected_two + 1
    }
    print(c(i,j))
  }
  detected_one_results[i] <- detected_one / ncol(y)
  detected_two_results[i] <- detected_two / ncol(y)
}

detected_one_results <- c(detected_one_results, 0)   #for start of ROC Curve
detected_two_results <- c(detected_two_results, 0)

df_sinusoidal <- data.frame(false_one = df_quantiles, false_two = df_quantiles,
                            detected_one = detected_one_results, detected_two = detected_two_results)

ggplot(df_sinusoidal) + geom_line(aes(x = false_one, y = detected_one, col = "Order One")) +
  geom_line(aes(x = false_two, y = detected_two, col = "Order Two")) +
  labs(x = "False Alarms", y = "Anomalies Detected") +
  ylim(0,1) +
  theme_idris()

#################################

#Noise Only

##################################

set.seed(10)

time <- 1:500
covar <- squared_exp_covar(1:length(time), 40, 0.3)
lt <- length(time)

n_train <- 100
n_test <- 1000

x_list <- vector("list", n_test)

x_underlying <- sin(2 * pi * time / 200) + cos(2*pi*time/200)

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
for(i in 1:ncol(y)){
  y[101:200, i] <- y[101:200,i] - x_underlying[101:200] + x_underlying[101]
}
bas <- create.bspline.basis(range(time), nbasis = 60, norder = 6)
yfd <- smooth.basis(time, y, bas)

quantiles <- c(1e-12, 1e-3, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05)

thresholds <- c(qnorm(1 - quantiles / (2 * (length(time) - 1))))

detected_one_results <- detected_two_results <- rep(0, length(thresholds))
delay_one_results <- delay_two_results <- rep(0, length(thresholds))

for(i in 1:length(thresholds)){
  detected_one <- detected_two <- 0

  for(j in 1:ncol(y)){

    test_one <- fast(time, x_list[[j]], 1, yfd$fd[j],
                     threshold = thresholds[i])
    test_two <- fast(time, x_list[[j]], 2, yfd$fd[j],
                     threshold = thresholds[i])

    if( !is.na(test_one$detection_time[1])){    #if detected
      detected_one <- detected_one + 1
    }
    if( !is.na(test_two$detection_time[1])){
      detected_two <- detected_two + 1
    }
    print(c(i,j))
  }
  detected_one_results[i] <- detected_one / ncol(y)
  detected_two_results[i] <- detected_two / ncol(y)
}

detected_one_results <- c(detected_one_results, 0)   #for start of ROC Curve
detected_two_results <- c(detected_two_results, 0)

df_quantiles <- c(quantiles, 0)

df_noise <- data.frame(false_one = df_quantiles, false_two = df_quantiles,
                            detected_one = detected_one_results, detected_two = detected_two_results)

ggplot(df_noise) + geom_line(aes(x = false_one, y = detected_one, col = "Order One")) +
  geom_line(aes(x = false_two, y = detected_two, col = "Order Two")) +
  labs(x = "False Alarms", y = "Anomalies Detected") +
  ylim(0,1) +
  theme_idris()

##############################

#ROC

###############################

g1 <- ggplot() +
  geom_line(aes(x = df_polynomial$false_one, y = df_polynomial$detected_one, col = "Polynomial")) +
  geom_line(aes(x = df_sinusoidal$false_one, y = df_sinusoidal$detected_one, col = "Sinusoidal")) +
  geom_line(aes(x = df_noise$false_one, y = df_noise$detected_one, col = "Noise")) +
  ylim(0,1) +
  labs(x = "Alpha", y = "True Positives") + theme_idris() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))
ggsave("./Results/Figure5a.png", g1)

g2 <- ggplot() +
  geom_line(aes(x = df_polynomial$false_two, y = df_polynomial$detected_two, col = "Polynomial")) +
  geom_line(aes(x = df_sinusoidal$false_two, y = df_sinusoidal$detected_two, col = "Sinusoidal")) +
  geom_line(aes(x = df_noise$false_two, y = df_noise$detected_two, col = "Noise")) +
  ylim(0,1) +
  labs(x = "Alpha", y = "True Positives") + theme_idris() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))
ggsave("./Results/Figure5b.png", g2)

