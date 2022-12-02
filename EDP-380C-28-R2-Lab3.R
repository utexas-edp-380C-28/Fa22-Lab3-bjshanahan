#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression
#        and Data Generation
#
# Name:  Brendan Shanahan
#---------------------------------------------------#

## SETUP -------------------------------------------
# Source rmvnorm()
source("scripts/rmvnorm.R")

## Begin Code --------------------------------------

## Generating data for simple regression
# set parameters for X and Y distributions
m_x <- -2
s_x <- 3
p_x <- list(m_x = m_x, s_x = s_x)
b0 <- 12
b1 <- 4
s_e <- 7
p_y <- list(b0 = b0, b1 = b1, s_e = s_e)

# create a function generating X and Y
# p_x = list of mean and st dev of x, respectively; p_y list of mean and st dev of y, respectively
generate_reg <- function(n, p_x, p_y) {
  x <- with(p_x, rnorm(n, m_x, s_x))  # generate X using predefined parameters
  y <- with(p_y, rnorm(n, b0 + b1 * x, s_e))  # generate y using equation 9
  return(as.matrix(cbind(x = x, y = y))) # output is a n x 2 column matrix
}

# create an analyze function to extract
bivariate_data_parameters <- function(df) {
  mod_summary = summary(lm(df[,2] ~ df[,1], ))
  m_x = mean(df[, 1])
  s_x = sd(df[, 1])
  m_y = mean(df[, 2])
  s_y = sd(df[, 2])
  cor_xy = cor(df[,1], df[,2])
  b0 = mod_summary$coefficients[1, 1]
  b1 = mod_summary$coefficients[2, 1]
  s_e = mod_summary$sigma
  c(m_x = m_x, m_y = m_y, s_x = s_x, s_y = s_y, cor = cor_xy, "b0" = b0, "b1" = b1, s_e = s_e)
}

set.seed(17290)
data_gen_reg <- replicate(500, generate_reg(100, p_x, p_y))

parameters <-  apply(data_gen_reg, 3, bivariate_data_parameters)

answers_1a <- list(mean = apply(parameters, 1, mean), sd = apply(parameters, 1, sd))

# answers_1a
# $mean
# m_x        m_y        s_x        s_y        cor         b0         b1        s_e
# -2.0061317  3.9810216  3.0034995 13.8596003  0.8631812 11.9780298  3.9873856  6.9749258

# $sd
# m_x        m_y        s_x        s_y        cor         b0         b1        s_e
# 0.31115637 1.41924855 0.20911389 1.03767609 0.02432744 0.86119359 0.24140973 0.49070350


produce_OLS_solution <- function(y, X) {
  # add column for intercept, if necessary
  ifelse(sum(X[, 1]) == nrow(X), X, X <- cbind(rep(1, nrow(X)), X) )

  # compute b0, b1...bp, SD(e), R2
  beta = solve(t(X) %*% X) %*% t(X) %*% y
  y_hat = X %*% beta
#  H = X %*% solve((t(X) %*% X)) %*% t(X)
  e = y - y_hat
  SD_e = sqrt( (t(e) %*% e) / (nrow(X) - ncol(X) - 1) )
  R2 = ( (t(beta) %*% cov(X) %*% beta) / ( (t(beta) %*% cov(X) %*% beta) + SD_e^2) )
  dfs_residual = nrow(X) - ncol(X) - 1
  OLS_solution = list(beta = beta, SD_e = SD_e, R2 = R2)

    # extract estimate, SE, t value and p value for b0, b1...bp
  Estimate = beta
    # compute standard errors for regression coefficients
  SE = sqrt((SD_e)^2 %*% diag(solve(t(X) %*%  X)))
    # compute t-value
  t_value = Estimate / t(SE)
    # compute p value
  p_value = (2 * pt(abs(t_value), dfs_residual, lower.tail = FALSE))

  output = data.frame(Estimate = c(b = beta, SD_e = SD_e, R2 = R2),
             SE = c(SE, "NA", "NA"),
             t_value = c(t_value, "NA", "NA"),
             p_value = c(p_value, "NA", "NA") )

  rownames(output) <- c(paste0(rep("b", nrow(beta)), 0:(nrow(beta) - 1)),
                        "SD(e)",
                        "R2")

  return(output)
}

output <- produce_OLS_solution(mtcars$mpg, cbind(mtcars$wt, mtcars$cyl, mtcars$gear))
output

summary(lm(mpg ~ wt + cyl + gear, data = mtcars))


## Generating data for a multiple regression by setting the total variance explained

# make function to translate parameter list to inputs expected by rmvnorm(n, mu, Sigma)
# input, p_x, is a list with the following elements:
# n = scalar; sample size
# cor = scalar; correlation between generated variabels
# mu_x = column vector of means
# sd_x = column vector of standard deviations
param_translate <- function(p_x) {
  n = n
  mu = mu_x
  diag_sigma = diag(c(sd_x))
  R_x = matrix(p_x$cor, dim(diag_sigma)[1], dim(diag_sigma)[1])
  diag(R_x) <- 1
  Sigma = diag_sigma %*% R_x %*% diag_sigma
  return(list(n = n, mu = mu, Sigma = Sigma))
}

# generate X1 and X2 using rmvnorm() function with the following parameters:

set.seed(21389)
n <- 100000

# set correlation between predictors to 0.3
r12 <- 0.3

# set mean vector
x <- c(5, 10)
mu_x <- matrix(x, nrow = length(x), ncol = 1)

# set st dev
sd_x <- 1:2
sd_x <- matrix(sd_x, nrow = length(sd_x), ncol = 1)

p_x <- list(n = n, cor = r12, mu_x = mu_x, sd_x = sd_x)

# use function to translate parameter list to inputs expected by rmvnorm(n, mu, Sigma)
param_rmvnorm <- param_translate(p_x)

# generate X given parameters
X <- rmvnorm(param_rmvnorm$n, param_rmvnorm$mu, param_rmvnorm$Sigma)

# set y parameters
set.seed(23921)

b1 = 1
b2 = 1
betas = matrix(c(b1, b2), nrow = 2, ncol = 1)
R2 = 0.6
mu_y = 10
sd_y = 5
p_y <- list(betas = betas, R2 = R2, mu_y = mu_y, sd_y = sd_y)

# generate y setting the total variance explained
generate_y_R2 <- function(X, p_y) {
  B = matrix(betas, ncol = 1)
  b0 = mu_y - t(mu_x) %*% B
  Sigma_x <- cov(X)
  sd_e <- sqrt( (t(B) %*% Sigma_x %*% B) %*% ((1 / R2) - 1) )
  y = rnorm(n, as.vector(b0) + X %*% B, sd_e)
}

y <- generate_y_R2(X, p_y)

solution3b <- produce_OLS_solution(y, X)

diff_y_param <- list(mean = mean(y) - mu_y, sd = sd(y) - sd_y)

#> solution3b <- produce_OLS_solution(y, X)
#> solution3b
#Estimate                  SE           t_value p_value
#b0    -4.9370527  0.0404683095833365 -121.997997586475       0
#b1     0.9891609  0.0067723597326932  146.058526650129       0
#b2     0.9997359 0.00337409509994351  296.297476432642       0
#SD(e)  2.0364187                  NA                NA      NA
#R2     0.5978293                  NA                NA      NA

#> diff_y_param <- list(mean = mean(y) - mu_y, sd = sd(y) - sd_y)
#> diff_y_param
#$mean
#[1] 0.001649459

#$sd
#[1] -1.78886


## Generating data setting the marginal correlations (rhoyx1, rhoyx2...)

set.seed(123782)
x1 = X[, 1]
x2 = X[, 2]
rho_yx1 =  0.3
rho_yx2 = -0.4
mu_y = 10
sd_y = 5
param_y_3c <- list(X = X,
                   mu_x = matrix(c(mean(X[,1]), mean(X[,2])), ncol = 1),
                   rho_yx1 =  0.3,
                   rho_yx2 = -0.4,
                   mu_y = 10,
                   sd_y = 5)

generate_y_setRhoxy <- function(p_y) {
  S_x = cov(X)
  S_xy = matrix(c((rho_yx1 * sd(X[,1]) * sd_y), (rho_yx2 * sd(X[,2]) * sd_y) ) )
  beta = solve(S_x) %*% S_xy
  var_e = (sd_y ^ 2) - crossprod(beta, S_xy)
  R_yx = matrix(c(rho_yx1 =  0.3, rho_yx2 = -0.4), nrow = 2, ncol = 1)
  R_x = as.matrix(cor(X))
  R_xy = t(R_yx)
  R2 = R_xy %*% solve(R_x) %*% R_yx
  b0 = mu_y - t(mu_x) %*% beta
  Sigma_x <- cov(p_y$X)
  sd_e <- sqrt( (t(beta) %*% Sigma_x %*% beta) %*% ((1 / R2) - 1) )
  y = rnorm(n, as.vector(b0) + X %*% beta, sd_e)
}

y_margin_corr <- generate_y_setRhoxy(param_y_3c)

solution3c <- produce_OLS_solution(y_margin_corr, X)

diff_y_param_3c <- list(mean = mean(y_margin_corr) - mu_y, sd = sd(y_margin_corr) - sd_y)
diff_y_param_3c

# > solution3c
# Estimate                  SE           t_value p_value
# b0    11.9151171  0.0795430601701363  149.794551894501       0
# b1     2.3209383  0.0133115077762783  174.355777041011       0
# b2    -1.3531807 0.00663200050404578 -204.038092166062       0
# SD(e)  4.0027116                  NA                NA      NA
# R2     0.3575643                  NA                NA      NA
#
# > diff_y_param_3c
# $mean
# [1] -0.01271435
#
# $sd
# [1] -0.006152551



# Data generation setting total variance explained

set.seed(6972)
# set parameters
n <- 100000
rho <- 0.15

# set mean vector
x <- c(2:6)
mu_x <- matrix(x, nrow = length(x), ncol = 1)

# set st dev
sd_x <- sqrt(1:5)
sd_x <- matrix(sd_x, nrow = length(sd_x), ncol = 1)

p_x <- list(n = n, cor = rho, mu_x = mu_x, sd_x = sd_x)

# create Covariance matrix, Sigma, from correlations and variances of x predictors
p_x <- param_translate(p_x)

# set p_y
b1 = 1
b2 = 1
b3 = 1
b4 = 1
b5 = 1
betas = matrix(c(b1, b2, b3, b4, b5), nrow = 5, ncol = 1)
R2 = 0.5
mu_y = 25
sd_y = 5
b0 = mu_y - t(mu_x) %*% betas
p_y <- list(intercept = b0, betas = betas, R2 = R2, mu_y = mu_y, sd_y = sd_y)

gen_data_setR2 <- function(n, p_x, p_y) {
  # generate X
  X = rmvnorm(p_x$n, p_x$mu, p_x$Sigma)
  y = generate_y_R2(X, p_y)
  output_data <- list(X = X, y = y)
  # add as attributes to the output the population parameters
  attr(output_data, 'pop param') <- p_y
  return(output_data)
}

data_method_1 <- gen_data_setR2(n, p_x, p_y)
data_method_1

solution3.d.2 <- produce_OLS_solution(data_method_1$y, data_method_1$X)
solution3.d.2

# > solution3.d.2
# Estimate                  SE          t_value p_value
# b0    5.0905785  0.0650750923269285 78.2262203281901       0
# b1    0.9812633  0.0156515010851895  62.694517412879       0
# b2    1.0199224  0.0110811906847618 92.0408646338139       0
# b3    0.9881525 0.00905975770179736 109.070524541459       0
# b4    0.9852174 0.00784277723624359 125.620985738752       0
# b5    1.0050576 0.00705536110150286 142.453033477147       0
# SD(e) 4.8135804                  NA               NA      NA
# R2    0.4991120                  NA               NA      NA

pop_param <- matrix(c(
  intercept = attr(data_method_1,"pop param")$intercept,
  b1 = attr(data_method_1,"pop param")$betas[1, 1],
  b2 = attr(data_method_1,"pop param")$betas[2, 1],
  b3 = attr(data_method_1,"pop param")$betas[3, 1],
  b4 = attr(data_method_1,"pop param")$betas[4, 1],
  b5 = attr(data_method_1,"pop param")$betas[5, 1],
  SD_e = attr(data_method_1,"pop param")$sd_y,
  R2 = attr(data_method_1,"pop param")$R2 ), 8, 1)

beta_diff <- matrix(c(
  intercept = solution3.d.2[1, 1] - attr(data_method_1,"pop param")$intercept,
  b1 = solution3.d.2[2, 1] - attr(data_method_1,"pop param")$betas[1, 1],
  b2 = solution3.d.2[3, 1] - attr(data_method_1,"pop param")$betas[2, 1],
  b3 = solution3.d.2[4, 1] - attr(data_method_1,"pop param")$betas[3, 1],
  b4 = solution3.d.2[5, 1] - attr(data_method_1,"pop param")$betas[4, 1],
  b5 = solution3.d.2[6, 1] - attr(data_method_1,"pop param")$betas[5, 1],
  SD_e = solution3.d.2[7, 1] - attr(data_method_1,"pop param")$sd_y,
  R2 = solution3.d.2[8, 1] - attr(data_method_1,"pop param")$R2 ), 8, 1)

Solution = data.frame(Estimate = solution3.d.2[, 1],
                    Paramater = pop_param,
                    Difference = beta_diff )

# > Solution
#    Estimate Paramater    Difference
# 1 5.0905785       5.0  0.0905785102
# 2 0.9812633       1.0 -0.0187366927
# 3 1.0199224       1.0  0.0199223718
# 4 0.9881525       1.0 -0.0118474752
# 5 0.9852174       1.0 -0.0147825927
# 6 1.0050576       1.0  0.0050575912
# 7 4.8135804       5.0 -0.1864196485
# 8 0.4991120       0.5 -0.0008880124


# Data generation setting marginal correlations

set.seed(1237)
# set parameters
n <- 100000
rho <- 0.15

# set mean vector
x <- c(2:6)
mu_x <- matrix(x, nrow = length(x), ncol = 1)

# set st dev
sd_x <- sqrt(1:5)
sd_x <- matrix(sd_x, nrow = length(sd_x), ncol = 1)

p_x <- list(n = n, cor = rho, mu_x = mu_x, sd_x = sd_x)

p_y <- list(rho_yx1 = -0.15,
            rho_yx2 = -0.50,
            rho_yx3 = 0.15,
            rho_yx4 = 0.30,
            rho_yx5 = .20,
            mu_y = 10,
            sd_y = 4)

# generate X covariance matrix

# unused function
# gen_covX <- function(X) {
# S_x = cov(X)
# sd_X = matrix(apply(X, 2, sd))
# rho_yx = c(rho_yx1 = -0.15,
#            rho_yx2 = -0.50,
#            rho_yx3 = 0.15,
#            rho_yx4 = 0.30,
#            rho_yx5 = .20)
# S_xy = matrix(c((rho_yx * sd_X * sd_y)))
# }

# computing the covariance matrix


gen_data_setRhoxy <- function(n, p_x, p_y) {
  #  # create Covariance matrix, Sigma, from correlations and variances of x predictors
  param_rmvnorm <- param_translate(p_x)
  # generate X
  X <- rmvnorm(param_rmvnorm$n, param_rmvnorm$mu, param_rmvnorm$Sigma)

  # generate Y
  y = rnorm(n, as.vector(b0) + X %*% beta, sd_e)
  cbind(X, y)
}

data_method_2 <- gen_data_setRhoxy(n, p_x, p_y)
data_method_2

solution3.d.2 <- produce_OLS_solution(y, X)
solution3.d.2

diff_y_param <- list(mean = mean(y) - mu_y, sd = sd(y) - sd_y)
diff_y_param



