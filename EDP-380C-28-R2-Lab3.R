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


## Generating data for a multiple regression

# make function to translate parameter list to inputs expected by rmvnorm(n, mu, Sigma)
# input, p_x, is a list with the following elements:
# n = scalar; sample size
# cor = scalar; correlation between generated variabels
# mu_x = column vector of means
# sd_x = column vector of standard deviations
param_translate <- function(p_x) {
  n = n
  mu = t(mu_x)
  diag_sigma = diag(c(sd_x[1], sd_x[2]), 2, 2)
  R_x = matrix(c(1, p_x$cor, p_x$cor, 1), 2, 2)
  Sigma = diag_sigma %*% R_x %*% diag_sigma
  return(list(n = n, mu = mu, Sigma = Sigma))
}

# generate X1 and X2 using rmvnorm() function with the following parameters;

set.seed(21389)
n <- 100000

# set correlation to 0.3
r12 <- 0.3

# set mean vector
mu_x <- matrix(c(5, 10), nrow = 1, ncol = 2)

# set st dev
sd_x <- matrix(1:2, nrow = 2, ncol = 1)

p_x <- list(n = n, cor = r12, mu_x = mu_x, sd_x = sigma_x)

param_rmvnorm <- param_translate(p_x)

data_3a <- rmvnorm(param_rmvnorm$n, param_rmvnorm$mu, param_rmvnorm$Sigma)

x1 = data_3a[, 1]
x2 = data_3a[, 2]
X = cbind(x1, x2)
b1 = 1
b2 = 1
R2 = 0.6
mu_y = 10
sd_y = 5
p_y <- list(X = X, b1 = b1, b2 = b2, R2 = R2, mu_y = mu_y, sd_y = sd_y)

generate_y <- function(p_y) {
  B = matrix(c(p_y$b1, p_y$b2), ncol = 1)
  Sigma_x <- cov(p_y$X)
  sd_e <- sqrt( (t(B) %*% Sigma_x %*% B) %*% ((1 / R2) - 1) )
  y = rnorm(n, X %*% B, sd_e)
}

set.seed(23921)

y <- generate_y(p_y)

solution3b <- produce_OLS_solution(y, X)

diff_y_param <- list(mean = mean(y) - mu_y, sd = sd(y) - sd_y)

#> produce_OLS_solution(data_3b, X)
#Estimate                  SE          t_value           p_value
#b0    0.06294727  0.0404683095833365 1.55547058359421 0.119837199875773
#b1    0.98916088 0.00677235973269319 146.058526650088                 0
#b2    0.99973586  0.0033740950999435 296.297476432644                 0
#SD(e) 2.03641868                  NA               NA                NA
#R2    0.59782931                  NA               NA                NA

#> diff_y_param <- list(mean = mean(y) - mu_y, sd = sd(y) - sd_y)
#> diff_y_param
#$mean
#[1] 4.99389

#$sd
#[1] -1.783684


set.seed(123782)
S_x <- cov(X)
S_xy <- cov(X, y)
R_x <- cor(X)
R_xy <- cor(X, y)
R_yx <- cor(y, X)
S_y <- var(y)

beta <- crossprod(S_x, S_xy)

resid_var <- S_y - crossprod(beta, S_xy)

R2 <- R_yx %*% R_x %*% R_xy

# shortcuts
tcrossprod(X) = X %*% t(X)
crossprod(X) = t(X) %*% X
crossprod(X, Y) = t(X) %*% Y

solve(crossprod(X) %*% crossprod(X, Y)) =
solve((crossprod(X)), crossprod(X, Y))



