#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression
#        and Data Generation
#
#' The function, rnvnorm(), provides an implementation for a function to generate data from a p-dimensional multivariate normal distribution
#'
#' @param n set sample size
#' @param mu set a column vector of means
#' @param Sigma set covariance matrix
#' @return function returns a matrix of data generated from p standard normal deviates
rmvnorm <- function(n, mu, Sigma) {
  if(nrow(mu) == nrow(Sigma)) {
    z = matrix(1, n, nrow(mu))
    for (i in 1:nrow(mu)) {
      z[, i] = rnorm(n)
    }
    unit_vector_nx1 <- matrix(1, n, 1)
    y <- unit_vector_nx1 %*% t(mu) + z %*% chol(Sigma)
    return(as.matrix(y))
  }

  if(nrow(mu) != nrow(Sigma)) {
    stop("input dimensions do not match")
  }
}
