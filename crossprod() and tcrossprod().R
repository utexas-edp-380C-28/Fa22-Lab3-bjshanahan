# matrix multiplication shortcuts
tcrossprod(X) = X %*% t(X)
crossprod(X) = t(X) %*% X
crossprod(X, Y) = t(X) %*% Y

solve(crossprod(X) %*% crossprod(X, Y)) =
solve((crossprod(X)), crossprod(X, Y))
