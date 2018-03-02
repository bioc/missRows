RVcoeff <- function(X, Y) {
  Y <- scale(Y, scale = FALSE)
  X <- scale(X, scale = FALSE)
  W1 <- X %*% t(X)
  W2 <- Y %*% t(Y)
  rv <- sum(diag(W1 %*% W2)) / (sum(diag(W1 %*% W1)) * sum(diag(W2 %*% W2)))^0.5
  return(rv)
}
