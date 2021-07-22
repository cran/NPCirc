matrixK_block <- function(x){

  n <- length(x)
  xx <- sort(x)
  xx1 <- c(diff(xx), xx[1] - xx[n])
  xx2 <- c(xx[2] - xx[n], diff(xx, lag=2), xx[1] - xx[n-1])
  a <- xx1 / xx2
  b <- c(xx1[n], xx1[-n]) / xx2
  a[is.na(a)] <- 0.5
  b[is.na(b)] <- 0.5
  c <- sqrt(a^2 + b^2 + 1)

  A <- diag(0, n)
  diag(A) <- -1 / c
  A[(row(A) - col(A)) == -1] <- (b / c)[-n]
  A[(row(A) - col(A)) == 1] <- (a / c)[-1]
  A[1, n] <- (a / c)[1]
  A[n, 1] <- (b / c)[n]

  A <- t(A) %*% A

  return(A)
}
