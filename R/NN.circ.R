NN.circ <- function(k, x, t){

  n <- length(x)
  x <- c(x - 2 * pi, x, x + 2 * pi)
  t <- c(t - 2 * pi, t, t + 2 * pi)
  distance <- abs(outer(x, t, "-"))
  distance <- apply(distance, 2, sort)
  range <- (n + 1):(2 * n)
  bw <- distance[k + 1, range]
  bw[(1 / bw^2) > 50] <- 0.1

  return(bw)
}
