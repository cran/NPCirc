NN <- function(k, x, t){

  n <- length(x)
  distance <- abs(outer(x, t, "-"))
  distance <- apply(distance, 2, sort)
  bw <- numeric(n)
  ind <- rep(k + 1, n)
  for (i in 1:n){
    while(distance[ind[i],i] == 0){ind[i] <- ind[i] + 1}
    bw[i] <- distance[ind[i],i]
  }

  return(bw)

}
