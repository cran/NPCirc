estimation_gamma <- function(g, mhat1, y, groupnumber){

  n <- length(g)
  g <- factor(g)
  levels <- levels(g)

  group <- levels[groupnumber]

  y_t <- y[g == group]
  mhat1_t <- mhat1[g == group]
  C_t <- sum(cos(y_t - mhat1_t))
  S_t <- sum(sin(y_t - mhat1_t))

  return(atan2(S_t, C_t))
}
