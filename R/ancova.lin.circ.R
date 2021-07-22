ancova.lin.circ <- function(x, y, g, bw, bw1, test = "eq", method = "LL", n_boot = 500){


  name1 <- deparse(substitute(x))
  name2 <- deparse(substitute(y))
  name3 <- deparse(substitute(g))
  DNAME <- paste(paste(name1, collapse = "\n"),
                 paste(name2, collapse = "\n"), "and", paste(name3, collapse = "\n"))

  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same number of observations")
  if (length(x) != length(g))
    stop("argument 'g' must have the same length as 'x' and 'y'")

  if (is.circular(y)) {
    datacircularp <- circularp(y)
  }else{
    datacircularp <- list(type = "angles", units = "radians",
                          template = "none", modulo = "2pi", zero = 0,
                          rotation = "counter")
  }

  y <- conversion.circular(y, units = "radians", zero = 0,
                           rotation = "counter", modulo = "2pi")

  nax <- is.na(x)
  nay <- is.na(y)
  x <- x[!nax & !nay]
  y <- y[!nax & !nay]
  if ((sum(nax) + sum(nay)) > 0)
    warning("Missing values were removed.", "\n")
  n <- length(x)
  if (n == 0)
    stop("No observations (at least after removing missing values)")

  if (missing(bw)){
    bw <- bw.reg.lin.circ(x, y, method = method)
  }else{
    if (is.numeric(bw)){
      if (bw < 0)
        stop("Argument 'bw' must be positive")
    }else{
      stop("Argument 'bw' must be numeric")
    }
  }
  if (missing(bw1)) {
    bw1 <- NN(8, x, x)
  }else if (is.numeric(bw1)){
    if (bw1 < 0)
      stop("Argument 'bw1' must be positive")
  }else{
    stop("Argument 'bw1' must be numeric")
  }


  if (!is.character(method))
    stop("method must be either ''LL'' or ''NW'' ")
  if (method!="LL"&method!="NW")
    stop("method must be either ''LL'' or ''NW'' ")
  if (!is.character(test))
    stop("test must be either ''eq'' or ''paral'' ")
  if (test!="eq"&test!="paral")
    stop("test must be either ''eq'' or ''paral'' ")
  if (!is.numeric(n_boot))
    stop(" 'n_boot' must be numeric")
  if (n_boot <= 0)
    stop(" 'n_boot' must be > 0")


  gf <- factor(g)
  levels <- levels(gf)
  ni <- table(gf)
  nlev <- length(levels)
  g <- as.character(g)

  istart <- 1
  xn <- numeric(n)
  yn <- numeric(n)
  gn <- character(n)
  for (i in 1:nlev){

    irange <- istart:(istart + ni[i] - 1)
    xn[irange] <- sort(x[g == levels[i]])
    ord <- order(x[g == levels[i]])
    yn[irange] <- y[g == levels[i]]
    yn[irange] <- (yn[irange])[ord]
    gn[irange] <- g[g==levels[i]]
    gn[irange] <- (gn[irange])[ord]
    istart <- istart + ni[i]

  }

  x <- xn
  y <- circular(yn)
  g <- gn


  if (test == "eq"){
    m_hat <- RegLinCirc(x, y, t = x, bw = bw, method = method)
    res <- (y - m_hat)
  }else if (test == "paral"){
    m_hat1 <- RegLinCirc(x, y, t = x, bw = bw1, method = method)
    gamma_hat <- sapply(1:nlev, function(i) estimation_gamma(g, m_hat1, y, i))
    gamma_hats <- rep(gamma_hat, times = ni)
    m_hat <- RegLinCirc(x, y - gamma_hats, t = x, bw = bw, method = method)
    res <- y - gamma_hats - m_hat
  }


  mi_hat <- numeric(n)
  istart <- 1
  for (i in 1:nlev){
    irange <- istart:(istart + ni[i] - 1)
    mi_hat[irange] <- RegLinCirc(x[g == levels[i]], y[g == levels[i]], t = x[g == levels[i]], bw = bw, method = method)
    istart <- istart + ni[i]
  }


  dispersion <- (1 / (n-nlev)) * sum(1 - cos(y - mi_hat))

  if (test == "eq"){
    obs <- (1 / dispersion) * sum(1 - cos(mi_hat - m_hat))
  }else if (test == "paral"){
    obs <- (1 / dispersion) * sum(1 - cos(mi_hat - m_hat - gamma_hats))
  }


  y_boot <- numeric(n_boot)
  stat_boot <- numeric(n_boot)

  for (b in 1:n_boot){
    res_boot <- sample(res, n, replace = T)
    if (test == "eq"){
      y_boot <- (m_hat + res_boot) %% (2 * pi)
      m_hat_boot <- RegLinCirc(x, y_boot, t = x, bw = bw, method = method)
    }else if (test == "paral"){
      y_boot <- (gamma_hats + m_hat + res_boot) %% (2 * pi)
      m_hat1_boot <- RegLinCirc(x, y_boot, t = x, bw = bw1, method = method)
      gamma_hat_boot <- sapply(1:nlev, function(i) estimation_gamma(g, m_hat1_boot, y_boot, i))
      gamma_hats_boot <- rep(gamma_hat_boot, times = ni)
      m_hat_boot <- RegLinCirc(x, y_boot - gamma_hats_boot, t = x, bw = bw, method = method)
    }

    mi_hat_boot <- numeric(n)
    istart <- 1
    for (i in 1:nlev){
      irange <- istart:(istart + ni[i] - 1)
      mi_hat_boot[irange] <- RegLinCirc(x[g == levels[i]], y_boot[g == levels[i]], t = x[g == levels[i]], bw = bw, method = method)
      istart <- istart + ni[i]
    }
    dispersion_boot <- (1 / (n - nlev)) * sum(1 - cos(y_boot - mi_hat_boot))
    if (test == "eq"){
      stat_boot[b] <- (1 / dispersion_boot) * sum(1 - cos(mi_hat_boot - m_hat_boot))
    }else if (test == "paral"){
      stat_boot[b] <- (1 / dispersion_boot) * sum(1 - cos(mi_hat_boot - m_hat_boot - gamma_hats_boot))
    }
  }

  p <- sum(ifelse(stat_boot >= obs, 1, 0)) / n_boot

  if (test == "eq"){
    meth <- "Equality test for a real-valued predictor and a circular response"
  }else if (test == "paral"){
    meth <- "Parallelism test for a real-valued predictor and a circular response"
  }

  STATISTIC <- obs
  names(STATISTIC) <- "C.obs"
  PARAMETER <- bw
  names(PARAMETER) <- "bw"
  structure(list( statistic = STATISTIC, alternative = "The curves are different",
                  p.value = p, method = meth,  parameter = PARAMETER, data.name = DNAME),
            class = "htest")
}
