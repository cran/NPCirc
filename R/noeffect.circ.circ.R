noeffect.circ.circ <- function (x, y, bw, method = "LL", n_boot = 500) {

  name1 <- deparse(substitute(x))
  name2 <- deparse(substitute(y))
  DNAME <- paste(paste(name1, collapse = "\n"), "and",
                 paste(name2, collapse = "\n"))

  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same number of observations")

  if (is.circular(x)) {
    datacircularp <- circularp(x)
  }else if (is.circular(y)) {
    datacircularp <- circularp(y)
  }else {
    datacircularp <- list(type = "angles", units = "radians",
                          template = "none", modulo = "2pi", zero = 0,
                          rotation = "counter")
  }

  x <- conversion.circular(x, units = "radians", zero = 0,
                           rotation = "counter", modulo = "2pi")
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
    bw<-bw.reg.circ.circ(x, y, method = method) * 4
  }else{
    if (is.numeric(bw)){
      if (bw < 0)
        stop("Argument 'bw' must be positive")
    }else{
      stop("Argument 'bw' must be numeric")
    }
  }

  if (!is.character(method))
    stop("method must be either ''LL'' or ''NW'' ")
  if (method!="LL"&method!="NW")
    stop("method must be either ''LL'' or ''NW'' ")
  if (!is.numeric(n_boot))
    stop(" 'n_boot' must be numeric")
  if (n_boot <= 0)
    stop(" 'n_boot' must be > 0")


  me <- (mean(y)) %% (2 * pi)
  gamma_hat <- rep(me, n)
  m_hat <- (RegCircCirc(x, y, x, method = method, bw = bw)) %% (2 * pi)

  rsd0 <- sum(1 - cos(y - gamma_hat))
  rsd1 <- sum(1 - cos(y - m_hat))

  obs <- (rsd0 - rsd1) / rsd1

  res <- y - circular(gamma_hat)



  y_boot <- numeric(n_boot)
  stat_boot <- numeric(n_boot)

  for (b in 1:n_boot){
    res_boot <- sample(res, n, replace=T)
    y_boot <- (gamma_hat + res_boot) %% (2 * pi)
    me_boot <- (mean(y_boot)) %% (2 * pi)
    gamma_hat_boot <- rep(me_boot, n)
    m_hat_boot <- (RegCircCirc(x, y_boot, t=x, bw=bw, method = method)) %% (2 * pi)

    rsd0_boot <- sum(1 - cos(y_boot - gamma_hat_boot))
    rsd1_boot <- sum(1 - cos(y_boot - m_hat_boot))


    stat_boot[b] <- (rsd0_boot - rsd1_boot) / rsd1_boot

  }


  p <- sum(ifelse(stat_boot >= obs, 1, 0)) / n_boot

  meth <- "No-effect test for a circular predictor and a circular response"
  STATISTIC <- obs
  names(STATISTIC) <- "C.obs"
  PARAMETER <- bw
  names(PARAMETER) <- "bw"
  structure(list( statistic = STATISTIC, alternative = "Significant effect",
                  p.value = p, method = meth,  parameter = PARAMETER, data.name = DNAME),
            class = "htest")

}
