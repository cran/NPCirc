noeffect.circ.lin <- function (x, y, bw, method = "LL", calib = "chisq", n_boot = 500) {

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
  }else{
    datacircularp <- list(type = "angles", units = "radians",
                          template = "none", modulo = "2pi", zero = 0,
                          rotation = "counter")
  }

  x <- conversion.circular(x, units = "radians", zero = 0,
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

  if (missing(bw)) {
    bw <- bw.reg.circ.lin(x, y, method = method) * 4
  } else{
    if (is.numeric(bw)){
      if (bw<0)
        stop("Argument 'bw' must be positive")
    }else{
      stop("Argument 'bw' must be numeric")
    }
  }
  if (!is.character(calib))
    stop("calib must be either ''chisq'' or ''boot'' ")
  if (calib!="chisq"&calib!="boot")
    stop("calib must be either ''chisq'' or ''boot'' ")
  if (!is.character(method))
    stop ("method must be either ''LL'' or ''NW'' ")
  if (method!="LL"&method!="NW")
    stop("method must be either ''LL'' or ''NW'' ")
  if (!is.numeric(n_boot))
    stop(" 'n_boot' must be numeric")
  if (n_boot <= 0)
    stop(" 'n_boot' must be > 0")

  me <- mean(y)
  S <- kernCL(x, x, bw, method = method, tol = 300)
  m_hat <- S %*% y

  L <- matrix(rep(1/n,n^2),nrow=n)
  I <- diag(n)
  A <- t(I - S) %*% (I - S)
  B <- I - L - A

  rss0 <- as.numeric(t(y) %*% (I-L) %*% y)
  rss1 <- as.numeric(t(y) %*% A %*% y)

  obs <- (rss0 - rss1) / rss1

  if (calib == "chisq") {
    T <- B - A * obs
    k1 <- sum(diag(T))
    C <- T %*% T
    k2 <- 2 * sum(diag(C))
    k3 <- 8 * sum(diag(C %*% T))
    aa <- abs(k3 / (4 * k2))
    bb <- (8 * k2^3) / k3^2
    cc <- k1 - aa * bb
    p <- 1 - pchisq(-cc / aa, bb)

  }else{

    res <- y - me

    y_boot <- numeric(n_boot)
    stat_boot <- numeric(n_boot)

    for (b in 1:n_boot){

      res_boot <- sample(res, n, replace=T)
      y_boot <- me + res_boot

      me_boot <- mean(y_boot)
      S_boot <- kernCL(x, x, bw, tol = 300)
      m_hat_boot <- S_boot %*% y_boot

      A_boot <- t(I - S_boot) %*% (I - S_boot)
      B_boot <- I - L - A_boot

      rss0_boot <- as.numeric(t(y_boot) %*% (I - L) %*% y_boot)
      rss1_boot <- as.numeric(t(y_boot) %*% A_boot %*% y_boot)

      stat_boot[b] <- (rss0_boot-rss1_boot) / rss1_boot

    }

    p <- sum(ifelse(stat_boot >= obs, 1, 0)) / n_boot

  }
  meth <- "No-effect test for a circular predictor and a real-valued response"
  STATISTIC <- obs
  names(STATISTIC) <- "C.obs"
  PARAMETER <- bw
  names(PARAMETER) <- "bw"

  structure(list(  statistic = STATISTIC,  alternative = "Significant effect",
                   p.value = p, method = meth,  parameter = PARAMETER, data.name = DNAME),
            class = "htest")

}
