ancova.circ.lin <- function(x, y, g, bw, bw1, test = "eq", method = "LL", calib = "chisq", n_boot = 500){

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

  if (is.circular(y))
    y <- as.numeric(y)
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

  if (!is.character(calib))
    stop("calib must be character")
  if (calib!="chisq"&calib!="boot")
    stop("calib must be either ''chisq'' or ''boot'' ")
  if (!is.character(test))
    stop("test must be either ''eq'' or ''paral'' ")
  if (test!="eq"&test!="paral")
    stop("test must be either ''eq'' or ''paral'' ")
  if (!is.character(method))
    stop("method must be either ''LL'' or ''NW'' ")
  if (method!="LL"&method!="NW")
    stop("method must be either ''LL'' or ''NW'' ")

  if (missing(bw)) {
    bw <- bw.reg.circ.lin(x, y, method=method)
  }else if (is.numeric(bw)){
    if (bw < 0)
      stop("Argument 'bw' must be positive")
  }else{
    stop("Argument 'bw' must be numeric")
  }

  if (missing(bw1)){
    bw1 <- 1 / (NN.circ(8, x, x))^2
  }else if (is.numeric(bw1)){
    if (bw1 < 0)
      stop("Argument 'bw1' must be positive")
  }else{
    stop("Argument 'bw1' must be numeric")
  }
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
    gn[irange] <- g[g == levels[i]]
    gn[irange] <- (gn[irange])[ord]
    istart <- istart + ni[i]
  }
  x <- circular(xn)
  y <- yn
  g <- gn


  S <- kernCL(x, x, bw, method = method, tol = 300)

  Sd <- diag(nrow = n, ncol = n)
  K<-matrix(0, nrow = n, ncol = n)
  istart <- 1
  for(j in 1:nlev){
    irange <- istart:(istart + ni[j] - 1)
    Sd[irange, irange] <- kernCL(x[irange],t = x[irange], bw = bw, tol = 300)
    K[irange, irange] <- matrixK_block(as.numeric(x[irange]))
    istart <- istart + ni[j]
  }
  K <- K / (n - nlev)


  if (test == "eq"){

    Q <- t(Sd - S) %*% (Sd - S)
    m_hat <- S %*% y

  }else if (test == "paral"){


    In <- diag(n)
    D <- matrix(0, nrow = n, ncol = nlev - 1)
    istart <- ni[1] + 1
    for(i in 2:nlev){

      D[istart:(istart + ni[i] - 1), i - 1]<-1
      istart <- istart + ni[i]
    }

    S1 <- kernCL(x, x, bw1, method = method, tol = 300)

    W <- solve(t(D) %*% t(In - S1) %*% (In - S1) %*% D) %*% t(D) %*% t(In - S1) %*% (In - S1)
    gam_hat <- W %*% y

    aux <- D %*% W
    P <- (S %*% (In - aux)) + aux - Sd
    Q <- t(P) %*% P

    gamaux <- c(0, as.vector(gam_hat))
    gamma_hat <- rep(gamaux, times = ni)
    m_hat <- S %*% (In - D %*% W) %*% y

  }

  est_sigma <- sqrt((y %*% K %*% y))[1,1]
  obs <- ((y %*% Q %*% y)[1,1] / est_sigma^2)


  if (calib == "chisq"){

    V <- Q - K * obs
    k1 <- sum(diag(V))
    C <- V %*% V
    k2 <- 2 * sum(diag(C))
    k3 <- 8 * sum(diag(C %*% V))
    aa <- abs(k3 / (4 * k2))
    bb <- (8 * k2^3) / k3^2
    cc <- k1 - aa * bb
    p <- 1 - pchisq(-cc / aa, bb)

  }else if (calib == "boot"){

    if (test == "eq"){
      res <- y - m_hat
    }else if (test == "paral"){
      res <- y - m_hat - gamma_hat
    }

    y_boot <- numeric(n_boot)
    stat_boot <- numeric(n_boot)

    for (b in 1:n_boot){
      res_boot <- sample(res, n, replace = T)
      if (test == "eq"){
        y_boot <- (m_hat + res_boot)
      }else if (test == "paral"){
        y_boot <- (gamma_hat + m_hat + res_boot)

        S1_boot <- kernCL(x, x, bw1, method = method, tol = 300)

        W_boot <- solve(t(D) %*% t(In - S1_boot) %*% (In - S1_boot) %*% D) %*% t(D) %*% t(In - S1_boot) %*% (In - S1_boot)
        gam_hat_boot <- W_boot %*% y_boot
      }


      S_boot <- kernCL(x, x, bw, method = method, tol = 300)
      Sd_boot <- diag(nrow = n, ncol = n)
      K_boot <- K
      istart <- 1
      for (j in 1:nlev){
        irange <- istart:(istart + ni[j] - 1)
        Sd_boot[irange,irange] <- kernCL(x[irange], t = x[irange], bw = bw, method = method, tol = 300)
        istart <- istart + ni[j]
      }


      est_sigma_boot <- sqrt((as.numeric(y_boot) %*% K_boot %*% as.numeric(y_boot)))[1,1]

      if (test == "eq") {
        Q_boot <- t(Sd_boot - S_boot) %*% (Sd_boot - S_boot)
        stat_boot[b] <- (as.numeric(y_boot) %*% Q_boot %*% as.numeric(y_boot))[1,1] / est_sigma_boot^2
      }else if (test == "paral"){
        aux_boot <- D %*% W_boot
        P_boot <- (S_boot %*% (In - aux_boot)) + aux_boot - Sd_boot
        Q_boot <- t(P_boot) %*% P_boot

        stat_boot[b] <- (as.numeric(y_boot) %*% Q_boot %*% as.numeric(y_boot))[1,1]/(est_sigma_boot^2)
      }
    }

    p <- sum(ifelse(stat_boot >= obs, 1, 0)) / n_boot

  }

  if (test == "eq"){
    meth <- "Equality test for a circular predictor and a real-valued response"
  }else if (test == "paral"){
    meth <- "Parallelism test for a circular predictor and a real-valued response"
  }
  STATISTIC <- obs
  names(STATISTIC) <- "C.obs"
  PARAMETER <- bw
  names(PARAMETER) <- "bw"
  structure(list( statistic = STATISTIC, alternative = "The curves are different",
                  p.value = p, method = meth,  parameter = PARAMETER, data.name = DNAME),
            class = "htest")
}
