dpreg.circ<-function(x, y, k = 2, ktilde = 1, startvmu = NULL, startvgam = NULL, tol= 0.000001, maxit = 300){

  name <- deparse(substitute(x))
  n <- length(data)
  datax <- x
  datay <- y
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same number of observations")
  if (!is.null(t) && is.circular(t)) {
    datacircularp <- circularp(t)
  }
  else if (is.circular(x)) {
    datacircularp <- circularp(x)
  }
  else {
    datacircularp <- list(type = "angles", units = "radians",
                          template = "none", modulo = "2pi", zero = 0,
                          rotation = "counter")
  }
  dc <- list()
  dc$type <- datacircularp$type
  dc$units <- datacircularp$units
  dc$template <- datacircularp$template
  dc$modulo <- datacircularp$modulo
  dc$zero <- datacircularp$zero
  dc$rotation <- datacircularp$rotation

  x <- conversion.circular(x, units = "radians", zero = 0,
                           rotation = "counter", modulo = "2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  nax <- is.na(x)
  nay <- is.na(y)
  x <- x[!nax & !nay]
  y <- y[!nax & !nay]
  if ((sum(nax) + sum(nay)) > 0)
    warning("Missing values were removed.", "\n")
  n <- length(x)
  if (n == 0)
    stop("No observations (at least after removing missing values)")
  if (!is.numeric(k))
    stop("argument 'k' must be numeric")
  if (k <= 0)
    stop("argument 'k' must be greater than 0")
  if (!is.numeric(ktilde))
    stop("argument 'ktilde' must be numeric")
  if (ktilde <= 0)
    stop("argument 'ktilde' must be greater than 0")
  if (is.null(startvmu)){
    startvmu<-c(log(mean(y)),rep(0,k-1))
  }else if(!is.numeric(startvmu)){
    warning("'startvmu' must be numeric. Default values were used")
    startvmu<-c(log(mean(y)),rep(0,k-1))
  }else if(length(startvmu)!=k){
    warning("'startvmu' must have the same length as 'k'. Default values were used")
    startvmu<-c(log(mean(y)),rep(0,k-1))
  }
  if (is.null(startvgam)){
    startvgam<-rep(0,ktilde)
  }else if(!is.numeric(startvgam)){
    warning("'startvgam' must be numeric. Default values were used")
    startvgam<-rep(0,ktilde)
  }else if(length(startvgam)!=ktilde){
    warning("'startvgam' must have the same length as 'ktilde'. Default values were used")
    startvgam<-rep(0,ktilde)
  }
  if (!is.numeric(tol))
    stop("argument 'tol' must be numeric")
  if (tol <= 0)
    stop("argument 'tol' must be greater than 0")
  if (!is.numeric(maxit))
    stop("argument 'maxit' must be numeric")
  if (maxit <= 0)
    stop("argument 'maxit' must be greater than 0")

  attr(x, "class") <- attr(x, "circularp") <- NULL


  est<-Param_DP(x, y, dim_mu = k, dim_gam = ktilde, startvmu = startvmu, startvgam = startvgam ,tol, maxit)
  structure(list(datax = datax, datay = datay, coefficients_mu = est[[1]], coefficients_gam = est[[2]],
                 numit = est[[3]], n = n, call = match.call(),
                 data.name = name, has.na = FALSE))

}

