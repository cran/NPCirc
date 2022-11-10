

modalreg.circ.lin<-function(x, y, t=NULL, bw=NULL, tol = 0.0001, maxit = 500, from = circular(0),
                            to = circular(2 * pi), len = 300){

  name <- deparse(substitute(x))
  datax <- x
  datay <- y
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same number of observations")
  if (tol<=0)
    stop("argument 'tol' must be greater than 0")
  if (maxit<=0)
    stop("argument 'maxit' must be greater than 0")
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
  if (is.null(t)) {
    if (!is.numeric(from))
      stop("argument 'from' must be numeric")
    if (!is.numeric(to))
      stop("argument 'to' must be numeric")
    if (!is.finite(from))
      stop("non-finite 'from'")
    if (!is.finite(to))
      stop("non-finite 'to'")
    if (!is.numeric(len))
      stop("argument 'len' must be numeric")
    if (len <= 0)
      stop("argument 'len' must be integer and positive")
    from <- conversion.circular(from, units = "radians",
                                zero = 0, rotation = "counter")
    attr(from, "class") <- attr(from, "circularp") <- NULL
    to <- conversion.circular(to, units = "radians",
                              zero = 0, rotation = "counter")
    attr(to, "class") <- attr(to, "circularp") <- NULL
    if (from > to)
      stop("argument 'from' must be smaller than argument 'to'")
    t <- circular(seq(from = from, to = to, length = len))
  }
  if (tol<=0)
    stop("argument 'tol' must be greater than 0")
  if (maxit<=0)
    stop("argument 'maxit' must be greater than 0")
  else {
    if (!is.numeric(t))
      stop("argument 't' must be numeric")
    t.na <- is.na(t)
    t <- t[!t.na]
    if (sum(t.na) > 0)
      warning("'t' contains missing values. They were removed")
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
  if(is.null(bw)){
    bw <- bw.modalreg.circ.lin(x, y, lower = NULL, upper = NULL, tol = 0.00001)
  }else{
    if(!is.numeric(bw) | length(bw)!=2)
      stop("'bw' must be a vector of length two")
    if(any(bw<0))
      stop("The components of 'bw' must be positive.  The values of 'bw' were computed by using modal cross-validation")
  }
  attr(x, "class") <- attr(x, "circularp") <- NULL
  tt <- conversion.circular(t, dc$units, dc$type, dc$template,
                            dc$modulo, dc$zero, dc$rotation)
  t <- conversion.circular(t, units = "radians", modulo = "2pi",
                           zero = 0, rotation = "counter")
  attr(t, "class") <- attr(t, "circularp") <- NULL

  fit <- R_modereg_CircLin(y,x,t,kappa=bw[1],h=bw[2],maxit,tol)
  structure(list(datax = datax, datay = datay, x = tt, estim = fit,
                 bw = bw, n = n, kernel = c("vonmises","normal"), call = match.call(),
                 data.name = name, has.na = FALSE))


}

