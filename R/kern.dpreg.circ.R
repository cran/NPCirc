kern.dpreg.circ<-function(x, y, t = NULL, bw, startvmu = NULL, startvgam = NULL, tol= 0.000001, maxit = 300, from = circular(0),
                       to = circular(2 * pi), len = 250){
  name <- deparse(substitute(x))
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
  if (is.null(startvmu)){
    startvmu<-c(log(mean(y)),0)
  }else if(!is.numeric(startvmu)){
    warning("'startvmu' must be numeric. Default values were used")
    startvmu<-c(log(mean(y)),0)
  }else if(length(startvmu)!=2){
    warning("'startvmu' must have length two. Default values were used")
    startvmu<-c(log(mean(y)),0)
  }
  if (is.null(startvgam)){
    startvgam<-c(0,0)
  }else if(!is.numeric(startvgam)){
    warning("'startvgam' must be numeric. Default values were used")
    startvgam<-c(0,0)
  }else if(length(startvgam)!=2){
    warning("'startvgam' must have length two. Default values were used")
    startvgam<-c(0,0)
  }
  if (!is.numeric(tol))
    stop("argument 'tol' must be numeric")
  if (tol <= 0)
    stop("argument 'tol' must be greater than 0")
  if (!is.numeric(maxit))
    stop("argument 'maxit' must be numeric")
  if (maxit <= 0)
    stop("argument 'maxit' must be greater than 0")
  if (!is.numeric(bw))
      stop("argument 'bw' must be numeric")
  if (any(bw<=0))
      stop("argument 'bw' must have components greater than 0")
  if(length(bw)!=2)
      stop("argument 'bw' must have length 2")

  attr(x, "class") <- attr(x, "circularp") <- NULL
  tt <- conversion.circular(t, dc$units, dc$type, dc$template,
                            dc$modulo, dc$zero, dc$rotation)
  t <- conversion.circular(t, units = "radians", modulo = "2pi",
                           zero = 0, rotation = "counter")
  attr(t, "class") <- attr(t, "circularp") <- NULL
  
  est<-my_fun_DoublePois(t,x,y,startvmu,startvgam,bw[1],bw[2],tol,maxit)
  structure(list(datax = datax, datay = datay, x = tt, estim = est,
                 bw = bw, n = n, kernel = c("vonmises","vonmises"), call = match.call(),
                 data.name = name, has.na = FALSE))
  
}