circ.local.lik<-function(x, y, t = NULL, bw = NULL, family, p = 1, startv = NULL, tol = 0.00001, maxit = 300, from = circular(0),
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
  if (p!=1 & p!=3)
    stop("argument 'p' must be either 1 or 3")
  if (family!="gaussian" & family!="poisson" & family!="bernoulli" & family!="gamma")
    stop("argument 'family' must be one of the following characters: 'gaussian', 'poisson', 'bernoulli' or 'gamma'" )
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
  if (is.null(startv)){
    if(family == "bernoulli"){
      startv <- c(log(mean(y)/(1-mean(y))),0,0,0)
    }else if(family == "poisson" | family == "gamma"){
      startv<-c(log(mean(y)),0,0,0)
    }
  }else{
    if (!is.numeric(startv))
      stop("argument 'startv' must be numeric")
  }
  if(is.null(bw)){
    if(p == 1){
      bw<-refined_rule(x,y,t,family=family,startv=startv,lower=0,upper=50,lower_ast=0,upper_ast=15,tol=tol,maxit=maxit)
      warning("the value of 'bw' was computed by using the refined rule")
    }else{
      bw<-CRSC(x,y,t,family=family,p=p,startv=startv,lower=0,upper=15,tol=tol,maxit=maxit)
      warning("the value of 'bw' was computed by using the CRSC rule")
    }
  }else{
    if (!is.numeric(bw))
      stop("argument 'bw' must be numeric")
    if (bw<=0)
      stop("argument 'bw' must be greater than 0")
  }
  attr(x, "class") <- attr(x, "circularp") <- NULL
  tt <- conversion.circular(t, dc$units, dc$type, dc$template,
                            dc$modulo, dc$zero, dc$rotation)
  t <- conversion.circular(t, units = "radians", modulo = "2pi",
                           zero = 0, rotation = "counter")
  attr(t, "class") <- attr(t, "circularp") <- NULL

  ghat <- loc_loglik(x = x, y = y, t = t, bw = bw, family = family, p = p, startv = startv, tol = tol, maxit = maxit)
  structure(list(datax = datax, datay = datay, x = tt, estim = ghat,
                 bw = bw, p = p, n = n, kernel = "vonmises", call = match.call(),
                 data.name = name, has.na = FALSE))

}
