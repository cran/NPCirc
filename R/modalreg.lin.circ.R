
modalreg.lin.circ<-function(x, y, t=NULL, bw=NULL, tol = 0.0001, maxit = 500, len=300){

  name <- deparse(substitute(x))
  datax <- x
  datay <- y
  if (!is.numeric(x)) stop("argument 'x' must be numeric")
  if (!is.numeric(y)) stop("argument 'y' must be numeric")
  if (length(x) != length(y)) stop("'x' and 'y' must have the same number of observations")
  if (is.circular(y)){
    datacircularp <- circularp(y)
  } else {
    datacircularp <- list(type = "angles", units = "radians", template = "none", modulo = "2pi", zero = 0, rotation = "counter")
  }
  dc <- list()
  dc$type <- datacircularp$type
  dc$units <- datacircularp$units
  dc$template <- datacircularp$template
  dc$modulo <- datacircularp$modulo
  dc$zero <- datacircularp$zero
  dc$rotation <- datacircularp$rotation
  if (is.null(t)){
    t <- seq(min(x), max(x), length=len)
  }else{
    if (!is.numeric(t)) stop("argument 't' must be numeric")
    t.na <- is.na(t)
    t <- t[!t.na]
    if (sum(t.na)>0) warning("'t' contains missing values. They were removed")
  }
  if (!is.numeric(tol))
    stop("argument 'tol' must be numeric")
  if (tol<=0)
    stop("argument 'tol' must be greater than 0")
  if (maxit<=0)
    stop("argument 'maxit' must be greater than 0")
  y <- conversion.circular(y, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
  nax <- is.na(x)
  nay <- is.na(y)
  x<-x[!nax & !nay]
  y<-y[!nax & !nay]
  if ((sum(nax)+sum(nay))>0) warning("Missing values were removed.", "\n")
  n <- length(x)
  if (n==0) stop("No observations (at least after removing missing values)")
  if(is.null(bw)){
    bw <- bw.modalreg.lin.circ(x, y, lower = NULL, upper = NULL, tol = 0.00001)
  }else{
    if(!is.numeric(bw) | length(bw)!=2)
      stop("'bw' must be a vector of length two")
    if(any(bw<0))
      stop("The components of 'bw' must be positive.")
  }

  attr(y, "class") <- attr(y, "circularp") <- NULL

  fit <- R_modereg_LinCirc(circular(y,units="radians"),x,t,bw[1],bw[2],maxit,tol)
  structure(list(datax = datax, datay = datay, x = t, estim = fit,
                 bw = bw, n = n, kernel = c("normal","vonmises"), call = match.call(),
                 data.name = name, has.na = FALSE))


}
