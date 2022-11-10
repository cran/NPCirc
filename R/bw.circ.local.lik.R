
bw.circ.local.lik<-function(x, y, t = NULL, rule = NULL, p, family, startv = NULL, lower = 0, upper = 50, lower_ast = 0, upper_ast = 15, tol = 0.00001, maxit = 300, from = circular(0), to = circular(2 * pi), len = 250){

  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same number of observations")
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
  if (!is.numeric(upper)) {
    warning("argument 'upper' must be numeric. Default upper boundary was used")
    upper <- 50
  }
  if (!is.numeric(lower)) {
    warning("argument 'lower' must be numeric. Default lower boundary was used")
    lower <- 0
  }
  if (lower < 0 | lower >= upper) {
    warning("The boundaries must be positive and 'lower' must be smaller that 'upper'. Default boundaries were used")
    upper <- 50
    lower <- 0
  }
  if (!is.numeric(upper_ast)) {
    warning("argument 'upper_ast' must be numeric. Default upper boundary was used")
    upper_ast <- 15
  }
  if (!is.numeric(lower_ast)) {
    warning("argument 'lower_ast' must be numeric. Default lower boundary was used")
    lower <- 0
  }
  if (lower_ast < 0 | lower_ast >= upper_ast) {
    warning("The boundaries must be positive and 'lower_ast' must be smaller that 'upper_ast'. Default boundaries were used")
    upper_ast <- 15
    lower_ast <- 0
  }
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
  if (p!=1 & p!=3)
    stop("argument 'p' must be either 1 or 3")
  if (family!="gaussian" & family!="poisson" & family!="bernoulli" & family!="gamma")
    stop("argument 'family' must be one of the following characters: 'gaussian', 'poisson', 'bernoulli' or 'gamma'" )
  if (tol <= 0)
    stop("argument 'tol' must be greater than 0")
  if (maxit <= 0)
    stop("argument 'maxit' must be greater than 0")
  if (p == 1){
    if(rule != "refined" & rule != "CRSC" & rule != "cv")
      stop("argument 'family' must be one of the following characters: 'refined', 'CRSC' or 'cv'" )
  }else{
    if(rule != "refined" & rule != "CRSC" & rule != "cv")
      stop("argument 'rule' must be one of the following characters: 'refined', 'CRSC' or 'cv'" )
    if(rule == "refined")
      stop("rule 'refined' is only available for p=1" )
  }
  if(rule == "refined"){
    bw <- refined_rule(x,y,t,family=family,startv=startv,lower=lower,upper=upper,lower_ast=lower_ast,upper_ast=upper_ast,tol=tol,maxit=maxit)
  }else if(rule == "CRSC"){
    bw <- CRSC(x,y,t,family=family,p=p,startv=startv,lower=lower,upper=upper,tol=tol,maxit=maxit)
  }else{
    bw <-cross_validation(x,y,family=family,startv=startv,lower=lower,upper=upper,tol=tol,maxit=maxit)
  }
  return(bw)
}
