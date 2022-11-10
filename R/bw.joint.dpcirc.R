bw.joint.dpcirc<-function(x, y, startvmu = NULL, startvgam = NULL, lower=c(0.05,0.05), upper=c(50,7), tol = 0.00001, maxit = 300){
  
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same number of observations")
  x <- conversion.circular(x, units = "radians", zero = 0,
                           rotation = "counter", modulo = "2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  nax <- is.na(x)
  nay <- is.na(y)
  x <- x[!nax & !nay]
  y <- y[!nax & !nay]
  if ((sum(nax) + sum(nay)) > 0)
    warning("Missing values were removed.", "\n")
  n <- length(x)
  if (n == 0)
    stop("No observations (at least after removing missing values)")
  if (!is.numeric(upper) ){
    warning("argument 'upper' must be numeric. Default upper boundaries were used")
    upper <- c(50,7)
  }else{
    if(length(upper)!=2){
      warning("argument 'upper' must have length 2. Default upper boundaries were used")
      upper <- c(50,7)
    }else if(any(upper<=0) ){
      warning("The boundaries must be positive. Default boundaries were used")
      upper <- c(50,7)
    }
  }
  if(!is.numeric(lower)){
    warning("argument 'lower' must be numeric. Default lower boundaries were used")
    lower <- c(0.05,0.05)
  }else{
    if(length(lower)!=2){
      warning("argument 'lower' must have length 2. Default lower boundaries were used")
      upper <- c(50,7)
    }else if(any(lower<=0)){
      warning("The boundaries must be positive. Default boundaries were used")
      lower <- c(0.05,0.05)
    }
  }
  if(any(lower>=upper)){
    warning("'lower' must be smaller that 'upper'. Default boundaries were used")
    upper <- c(50,7)
    lower <- c(0.05,0.05)
  }
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
 
  kappa_cv<-optimize(function(kappa)loglik_cv(x,y,kappa,startvmu),lower=lower[1],upper=upper[1],maximum=TRUE)$maximum
  nucv_seq<-seq(lower[2],upper[2],length=50)
  n_cv<-sapply(1:length(nucv_seq),function(i){my_fun_loglik_cv_nu(x,y,startvmu,startvgam,kappa_cv,nucv_seq[i])})
  nu_cv<-nucv_seq[which.max(n_cv)]	 
  
  sm_param<-c(kappa_cv,nu_cv)
  return(sm_param)
  
}

