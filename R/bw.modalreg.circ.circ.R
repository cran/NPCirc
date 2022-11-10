
bw.modalreg.circ.circ<-function(x, y, lower = NULL, upper = NULL, maxit = 500, tol = 0.00001){

  if (!is.numeric(x)) stop("argument 'x' must be numeric")
  if (!is.numeric(y)) stop("argument 'y' must be numeric")
  if (length(x) != length(y)) stop("'x' and 'y' must have the same number of observations")
  x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  y <- conversion.circular(y, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
  attr(y, "class") <- attr(y, "circularp") <- NULL
  nax <- is.na(x)
  nay <- is.na(y)
  x<-x[!nax & !nay]
  y<-y[!nax & !nay]
  if ((sum(nax)+sum(nay))>0) warning("Missing values were removed.", "\n")
  n <- length(x)
  if (n==0) stop("No observations (at least after removing missing values)")
  if (!is.numeric(upper) & !is.null(upper)){
    warning("argument 'upper' must be numeric. Default upper boundaries were used")
    upper <- NULL
  }else if(is.numeric(upper) & !is.null(upper)){
    if(length(upper)!=2){
      warning("argument 'upper' must have length 2. Default upper boundaries were used")
      upper <- NULL
    }else if(any(upper<=0) ){
      warning("The boundaries must be positive. Default boundaries were used")
      upper <- NULL
    }
  }
  if(!is.numeric(lower) & !is.null(lower)){
    warning("argument 'lower' must be numeric. Default lower boundaries were used")
    lower <- NULL
  }else if(is.numeric(lower) & !is.null(lower)){
    if(length(lower)!=2){
      warning("argument 'lower' must have length 2. Default lower boundaries were used")
      upper <- NULL
    }else if(any(lower<=0)){
      warning("The boundaries must be positive. Default boundaries were used")
      lower <- NULL
    }
  }
  if(any(lower>=upper)){
    warning("'lower' must be smaller that 'upper'. Default boundaries were used")
    upper <- NULL
    lower <- NULL
  }
  if (!is.numeric(tol))
    stop("argument 'tol' must be numeric")
  if (tol<=0)
    stop("argument 'tol' must be greater than 0")
  if (maxit<=0)
    stop("argument 'maxit' must be greater than 0")



  if(is.null(lower)){
    lower_ast <- c(0.135, 0.135)
    lower <- c(1/2^2,1/2^2)
  }else{
    lower_ast <- c(1/sqrt(lower[2]),1/sqrt(lower[1]))
  }
  if(is.null(upper)){
    upper_ast <- c(2,2)
    upper <- c(1/0.135^2,1/0.135^2)
  }else{
    upper_ast <- c(1/sqrt(upper[2]),1/sqrt(upper[1]))
  }

  e<-R_EM_vMvM_select(y,x,10,0.0001,1000,15,50)
  param<-optim(c(0.5,0.5),MISECC,method="L-BFGS-B",mu=as.numeric(e[[1]][[4]]),kappa=as.numeric(e[[1]][[5]]),m=as.numeric(e[[1]][[6]]),nu=as.numeric(e[[1]][[7]]),pip=as.numeric(e[[1]][[3]]),n=n,lower=lower_ast,upper=upper_ast)$par
  kappa<-1/param[1]^2
  nu<-1/param[2]^2
  kappaseq<-seq(lower[2],kappa/2,length=30)
  cv1<-R_CV_modereg_CircCirc(circular(y,units="radians"), x, nu, kappaseq, 1000, 0.00001)
  kappa<-kappaseq[which.min(cv1)]
  nuseq<-seq(lower[1],upper[1],length=30)
  cv2<- R_cv_re_modereg_CircCirc(circular(y,units="radians"),x,nuseq,kappa,1000,0.00001)
  nu<-nuseq[which.min(cv2)]


  return(c(nu,kappa))

}
