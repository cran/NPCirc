
bw.modalreg.circ.lin<-function(x, y, lower = NULL, upper = NULL, maxit = 500, tol = 0.00001){

  if (!is.numeric(x)) stop("argument 'x' must be numeric")
  if (!is.numeric(y)) stop("argument 'y' must be numeric")
  if (length(x) != length(y)) stop("'x' and 'y' must have the same number of observations")
  x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (!is.numeric(x)) stop("argument 'x' must be numeric")
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

  range2<-(max(y)-min(y))/2

  if(is.null(lower)){
    lower = c(0.05, 0.05)
  }else{
    lower[1] = 1/lower[1]^2
  }
  if(is.null(upper)){
    upper = c(3,range2)
  }else{
    upper[1] = 1/upper[1]^2
  }

  e<-R_EM_vMN_select(x,y,10,0.0001,1000,15,50)
  param<-optim(c(0.5,1),MISECL,method="L-BFGS-B",mu=as.numeric(e[[1]][[4]]),kappa=as.numeric(e[[1]][[5]]),m=as.numeric(e[[1]][[6]]),s=as.numeric(e[[1]][[7]]),pip=as.numeric(e[[1]][[3]]),n=n,lower=lower,upper=upper)$par
  kappa<-1/param[1]^2
  hseq<-seq(param[2]*0.75,range2,length=30)
  cv<-R_CV_modereg_CircLin2(y, x, kappa, hseq, maxit, tol)
  h<-hseq[which.min(cv)]
  kappaseq<-seq(0.5,50,length=30)
  cv2<-R_re_CV_modereg_CircLin2(y, x, kappaseq, h, maxit, tol)
  kappa<-kappaseq[which.min(cv2)]

  return(c(kappa,h))

}
