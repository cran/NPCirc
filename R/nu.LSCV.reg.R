nu.LSCV.reg<-function(x,y,method="LL",lower=0,upper=100,tol=0.1){
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	if (!is.numeric(y)) stop("argument 'y' must be numeric")
	if (length(x) != length(y)) stop("'x' and 'y' must have the same number of observations")
	nax <- is.na(x)
	nay <- is.na(y)
	x<-x[!nax & !nay]
	y<-y[!nax & !nay]
	if ((sum(nax)+sum(nay))>0) warning("Missing values were removed.", "\n")
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
	if (any(x<0) | any(x>2*pi)) stop("The sample of angles 'x' must be in radians between 0 and 2*pi")
	if (!is.numeric(upper)){ 
		warning("argument 'upper' must be numeric. Default upper boundary was used")
		upper=100
	}
	if (!is.numeric(lower)){
		warning("argument 'lower' must be numeric. Default lower boundary was used")
		lower=0
	}
	if (lower<0 | lower>=upper){
      	warning("The boundaries must be positive and 'lower' must be smaller that 'upper'. Default boundaries were used")
		upper=100
		lower=0
	}
	if (!is.numeric(tol)) stop("argument 'tol' must be numeric")
	lscv<-function(x,nu){
		error<-numeric(n)
		for (j in 1:n){
			error[j]<-(y[j]-kern.reg.circ(x[-j],y[-j],x[j],nu,method=method))^2
		}
		return(mean(error))
	}
	nu.lscv <- optimize(function(h)lscv(x,h),interval=c(lower,upper),tol=tol)$minimum
 	if (nu.lscv < lower + tol | nu.lscv > upper - tol) 
      warning("minimum occurred at one end of the range")
	return(nu.lscv)
}
