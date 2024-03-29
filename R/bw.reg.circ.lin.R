bw.reg.circ.lin<-function(x,y,method="LL",lower=0,upper=50,tol=1e-2){
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
	if (!is.numeric(upper)){
		warning("argument 'upper' must be numeric. Default upper boundary was used")
		upper <- 50
	}
	if (!is.numeric(lower)){
		warning("argument 'lower' must be numeric. Default lower boundary was used")
		lower <- 0
	}
	if (lower<0 | lower>=upper){
      	warning("The boundaries must be positive and 'lower' must be smaller that 'upper'. Default boundaries were used")
		upper <- 50
		lower <- 0
	}
	if (!is.numeric(tol)) stop("argument 'tol' must be numeric")

	lscv<-function(x,bw){
		error<-numeric(n)
		for (j in 1:n){
			error[j]<-(y[j]-RegCircLin(x[-j],y[-j],x[j],bw,method=method,tol=tol))^2
		}
		return(mean(error))
	}
	bw.lscv <- optimize(function(h)lscv(x,h),interval=c(lower,upper),tol=tol)$minimum
 	if (bw.lscv < lower + tol | bw.lscv > upper - tol)
      warning("minimum occurred at one end of the range")
	return(bw.lscv)
}
