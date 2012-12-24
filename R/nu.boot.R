nu.boot <- function(x,lower=0,upper=100,np=500,tol=0.1){
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	x <- x[!is.na(x)]
	n <- length(x)
	if (sum(is.na(x))>0) warning("Missing values were removed")
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
	if (!is.numeric(np)) stop("argument 'np' must be numeric")
	if (np<=0){ 
		warning("'np' must be positive. Default value of 'np' was used")
		np=500
	}
	t <- seq(0,2*pi,length=np)
	if (!is.numeric(tol)) stop("argument 'tol' must be numeric")
	MISEboot<-function(t,x,nu){
		n <- length(x)
		costx <- cos(outer(t,x,"-"))
		exptx <- exp(costx*nu)
		fhat <- apply(exptx,1,sum)/(n*2*pi*besselI(nu,0))
		EBf <- apply((besselI(abs(2*nu*cos(outer(t,x,"-")/2)),0)/besselI(nu,0)),1,sum)/(2*pi*n*besselI(nu,0))
		EBff2 <-(apply(besselI(nu*sqrt(5+4*costx),0),1,sum)/besselI(nu,0))/((2*pi*n*besselI(nu,0))^2)+(EBf-fhat)^2 - EBf^2/n
		return(int.Simp(t,EBff2))
	}
	nuboot<-optimize(function(h)MISEboot(t,x,h),interval=c(lower,upper),tol=tol)$minimum
	return(nuboot)
}
