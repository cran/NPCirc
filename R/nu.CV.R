nu.CV<-function(x,method="LCV",lower=0,upper=100,tol=0.1,np=500){
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
	if (!is.numeric(tol)) stop("argument 'tol' must be numeric")
	if (method=="LSCV"){
		if (np<=0){ 
			warning("'np' must be positive. Default value of 'np' was used")
			np=500
		}
		t<-seq(0,2*pi,length=np)
		lscv<-function(x,nu){
			est<-kern.den.circ(x,t,nu)
			n<-length(x)
			cv<-numeric(n)
			for (j in 1:n){
			cv[j]<-kern.den.circ(x[-j],x[j],nu)
			}
			return(int.Simp(t,est^2)-2*mean(cv))	
		}
		nu <- optimize(function(h)lscv(x,h),interval=c(lower,upper),tol=tol)$minimum
	}else{
		cv <- function(x,nu){ 
			n<-length(x)
			logdens<-numeric(n)
			for (j in 1:n){
				logdens[j]<-log(kern.den.circ(x[-j],x[j],nu))
			}
			return(sum(logdens))
		}
		nu <- optimize(function(h)cv(x,h),interval=c(lower,upper),tol=tol,maximum=TRUE)$maximum
	}
	if (nu < lower + tol | nu > upper - tol) 
      warning("minimum/maximum occurred at one end of the range")
	return(nu)
}

