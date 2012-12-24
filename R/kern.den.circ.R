kern.den.circ<-function(x,t=NULL,nu,from=0,to=2*pi,len=250){
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	x <- x[!is.na(x)]
	n <- length(x)
	if (sum(is.na(x))>0) warning("Missing values were removed")
	if (n==0) stop("No observations (at least after removing missing values)")
	if (any(x<0) | any(x>2*pi)) stop("The sample of angles 'x' must be in radians between 0 and 2*pi")
	if (!is.numeric(nu)) stop("argument 'nu' must be numeric")
	if (nu<0) stop("argument 'nu' must be integer and positive")
	if (!is.numeric(from)) stop("argument 'from' must be numeric")
	if (!is.numeric(to)) stop("argument 'to' must be numeric")
	if (!is.finite(from)) stop("non-finite 'from'")
	if (!is.finite(to)) stop("non-finite 'to'")
	if (from>to) stop("argument 'from' must be smaller than argument 'to'")
	if (!is.numeric(len)) stop("argument 'len' must be numeric")
	if (len <= 0) stop("argument 'len' must be integer and positive")
	if (is.null(t)){
		t<-seq(from=from,to=to,length=len)
	}else{
		if (any(t<0) | any(t>2*pi)){
			if (!is.numeric(t)) stop("argument 't' must be numeric")
			warning("The estimator was computed for values of 't' between 0 and 2*pi")
			t<-t[t>0 & t<2*pi]
		}
	}
	exptx<-exp(cos(outer(t,x,"-"))*nu)
	return(rowMeans(exptx)/(2*pi*besselI(nu,0)))
}