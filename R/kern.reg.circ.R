kern.reg.circ<-function(x,y,t=NULL,nu,method="LL",tol=300,from=0,to=2*pi,len=250){
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
		if (!is.numeric(t)) stop("argument 't' must be numeric")
		if (any(t<0) | any(t>2*pi)){
			warning("The estimator was computed for values of 't' between 0 and 2*pi")
			t<-t[t>0 & t<2*pi]
		}
	}
	if (method=="NW"){
		t_x<-outer(t,x,"-")
		if (max(nu)>tol){
			z<-apply(nu*cos(t_x),1,max)-tol
			B<-exp(nu*cos(t_x)-z) 
		}else {
			B<-exp(nu*cos(t_x))
		}
		L<-B/apply(B,1,sum)
	}else{
		x_t<-outer(x,t,"-")
		if (max(nu)>tol){
			z<-apply(nu*cos(x_t),2,max)-tol
			kxt<-exp(t(t(nu*cos(x_t))-z)) 
		}else {
			kxt<-exp(nu*cos(x_t))
		}
		m<-sin(x_t)
		sn1<-apply(kxt*m,2,sum)
		sn2<-apply(kxt*m^2,2,sum)
		bjti<-t(kxt)*(sn2-t(m)*sn1)
		L<-bjti/apply(bjti,1,sum)
	}
	return(L%*%y)
}