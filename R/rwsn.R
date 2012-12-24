rwsn<-function(n,xi,eta,lambda){
	if (!is.numeric(n)) stop("argument 'n' must be numeric")
	if (!is.numeric(xi)) stop("argument 'xi' must be numeric")
	if (!is.numeric(eta)) stop("argument 'eta' must be numeric")
	if (!is.numeric(lambda)) stop("argument 'lambda' must be numeric")
	modulo<-function(x,m){
		t1<-floor(x/m)
		return(x-t1*m)
	}
	mas<-modulo(rsn(n,xi,eta,lambda),2*pi)
	return(mas)
}