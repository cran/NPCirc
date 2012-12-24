dwsn<- function(x,xi,eta,lambda,l=20){
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	if (!is.numeric(xi)) stop("argument 'xi' must be numeric")
	if (!is.numeric(eta)) stop("argument 'eta' must be numeric")
	if (!is.numeric(lambda)) stop("argument 'lambda' must be numeric")
	x <- x[!is.na(x)]
	n <- length(x)
	if (sum(is.na(x))>0) warning("Missing values were removed")
	if (n==0) stop("No observations (at least after removing missing values)")
	if (!is.numeric(l)){
		warning("argument 'l' must be numeric. Default value of 'l' was used")	
		l=20
	}
	if (l<=0){
		warning("argument 'l' must be a positive integer. Default value of 'l' was used")
		l=20
	}
	fx<-numeric(n)
	for (i in 1:n){
		suma<-0
		for(r in -l:l){
			val<-(x[i]+2*pi*r-xi)/eta
			suma<-suma+dnorm(val)*pnorm(lambda*val)
		}
		fx[i]<-2/eta*suma
	}
	return(fx)
}


