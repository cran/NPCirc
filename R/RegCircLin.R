RegCircLin<-function(x,y,t,bw,method,tol){
	L<-kernCL(x, t, bw, method = method , tol = tol)
	fhat<-as.vector(L%*%y)
	return(fhat)
}
