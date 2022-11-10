
pois_p1<-function(xseq,x,y,kappa,startv){
  
  k<-length(xseq)
  beta0<-rep(startv[1],k)
  beta1<-rep(startv[2],k)
  
  x_t <- outer(x, xseq, "-")
  kxt <- exp(kappa * cos(x_t))/(besselI(kappa,0)*2*pi)
  s <- sin(x_t)
  ts <- t(beta0+beta1*t(s))
  ets <- exp(ts)
  bn0 <- apply(kxt * ets, 2, sum)
  bn1 <- apply(kxt * ets * s, 2, sum)
  bn2 <- apply(kxt * ets * s^2, 2, sum)
  bjti0 <-  t(kxt)*(bn2 - t(s) * bn1)
  
  
  bjti1 <-  t(kxt)*(t(s)*bn0 - bn1)
  
  K <- (y-ets)
  
  
  beta0<-beta0 + colSums(t(bjti0) * K)/(bn2*bn0-bn1^2)
  beta1<-beta1 + colSums(t(bjti1) * K)/(bn2*bn0-bn1^2)
  
  res<-5000
  tol<-0.001
  b<-0
  while(res>tol){
    
    
    beta0_old<-beta0
    beta1_old<-beta1
    
    ts <- t(beta0+beta1*t(s))
    ets <- exp(ts)
    K <- (y-ets)
    bn0 <- apply(kxt * ets, 2, sum)
    bn1 <- apply(kxt * ets * s, 2, sum)
    bn2 <- apply(kxt * ets * s^2, 2, sum)
    bjti0 <- t(kxt)*(bn2 - t(s) * bn1)
    bjti1 <- t(kxt)*(t(s)*bn0 - bn1)
    
    
    beta0<-beta0 + colSums(t(bjti0) * K)/(bn2*bn0-bn1^2)
    beta1<-beta1 + colSums(t(bjti1) * K)/(bn2*bn0-bn1^2)
    
    res<-max(c(abs(beta0-beta0_old),abs(beta1-beta1_old)))
    b<-b+1

  }
  return(list(beta0,beta1))
}



loglik_cv <- function(x,y,kappa,startv) {	
  n<-length(x)
  error <- numeric(n)
  for (j in 1:n) {
    
    fitj<-pois_p1(x[j],x[-j],y[-j],kappa,startv)[[1]]
    error[j] <- y[j]*fitj - exp(fitj)
  }
  me<-mean(error)
  
  return(me)
  
  
}