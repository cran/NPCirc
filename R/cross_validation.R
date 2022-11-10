
aux_cross_validation<-function(x,y,bw,family,startv,tol,maxit){

  n<-length(x)
  error <- numeric(n)
  for (j in 1:n) {

    fitj<-loc_loglik(x[-j],y[-j],x[j],bw,family,p=1,startv,tol,maxit)[[1]]

    if(family=="gaussian"){
      error[j] <- (y[j] - RegCircLin(x[-j], y[-j], x[j], bw, method = "LL", tol = 0.01))^2 #!!!!!!!!!!!!!!!!!!!!!!
    }else if(family=="bernoulli"){
      error[j] <- (y[j]*fitj - log(1+exp(fitj)))
    }else if(family=="poisson"){
      error[j] <- y[j]*fitj - exp(fitj)
    }else if(family=="gamma"){
      error[j] <- (y[j] - exp(fitj))^2
    }

  }

  me<-mean(error)

  return(me)
}



cross_validation<-function(x,y,family,startv,lower,upper,tol,maxit){

  if(family=="bernoulli"){

    bw_seq<-seq(lower,upper,length=200)
    k_cv<-sapply(1:length(bw_seq),function(i){aux_cross_validation(x,y,bw_seq[i],family,startv,tol,maxit)})
    bw_select<-bw_seq[which.max(k_cv)]

  }else if(family=="poisson"){

    bw_select<-optimize(function(bw) aux_cross_validation(x,y,bw,family,startv,tol,maxit),lower=lower,upper=upper,maximum=TRUE)$maximum

  }else{

    bw_select<-optimize(function(bw) aux_cross_validation(x,y,bw,family,startv,tol,maxit),lower=lower,upper=upper)$minimum

  }

  return(bw_select)

}
