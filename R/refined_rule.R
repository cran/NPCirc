
aux_refined<-function(x,y,t,bw,family,startv,lower_ast,upper_ast,tol,maxit){

  bw_ast <- CRSC(x,y,t,family,p=3,startv,lower_ast,upper_ast,tol,maxit)
  bw_ast <- bw_ast * (192/167)^(2/9)

  bias<-est_bias(x,y,t,bw,bw_ast,family,startv,tol,maxit)
  variance<-est_variance(x,y,t,bw,bw_ast,family,startv,tol,maxit)
  MSE<-bias^2+variance

  MISE<-Bolstad2::sintegral(t,MSE)$int

  return(MISE)
}






refined_rule<-function(x,y,t,family,startv,lower=0,upper=50,lower_ast=0,upper_ast=15,tol,maxit){

  if(family=="bernoulli"){

    bw_seq<-seq(lower,upper,length=200)
    k_cv<-sapply(1:length(bw_seq),function(i){aux_refined(x,y,t,bw_seq[i],family,startv,lower_ast,upper_ast,tol,maxit)})
    bw_select<-bw_seq[which.min(k_cv)]

  }else{

    bw_select<-optimize(function(bw)aux_refined(x,y,t,bw,family,startv,lower_ast,upper_ast,tol,maxit),lower=lower,upper=upper)$minimum

  }

  return(bw_select)

}
