# Internal functions PDP


R_myfun_DPoisM<-function(x,y,startv,gammai){

  n<-length(x)
  aux<-1:((length(startv)-1)/2)
  aux2<-outer(x,aux,"*")
  T<-cbind(rep(1,n),cos(aux2),sin(aux2))


  W<-diag(exp(as.numeric(t(T%*%startv)))/gammai)
  F<-t(T)%*%W%*%T
  F_1<-solve(F)


  v<-t(T)%*%(as.numeric((y-exp(t(T%*%startv)))/gammai))

  beta<-startv+F_1%*%v
  return(beta)
}



d_fun<-function(y,g){

  res<-2*(y*log(y)-y-y*g+exp(g))
  ind<-y==0
  res[ind]<-2*g[ind]

  return(res)
}


R_myfun_DPoisG<-function(x,y,startv,gi){

  n<-length(x)
  aux<-1:((length(startv)-1)/2)
  aux2<-outer(x,aux,"*")
  T<-cbind(rep(1,n),cos(aux2),sin(aux2))

  W<-diag(0.5,n)
  F<-t(T)%*%W%*%T
  F_1<-solve(F)

  dd<-d_fun(y,gi)
  ets<-exp(t(T%*%startv))


  v<- -0.5 * t(T)%*%as.numeric((1-dd/ets))


  beta<-startv+F_1%*%v
  return(beta)
}




Param_DP<-function(x,y,dim_mu,dim_gam,startvmu,startvgam,tol,maxit){


  a1<-startvmu
  a2<-startvgam
  n<-length(x)
  gammai<-rep(1,n)

  aux_1<-1:((dim_mu-1)/2)
  aux2_1<-outer(x,aux_1,"*")
  T_1<-cbind(rep(1,n),cos(aux2_1),sin(aux2_1))

  aux_2<-1:((dim_gam-1)/2)
  aux2_2<-outer(x,aux_2,"*")
  T_2<-cbind(rep(1,n),cos(aux2_2),sin(aux2_2))



  res<-10000
  it<-0
  while((res>=1+tol | res<=1-tol )& it<maxit){

    a1_old<-a1
    gammai<-as.numeric(exp(t(T_2%*%a2)))
    startv<-a1
    a1<-R_myfun_DPoisM(x,y,startv,gammai)

    gi<-t(T_1%*%a1)
    a2_old<-a2
    startv2<-a2
    a2<-R_myfun_DPoisG(x,y,startv2,gi)

    res<-max(c(abs(exp(a2)/exp(a2_old)),abs(exp(a1)/exp(a1_old))))
    it<-it+1

  }

  return(list(beta_mu=a1,beta_gam=a2,iterations=it,res=res))
}
