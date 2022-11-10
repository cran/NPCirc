
est_variance<-function(x,y,t,bw,bw_ast,family,startv,tol,maxit){

  x_t <- outer(x, t, "-")
  kxt <- exp(bw * cos(x_t))
  kxt2<-kxt^2
  s <- sin(x_t)

  if(family=="gaussian"){
    y_hat_p3<-loc_loglik(x,y,x,bw_ast,family,p=3,startv,tol,maxit)[[1]]
    yY2<-(y-y_hat_p3)^2
    x_t <- outer(x, t, "-")
    kxt_ast <- exp(bw_ast * cos(x_t))
    y_kxt_ast<- yY2*kxt_ast
    sigma2<-apply(y_kxt_ast,2,sum)/apply(kxt_ast,2,sum)

    sn0 <- apply(kxt , 2, sum)
    sn1 <- apply(kxt * s, 2, sum)
    sn2 <- apply(kxt * s^2, 2, sum)
    gn0 <- apply(kxt2 , 2, sum)
    gn1 <- apply(kxt2 * s, 2, sum)
    gn2 <- apply(kxt2 * s^2, 2, sum)
    eSnGnSne<-(sn2^2*gn0-2*sn1*sn2*gn1+sn1^2*gn2)/(sn0^2*sn2^2-2*sn1^2*sn0*sn2+sn1^4)

  }else{

    if(family=="gamma"){

      ggp3<-loc_loglik(x,y,x,bw_ast,family,p=3,startv,tol,maxit)[[1]]
      e_gY_1_2<-(y*exp(-ggp3)-1)^2


      kxt_ast <- exp(bw_ast * cos(x_t))
      y_kxt_ast<- e_gY_1_2*kxt_ast
      sigma2<-apply(y_kxt_ast,2,sum)/apply(kxt_ast,2,sum)

    }

    estp1<-loc_loglik(x,y,t,bw,family,p=1,startv,tol,maxit)
    ts <- t(estp1[[1]]+estp1[[2]]*t(s))

    if(family=="poisson"){
      ets <- exp(ts)
      K <- y-ets
      sigma2 <- exp(estp1[[1]])
    }else if(family=="bernoulli"){
      ets1 <- exp(ts)/(1+exp(ts))
      ets <- ets1/(1+exp(ts))
      K <- y-ets1
      sigma2 <- exp(estp1[[1]])/(1+exp(estp1[[1]]))^2
    }else{
      ets <- 1*y*exp(-ts)
      K <- ets-1
    }

    bn0 <- apply(kxt * ets, 2, sum)
    bn1 <- apply(kxt * ets * s, 2, sum)
    bn2 <- apply(kxt * ets * s^2, 2, sum)
    gn0 <- apply(kxt2 , 2, sum)
    gn1 <- apply(kxt2 * s, 2, sum)
    gn2 <- apply(kxt2 * s^2, 2, sum)
    eSnGnSne<-(bn2^2*gn0-2*bn1*bn2*gn1+bn1^2*gn2)/(bn0^2*bn2^2-2*bn1^2*bn0*bn2+bn1^4)

  }

  return(sigma2*eSnGnSne)

}
