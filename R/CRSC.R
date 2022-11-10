

aux_CRSC<-function(x,y,t,bw,family,p,startv,tol,maxit){

  x_t <- outer(x, t, "-")
  s <- sin(x_t)
  kxt <- exp(bw * cos(x_t))/(besselI(bw,0)*2*pi)

  if(family=="gaussian"){

    y_hat<-loc_loglik(x,y,x,bw,family,p,startv,tol,maxit)[[1]]
    yY2<-(y-y_hat)^2
    y_kxt<-yY2*kxt
    num<-apply(y_kxt,2,sum)


  }else{

    est<-loc_loglik(x,y,t,bw,family,p,startv,tol,maxit)

    if(family=="poisson"){

      e_est_1<-exp(est[[1]])

      if(p==1){
        ts <- t(est[[1]]+est[[2]]*t(s))
      }else if(p==3){
        ts <- t(est[[1]]+est[[2]]*t(s)+est[[3]]*t(s^2)+est[[4]]*t(s^3))
      }

      ets <- exp(ts)
      zZ2<-t(t(y-ets)/e_est_1)^2


    }else if(family=="bernoulli"){

      e_est_1<-exp(est[[1]])/(1+exp(est[[1]]))

      if(p==1){
        ts <- t(est[[1]]+est[[2]]*t(s))
      }else if(p==3){
        ts <- t(est[[1]]+est[[2]]*t(s)+est[[3]]*t(s^2)+est[[4]]*t(s^3))
      }

      ets1 <- exp(ts)/(1+exp(ts))
      zZ2<-t(t(y-ets1)*(e_est_1/(1+e_est_1)^2))^2

    }else{

      e_est_1<-exp(est[[1]])

      if(p==1){
        ts <- t(est[[1]]+est[[2]]*t(s))
      }else if(p==3){
        ts <- t(est[[1]]+est[[2]]*t(s)+est[[3]]*t(s^2)+est[[4]]*t(s^3))
      }

      ets <- 1*y*exp(-ts)
      zZ2<-t(t(ets/1-1))^2

    }

    z_kxt<- zZ2*kxt
    num<-apply(z_kxt,2,sum)

  }

  if(p==1){

    sn0 <- apply(kxt , 2, sum)
    sn1 <- apply(kxt * s, 2, sum)
    sn2 <- apply(kxt * s^2, 2, sum)
    gn0 <- apply(kxt^2 , 2, sum)
    gn1 <- apply(kxt^2 * s, 2, sum)
    gn2 <- apply(kxt^2 * s^2, 2, sum)

    bjti <-  t(kxt) * (t(s)*sn0 - sn1)
    determinante <- apply(bjti, 1, sum)
    denom <- sn0 - (sn2*gn0 + sn0*gn2 - 2*sn1*gn1)/(sn0*sn2-sn1^2)
    sigma2<-num/denom
    eSnGnSne <- (sn2^2*gn0-2*sn1*sn2*gn1+sn1^2*gn2)/((sn0*sn2-sn1^2)^2)

  }else if(p==3){

    sn0 <- apply(kxt , 2, sum)
    sn1 <- apply(kxt * s, 2, sum)
    sn2 <- apply(kxt * s^2, 2, sum)
    sn3 <- apply(kxt * s^3, 2, sum)
    sn4 <- apply(kxt * s^4, 2, sum)
    sn5 <- apply(kxt * s^5, 2, sum)
    sn6 <- apply(kxt * s^6, 2, sum)
    gn0 <- apply(kxt^2 , 2, sum)
    gn1 <- apply(kxt^2 * s, 2, sum)
    gn2 <- apply(kxt^2 * s^2, 2, sum)
    gn3 <- apply(kxt^2 * s^3, 2, sum)
    gn4 <- apply(kxt^2 * s^4, 2, sum)
    gn5 <- apply(kxt^2 * s^5, 2, sum)
    gn6 <- apply(kxt^2 * s^6, 2, sum)

    A11 <- sn2*sn4*sn6 +2*sn3*sn4*sn5 - sn4^3 - sn3^2*sn6 -sn2*sn5^2
    A12 <- -sn1*sn4*sn6 -sn2*sn4*sn5 -sn3^2*sn5 + sn3*sn4^2 +sn2*sn3*sn6 + sn1*sn5^2
    A13 <- sn1*sn3*sn6 + sn2*sn4^2 + sn2*sn3*sn5 - sn3^2*sn4 - sn2^2*sn6 - sn1*sn4*sn5
    A14 <- -sn1*sn3*sn5 - 2*sn2*sn3*sn4+sn3^3 +sn2^2*sn5 +sn1*sn4^2
    bjti <-  t(kxt) * (A11 + t(s) * A12 + t(s^2) * A13  + t(s^3) * A14)
    determinante <- apply(bjti, 1, sum)

    A21 <- - sn1*sn4*sn6 -sn3^2*sn5 -sn4*sn2*sn5 +sn4^2*sn3 +sn3*sn2*sn6 +sn1*sn5^2
    A22 <- sn0*sn4*sn6 +sn2*sn5*sn3 + sn3*sn5^2 -sn3^2*sn4 - sn2^2*sn6 -sn0*sn5^2
    A23 <- -sn0*sn3*sn6 -sn2*sn4*sn3 - sn3*sn1*sn5 + sn3^3 +sn2*sn1*sn6 +sn0*sn4*sn5
    A24 <- sn0*sn3*sn5 + sn2^2*sn4 +sn3*sn1*sn4 - sn3^2*sn2 -sn2*sn1*sn5 - sn0*sn4^2

    A31 <- sn1*sn3*sn6 + sn2*sn5*sn3 + sn4^2*sn2 - sn4*sn3^2 - sn2^2*sn6 - sn1*sn5*sn4
    A32 <- -sn0*sn3*sn6 - sn1*sn5*sn3 - sn3*sn2*sn4 +sn3^3 + sn1*sn2*sn6 + sn0*sn5*sn4
    A33 <- sn0*sn2*sn6 + 2*sn1*sn3*sn4 - sn3^2*sn2 - sn1^2*sn6 - sn0*sn4^2
    A34 <- -sn0*sn2*sn5 - sn1*sn4*sn2 - sn3^2*sn1 + sn3*sn2^2 + sn1^2*sn5 + sn0*sn3*sn4

    A41 <- -sn1*sn3*sn5 - 2*sn2*sn3*sn4 + sn3^3 + sn2^2*sn5 + sn1*sn4^2
    A42 <- sn0*sn3*sn5 + sn1*sn4*sn3 + sn2^2*sn4 - sn2*sn3^2 -sn1*sn2*sn5 -sn0*sn4^2
    A43 <- -sn0*sn2*sn5 - sn1*sn3^2 - sn2*sn1*sn4 + sn2^2*sn3 + sn1^2*sn5 +sn0*sn3*sn4
    A44 <- sn0*sn2*sn4 + 2*sn1*sn2*sn3 - sn2^3 - sn1^2*sn4 - sn0*sn3^2

    traceSnGn_num <- (A11*gn0 + A12*gn1 + A13*gn2 + A14 *gn3 +
                        A21*gn1 + A22*gn2 + A23*gn3 + A24*gn4 +
                        A31*gn2 + A32*gn3 + A33*gn4 + A34*gn5 +
                        A41*gn3 + A42*gn4 + A43*gn5 + A44*gn6)

    traceSnGn <- traceSnGn_num / (determinante)
    denom <- sn0 - traceSnGn
    sigma2<-num/denom

    b11 <- A11*gn0+A12*gn1+A13*gn2+A14*gn3
    b12 <- A11*gn1+A12*gn2+A13*gn3+A14*gn4
    b13 <- A11*gn2+A12*gn3+A13*gn4+A14*gn5
    b14 <- A11*gn3+A12*gn4+A13*gn5+A14*gn6
    eSnGnSne <- (b11*A11+b12*A21 + b13*A31 + b14*A41)/(determinante^2)

  }

  result <- sigma2*(1+(p+1)*eSnGnSne)


  int<-Bolstad2::sintegral(t,result)$int
  return(int)

}



CRSC<-function(x,y,t,family,p,startv,lower,upper,tol,maxit){

  bw_select<-optimize(function(bw)aux_CRSC(x,y,t,bw,family,p,startv,tol,maxit),lower=lower,upper=upper)$minimum

  return(bw_select)

}
