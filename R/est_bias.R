
est_bias<-function(x,y,t,bw,bw_ast,family,startv,tol,maxit){

  x_t <- outer(x, t, "-")
  s <- sin(x_t)
  kxt <- exp(bw * cos(x_t))

  estp3<-loc_loglik(x,y,t,bw_ast,family,p=3,startv,tol,maxit)
  r<-estp3[[3]]*t(s^2)+estp3[[4]]*t(s^3)

  if(family=="gaussian"){

    sn1 <- apply(kxt * s, 2, sum)
    sn2 <- apply(kxt * s^2, 2, sum)
    bjti0 <- t(kxt) * (sn2 - t(s) * sn1)
    L_r<- bjti0/apply(bjti0, 1, sum)
    resa<-L_r%*%t(r)
    res<-diag(resa)

  }else{

    k<-length(t)
    beta0<-rep(startv[1],k)
    beta1<-rep(startv[2],k)

    ts <- t(beta0+beta1*t(s))

    if(family=="poisson"){
      ets <- exp(ts)
      K <- y-ets
    }else if(family=="bernoulli"){
      ets1 <- exp(ts)/(1+exp(ts))
      ets <- ets1/(1+exp(ts))
      K <- y-ets1
    }else{
      ets <- 1*y*exp(-ts)
      K <- ets-1
    }

    bn0 <- apply(kxt * ets, 2, sum)
    bn1 <- apply(kxt * ets * s, 2, sum)
    bn2 <- apply(kxt * ets * s^2, 2, sum)
    bjti0 <-  t(kxt)*(bn2 - t(s) * bn1)
    bjti1 <-  t(kxt)*(t(s)*bn0 - bn1)

    beta0<-beta0 + colSums(t(bjti0) * K)/(bn2*bn0-bn1^2)
    beta1<-beta1 + colSums(t(bjti1) * K)/(bn2*bn0-bn1^2)

    res<-5000
    b<-0
    while(res>tol){

      beta0_old<-beta0
      beta1_old<-beta1

      ts <- t(beta0+beta1*t(s))

      if(family=="poisson"){
        ets <- exp(ts)
        K <- y-ets
      }else if(family=="bernoulli"){
        ets1 <- exp(ts)/(1+exp(ts))
        ets <- ets1/(1+exp(ts))
        K <- y-ets1
      }else{
        ets <- 1*y*exp(-ts)
        K <- ets-1
      }

      bn0 <- apply(kxt * ets, 2, sum)
      bn1 <- apply(kxt * ets * s, 2, sum)
      bn2 <- apply(kxt * ets * s^2, 2, sum)
      bjti0 <- t(kxt)*(bn2 - t(s) * bn1)
      bjti1 <- t(kxt)*(t(s)*bn0 - bn1)

      beta0<-beta0 + colSums(t(bjti0) * K)/(bn2*bn0-bn1^2)
      beta1<-beta1 + colSums(t(bjti1) * K)/(bn2*bn0-bn1^2)

      res<-max(c(abs(beta0-beta0_old),abs(beta1-beta1_old)))
    }

    ts <- t(beta0+beta1*t(s)+r)

    if(family=="poisson"){
      ets <- exp(ts)
      K <- y-ets
    }else if(family=="bernoulli"){
      ets1 <- exp(ts)/(1+exp(ts))
      ets <- ets1/(1+exp(ts))
      K <- y-ets1
    }else{
      ets <- 1*y*exp(-ts)
      K <- ets-1
    }

    bn0 <- apply(kxt * ets, 2, sum)
    bn1 <- apply(kxt * ets * s, 2, sum)
    bn2 <- apply(kxt * ets * s^2, 2, sum)
    bjti0 <- t(kxt)*(bn2 - t(s) * bn1)

    res <- - colSums(t(bjti0) * K)/(bn2*bn0-bn1^2)

  }

  return(res)

}
