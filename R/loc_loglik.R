
loc_loglik<-function(x,y,t,bw,family,p,startv,tol,maxit){

  x_t <- outer(x, t, "-")
  kxt <- exp(bw * cos(x_t))
  s <- sin(x_t)

  if(p==1){

    if(family=="gaussian"){

      sn1 <- apply(kxt * s, 2, sum)
      sn2 <- apply(kxt * s^2, 2, sum)
      bjti0 <- t(kxt) * (sn2 - t(s) * sn1)
      L0 <- bjti0/apply(bjti0, 1, sum)
      beta0<-as.vector(L0 %*% y)

      sn0 <- apply(kxt, 2, sum)
      bjti1 <- t(kxt) * (t(s)*sn0 - sn1)
      L1 <- bjti1/apply(bjti0, 1, sum)
      beta1<-as.vector(L1 %*% y)

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
        K <- y-ets
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
      while(res>tol & b<maxit){

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
        bjti0 <-  t(kxt)*(bn2 - t(s) * bn1)
        bjti1 <-  t(kxt)*(t(s)*bn0 - bn1)

        beta0<-beta0 + colSums(t(bjti0) * K)/(bn2*bn0-bn1^2)
        beta1<-beta1 + colSums(t(bjti1) * K)/(bn2*bn0-bn1^2)

        indexes<-c(is.nan(beta0),is.nan(beta1))
        if(any(indexes))
          stop("convergence not reached, try other values of the smoothing parameter")

        res<-max(c(abs(beta0-beta0_old),abs(beta1-beta1_old)))
        b<-b+1
      }
      if (b==maxit)
        stop("convergence not reached, try other values of the smoothing parameter")
    }
    return(list(beta0,beta1))

  }else{

    if(family=="gaussian"){

      sn0 <- apply(kxt, 2, sum)
      sn1 <- apply(kxt * s, 2, sum)
      sn2 <- apply(kxt * s^2, 2, sum)
      sn3 <- apply(kxt * s^3, 2, sum)
      sn4 <- apply(kxt * s^4, 2, sum)
      sn5 <- apply(kxt * s^5, 2, sum)
      sn6 <- apply(kxt * s^6, 2, sum)

      A11 <- sn2*sn4*sn6 +2*sn3*sn4*sn5 - sn4^3 - sn3^2*sn6 -sn2*sn5^2
      A12 <- -sn1*sn4*sn6 -sn2*sn4*sn5 -sn3^2*sn5 + sn3*sn4^2 +sn2*sn3*sn6 + sn1*sn5^2
      A13 <- sn1*sn3*sn6 + sn2*sn4^2 + sn2*sn3*sn5 - sn3^2*sn4 - sn2^2*sn6 - sn1*sn4*sn5
      A14 <- -sn1*sn3*sn5 - 2*sn2*sn3*sn4+sn3^3 +sn2^2*sn5 +sn1*sn4^2
      bjti0 <- t(kxt) * (A11 + t(s) * A12 + t(s^2) * A13  + t(s^3) * A14)
      L0 <- bjti0/apply(bjti0, 1, sum)

      A21 <- - sn1*sn4*sn6 -sn3^2*sn5 -sn4*sn2*sn5 +sn4^2*sn3 +sn3*sn2*sn6 +sn1*sn5^2
      A22 <- sn0*sn4*sn6 + 2*sn2*sn5*sn3  -sn3^2*sn4 - sn2^2*sn6 -sn0*sn5^2
      A23 <- -sn0*sn3*sn6 -sn2*sn4*sn3 - sn3*sn1*sn5 + sn3^3 +sn2*sn1*sn6 +sn0*sn4*sn5
      A24 <- sn0*sn3*sn5 + sn2^2*sn4 +sn3*sn1*sn4 - sn3^2*sn2 -sn2*sn1*sn5 - sn0*sn4^2
      bjti1 <- t(kxt) * (A21 + t(s) * A22 + t(s^2) * A23  + t(s^3) * A24)
      L1 <- bjti1/apply(bjti0, 1, sum)

      A31 <- sn1*sn3*sn6 + sn2*sn5*sn3 + sn4^2*sn2 - sn4*sn3^2 - sn2^2*sn6 - sn1*sn5*sn4
      A32 <- -sn0*sn3*sn6 - sn1*sn5*sn3 - sn3*sn2*sn4 +sn3^3 + sn1*sn2*sn6 + sn0*sn5*sn4
      A33 <- sn0*sn2*sn6 + 2*sn1*sn3*sn4 - sn3^2*sn2 - sn1^2*sn6 - sn0*sn4^2
      A34 <- -sn0*sn2*sn5 - sn1*sn4*sn2 - sn3^2*sn1 + sn3*sn2^2 + sn1^2*sn5 + sn0*sn3*sn4
      bjti2 <- t(kxt) * (A31 + t(s) * A32 + t(s^2) * A33  + t(s^3) * A34)
      L2 <- bjti2/apply(bjti0, 1, sum)

      A41 <- -sn1*sn3*sn5 - 2*sn2*sn3*sn4 + sn3^3 + sn2^2*sn5 + sn1*sn4^2
      A42 <- sn0*sn3*sn5 + sn1*sn4*sn3 + sn2^2*sn4 - sn2*sn3^2 -sn1*sn2*sn5 -sn0*sn4^2
      A43 <- -sn0*sn2*sn5 - sn1*sn3^2 - sn2*sn1*sn4 + sn2^2*sn3 + sn1^2*sn5 +sn0*sn3*sn4
      A44 <- sn0*sn2*sn4 + 2*sn1*sn2*sn3 - sn2^3 - sn1^2*sn4 - sn0*sn3^2
      bjti3 <- t(kxt) * (A41 + t(s) * A42 + t(s^2) * A43  + t(s^3) * A44)
      L3 <- bjti3/apply(bjti0, 1, sum)

      beta0<-as.vector(L0%*%y)
      beta1<-as.vector(L1%*%y)
      beta2<-as.vector(L2%*%y)
      beta3<-as.vector(L3%*%y)

    }else{
      k<-length(t)
      beta0<-rep(startv[1],k)
      beta1<-rep(startv[2],k)
      beta2<-rep(startv[3],k)
      beta3<-rep(startv[4],k)

      beta0_old<-beta0
      beta1_old<-beta1
      beta2_old<-beta2
      beta3_old<-beta3

      ts <- t(beta0+beta1*t(s)+beta2*t(s^2)+beta3*t(s^3))

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
      bn3 <- apply(kxt * ets * s^3, 2, sum)
      bn4 <- apply(kxt * ets * s^4, 2, sum)
      bn5 <- apply(kxt * ets * s^5, 2, sum)
      bn6 <- apply(kxt * ets * s^6, 2, sum)

      A11 <- bn2*bn4*bn6 +2*bn3*bn4*bn5 - bn4^3 - bn3^2*bn6 -bn2*bn5^2
      A12 <- -bn1*bn4*bn6 -bn2*bn4*bn5 -bn3^2*bn5 + bn3*bn4^2 +bn2*bn3*bn6 + bn1*bn5^2
      A13 <- bn1*bn3*bn6 + bn2*bn4^2 + bn2*bn3*bn5 - bn3^2*bn4 - bn2^2*bn6 - bn1*bn4*bn5
      A14 <- -bn1*bn3*bn5 - 2*bn2*bn3*bn4+bn3^3 +bn2^2*bn5 +bn1*bn4^2
      bjti0 <-   t(kxt) * (A11 + t(s) * A12 + t(s^2) * A13  + t(s^3) * A14)

      A21 <- - bn1*bn4*bn6 -bn3^2*bn5 -bn4*bn2*bn5 +bn4^2*bn3 +bn3*bn2*bn6 +bn1*bn5^2
      A22 <- bn0*bn4*bn6 + 2*bn2*bn5*bn3 + -bn3^2*bn4 - bn2^2*bn6 -bn0*bn5^2
      A23 <- -bn0*bn3*bn6 -bn2*bn4*bn3 - bn3*bn1*bn5 + bn3^3 +bn2*bn1*bn6 +bn0*bn4*bn5
      A24 <- bn0*bn3*bn5 + bn2^2*bn4 +bn3*bn1*bn4 - bn3^2*bn2 -bn2*bn1*bn5 - bn0*bn4^2
      bjti1 <- t(kxt) * (A21 + t(s) * A22 + t(s^2) * A23  + t(s^3) * A24)

      A31 <- bn1*bn3*bn6 + bn2*bn5*bn3 + bn4^2*bn2 - bn4*bn3^2 - bn2^2*bn6 - bn1*bn5*bn4
      A32 <- -bn0*bn3*bn6 - bn1*bn5*bn3 - bn3*bn2*bn4 +bn3^3 + bn1*bn2*bn6 + bn0*bn5*bn4
      A33 <- bn0*bn2*bn6 + 2*bn1*bn3*bn4 - bn3^2*bn2 - bn1^2*bn6 - bn0*bn4^2
      A34 <- -bn0*bn2*bn5 - bn1*bn4*bn2 - bn3^2*bn1 + bn3*bn2^2 + bn1^2*bn5 + bn0*bn3*bn4
      bjti2 <- t(kxt) * (A31 + t(s) * A32 + t(s^2) * A33  + t(s^3) * A34)

      A41 <- -bn1*bn3*bn5 - 2*bn2*bn3*bn4 + bn3^3 + bn2^2*bn5 + bn1*bn4^2
      A42 <- bn0*bn3*bn5 + bn1*bn4*bn3 + bn2^2*bn4 - bn2*bn3^2 -bn1*bn2*bn5 -bn0*bn4^2
      A43 <- -bn0*bn2*bn5 - bn1*bn3^2 - bn2*bn1*bn4 + bn2^2*bn3 + bn1^2*bn5 +bn0*bn3*bn4
      A44 <- bn0*bn2*bn4 + 2*bn1*bn2*bn3 - bn2^3 - bn1^2*bn4 - bn0*bn3^2
      bjti3 <- t(kxt) * (A41 + t(s) * A42 + t(s^2) * A43  + t(s^3) * A44)

      determinante<-bn3^4 - (3* bn2 * bn4 + 2 *bn1* bn5 + bn0 *bn6)* bn3^2 + 2 *((bn2^2 + bn0 *bn4)* bn5 + bn1 *(bn4^2 + bn2 *bn6))* bn3 - bn0* bn4^3 + bn2^2* bn4^2 + bn1^2 *bn5^2 - bn2^3* bn6 - bn1^2* bn4* bn6 + bn2* (bn0* (bn4* bn6 - bn5^2) - 2* bn1 *bn4* bn5)

      beta0<-beta0 + colSums(t(bjti0) * K)/(determinante)
      beta1<-beta1 + colSums(t(bjti1) * K)/(determinante)
      beta2<-beta2 + colSums(t(bjti2) * K)/(determinante)
      beta3<-beta3 + colSums(t(bjti3) * K)/(determinante)

      b<-0
      res<-5000
      while(res>tol & b<maxit){
        beta0_old<-beta0
        beta1_old<-beta1
        beta2_old<-beta2
        beta3_old<-beta3

        ts <- t(beta0+beta1*t(s)+beta2*t(s^2)+beta3*t(s^3))

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
        bn3 <- apply(kxt * ets * s^3, 2, sum)
        bn4 <- apply(kxt * ets * s^4, 2, sum)
        bn5 <- apply(kxt * ets * s^5, 2, sum)
        bn6 <- apply(kxt * ets * s^6, 2, sum)

        A11 <- bn2*bn4*bn6 +2*bn3*bn4*bn5 - bn4^3 - bn3^2*bn6 -bn2*bn5^2
        A12 <- -bn1*bn4*bn6 -bn2*bn4*bn5 -bn3^2*bn5 + bn3*bn4^2 +bn2*bn3*bn6 + bn1*bn5^2
        A13 <- bn1*bn3*bn6 + bn2*bn4^2 + bn2*bn3*bn5 - bn3^2*bn4 - bn2^2*bn6 - bn1*bn4*bn5
        A14 <- -bn1*bn3*bn5 - 2*bn2*bn3*bn4 + bn3^3 +bn2^2*bn5 +bn1*bn4^2
        bjti0 <-   t(kxt) * (A11 + t(s) * A12 + t(s^2) * A13  + t(s^3) * A14)

        A21 <- - bn1*bn4*bn6 -bn3^2*bn5 -bn4*bn2*bn5 +bn4^2*bn3 +bn3*bn2*bn6 +bn1*bn5^2
        A22 <- bn0*bn4*bn6 + 2*bn2*bn5*bn3 + -bn3^2*bn4 - bn2^2*bn6 -bn0*bn5^2
        A23 <- -bn0*bn3*bn6 -bn2*bn4*bn3 - bn3*bn1*bn5 + bn3^3 +bn2*bn1*bn6 +bn0*bn4*bn5
        A24 <- bn0*bn3*bn5 + bn2^2*bn4 +bn3*bn1*bn4 - bn3^2*bn2 -bn2*bn1*bn5 - bn0*bn4^2
        bjti1 <- t(kxt) * (A21 + t(s)*A22 + t(s^2)*A23  + t(s^3)*A24)

        A31 <- bn1*bn3*bn6 + bn2*bn5*bn3 + bn4^2*bn2 - bn4*bn3^2 - bn2^2*bn6 - bn1*bn5*bn4
        A32 <- -bn0*bn3*bn6 - bn1*bn5*bn3 - bn3*bn2*bn4 +bn3^3 + bn1*bn2*bn6 + bn0*bn5*bn4
        A33 <- bn0*bn2*bn6 + 2*bn1*bn3*bn4 - bn3^2*bn2 - bn1^2*bn6 - bn0*bn4^2
        A34 <- -bn0*bn2*bn5 - bn1*bn4*bn2 - bn3^2*bn1 + bn3*bn2^2 + bn1^2*bn5 + bn0*bn3*bn4
        bjti2 <- t(kxt) * (A31 + t(s) * A32 + t(s^2) * A33  + t(s^3) * A34)

        A41 <- -bn1*bn3*bn5 - 2*bn2*bn3*bn4 + bn3^3 + bn2^2*bn5 + bn1*bn4^2
        A42 <- bn0*bn3*bn5 + bn1*bn4*bn3 + bn2^2*bn4 - bn2*bn3^2 -bn1*bn2*bn5 -bn0*bn4^2
        A43 <- -bn0*bn2*bn5 - bn1*bn3^2 - bn2*bn1*bn4 + bn2^2*bn3 + bn1^2*bn5 +bn0*bn3*bn4
        A44 <- bn0*bn2*bn4 + 2*bn1*bn2*bn3 - bn2^3 - bn1^2*bn4 - bn0*bn3^2
        bjti3 <- t(kxt) * (A41 + t(s) * A42 + t(s^2) * A43  + t(s^3) * A44)

        determinante<-bn3^4 - (3* bn2 * bn4 + 2 *bn1* bn5 + bn0 *bn6)* bn3^2 + 2 *((bn2^2 + bn0 *bn4)* bn5 + bn1 *(bn4^2 + bn2 *bn6))* bn3 - bn0* bn4^3 + bn2^2* bn4^2 + bn1^2 *bn5^2 - bn2^3* bn6 - bn1^2* bn4* bn6 + bn2* (bn0* (bn4* bn6 - bn5^2) - 2* bn1 *bn4* bn5)

        beta0<-beta0 + colSums(t(bjti0) * K)/(determinante)
        beta1<-beta1 + colSums(t(bjti1) * K)/(determinante)
        beta2<-beta2 + colSums(t(bjti2) * K)/(determinante)
        beta3<-beta3 + colSums(t(bjti3) * K)/(determinante)

        indexes<-c(is.nan(beta0),is.nan(beta1),is.nan(beta2),is.nan(beta3))
        if(any(indexes))
          stop("convergence not reached, try other values of the smoothing parameter")

        res<-max(c(abs(beta0-beta0_old),abs(beta1-beta1_old),abs(beta2-beta2_old),abs(beta3-beta3_old)))
        b<-b+1
      }
      if (b==maxit)
        stop("convergence not reached, try other values of the smoothing parameter")
    }

    return(list(beta0,beta1,beta2,beta3))

  }

}
