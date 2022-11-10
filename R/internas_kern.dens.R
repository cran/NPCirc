
bw.dpi=function(data,deriv.order=0,nstage=2,kernel="vonmises",M=NULL,commonkappa=TRUE,approximate=TRUE,Q1=NULL,Q2=NULL){



  if (approximate != T & approximate != F) {
    warning("Argument 'approximate' must be T or F. Default value of 'approximate' was used")
    approximate=T
  }

  n=length(data)

  # Fitting with a mixture of von Mises (Step 1.a-b)
  mlemixVM=mle.mixvm(data,M,commonkappa)

  if(length(mlemixVM[[1]])==0){
    mlemixVM=mle.mixvm(data,1,commonkappa)
    warning("The smoothing parameter was computed by using a von Mises at stage 0 (M=1)")
  }

  # Rule of thumb at stage 0 (Step 1.c)
  intf2si=intmixddvonmises(mlemixVM[[1]],mlemixVM[[2]],mlemixVM[[3]],deriv.order+nstage+2)

  # Plug-in concentration parameter of phi (Steps 2-3)
  if(nstage>0){
    intf2si=(-1)^(deriv.order+nstage+2)*intf2si

    for(ni in nstage:1){
      bwsi=amsephi(2*(deriv.order+ni)+2,intf2si,n,kernel,Q1[2*(deriv.order+ni)+2],approximate)
      intf2si=sum(kern.den.circ(data,z=data,bw=bwsi,deriv.order=2*(deriv.order+ni)+2)$y)/n
    }

    intf2si=intf2si*(-1)^(deriv.order+2)
  }

  # Plug-in AMISE smoothing parameter (Step 4)
  bwpi=amisefr(2*deriv.order,intf2si,n,kernel,Q2,approximate)

  return(bwpi)

}





###################################################################
###################################################################




bw.ste=function(data,deriv.order=0,nstage=2,kernel="vonmises",M=NULL,commonkappa=TRUE,lower=10^(-10),upper=0.5,tol=.Machine$double.eps^0.25,Q1=NULL,Q2=NULL){

  if ((length(lower) != 1) | (lower < 0)) {
    "Argument 'lower' must be a positive real number. Default value of 'lower' was used"
    lower=10^(-10)
    if(kernel=="vonmises"){lower=10^(-3)}
  }

  if ((length(upper) != 1) | (upper < 0)) {
    "Argument 'upper' must be a positive real number. Default value of 'lower' was used"
    upper=0.5
    if(kernel=="vonmises"){upper=pi^2/3}
  }

  n=length(data)

  # Fitting with a mixture of von Mises (Step 1)
  mlemixVM=mle.mixvm(data,M,commonkappa)

  if(length(mlemixVM[[1]])==0){
    mlemixVM=mle.mixvm(data,1,commonkappa)
    warning("The smoothing parameter was computed by using a von Mises at stage 0 (M=1)")
  }

  # Rule of thumb to estimate the two needed functionals
  if(nstage==2){
    intf2s1=intmixddvonmises(mlemixVM[[1]],mlemixVM[[2]],mlemixVM[[3]],deriv.order+3)
    intf2s2=intmixddvonmises(mlemixVM[[1]],mlemixVM[[2]],mlemixVM[[3]],deriv.order+4)
    intf2s1=(-1)^(deriv.order+3)*intf2s1
    intf2s2=(-1)^(deriv.order+4)*intf2s2

  }

  # If one prefers to estimate with the rule of thumb other functional

  if(nstage>2){
    intf2si=intmixddvonmises(mlemixVM[[1]],mlemixVM[[2]],mlemixVM[[3]],deriv.order+nstage+2)
    intf2si=(-1)^(deriv.order+nstage+2)*intf2si
    intf2sin=intf2si

    for(ni in nstage:2){
      intf2sio=intf2sin
      bwsi=amsephi(2*(deriv.order+ni)+2,intf2si,n,kernel,Q1[2*(deriv.order+ni)+2])
      intf2si=sum(kern.den.circ(data,z=data,bw=bwsi,deriv.order=2*(deriv.order+ni)+2)$y)/n
      intf2sin=intf2si
    }

    intf2s2=intf2sio
    intf2s1=intf2sin

  }



  # Estimate the two functionals phi_{2r+4} and  phi_{2r+6}
  # needed to stablish the relation between the two smoothing parameters (Step 2)


  bws1=amsephi(2*deriv.order+4,intf2s1,n,kernel,Q1[2*deriv.order+4],approximate=TRUE)
  intf2s3=sum(kern.den.circ(data,z=data,bw=bws1,deriv.order=2*deriv.order+4)$y)/n
  bws2=amsephi(2*deriv.order+6,intf2s2,n,kernel,Q1[2*deriv.order+6],approximate=TRUE)
  intf2s4=sum(kern.den.circ(data,z=data,bw=bws2,deriv.order=2*deriv.order+6)$y)/n

  # Estimate the functional phi_{2r+4} with the smoothing parameter gamma(h) (Step 3)

  gammabwt=function(h){
    ht=gammabw(h,deriv.order,intf2s3,intf2s4,kernel,Q1[2*deriv.order+4],Q2)

    if(kernel=="vonmises"){kappat=1/ht}
    if(kernel=="wrappednormal"){
      if(ht<1){kappat=(1-ht)^(1/4)}else{kappat=0}
    }
    return(kappat)
  }

  intf2t=function(h){
    solt=try(sum(kern.den.circ(data,z=data,bw=gammabwt(h),deriv.order=2*deriv.order+4)$y)/n,silent=T)
    if(inherits(solt, "try-error")){solt=NA}
    return(solt)
  }


  lowero=lower
  uppero=upper
  while(is.na(intf2t(lower))){lower=min(lower/0.9,lower+0.01)}
  while((intf2t(upper)==0)|(is.na(intf2t(upper)))){upper=max(0.99*upper,upper-0.01)}
  if((lower>lowero)|(upper<uppero))
    warning(paste("modifying bw.ste() search interval to (",lower,",",upper,")",sep=""))

  # Equation that needs to be solved to find the optimal h (Step 4)

  hsolve=function(h) amisefr2(deriv.order,(-1)^(deriv.order)*intf2t(h),n,kernel="vonmises",Q2,approximate=TRUE)


  solrot=try(uniroot(function(h) hsolve(h)-h,lower=lower,upper=upper,tol=tol)$root,silent=T)

  if(inherits(solrot, "try-error")){
    hopt=1
  }else{
    hopt=solrot
  }

  if(kernel=="vonmises"){kappaopt=1/hopt}
  if(kernel=="wrappednormal"){
    if(hopt<1){kappaopt=(1-hopt)^(1/4)}else{kappaopt=0}
  }


  if((kernel!="vonmises")&(kernel!="wrappednormal")){

    Aj=function(kappa,j) cosmoments(kappa,j,kernel)

    hfunc=function(kappa,seqj=1:100){
      Ajtot=outer(kappa,seqj,Aj)
      hkappa=pi^2/3+4*colSums((-1)^(seqj)*t(Ajtot)/(seqj)^2)
      return(hkappa)
    }

    kappaopt=uniroot(function(kappa) hfunc(kappa)-hopt,lower=10^(-5),upper=1-10^(-5))$root

  }

  return(kappaopt)

}



###################################################################
###################################################################



### Auxiliary functions


# This function computes the MLE for a mixture of von Mises
# of M components (it choses the 'best' M using AIC)

# If commonkappa=TRUE, the concentration parameter is the same in all components

# kappamax returns the maximum value of the concentration parameter
# in one of the components (when considering M>1)

mle.mixvm=function(x,M=NULL,commonkappa=TRUE,kappamax=200){

  if(is.null(M)){
    M=1:5
  }




  meanf=numeric()
  kappaf=numeric()
  AICf=Inf
  pf=0

  for(indM in M){

    if(indM==1){

      mlevm=mle.vonmises(x)
      mean=mlevm$mu
      kappa=mlevm$kappa
      AIC <-dvonmises(x, circular(mean),kappa)
      AIC <- -2 * sum(log(AIC)) + 4*indM + 2*(indM-1)

      if(AIC<AICf){
        AICf=AIC
        meanf=mean
        kappaf=kappa
        pf=1
      }

    }else{

      z <- cbind(cos(x), sin(x))
      if(commonkappa==F){
        y2 <- try(movMF(z, indM, start = "S"), TRUE)
      }else{
        y2 <- try(movMF(z, indM, start = "S",kappa = list(common = TRUE)), TRUE)
      }
      norm <- mean <- kappa <- p <- numeric(indM)
      mu <- matrix(NA, indM, 2)
      AIC <- 0

      if(inherits(y2, "try-error")){
        for (i in 1:indM) {
          norm[i] <- sqrt(sum(y2$theta[i, ]^2))
          mu[i, ] <- y2$theta[i, ]/norm[i]
          mean[i] <- atan2(mu[i, 2], mu[i, 1])
          kappa[i] <- y2$theta[i, 1]/mu[i, 1]
          p[i] <- y2$alpha[i]
          AIC <- AIC + p[i] * dvonmises(x, circular(mean[i]),kappa[i])
        }
        AIC <- -2 * sum(log(AIC)) + 4*indM + 2*(indM-1)

        # The integral cannot be computed if the concentration parameter is too large
        # This is controled with the argument kappamax
        if(sum(kappa>kappamax)>0){AIC=Inf}
        if(AIC<AICf){
          AICf=AIC
          meanf=mean
          kappaf=kappa
          pf=p
        }

      }

    }
  }

  return(list(meanf,kappaf,pf,AICf))

}


###################################################################
###################################################################

# This function computes the integral of the (r)th derivative squared

intmixddvonmises=function(mu,kappa,p,r){
  funmixvmdd=function(x) dmixdvonmises(x,mu,kappa,p,r)^2
  valint=integrate(funmixvmdd,-pi,pi)$val
  return(valint)
}


###################################################################
###################################################################

# Density function of the derivative of the von Mises mixture

dmixdvonmises=function(x,mu,kappa,p,r){

  res=rep(0,length(x))
  for(k in 1:length(mu)){
    res=res+p[k]*ddvonmises(x,mu[k],kappa[k],r)
  }
  return(res)

}


# Density function of the derivative of the von Mises density

ddvonmises=function(x,mu,kappa,r){

  if(r==0){valdd=dvonmisesd0(x,mu,kappa)}
  if(r==1){valdd=dvonmisesd1(x,mu,kappa)}
  if(r==2){valdd=dvonmisesd2(x,mu,kappa)}
  if(r==3){valdd=dvonmisesd3(x,mu,kappa)}
  if(r==4){valdd=dvonmisesd4(x,mu,kappa)}
  if(r==5){valdd=dvonmisesd5(x,mu,kappa)}
  if(r==6){valdd=dvonmisesd6(x,mu,kappa)}
  if(r==8){valdd=dvonmisesd8(x,mu,kappa)}

  if(r==7|r>8){valdd=dcirculardr(x,mu,kappa,r,"vonmises")}

  return(valdd)
}


# Trigonometric moments of the von Mises/wrapped normal
# If one needs to do KDE with other kernels, they should be added here

cosmoments=function(kappa,j,density="vonmises"){

  if (density != "vonmises" & density != "wrappednormal") {
    stop("The employed kernel should be one of 'vonmises' or 'wrappednormal'")
  }

  if(density=="vonmises"){
    result=besselI(kappa, nu = j, expon.scaled = TRUE)/besselI(kappa,nu = 0,expon.scaled = TRUE)
  }
  if(density=="wrappednormal"){
    result=kappa^(j^2)
  }
  return(result)
}


# Generic derivative of any circular density
# (obtained from the Fourier series representation)

dcirculardr=function(x,mu,kappa,r,density="vonmises",tol=.Machine$double.eps^(0.6)){

  minusxmu=outer(x,mu,"-")
  seqj=1:10
  multterm=seqj^r*cosmoments(kappa,seqj,density)

  while(multterm[length(multterm)]>tol){
    seqj2=seqj[length(seqj)]+1:10
    seqj=c(seqj,seqj2)
    multterm=c(multterm,seqj2^r*cosmoments(kappa,seqj2,density))
  }

  if(r%%4==1|r%%4==2){
    multterm=-multterm
  }
  if(r==0){sumf=1}else{sumf=0}
  if(r%%2==0){
    for(j in seqj){
      sumf=sumf+2*multterm[j]*cos(j*minusxmu)
    }
  }else{
    for(j in seqj){
      sumf=sumf+2*multterm[j]*sin(j*minusxmu)
    }
  }

  return(sumf/(2*pi))

}


# Explicit expressions of the derivative (for efficiency and accuracy)

dvonmisesd0=function(x,mu,kappa){

  cosf=function(x,mu) cos(x-mu)

  cosxmu=outer(x,mu,cosf)

  val1mat=exp(kappa*cosxmu)
  val1=val1mat/(2*pi*besselI(kappa,0))
  return(val1)
}

dvonmisesd1=function(x,mu,kappa){

  cosf=function(x,mu) cos(x-mu)
  sinf=function(x,mu) sin(x-mu)

  cosxmu=outer(x,mu,cosf)
  sinxmu=outer(x,mu,sinf)

  val1mat=exp(kappa*cosxmu)*sinxmu
  val1=-kappa*val1mat/(2*pi*besselI(kappa,0))
  return(val1)
}


dvonmisesd2=function(x,mu,kappa){

  n=length(x)
  cosf=function(x,mu) cos(x-mu)
  sinf=function(x,mu) sin(x-mu)

  cosxmu=outer(x,mu,cosf)
  sinxmu=outer(x,mu,sinf)

  val1mat=exp(kappa*cosxmu)*(kappa*sinxmu^2-cosxmu)
  val1=kappa*val1mat/(2*pi*besselI(kappa,0))
  return(val1)
}

dvonmisesd3=function(x,mu,kappa){

  n=length(x)
  cosf=function(x,mu) cos(x-mu)
  sinf=function(x,mu) sin(x-mu)

  cosxmu=outer(x,mu,cosf)
  sinxmu=outer(x,mu,sinf)

  val1mat=-exp(kappa*cosxmu)*sinxmu*(kappa^2*sinxmu^2-3*kappa*cosxmu-1)
  val1=kappa*val1mat/(2*pi*besselI(kappa,0))
  return(val1)
}

dvonmisesd4=function(x,mu,kappa){

  cosf=function(x,mu) cos(x-mu)
  sinf=function(x,mu) sin(x-mu)


  cosxmu=outer(x,mu,cosf)
  sinxmu=outer(x,mu,sinf)
  ekcosxmu=exp(kappa*cosxmu)

  val1mat=kappa^3*sinxmu^4-6*kappa^2*sinxmu^2*cosxmu+3*kappa*cosxmu^2-
    4*kappa*sinxmu^2+cosxmu
  val1=(kappa*ekcosxmu*val1mat)/(2*pi*besselI(kappa,0))
  return(val1)
}


dvonmisesd5=function(x,mu,kappa){

  cosf=function(x,mu) cos(x-mu)
  sinf=function(x,mu) sin(x-mu)


  cosxmu=outer(x,mu,cosf)
  sinxmu=outer(x,mu,sinf)
  ekcosxmu=exp(kappa*cosxmu)
  sinxmu2=sinxmu^2
  sinxmu4=sinxmu^4
  cosxmu2=cosxmu^2

  val1mat=-10*kappa^2+kappa^4*sinxmu4+25*kappa^2*cosxmu2-
    5*kappa*cosxmu*(2*kappa^2*sinxmu2-3)+1

  val1=-(kappa*sinxmu*ekcosxmu*val1mat)/(2*pi*besselI(kappa,0))
  return(val1)
}

dvonmisesd6=function(x,mu,kappa){

  cosf=function(x,mu) cos(x-mu)
  sinf=function(x,mu) sin(x-mu)


  cosxmu=outer(x,mu,cosf)
  sinxmu=outer(x,mu,sinf)
  ekcosxmu=exp(kappa*cosxmu)
  sinxmu2=sinxmu^2
  sinxmu4=sinxmu^4
  cosxmu2=cosxmu^2
  cosxmu3=cosxmu^3

  val1mat=-15*kappa^2*cosxmu3+15*kappa*cosxmu2*(3*kappa^2*sinxmu2-1)+
    kappa*sinxmu2*(kappa^4*sinxmu4-20*kappa^2*sinxmu2+16)+cosxmu*
    (-15*kappa^4*sinxmu4+75*kappa^2*sinxmu2-1)

  val1=(kappa*ekcosxmu*val1mat)/(2*pi*besselI(kappa,0))
  return(val1)
}



dvonmisesd8=function(x,mu,kappa){

  cosf=function(x,mu) cos(x-mu)
  sinf=function(x,mu) sin(x-mu)


  cosxmu=outer(x,mu,cosf)
  sinxmu=outer(x,mu,sinf)
  ekcosxmu=exp(kappa*cosxmu)
  sinxmu2=sinxmu^2
  sinxmu4=sinxmu^4
  sinxmu6=sinxmu^6
  cosxmu2=cosxmu^2
  cosxmu3=cosxmu^3
  cosxmu4=cosxmu^4


  val1mat=105*kappa^3*cosxmu4-210*kappa^2*cosxmu3*(2*kappa^2*sinxmu2-1)+
    21*kappa*cosxmu2*(10*kappa^4*sinxmu4-60*kappa^2*sinxmu2+3)+
    kappa*sinxmu2*(kappa^6*sinxmu6-56*kappa^4*sinxmu4+336*kappa^2*sinxmu2-64)+
    cosxmu*(-28*kappa^6*sinxmu6+630*kappa^4*sinxmu4-756*kappa^2*sinxmu2+1)

  val1=(kappa*ekcosxmu*val1mat)/(2*pi*besselI(kappa,0))
  return(val1)
}














###################################################################
###################################################################

# This function computes the optimal concentration parameter
# for the density functional

amsephi=function(r,intf2,n,kernel="vonmises",Q1,approximate=TRUE){

  hopt=(-2*Q1/(n*intf2))^(2/(r+3))


  if(kernel=="vonmises"){kappaopt=1/hopt}
  if(kernel=="wrappednormal"){
    if(hopt<1){kappaopt=(1-hopt)^(1/4)}else{kappaopt=0}
  }

  if((kernel!="vonmises")&(kernel!="wrappednormal")){

    Aj=function(kappa,j) cosmoments(kappa,j,kernel)

    hfunc=function(kappa,seqj=1:100){
      Ajtot=outer(kappa,seqj,Aj)
      hkappa=pi^2/3+4*colSums((-1)^(seqj)*t(Ajtot)/(seqj)^2)
      return(hkappa)
    }

    kappaopt=uniroot(function(kappa) hfunc(kappa)-hopt,lower=10^(-5),upper=1-10^(-5))$root

  }

  # A slightly more accurate concentration parameter (see Section 10)

  if(approximate==FALSE){

    kappaini=kappaopt


    if(kernel=="vonmises"){

      Aj=function(kappa,j) cosmoments(kappa,j,"vonmises")

      hfunc=function(kappa,seqj=1:100){
        Ajtot=outer(kappa,seqj,Aj)
        hkappa=pi^2/3+4*colSums((-1)^(seqj)*t(Ajtot)/(seqj)^2)
        return(hkappa)
      }

      maxseqj=c(10*(1:10),100*(2:10))
      maxseqj=maxseqj[which((Aj(2*kappaini,maxseqj)/Aj(2*kappaini,1))<10^(-10))[1]]
      seqj=1:maxseqj


      amsef=function(kappa) (ddvonmises(0,0,kappa,r)/n+hfunc(kappa,seqj)*intf2/2)^2
      solopt=try(optim(kappaini,amsef,lower=0,method="L-BFGS-B"),silent=T)
    }else{

      Aj=function(kappa,j) cosmoments(kappa,j,kernel)

      hfunc=function(kappa,seqj=1:100){
        Ajtot=outer(kappa,seqj,Aj)
        hkappa=pi^2/3+4*colSums((-1)^(seqj)*t(Ajtot)/(seqj)^2)
        return(hkappa)
      }

      maxseqj=c(10*(1:10),100*(2:10))
      maxseqj=maxseqj[which((Aj(10^(log10(kappaini)/10),maxseqj)/Aj(10^(log10(kappaini)/10),1))<10^(-10))[1]]
      seqj=1:maxseqj

      amsef=function(kappa) (dcirculardr(0,0,kappa,r,density=kernel)/n+hfunc(kappa,seqj)*intf2/2)^2
      solopt=try(optim(kappaini,amsef,lower=0.001,upper=max(0.99,10^(log10(kappaini)/100)),method="L-BFGS-B"),silent=T)
    }

    if(inherits(solopt, "try-error")){
      kappaopt=kappaini
    }else{
      kappaopt=solopt$par
    }

  }


  return(kappaopt)
}







###################################################################
###################################################################

# This function computes the optimal concentration parameter
# for the density derivative

amisefr=function(r,intf2,n,kernel="vonmises",Q2,approximate=TRUE){


  hopt=((2*r+1)*Q2/(n*intf2))^(2/(2*r+5))
  if(kernel=="vonmises"){kappaopt=1/hopt}
  if(kernel=="wrappednormal"){
    if(hopt<1){kappaopt=(1-hopt)^(1/4)}
    else{kappaopt=0}
  }

  if((kernel!="vonmises")&(kernel!="wrappednormal")){

    Aj=function(kappa,j) cosmoments(kappa,j,kernel)

    hfunc=function(kappa,seqj=1:100){
      Ajtot=outer(kappa,seqj,Aj)
      hkappa=pi^2/3+4*colSums((-1)^(seqj)*t(Ajtot)/(seqj)^2)
      return(hkappa)
    }

    kappaopt=uniroot(function(kappa) hfunc(kappa)-hopt,lower=10^(-5),upper=1-10^(-5))$root

  }

  # A slightly more accurate concentration parameter (see Remark 1)

  if(approximate==FALSE){

    kappaini=kappaopt


    if(kernel=="vonmises"){

      Aj=function(kappa,j) cosmoments(kappa,j,"vonmises")

      hfunc=function(kappa,seqj=1:100){
        Ajtot=outer(kappa,seqj,Aj)
        hkappa=pi^2/3+4*colSums((-1)^(seqj)*t(Ajtot)/(seqj)^2)
        return(hkappa)
      }

      maxseqj=c(10*(1:10),100*(2:10))
      maxseqj=maxseqj[which((Aj(2*kappaini,maxseqj)/Aj(2*kappaini,1))<10^(-10))[1]]
      seqj=1:maxseqj

      R2=function(r,kappa,seqj=1:100){

        if(r==0){
          val=besselI(2*kappa, nu = 0, expon.scaled = TRUE)/((besselI(kappa,nu = 0,expon.scaled = TRUE))^2)
          val=val/(2*pi)
        }else{
          Ajtot=outer(kappa,seqj,Aj)
          val=colSums(t(Ajtot^2)*(seqj)^(2*r))/pi
        }

        return(val)
      }

      amisef=function(kappa) hfunc(kappa,seqj)^2*intf2/4+R2(r,kappa,seqj)/n
      solopt=try(optim(kappaini,amisef,lower=0,method="L-BFGS-B"),silent=T)
    }else{

      Aj=function(kappa,j) cosmoments(kappa,j,kernel)

      hfunc=function(kappa,seqj=1:100){
        Ajtot=outer(kappa,seqj,Aj)
        hkappa=pi^2/3+4*colSums((-1)^(seqj)*t(Ajtot)/(seqj)^2)
        return(hkappa)
      }

      maxseqj=c(10*(1:10),100*(2:10))
      maxseqj=maxseqj[which((Aj(10^(log10(kappaini)/10),maxseqj)/Aj(10^(log10(kappaini)/10),1))<10^(-10))[1]]
      seqj=1:maxseqj

      R2=function(r,kappa,seqj=1:100){

        Ajtot=outer(kappa,seqj,Aj)

        if(r==0){
          val=colSums(t(Ajtot^2))
          val=(1+2*val)/(2*pi)
        }else{
          val=colSums(t(Ajtot^2)*(seqj)^(2*r))/pi
        }

        return(val)
      }

      amisef=function(kappa) hfunc(kappa,seqj)^2*intf2/4+R2(r,kappa,seqj)/n
      solopt=try(optim(kappaini,amisef,lower=0.001,upper=max(0.99,10^(log10(kappaini)/100)),method="L-BFGS-B"),silent=T)
    }

    if(inherits(solopt, "try-error")){
      kappaopt=kappaini
    }else{
      kappaopt=solopt$par
    }

  }


  return(kappaopt)


}





###################################################################
###################################################################

# This function settles the relation between the two smoothing parameters
# the one of the density derivative and the one of the functional

gammabw=function(h,r,intf2s1,intf2s2,kernel="vonmises",Q1,Q2){

  hopt=((-1)^(r+1)*2*Q1*intf2s1/((2*r+1)*Q2*intf2s2))^(2/(2*r+7))*h^((2*r+5)/(2*r+7))

  return(hopt)


}



###################################################################
###################################################################


# Equation that needs to be solved to find the optimal h (Step 4 in Algortihm 2)

amisefr2=function(r,intf2,n,kernel="vonmises",Q2,approximate=TRUE){

  hopt=((2*r+1)*Q2/(n*intf2))^(2/(2*r+5))

  return(hopt)


}



