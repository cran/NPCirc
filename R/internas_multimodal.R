
Omega_a<-function(a,m1,m2,s1,s2,g){
  res<-dnorm(m1-m2,mean=0,sd=sqrt(a*g^2+s1^2+s2^2))
  return(res)
}

intP1<-function(x,h,mu1,mu2,kappa1,kappa2){
  tmp<-sqrt(h^(-4)+kappa1^2+2*h^(-2)*kappa1*cos(x-mu1))
  Bx<-exp(kappa2*cos(x-mu2))*2*pi*besselI(tmp,0,expon.scaled=TRUE)
  #logBx<-kappa2*cos(x-mu2)+log(2*pi)+log(besselI(tmp,0,expon.scaled=TRUE))
  #logAx<-logBx+tmp
  res<-Bx*exp(tmp)
  #res<-exp(kappa2*cos(x-mu2))*2*pi*besselI(sqrt(h^(-4)+kappa1^2+2*h^(-2)*kappa1*cos(x-mu1)),0)
  return(res)
}

intP2<-function(x,h,mu1,mu2,kappa1,kappa2){
  res<-(2*pi)^2*besselI(sqrt(h^(-4)+kappa1^2+2*h^(-2)*kappa1*cos(x-mu1)),0)*besselI(sqrt(h^(-4)+kappa2^2+2*h^(-2)*kappa2*cos(x-mu2)),0)
}


MISECL<-function(par,mu,kappa,m,s,pip,n){
  
  h<-par[1]
  g<-par[2]
  k<-length(pip)
  pip2<-pip/(2*pi*besselI(kappa,0))
  Psi0<-matrix(nrow=k,ncol=k)
  Psi1<-matrix(nrow=k,ncol=k)
  Psi2<-matrix(nrow=k,ncol=k)
  Omega0<-matrix(nrow=k,ncol=k)
  Omega1<-matrix(nrow=k,ncol=k)
  Omega2<-matrix(nrow=k,ncol=k)
  for (i in 1:k){
    for (j in 1:k){
      Omega0[i,j]<-Omega_a(0,m[i],m[j],s[i],s[j],g)
      Omega1[i,j]<-Omega_a(1,m[i],m[j],s[i],s[j],g)
      Omega2[i,j]<-Omega_a(2,m[i],m[j],s[i],s[j],g)
      
      Int1<-integrate(intP1,lower=0,upper=2*pi,h=h,mu1=mu[i],mu2=mu[j],kappa1=kappa[i],kappa2=kappa[j])$value
      Int2<-integrate(intP2,lower=0,upper=2*pi,h=h,mu1=mu[i],mu2=mu[j],kappa1=kappa[i],kappa2=kappa[j])$value
      
      Psi0[i,j]<- 2*pi*besselI(sqrt(kappa[i]^2+kappa[j]^2+2*kappa[i]*kappa[j]*cos(mu[i]-mu[j])),0)
      Psi1[i,j]<- Int1
      Psi2[i,j]<- Int2
    }
  }
  Psi1<-((2*pi)*besselI(h^(-2),0))^(-1)*Psi1
  Psi2<-((2*pi)^2*besselI(h^(-2),0)^2)^(-1)*Psi2
  A<-(1-1/n)*Psi2*Omega2 -2*Psi1*Omega1 + Psi0*Omega0
  B<-2*pi*2*sqrt(pi)*g*n*besselI(h^(-2),0)^2/(besselI(2*h^(-2),0))
  res<-(1/B)+t(pip2)%*%A%*%pip2
  return(res)
}


# Functions for selection of the reference parameters


intP1<-function(x,h,mu1,mu2,kappa1,kappa2){
  res<-exp(kappa2*cos(x-mu2))*2*pi*besselI(sqrt(h^(-4)+kappa1^2+2*h^(-2)*kappa1*cos(x-mu1)),0)
}

intP2<-function(x,h,mu1,mu2,kappa1,kappa2){
  res<-(2*pi)^2*besselI(sqrt(h^(-4)+kappa1^2+2*h^(-2)*kappa1*cos(x-mu1)),0)*besselI(sqrt(h^(-4)+kappa2^2+2*h^(-2)*kappa2*cos(x-mu2)),0)
}


MISECC<-function(par,mu,kappa,m,nu,pip,n){
  
  h<-par[1]
  g<-par[2]
  k<-length(pip)
  pip2<-pip/(2*pi*besselI(kappa,0)*2*pi*besselI(nu,0))
  Psi0h<-matrix(nrow=k,ncol=k)
  Psi1h<-matrix(nrow=k,ncol=k)
  Psi2h<-matrix(nrow=k,ncol=k)
  Psi0g<-matrix(nrow=k,ncol=k)
  Psi1g<-matrix(nrow=k,ncol=k)
  Psi2g<-matrix(nrow=k,ncol=k)
  for (i in 1:k){
    for (j in 1:k){
      
      Int1h<-integrate(intP1,lower=0,upper=2*pi,h=h,mu1=mu[i],mu2=mu[j],kappa1=kappa[i],kappa2=kappa[j])$value
      Int2h<-integrate(intP2,lower=0,upper=2*pi,h=h,mu1=mu[i],mu2=mu[j],kappa1=kappa[i],kappa2=kappa[j])$value
      
      Int1g<-integrate(intP1,lower=0,upper=2*pi,h=g,mu1=m[i],mu2=m[j],kappa1=nu[i],kappa2=nu[j])$value
      Int2g<-integrate(intP2,lower=0,upper=2*pi,h=g,mu1=m[i],mu2=m[j],kappa1=nu[i],kappa2=nu[j])$value
      
      
      Psi0h[i,j]<- 2*pi*besselI(sqrt(kappa[i]^2+kappa[j]^2+2*kappa[i]*kappa[j]*cos(mu[i]-mu[j])),0)
      Psi1h[i,j]<- Int1h
      Psi2h[i,j]<- Int2h
      
      Psi0g[i,j]<- 2*pi*besselI(sqrt(nu[i]^2+nu[j]^2+2*nu[i]*nu[j]*cos(m[i]-m[j])),0)
      Psi1g[i,j]<- Int1g
      Psi2g[i,j]<- Int2g
    }
  }
  Psi1h<-((2*pi)*besselI(h^(-2),0))^(-1)*Psi1h
  Psi2h<-((2*pi)^2*besselI(h^(-2),0)^2)^(-1)*Psi2h
  Psi1g<-((2*pi)*besselI(g^(-2),0))^(-1)*Psi1g
  Psi2g<-((2*pi)^2*besselI(g^(-2),0)^2)^(-1)*Psi2g
  A<-(1-1/n)*Psi2h*Psi2g -2*Psi1h*Psi1g + Psi0h*Psi0g
  B<-(2*pi)^2*n*besselI(h^(-2),0)^2*besselI(g^(-2),0)^2/(besselI(2*h^(-2),0)*besselI(2*g^(-2),0))
  res<-(1/B)+t(pip2)%*%A%*%pip2
  return(res)
}
