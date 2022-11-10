#include <RcppArmadillo.h>
#include <cmath>



using namespace Rcpp;
using namespace arma;
using namespace std;


  double BESSI0(double X) {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067492;
      P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
      Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2;
      Q4=-0.157565e-2; Q5=0.916281e-2; Q6=-0.2057706e-1;
      Q7=0.2635537e-1; Q8=-0.1647633e-1; Q9=0.392377e-2;
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }



double A1invCpp(double x){
	double res;
	if( (0 <= x) & (x < 0.53) ){
		res =  2 * x + x*x*x + (5 * x*x*x*x*x)/6;
	} else{
		if(x < 0.85){
			res = -0.4 + 1.39 * x + 0.43/(1 - x);
		}else{
			res = 1/(x*x*x - 4 * x*x + 3 * x);
		}
	}
	return res;
}		


// [[Rcpp::depends(RcppArmadillo)]]
List mixVMN_EM(arma::vec T, arma::vec X, arma::vec pi, arma::vec mu, arma::vec kappa, arma::vec m, arma::vec sigma, const double tol, const int maxiters){
  int n = X.size();
  int J = pi.size();
  List L(9);

  double con = 2.506628274631;

 // Initial loglikelihood
 
 arma::mat loglik01(n,J);
 
  for (int i=0; i<n; ++i){
	  for (int j=0; j<J; ++j){
		  if(pi(j)<0.001){	
			  loglik01(i,j) = 0;	
		  }else{
			  loglik01(i,j) = pi(j)* exp(kappa(j)*cos(T(i)-mu(j)))*exp(-(X(i)-m(j))*(X(i)-m(j))/(2*sigma(j)*sigma(j)))/(con*con*con*sigma(j)*BESSI0(kappa(j)));
		  }
	  }
  }
  
  arma::vec loglik02 = log(sum(loglik01,1)) ;
  
  double loglik0 = sum(loglik02);
  
  arma::mat taus(n,J);
  taus = loglik01.each_col()/sum(loglik01,1) ;

  
  // estimation of the weights
  
  arma::vec pi0 = mean(taus.t(),1);
  
  // estimation of the parameters: mu0, kappa0, m0, sigma0
  
  arma::vec mu0(J), kappa0(J), m0(J), sigma0(J);
  for (int j=0; j<J; ++j){
	  mu0(j) = atan2(sum(taus.col(j)%sin(T)),sum(taus.col(j)%cos(T)));
	  m0(j) = sum(taus.col(j)%X)/sum(taus.col(j));
	  sigma0(j) = sqrt(sum(taus.col(j)%(X-m0(j))%(X-m0(j)))/sum(taus.col(j)));
	  double tmp =sum(taus.col(j)%cos(T-mu0(j)))/sum(taus.col(j));
	  if (tmp > 0){		
			double tmp2 = A1invCpp( tmp );	
			if (tmp2>100){	
				kappa0(j)=100;	
			}else{
				kappa0(j)=tmp2;
			}
		}else{	
			kappa0(j) = 0;	
		}
  }
  
  // loglikelihood with new parameters
  
  arma::mat loglik11(n,J);
 
  for (int i=0; i<n; ++i){
	  for (int j=0; j<J; ++j){
		  if(pi0(j)<0.001){
			  loglik11(i,j) = 0;
		  }else{
			  loglik11(i,j) = pi0(j)* exp(kappa0(j)*cos(T(i)-mu0(j)))*exp(-(X(i)-m0(j))*(X(i)-m0(j))/(2*sigma0(j)*sigma0(j)))/(con*con*con*sigma0(j)*BESSI0(kappa0(j)));
		  }
	  }
  }
  
  arma::vec loglik12 = log(sum(loglik11,1)) ;
  
  double loglik1 = sum(loglik12);
  
  
   double res1 = loglik0;
   double res2 = loglik1;
  // // while loglik1>loglik0 update taus, pi0 and parameters 
   int i = 2;
  while(i++<maxiters && loglik1 - loglik0 > tol){
	  loglik0 = loglik1;
	  taus = loglik11.each_col()/sum(loglik11,1) ;
	  pi0 = mean(taus.t(),1);
	  for (int j=0; j<J; ++j){
		mu0(j) = atan2(sum(taus.col(j)%sin(T)),sum(taus.col(j)%cos(T)));
		m0(j) = sum(taus.col(j)%X)/sum(taus.col(j));
		sigma0(j) = sqrt(sum(taus.col(j)%(X-m0(j))%(X-m0(j)))/sum(taus.col(j)));
		double tmp = sum(taus.col(j)%cos(T-mu0(j))) / sum(taus.col(j));
		if (tmp > 0){		
			double tmp2 = A1invCpp( tmp );	
			if (tmp2>100){		
				kappa0(j)=100;		
			}else{
				kappa0(j)=tmp2;
			}
		}else{
			kappa0(j) = 0;		
		}
		
	  }
	  for (int i=0; i<n; ++i){
		for (int j=0; j<J; ++j){
		  if(pi0(j)<0.001){		
			  loglik11(i,j) = 0;
		  }else{
			  loglik11(i,j) = pi0(j)* exp(kappa0(j)*cos(T(i)-mu0(j)))*exp(-(X(i)-m0(j))*(X(i)-m0(j))/(2*sigma0(j)*sigma0(j)))/(con*con*con*sigma0(j)*BESSI0(kappa0(j)));
		  }
		}
	  }
	  loglik12 = log(sum(loglik11,1)) ;
	  loglik1 = sum(loglik12);
	  
  }
  
	L[0] = res1;
	L[1] = res2;
	L[2] = pi0;
	L[3] = mu0;
	L[4] = kappa0;
	L[5] = m0;
	L[6] = sigma0;
	L[7] = loglik1;
	L[8] = i;
	

  
  return L;
}



// [[Rcpp::depends(RcppArmadillo)]]
int vecmaxInd(arma::vec x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference 
  return it - x.begin();
}

// [[Rcpp::depends(RcppArmadillo)]]
int vecminInd(arma::vec x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return it - x.begin();
}

// [[Rcpp::depends(RcppArmadillo)]]
List EM_vMN_init(arma::vec T, arma::vec X, int J, const double tolEM, int maxitersEM, int maxitersInit, int rep){
  List Lpi(rep);
  List Lmu(rep);
  List Lkappa(rep);
  List Lm(rep);
  List Lsigma(rep);
  arma::vec Llik(rep);

  double con = 2.506628274631;
 
  for (int i=0; i<rep; ++i){
	arma::vec tmp = runif(J); 	
	Lpi[i] = tmp/sum(tmp);	
	Lmu[i] = runif(J,0,con*con);	
	Lkappa[i] = runif(J,0,150);		
	Lm[i] = runif(J,min(X),max(X));		
	Lsigma[i] = runif(J,0,(max(X)-min(X))/4);		
	
	List tmp4 = mixVMN_EM(T,X,Lpi[i],Lmu[i],Lkappa[i],Lm[i],Lsigma[i],tolEM,maxitersInit);
	Llik(i) = tmp4[7];
  }
  int ind = vecmaxInd(Llik);
  List res = mixVMN_EM(T,X,Lpi[ind],Lmu[ind],Lkappa[ind],Lm[ind],Lsigma[ind],tolEM,maxitersEM);
  
  return res;
}




// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List R_EM_vMN_select(arma::vec T, arma::vec X, int M, const double tolEM, int maxitersEM, int maxitersInit, int rep){
  
  int n = X.size();
  List models(M-1);
  arma::vec nparam(M-1);
  arma::vec loglik(M-1);
  arma::vec BIC(M-1);
  List total(3);
  
  for (int i=0; i<M-1; ++i){
	  List tmp = EM_vMN_init(T,X,i+2,tolEM,maxitersEM,maxitersInit,rep);
	  models[i] = tmp;
	  nparam(i) = 4*(i+2)+i-1;
	  loglik(i) = tmp[7];
	  BIC(i) = nparam(i)*log(n)-2*loglik(i);
  }
  
  int ind = vecminInd(BIC);
  total[0] = models[ind];
  total[1] = BIC;
  total[2] = ind+2;
  
  return total;
}


// [[Rcpp::depends(RcppArmadillo)]]
List mixVMVM_EM(arma::vec T, arma::vec X, arma::vec pi, arma::vec mu, arma::vec kappa, arma::vec m, arma::vec nu, const double tol, const int maxiters){
  int n = X.size();
  int J = pi.size();
  List L(9);

  double con = 2.506628274631;

 // Initial loglikelihood
 
 arma::mat loglik01(n,J);
 
  for (int i=0; i<n; ++i){
	  for (int j=0; j<J; ++j){
		  if(pi(j)<0.001){	// !!!!!!!!!!!!!!!!!!
			  loglik01(i,j) = 0;	// !!!!!!!!!!!!!!!!!!
		  }else{
			  loglik01(i,j) = pi(j)* exp(kappa(j)*cos(T(i)-mu(j)))*exp(nu(j)*cos(X(i)-m(j)))/(con*con*con*con*BESSI0(kappa(j))*BESSI0(nu(j)));
		  }
	  }
  }
  
  arma::vec loglik02 = log(sum(loglik01,1)) ;
  
  double loglik0 = sum(loglik02);
  
  arma::mat taus(n,J);
  taus = loglik01.each_col()/sum(loglik01,1) ;

  
  // estimation of the weights
  
  arma::vec pi0 = mean(taus.t(),1);
  
  // estimation of the parameters: mu0, kappa0, m0, sigma0
  
  arma::vec mu0(J), kappa0(J), m0(J), nu0(J);
  for (int j=0; j<J; ++j){
	  mu0(j) = atan2(sum(taus.col(j)%sin(T)),sum(taus.col(j)%cos(T)));
	  m0(j) = atan2(sum(taus.col(j)%sin(X)),sum(taus.col(j)%cos(X)));
	  double tmp =sum(taus.col(j)%cos(T-mu0(j)))/sum(taus.col(j));
	  if (tmp > 0){		// !!!!!!!!!!!!!!!!!!
			double tmp2 = A1invCpp( tmp );	// !!!!!!!!!!!!!!!!!!
			if (tmp2>100){	// !!!!!!!!!!!!!!!!!!
				kappa0(j)=100;	// !!!!!!!!!!!!!!!!!!
			}else{
				kappa0(j)=tmp2;
			}
		}else{	// !!!!!!!!!!!!!!!!!!
			kappa0(j) = 0;	// !!!!!!!!!!!!!!!!!!
		}
	  double tmpnu =sum(taus.col(j)%cos(X-m0(j)))/sum(taus.col(j));
	  if (tmpnu > 0){		// !!!!!!!!!!!!!!!!!!
			double tmpnu2 = A1invCpp( tmpnu );	// !!!!!!!!!!!!!!!!!!
			if (tmpnu2>100){	// !!!!!!!!!!!!!!!!!!
				nu0(j)=100;	// !!!!!!!!!!!!!!!!!!
			}else{
				nu0(j)=tmpnu2;
			}
		}else{	// !!!!!!!!!!!!!!!!!!
			nu0(j) = 0;	// !!!!!!!!!!!!!!!!!!
		}
  }
  
  // loglikelihood with new parameters
  
  arma::mat loglik11(n,J);
 
  for (int i=0; i<n; ++i){
	  for (int j=0; j<J; ++j){
		  if(pi0(j)<0.001){
			  loglik11(i,j) = 0;
		  }else{
			  loglik11(i,j) = pi0(j)* exp(kappa0(j)*cos(T(i)-mu0(j)))*exp(nu0(j)*cos(X(i)-m0(j)))/(con*con*con*con*BESSI0(kappa0(j))*BESSI0(nu0(j)));
		  }
	  }
  }
  
  arma::vec loglik12 = log(sum(loglik11,1)) ;
  
  double loglik1 = sum(loglik12);
  
  
   double res1 = loglik0;
   double res2 = loglik1;
  // // while loglik1>loglik0 update taus, pi0 and parameters 
   int i = 2;
  while(i++<maxiters && loglik1 - loglik0 > tol){
	  loglik0 = loglik1;
	  taus = loglik11.each_col()/sum(loglik11,1) ;
	  pi0 = mean(taus.t(),1);
	  for (int j=0; j<J; ++j){
		mu0(j) = atan2(sum(taus.col(j)%sin(T)),sum(taus.col(j)%cos(T)));
		m0(j) = atan2(sum(taus.col(j)%sin(X)),sum(taus.col(j)%cos(X)));
		double tmp = sum(taus.col(j)%cos(T-mu0(j))) / sum(taus.col(j));
		if (tmp > 0){		// !!!!!!!!!!!!!!!!!!
			double tmp2 = A1invCpp( tmp );	// !!!!!!!!!!!!!!!!!!
			if (tmp2>100){		// !!!!!!!!!!!!!!!!!!
				kappa0(j)=100;		// !!!!!!!!!!!!!!!!!!
			}else{
				kappa0(j)=tmp2;
			}
		}else{
			kappa0(j) = 0;		// !!!!!!!!!!!!!!!!!!
		}
		double tmpnu =sum(taus.col(j)%cos(X-m0(j)))/sum(taus.col(j));
	    if (tmpnu > 0){		// !!!!!!!!!!!!!!!!!!
			double tmpnu2 = A1invCpp( tmpnu );	// !!!!!!!!!!!!!!!!!!
			if (tmpnu2>100){	// !!!!!!!!!!!!!!!!!!
				nu0(j)=100;	// !!!!!!!!!!!!!!!!!!
			}else{
				nu0(j)=tmpnu2;
			}
		}else{	// !!!!!!!!!!!!!!!!!!
			nu0(j) = 0;	// !!!!!!!!!!!!!!!!!!
		}
		
	  }
	  
	  for (int i=0; i<n; ++i){
		for (int j=0; j<J; ++j){
		  if(pi0(j)<0.001){		// !!!!!!!!!!!!!!!!!!
			  loglik11(i,j) = 0;	// !!!!!!!!!!!!!!!!!!
		  }else{
			  loglik11(i,j) = pi0(j)* exp(kappa0(j)*cos(T(i)-mu0(j)))*exp(nu0(j)*cos(X(i)-m0(j)))/(con*con*con*con*BESSI0(kappa0(j))*BESSI0(nu0(j)));
		  }
		}
	  }
	  loglik12 = log(sum(loglik11,1)) ;
	  loglik1 = sum(loglik12);
	  
  }
  
	L[0] = res1;
	L[1] = res2;
	L[2] = pi0;
	L[3] = mu0;
	L[4] = kappa0;
	L[5] = m0;
	L[6] = nu0;
	L[7] = loglik1;
	L[8] = i;
	

  
  return L;
}


// [[Rcpp::depends(RcppArmadillo)]]
List EM_vMvM_init(arma::vec T, arma::vec X, int J, const double tolEM, int maxitersEM, int maxitersInit, int rep){
  List Lpi(rep);
  List Lmu(rep);
  List Lkappa(rep);
  List Lm(rep);
  List Lnu(rep);
  arma::vec Llik(rep);

  double con = 2.506628274631;
 
  for (int i=0; i<rep; ++i){
	arma::vec tmp = runif(J); 	// !!!!!!!!!!!!!!!!!!
	Lpi[i] = tmp/sum(tmp);	// !!!!!!!!!!!!!!!!!!
	Lmu[i] = runif(J,0,con*con);	// !!!!!!!!!!!!!!!!!!
	Lkappa[i] = runif(J,0,100);		// !!!!!!!!!!!!!!!!!!
	Lm[i] = runif(J,0,con*con);		// !!!!!!!!!!!!!!!!!!
	Lnu[i] = runif(J,0,100);		// !!!!!!!!!!!!!!!!!!
	
	List tmp4 = mixVMVM_EM(T,X,Lpi[i],Lmu[i],Lkappa[i],Lm[i],Lnu[i],tolEM,maxitersInit);
	Llik(i) = tmp4[7];
  }
  int ind = vecmaxInd(Llik);
  List res = mixVMVM_EM(T,X,Lpi[ind],Lmu[ind],Lkappa[ind],Lm[ind],Lnu[ind],tolEM,maxitersEM);
  
  return res;
}




// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List R_EM_vMvM_select(arma::vec T, arma::vec X, int M, const double tolEM, int maxitersEM, int maxitersInit, int rep){
  
  int n = X.size();
  List models(M-1);
  arma::vec nparam(M-1);
  arma::vec loglik(M-1);
  arma::vec BIC(M-1);
  List total(3);
  
  for (int i=0; i<M-1; ++i){
	  List tmp = EM_vMvM_init(T,X,i+2,tolEM,maxitersEM,maxitersInit,rep);
	  models[i] = tmp;
	  nparam(i) = 4*(i+2)+i-1;
	  loglik(i) = tmp[7];
	  BIC(i) = nparam(i)*log(n)-2*loglik(i);
  }
  
  int ind = vecminInd(BIC);
  total[0] = models[ind];
  total[1] = BIC;
  total[2] = ind+2;
  
  return total;
}






