#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;
using namespace std;


NumericVector my_fun(NumericVector x){
    // calling order()
    Function f("order");   

    // Next code is interpreted as order(x)
    return f(x);
}


NumericVector my_fun2(int s, int fin){
    // calling seq()
    Function f("seq");   

    // Next code is interpreted as seq(s,fin)
    return f(s,fin);
}


NumericVector quantileCpp(NumericVector x, NumericVector probs) {
  Environment stats("package:stats");
  Function quantile = stats["quantile"];
  int npr = probs.size();
  NumericVector ans(npr);
  for(int i=0; i<npr; i++){
    ans[i] = as<double>(quantile(x, probs[i]));
  }
return ans;
}


NumericVector rcpp_package_circ_quantile(NumericVector x, NumericVector p){

    // Obtaining namespace of circular package
    Environment pkg = Environment::namespace_env("circular");

    // Picking up quantile.circular() function from circular package
    Function f = pkg["quantile.circular"];

    // Executing quantile.circular( x, p )
    return f( x, p);
}

// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List R_modereg_CircLin(NumericVector Y, NumericVector X, NumericVector T, const double kappa, const double h, int max_it, double tol){
	
	int n = X.size();
	int nt = T.size();
	//List L(nt);
	int kk = round(n/5);
	// temp variables
	mat Ku0ij(n,nt);
	NumericMatrix dist(n,nt);
	int l = 1;
	double KGj = 0;
	double KG_tot=0;
	double YKG_tot=0;
	double oldy=0;
	double newy=0;
	int iter_now;
	double err_now;
	

	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			Ku0ij(i,j) = exp( kappa*cos(X[i]-T[j]));
			dist(i,j) = 1 - cos(X[i]-T[j]);
		}
	}
	
	List L(nt);
	for(int j=0; j<nt; ++j){
		NumericVector index2 = my_fun2(l,kk);
		NumericVector index = my_fun(dist(_,j));
		index = index[index2];
		NumericVector yneig = Y[index-1];
		NumericVector ytemp(5);
		ytemp(0)=min(yneig);
		ytemp(1)=max(yneig);
		NumericVector aux1 = {0.25,0.5,0.75};
		NumericVector aux2 = quantileCpp(yneig,aux1);
		ytemp(2)=aux2(0);
		ytemp(3)=aux2(1);
		ytemp(4)=aux2(2);
		int ytemps = ytemp.size();
		NumericVector ym(ytemps);
		for (int jj=0; jj<ytemps; ++jj){
			newy=ytemp[jj];
			err_now=1e10;
			iter_now=0;
			while((iter_now < max_it)&&(err_now > tol)){
				//R_CheckUserInterrupt();
				YKG_tot=0;
				KG_tot=0;
				oldy = newy;

				for(int i=0; i<n; ++i){
					KGj = Ku0ij(i,j)*exp( -0.5*std::pow(((newy-Y[i])/h), 2) );
					KG_tot += KGj;
					YKG_tot += Y[i]*KGj;
				}
				if((KG_tot)<1e-10){
					newy=NA_REAL;
					break;
				}else{
					newy = YKG_tot/KG_tot;
					err_now = std::abs(newy-oldy);
					iter_now++;
				}
			}
			if((iter_now==max_it)&&(err_now > (tol*10))) newy=NA_REAL;
			ym[jj]=newy;
		}
		L[j]=ym;
		
	}
	
	
	return L;
}



// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List R_modereg_CircCirc(NumericVector Y, NumericVector X, NumericVector T, const double kappa1, const double kappa2, int max_it, double tol){
	
	int n = X.size();
	int nt = T.size();
	int kk = round(n/5);
	//List L(nt);
	
	// temp variables
	mat Ku0ij(n,nt);
	NumericMatrix dist(n,nt);
	double KGj = 0;
	double KG_tot=0;
	double YKG_tot1=0;
	double YKG_tot2=0;
	double oldy=0;
	double newy=0;
	int iter_now;
	double err_now;
	
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			Ku0ij(i,j) = exp(kappa1*cos(X[i]-T[j]));
			dist(i,j) = 1-cos(X[i]-T[j]);
		}
	}
	
	List L(nt);
	for(int j=0; j<nt; ++j){
		NumericVector index2 = my_fun2(1,kk);
		NumericVector index = my_fun(dist(_,j));
		index = index[index2];
		NumericVector yneig = Y[index-1];
		NumericVector aux1 = {0.05,0.25,0.5,0.75,0.95};
		NumericVector ytemp = rcpp_package_circ_quantile(yneig,aux1);
		int ytemps = ytemp.size();
		NumericVector ym(ytemps);
		for (int jj=0; jj<ytemps; ++jj){
			newy=ytemp[jj];
			err_now=2;
			iter_now=0;
			while((iter_now < max_it)&&(err_now > tol)){
				//R_CheckUserInterrupt();
				YKG_tot1=0;
				YKG_tot2=0;
				KG_tot=0;
				oldy = newy;

				for(int i=0; i<n; ++i){
					KGj = Ku0ij(i,j)*exp(kappa2*cos(newy-Y[i]));
					KG_tot += KGj;
					YKG_tot1 += sin(Y[i])*KGj;
					YKG_tot2 += cos(Y[i])*KGj;
				}
				if((KG_tot)<1e-10){
					newy=NA_REAL;
					break;
				}else{
					newy = atan2(YKG_tot1,YKG_tot2);
					err_now = 1-cos(newy-oldy);
					iter_now++;
				}
			}
			if((iter_now==max_it)&&(err_now > (tol*10))) newy=NA_REAL;
			ym[jj]=newy;
		}
		L[j]=ym;
	}

	
	return L;
}


// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List R_modereg_LinCirc(NumericVector Y, NumericVector X, NumericVector T, const double h, const double kappa, int max_it, double tol){
	
	int n = X.size();
	int nt = T.size();
	//List L(nt);
	int kk = round(n/5);
	
	// temp variables
	mat Ku0ij(n,nt);
	NumericMatrix dist(n,nt);
	double KGj = 0;
	double KG_tot=0;
	double YKG_tot1=0;
	double YKG_tot2=0;
	double oldy=0;
	double newy=0;
	int iter_now;
	double err_now;
	
	for(int i=0; i<n; ++i){
		for(int j=0; j<nt; ++j){
			Ku0ij(i,j) = exp( -0.5*std::pow(((T[j]-X[i])/h), 2) );
			dist(i,j) = std::abs(X[i]-T[j]);
		}
	}
	
	List L(nt);
	for(int j=0; j<nt; ++j){
		NumericVector index2 = my_fun2(1,kk);
		NumericVector index = my_fun(dist(_,j));
		index = index[index2];
		NumericVector yneig = Y[index-1];
		NumericVector aux1 = {0.05,0.25,0.5,0.75,0.95};
		NumericVector ytemp = rcpp_package_circ_quantile(yneig,aux1);
		int ytemps = ytemp.size();
		NumericVector ym(ytemps);
		for (int jj=0; jj<ytemps; ++jj){
			newy=ytemp[jj];
			err_now=2;
			iter_now=0;
			while((iter_now < max_it)&&(err_now > tol)){
				//R_CheckUserInterrupt();
				YKG_tot1=0;
				YKG_tot2=0;
				KG_tot=0;
				oldy = newy;

				for(int i=0; i<n; ++i){
					KGj = Ku0ij(i,j)*exp(kappa*cos(newy-Y[i]));
					KG_tot += KGj;
					YKG_tot1 += sin(Y[i])*KGj;
					YKG_tot2 += cos(Y[i])*KGj;
				}
				if((KG_tot)<1e-10){
					newy=NA_REAL;
					break;
				}else{
					newy = atan2(YKG_tot1,YKG_tot2);
					err_now = 1-cos(newy-oldy);
					iter_now++;
				}
			}
			if((iter_now==max_it)&&(err_now > (tol*10))) newy=NA_REAL;
			ym[jj]=newy;
		}
		L[j]=ym;
	}

	
	return L;
}



// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericVector R_CV_modereg_CircLin2(NumericVector Y, NumericVector X, double h1, NumericVector h2, int max_it, double tol){


  int n = X.size();
  int nh2 = h2.size();
  int kk = round(n/5);
  
  // temp variable
  NumericVector CVh(nh2);
  mat Ku0ij(n,n);
  NumericMatrix dist(n,n);
  double KGj = 0;
  double KG_tot=0;
  double YKG_tot=0;
  double oldy=0;
  double newy=0;
  int iter_now;
  double err_now;
  
 

	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
			Ku0ij(i,j) = exp( h1*cos(X[i]-X[j]));
			dist(i,j) = 1 - cos(X[i]-X[j]);
		}
	}
  
	
    for(int k2=0; k2<nh2; ++k2){
      double cvj=0;
      for(int j=0; j<n; ++j){
		NumericVector index2 = my_fun2(1,kk);
		NumericVector index = my_fun(dist(_,j));
		index = index[index2];
		NumericVector yneig = Y[index-1];
		NumericVector ytemp(5);
		ytemp(0)=min(yneig);
		ytemp(1)=max(yneig);
		NumericVector aux1 = {0.25,0.5,0.75};
		NumericVector aux2 = quantileCpp(yneig,aux1);
		ytemp(2)=aux2(0);
		ytemp(3)=aux2(1);
		ytemp(4)=aux2(2);
		NumericVector ym(5);
		for (int jj=0; jj<5; ++jj){
			newy=ytemp[jj];
            err_now=1e10;
            iter_now=0;
            while((iter_now < max_it)&&(err_now > tol)){
              //R_CheckUserInterrupt();
              YKG_tot=0;
              KG_tot=0;
              oldy = newy;
              for(int i=0; i<n; ++i){
                if(i!=j){
                  KGj = Ku0ij(i,j)*exp( -0.5*std::pow(((newy-Y[i])/h2[k2]), 2) );
                  KG_tot += KGj;
                  YKG_tot += Y[i]*KGj;
                }
              }
              if((KG_tot)<1e-10){
                newy=99999;
                break;
              }else{
                newy = YKG_tot/KG_tot;
                err_now = std::abs(newy-oldy);
                iter_now++;
              }
            }
            if((iter_now==max_it)&&(err_now > (tol*10))) newy=99999;
            ym[jj]=newy;
          }
          ym = round(ym*100,0);
          NumericVector ymnew = unique(ym)/100;
          cvj += std::pow(min(abs(ymnew-Y[j]))*(ymnew.size()), 2);
        }
        CVh(k2) = cvj/(n+0.0);
    }
  
  
  return CVh;

}



double CV_re_aux_modereg_CircLin2(NumericVector Y, NumericVector X, double h1, double h2, int max_it, double tol){


  int n = X.size();
  int kk = round(n/5);
  
  // temp variable
  mat Ku0ij(n,n);
  NumericMatrix dist(n,n);
  double KGj = 0;
  double KG_tot=0;
  double YKG_tot=0;
  double oldy=0;
  double newy=0;
  int iter_now;
  double err_now;
  
 

	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
			Ku0ij(i,j) = exp( h1*cos(X[i]-X[j]));
			dist(i,j) = 1 - cos(X[i]-X[j]);
		}
	}
  
	
	  double cvj = 0;
      for(int j=0; j<n; ++j){
		NumericVector ym(5);
		for (int jj=0; jj<5; ++jj){
			NumericVector index2 = my_fun2(1,kk);
			NumericVector index = my_fun(dist(_,j));
			index = index[index2];
			NumericVector yneig = Y[index-1];
			NumericVector ytemp(5);
			ytemp(0)=min(yneig);
			ytemp(1)=max(yneig);
			NumericVector aux1 = {0.25,0.5,0.75};
			NumericVector aux2 = quantileCpp(yneig,aux1);
			ytemp(2)=aux2(0);
			ytemp(3)=aux2(1);
			ytemp(4)=aux2(2);
			newy=ytemp[jj];
            err_now=1e10;
            iter_now=0;
            while((iter_now < max_it)&&(err_now > tol)){
              //R_CheckUserInterrupt();
              YKG_tot=0;
              KG_tot=0;
              oldy = newy;
              for(int i=0; i<n; ++i){
                if(i!=j){
                  KGj = Ku0ij(i,j)*exp( -0.5*std::pow(((newy-Y[i])/h2), 2) );
                  KG_tot += KGj;
                  YKG_tot += Y[i]*KGj;
                }
              }
              if((KG_tot)<1e-10){
                newy=99999;
                break;
              }else{
                newy = YKG_tot/KG_tot;
                err_now = std::abs(newy-oldy);
                iter_now++;
              }
            }
            if((iter_now==max_it)&&(err_now > (tol*10))) newy=99999;
            ym[jj]=newy;
          }
          ym = round(ym*100,0);
          NumericVector ymnew = unique(ym)/100;
          cvj += std::pow(min(abs(ymnew-Y[j]))*(ymnew.size()), 2);
        }

  
  
  return cvj;

}



// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericVector R_re_CV_modereg_CircLin2(NumericVector Y, NumericVector X, NumericVector h1, double h2, int max_it, double tol){
	
	int nh1 = h1.size();
	
	NumericVector CV(nh1);
	for (int i=0; i<nh1; ++i){
		CV(i) = CV_re_aux_modereg_CircLin2(Y,X,h1(i),h2,max_it,tol);
	}
	
	return CV;
}



// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericVector R_CV_modereg_LinCirc(NumericVector Y, NumericVector X, const double h, NumericVector kappa, int max_it, double tol){
	
	int n = X.size();
	int nh2 = kappa.size();
	int kk = round(n/5);
	//List L(nt);
	
	// temp variables
	NumericVector CVh(nh2);
	mat Ku0ij(n,n);
	NumericMatrix dist(n,n);
	double KGj = 0;
	double KG_tot=0;
	double YKG_tot1=0;
	double YKG_tot2=0;
	double oldy=0;
	double newy=0;
	int iter_now;
	double err_now;
	
	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
			Ku0ij(i,j) = exp(-(X[i]-X[j])*(X[i]-X[j])/(2*h*h));
			dist(i,j) = std::abs(X[i]-X[j]);
			
		}
	}
	
	for(int k2=0; k2<nh2; ++k2){
    double cvj=0;
	for(int j=0; j<n; ++j){
		NumericVector index2 = my_fun2(1,kk);
		NumericVector index = my_fun(dist(_,j));
		index = index[index2];
		NumericVector yneig = Y[index-1];
		NumericVector aux1 = {0.05,0.25,0.5,0.75,0.95};
		NumericVector ytemp = rcpp_package_circ_quantile(yneig,aux1);
		int ytemps = ytemp.size();
		NumericVector ym(ytemps);
		for (int jj=0; jj<ytemp.size(); ++jj){
			newy=ytemp[jj];
			err_now=2;
			iter_now=0;
			while((iter_now < max_it)&&(err_now > tol)){
				//R_CheckUserInterrupt();
				YKG_tot1=0;
				YKG_tot2=0;
				KG_tot=0;
				oldy = newy;

				for(int i=0; i<n; ++i){
					if(i!=j){
						KGj = Ku0ij(i,j)*exp(kappa[k2]*cos(newy-Y[i]));
						KG_tot += KGj;
						YKG_tot1 += sin(Y[i])*KGj;
						YKG_tot2 += cos(Y[i])*KGj;
					}
				}
				if((KG_tot)<1e-10){
					newy=NA_REAL;
					break;
				}else{
					newy = atan2(YKG_tot1,YKG_tot2);
					err_now = 1-cos(newy-oldy);
					iter_now++;
				}
			}
			if((iter_now==max_it)&&(err_now > (tol*10))) newy=NA_REAL;
			ym[jj]=newy;
		}
		ym = round(ym*100,0);
        NumericVector ymnew = unique(ym)/100;
        cvj += min(1-cos(ymnew-Y[j]))*std::pow((ymnew.size()), 1);
	}
    CVh(k2) = cvj/(n+0.0);
    }
	
	return CVh;
}




double CV_modereg_LinCirc_single(NumericVector Y, NumericVector X, const double h, const double kappa, int max_it, double tol){
	
	int n = X.size();
	int kk = round(n/5);
	//List L(nt);
	
	// temp variables
	mat Ku0ij(n,n);
	NumericMatrix dist(n,n);
	double KGj = 0;
	double KG_tot=0;
	double YKG_tot1=0;
	double YKG_tot2=0;
	double oldy=0;
	double newy=0;
	int iter_now;
	double err_now;
	
	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
			Ku0ij(i,j) = exp(-(X[i]-X[j])*(X[i]-X[j])/(2*h*h));
			dist(i,j) = std::abs(X[i]-X[j]);
			
		}
	}
	

    double cvj=0;
	for(int j=0; j<n; ++j){
		NumericVector index2 = my_fun2(1,kk);
		NumericVector index = my_fun(dist(_,j));
		index = index[index2];
		NumericVector yneig = Y[index-1];
		NumericVector aux1 = {0.05,0.25,0.5,0.75,0.95};
		NumericVector ytemp = rcpp_package_circ_quantile(yneig,aux1);
		int ytemps = ytemp.size();
		NumericVector ym(ytemps);
		for (int jj=0; jj<ytemp.size(); ++jj){
			newy=ytemp[jj];
			err_now=2;
			iter_now=0;
			while((iter_now < max_it)&&(err_now > tol)){
				//R_CheckUserInterrupt();
				YKG_tot1=0;
				YKG_tot2=0;
				KG_tot=0;
				oldy = newy;

				for(int i=0; i<n; ++i){
					if(i!=j){
						KGj = Ku0ij(i,j)*exp(kappa*cos(newy-Y[i]));
						KG_tot += KGj;
						YKG_tot1 += sin(Y[i])*KGj;
						YKG_tot2 += cos(Y[i])*KGj;
					}
				}
				if((KG_tot)<1e-10){
					newy=NA_REAL;
					break;
				}else{
					newy = atan2(YKG_tot1,YKG_tot2);
					err_now = 1-cos(newy-oldy);
					iter_now++;
				}
			}
			if((iter_now==max_it)&&(err_now > (tol*10))) newy=NA_REAL;
			ym[jj]=newy;
		}
		ym = round(ym*100,0);
        NumericVector ymnew = unique(ym)/100;
        cvj += min(1-cos(ymnew-Y[j]))*std::pow((ymnew.size()), 1);
	}
	
	return cvj;
}

// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericVector R_cv_re_modereg_LinCirc(NumericVector Y, NumericVector X, NumericVector h1, double h2, int max_it, double tol){
	
	int nh1 = h1.size();
	
	NumericVector CV(nh1);
	for (int i=0; i<nh1; ++i){
		CV(i) = CV_modereg_LinCirc_single(Y,X,h1(i),h2,max_it,tol);
	}
	
	return CV;
}


// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericMatrix R_CV_grid_modereg_LinCirc(NumericVector Y, NumericVector X, NumericVector h1, NumericVector h2, int max_it, double tol){
	
	int nh1 = h1.size();
	int nh2 = h2.size();
	
	//output variable
	NumericMatrix CV(nh1,nh2);
	
	for(int i=0; i<nh1; ++i){
		for (int j=0; j<nh2; ++j){
			CV(i,j) = CV_modereg_LinCirc_single(Y,  X, h1(i), h2(j),  max_it,  tol);
		}
	}
	
	return CV;
}




// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericVector R_CV_modereg_CircCirc(NumericVector Y, NumericVector X, const double kappa1, NumericVector kappa2, int max_it, double tol){
	
	int n = X.size();
	int nh2 = kappa2.size();
	int kk = round(n/5);
	//List L(nt);
	
	// temp variables
	NumericVector CVh(nh2);
	mat Ku0ij(n,n);
	NumericMatrix dist(n,n);
	double KGj = 0;
	double KG_tot=0;
	double YKG_tot1=0;
	double YKG_tot2=0;
	double oldy=0;
	double newy=0;
	int iter_now;
	double err_now;
	
	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
			Ku0ij(i,j) = exp(kappa1*cos(X[i]-X[j]));
			dist(i,j) = 1-cos(X[i]-X[j]);
		}
	}
	
	for(int k2=0; k2<nh2; ++k2){
    double cvj=0;
	for(int j=0; j<n; ++j){
		NumericVector index2 = my_fun2(1,kk);
		NumericVector index = my_fun(dist(_,j));
		index = index[index2];
		NumericVector yneig = Y[index-1];
		NumericVector aux1 = {0.05,0.25,0.5,0.75,0.95};
		NumericVector ytemp = rcpp_package_circ_quantile(yneig,aux1);
		int ytemps = ytemp.size();
		NumericVector ym(ytemps);
		for (int jj=0; jj<ytemps; ++jj){
			newy=ytemp[jj];
			err_now=2;
			iter_now=0;
			while((iter_now < max_it)&&(err_now > tol)){
				//R_CheckUserInterrupt();
				YKG_tot1=0;
				YKG_tot2=0;
				KG_tot=0;
				oldy = newy;

				for(int i=0; i<n; ++i){
					if(i!=j){
						KGj = Ku0ij(i,j)*exp(kappa2[k2]*cos(newy-Y[i]));
						KG_tot += KGj;
						YKG_tot1 += sin(Y[i])*KGj;
						YKG_tot2 += cos(Y[i])*KGj;
					}
				}
				if((KG_tot)<1e-10){
					newy=NA_REAL;
					break;
				}else{
					newy = atan2(YKG_tot1,YKG_tot2);
					err_now = 1-cos(newy-oldy);
					iter_now++;
				}
			}
			if((iter_now==max_it)&&(err_now > (tol*10))) newy=NA_REAL;
			ym[jj]=newy;
		}
		ym = round(ym*100,0);
        NumericVector ymnew = unique(ym)/100;
        cvj += std::pow(min(1-cos(ymnew-Y[j]))*(ymnew.size()), 2);
	}
    CVh(k2) = cvj/(n+0.0);
    }
	
	return CVh;
}



double CV_modereg_CircCirc_single(NumericVector Y, NumericVector X, const double kappa1, const double kappa2, int max_it, double tol){
	
	int n = X.size();
	int kk = round(n/5);
	//List L(nt);
	
	// temp variables
	mat Ku0ij(n,n);
	NumericMatrix dist(n,n);
	double KGj = 0;
	double KG_tot=0;
	double YKG_tot1=0;
	double YKG_tot2=0;
	double oldy=0;
	double newy=0;
	int iter_now;
	double err_now;
	
	for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
			Ku0ij(i,j) = exp(kappa1*cos(X[i]-X[j]));
			dist(i,j) = 1-cos(X[i]-X[j]);
		}
	}
	

    double cvj=0;
	for(int j=0; j<n; ++j){
		NumericVector index2 = my_fun2(1,kk);
		NumericVector index = my_fun(dist(_,j));
		index = index[index2];
		NumericVector yneig = Y[index-1];
		NumericVector aux1 = {0.05,0.25,0.5,0.75,0.95};
		NumericVector ytemp = rcpp_package_circ_quantile(yneig,aux1);
		int ytemps = ytemp.size();
		NumericVector ym(ytemps);
		for (int jj=0; jj<ytemp.size(); ++jj){
			newy=ytemp[jj];
			err_now=2;
			iter_now=0;
			while((iter_now < max_it)&&(err_now > tol)){
				//R_CheckUserInterrupt();
				YKG_tot1=0;
				YKG_tot2=0;
				KG_tot=0;
				oldy = newy;

				for(int i=0; i<n; ++i){
					if(i!=j){
						KGj = Ku0ij(i,j)*exp(kappa2*cos(newy-Y[i]));
						KG_tot += KGj;
						YKG_tot1 += sin(Y[i])*KGj;
						YKG_tot2 += cos(Y[i])*KGj;
					}
				}
				if((KG_tot)<1e-10){
					newy=NA_REAL;
					break;
				}else{
					newy = atan2(YKG_tot1,YKG_tot2);
					err_now = 1-cos(newy-oldy);
					iter_now++;
				}
			}
			if((iter_now==max_it)&&(err_now > (tol*10))) newy=NA_REAL;
			ym[jj]=newy;
		}
		ym = round(ym*100,0);
        NumericVector ymnew = unique(ym)/100;
        cvj += min(1-cos(ymnew-Y[j]))*std::pow((ymnew.size()), 1);
	}
	
	return cvj;
}

// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericVector R_cv_re_modereg_CircCirc(NumericVector Y, NumericVector X, NumericVector kappa1, double kappa2, int max_it, double tol){
	
	int nh1 = kappa1.size();
	
	NumericVector CV(nh1);
	for (int i=0; i<nh1; ++i){
		CV(i) = CV_modereg_CircCirc_single(Y,X,kappa1(i),kappa2,max_it,tol);
	}
	
	return CV;
}


// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
NumericMatrix R_CV_grid_modereg_CircCirc(NumericVector Y, NumericVector X, NumericVector h1, NumericVector h2, int max_it, double tol){
	
	int nh1 = h1.size();
	int nh2 = h2.size();
	
	//output variable
	NumericMatrix CV(nh1,nh2);
	
	for(int i=0; i<nh1; ++i){
		for (int j=0; j<nh2; ++j){
			CV(i,j) = CV_modereg_CircCirc_single(Y,  X, h1(i), h2(j),  max_it,  tol);
		}
	}
	
	return CV;
}


