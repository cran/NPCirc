// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// my_fun_DoublePois
List my_fun_DoublePois(NumericVector xseq, NumericVector x, NumericVector y, NumericVector startv, NumericVector startv2, const double kappa, const double nu, double tol, int maxit);
RcppExport SEXP _NPCirc_my_fun_DoublePois(SEXP xseqSEXP, SEXP xSEXP, SEXP ySEXP, SEXP startvSEXP, SEXP startv2SEXP, SEXP kappaSEXP, SEXP nuSEXP, SEXP tolSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xseq(xseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type startv(startvSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type startv2(startv2SEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(my_fun_DoublePois(xseq, x, y, startv, startv2, kappa, nu, tol, maxit));
    return rcpp_result_gen;
END_RCPP
}
// my_fun_loglik_cv_nu
double my_fun_loglik_cv_nu(NumericVector x, NumericVector y, NumericVector startv, NumericVector startv2, const double kappa, const double nu);
RcppExport SEXP _NPCirc_my_fun_loglik_cv_nu(SEXP xSEXP, SEXP ySEXP, SEXP startvSEXP, SEXP startv2SEXP, SEXP kappaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type startv(startvSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type startv2(startv2SEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(my_fun_loglik_cv_nu(x, y, startv, startv2, kappa, nu));
    return rcpp_result_gen;
END_RCPP
}
// R_modereg_CircLin
List R_modereg_CircLin(NumericVector Y, NumericVector X, NumericVector T, const double kappa, const double h, int max_it, double tol);
RcppExport SEXP _NPCirc_R_modereg_CircLin(SEXP YSEXP, SEXP XSEXP, SEXP TSEXP, SEXP kappaSEXP, SEXP hSEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_modereg_CircLin(Y, X, T, kappa, h, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_modereg_CircCirc
List R_modereg_CircCirc(NumericVector Y, NumericVector X, NumericVector T, const double kappa1, const double kappa2, int max_it, double tol);
RcppExport SEXP _NPCirc_R_modereg_CircCirc(SEXP YSEXP, SEXP XSEXP, SEXP TSEXP, SEXP kappa1SEXP, SEXP kappa2SEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa1(kappa1SEXP);
    Rcpp::traits::input_parameter< const double >::type kappa2(kappa2SEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_modereg_CircCirc(Y, X, T, kappa1, kappa2, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_modereg_LinCirc
List R_modereg_LinCirc(NumericVector Y, NumericVector X, NumericVector T, const double h, const double kappa, int max_it, double tol);
RcppExport SEXP _NPCirc_R_modereg_LinCirc(SEXP YSEXP, SEXP XSEXP, SEXP TSEXP, SEXP hSEXP, SEXP kappaSEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type T(TSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_modereg_LinCirc(Y, X, T, h, kappa, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_CV_modereg_CircLin2
NumericVector R_CV_modereg_CircLin2(NumericVector Y, NumericVector X, double h1, NumericVector h2, int max_it, double tol);
RcppExport SEXP _NPCirc_R_CV_modereg_CircLin2(SEXP YSEXP, SEXP XSEXP, SEXP h1SEXP, SEXP h2SEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_CV_modereg_CircLin2(Y, X, h1, h2, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_re_CV_modereg_CircLin2
NumericVector R_re_CV_modereg_CircLin2(NumericVector Y, NumericVector X, NumericVector h1, double h2, int max_it, double tol);
RcppExport SEXP _NPCirc_R_re_CV_modereg_CircLin2(SEXP YSEXP, SEXP XSEXP, SEXP h1SEXP, SEXP h2SEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< double >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_re_CV_modereg_CircLin2(Y, X, h1, h2, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_CV_modereg_LinCirc
NumericVector R_CV_modereg_LinCirc(NumericVector Y, NumericVector X, const double h, NumericVector kappa, int max_it, double tol);
RcppExport SEXP _NPCirc_R_CV_modereg_LinCirc(SEXP YSEXP, SEXP XSEXP, SEXP hSEXP, SEXP kappaSEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_CV_modereg_LinCirc(Y, X, h, kappa, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_cv_re_modereg_LinCirc
NumericVector R_cv_re_modereg_LinCirc(NumericVector Y, NumericVector X, NumericVector h1, double h2, int max_it, double tol);
RcppExport SEXP _NPCirc_R_cv_re_modereg_LinCirc(SEXP YSEXP, SEXP XSEXP, SEXP h1SEXP, SEXP h2SEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< double >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_cv_re_modereg_LinCirc(Y, X, h1, h2, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_CV_grid_modereg_LinCirc
NumericMatrix R_CV_grid_modereg_LinCirc(NumericVector Y, NumericVector X, NumericVector h1, NumericVector h2, int max_it, double tol);
RcppExport SEXP _NPCirc_R_CV_grid_modereg_LinCirc(SEXP YSEXP, SEXP XSEXP, SEXP h1SEXP, SEXP h2SEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_CV_grid_modereg_LinCirc(Y, X, h1, h2, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_CV_modereg_CircCirc
NumericVector R_CV_modereg_CircCirc(NumericVector Y, NumericVector X, const double kappa1, NumericVector kappa2, int max_it, double tol);
RcppExport SEXP _NPCirc_R_CV_modereg_CircCirc(SEXP YSEXP, SEXP XSEXP, SEXP kappa1SEXP, SEXP kappa2SEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa1(kappa1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa2(kappa2SEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_CV_modereg_CircCirc(Y, X, kappa1, kappa2, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_cv_re_modereg_CircCirc
NumericVector R_cv_re_modereg_CircCirc(NumericVector Y, NumericVector X, NumericVector kappa1, double kappa2, int max_it, double tol);
RcppExport SEXP _NPCirc_R_cv_re_modereg_CircCirc(SEXP YSEXP, SEXP XSEXP, SEXP kappa1SEXP, SEXP kappa2SEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa1(kappa1SEXP);
    Rcpp::traits::input_parameter< double >::type kappa2(kappa2SEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_cv_re_modereg_CircCirc(Y, X, kappa1, kappa2, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_CV_grid_modereg_CircCirc
NumericMatrix R_CV_grid_modereg_CircCirc(NumericVector Y, NumericVector X, NumericVector h1, NumericVector h2, int max_it, double tol);
RcppExport SEXP _NPCirc_R_CV_grid_modereg_CircCirc(SEXP YSEXP, SEXP XSEXP, SEXP h1SEXP, SEXP h2SEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(R_CV_grid_modereg_CircCirc(Y, X, h1, h2, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// R_EM_vMN_select
List R_EM_vMN_select(arma::vec T, arma::vec X, int M, const double tolEM, int maxitersEM, int maxitersInit, int rep);
RcppExport SEXP _NPCirc_R_EM_vMN_select(SEXP TSEXP, SEXP XSEXP, SEXP MSEXP, SEXP tolEMSEXP, SEXP maxitersEMSEXP, SEXP maxitersInitSEXP, SEXP repSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const double >::type tolEM(tolEMSEXP);
    Rcpp::traits::input_parameter< int >::type maxitersEM(maxitersEMSEXP);
    Rcpp::traits::input_parameter< int >::type maxitersInit(maxitersInitSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    rcpp_result_gen = Rcpp::wrap(R_EM_vMN_select(T, X, M, tolEM, maxitersEM, maxitersInit, rep));
    return rcpp_result_gen;
END_RCPP
}
// R_EM_vMvM_select
List R_EM_vMvM_select(arma::vec T, arma::vec X, int M, const double tolEM, int maxitersEM, int maxitersInit, int rep);
RcppExport SEXP _NPCirc_R_EM_vMvM_select(SEXP TSEXP, SEXP XSEXP, SEXP MSEXP, SEXP tolEMSEXP, SEXP maxitersEMSEXP, SEXP maxitersInitSEXP, SEXP repSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const double >::type tolEM(tolEMSEXP);
    Rcpp::traits::input_parameter< int >::type maxitersEM(maxitersEMSEXP);
    Rcpp::traits::input_parameter< int >::type maxitersInit(maxitersInitSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    rcpp_result_gen = Rcpp::wrap(R_EM_vMvM_select(T, X, M, tolEM, maxitersEM, maxitersInit, rep));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NPCirc_my_fun_DoublePois", (DL_FUNC) &_NPCirc_my_fun_DoublePois, 9},
    {"_NPCirc_my_fun_loglik_cv_nu", (DL_FUNC) &_NPCirc_my_fun_loglik_cv_nu, 6},
    {"_NPCirc_R_modereg_CircLin", (DL_FUNC) &_NPCirc_R_modereg_CircLin, 7},
    {"_NPCirc_R_modereg_CircCirc", (DL_FUNC) &_NPCirc_R_modereg_CircCirc, 7},
    {"_NPCirc_R_modereg_LinCirc", (DL_FUNC) &_NPCirc_R_modereg_LinCirc, 7},
    {"_NPCirc_R_CV_modereg_CircLin2", (DL_FUNC) &_NPCirc_R_CV_modereg_CircLin2, 6},
    {"_NPCirc_R_re_CV_modereg_CircLin2", (DL_FUNC) &_NPCirc_R_re_CV_modereg_CircLin2, 6},
    {"_NPCirc_R_CV_modereg_LinCirc", (DL_FUNC) &_NPCirc_R_CV_modereg_LinCirc, 6},
    {"_NPCirc_R_cv_re_modereg_LinCirc", (DL_FUNC) &_NPCirc_R_cv_re_modereg_LinCirc, 6},
    {"_NPCirc_R_CV_grid_modereg_LinCirc", (DL_FUNC) &_NPCirc_R_CV_grid_modereg_LinCirc, 6},
    {"_NPCirc_R_CV_modereg_CircCirc", (DL_FUNC) &_NPCirc_R_CV_modereg_CircCirc, 6},
    {"_NPCirc_R_cv_re_modereg_CircCirc", (DL_FUNC) &_NPCirc_R_cv_re_modereg_CircCirc, 6},
    {"_NPCirc_R_CV_grid_modereg_CircCirc", (DL_FUNC) &_NPCirc_R_CV_grid_modereg_CircCirc, 6},
    {"_NPCirc_R_EM_vMN_select", (DL_FUNC) &_NPCirc_R_EM_vMN_select, 7},
    {"_NPCirc_R_EM_vMvM_select", (DL_FUNC) &_NPCirc_R_EM_vMvM_select, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_NPCirc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}