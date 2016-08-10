// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// pn
double pn(mat y, mat mu, mat sigma);
RcppExport SEXP mcif_pn(SEXP ySEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< mat >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(pn(y, mu, sigma));
    return __result;
END_RCPP
}
// loglikfull
vec loglikfull(mat y, mat b, mat u, mat sigma, mat alph, mat dalph, mat tau, bool full);
RcppExport SEXP mcif_loglikfull(SEXP ySEXP, SEXP bSEXP, SEXP uSEXP, SEXP sigmaSEXP, SEXP alphSEXP, SEXP dalphSEXP, SEXP tauSEXP, SEXP fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< mat >::type u(uSEXP);
    Rcpp::traits::input_parameter< mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< mat >::type alph(alphSEXP);
    Rcpp::traits::input_parameter< mat >::type dalph(dalphSEXP);
    Rcpp::traits::input_parameter< mat >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< bool >::type full(fullSEXP);
    __result = Rcpp::wrap(loglikfull(y, b, u, sigma, alph, dalph, tau, full));
    return __result;
END_RCPP
}
// Dloglikfull
mat Dloglikfull(mat y, mat b, mat u, mat sigma, mat alph, mat dalph, mat tau, bool full);
RcppExport SEXP mcif_Dloglikfull(SEXP ySEXP, SEXP bSEXP, SEXP uSEXP, SEXP sigmaSEXP, SEXP alphSEXP, SEXP dalphSEXP, SEXP tauSEXP, SEXP fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< mat >::type u(uSEXP);
    Rcpp::traits::input_parameter< mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< mat >::type alph(alphSEXP);
    Rcpp::traits::input_parameter< mat >::type dalph(dalphSEXP);
    Rcpp::traits::input_parameter< mat >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< bool >::type full(fullSEXP);
    __result = Rcpp::wrap(Dloglikfull(y, b, u, sigma, alph, dalph, tau, full));
    return __result;
END_RCPP
}
// D2loglikfull
mat D2loglikfull(mat y, mat b, mat u, mat sigma, mat alph, mat dalph, mat tau, bool full);
RcppExport SEXP mcif_D2loglikfull(SEXP ySEXP, SEXP bSEXP, SEXP uSEXP, SEXP sigmaSEXP, SEXP alphSEXP, SEXP dalphSEXP, SEXP tauSEXP, SEXP fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< mat >::type u(uSEXP);
    Rcpp::traits::input_parameter< mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< mat >::type alph(alphSEXP);
    Rcpp::traits::input_parameter< mat >::type dalph(dalphSEXP);
    Rcpp::traits::input_parameter< mat >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< bool >::type full(fullSEXP);
    __result = Rcpp::wrap(D2loglikfull(y, b, u, sigma, alph, dalph, tau, full));
    return __result;
END_RCPP
}
// loglik
vec loglik(mat y, mat b, mat sigma, mat alph, mat dalph, mat tau, mat eb0, int nq, double stepsize, unsigned iter, bool debug);
RcppExport SEXP mcif_loglik(SEXP ySEXP, SEXP bSEXP, SEXP sigmaSEXP, SEXP alphSEXP, SEXP dalphSEXP, SEXP tauSEXP, SEXP eb0SEXP, SEXP nqSEXP, SEXP stepsizeSEXP, SEXP iterSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< mat >::type alph(alphSEXP);
    Rcpp::traits::input_parameter< mat >::type dalph(dalphSEXP);
    Rcpp::traits::input_parameter< mat >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< mat >::type eb0(eb0SEXP);
    Rcpp::traits::input_parameter< int >::type nq(nqSEXP);
    Rcpp::traits::input_parameter< double >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< unsigned >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    __result = Rcpp::wrap(loglik(y, b, sigma, alph, dalph, tau, eb0, nq, stepsize, iter, debug));
    return __result;
END_RCPP
}
// EB
mat EB(mat y, mat b, mat sigma, mat alph, mat dalph, mat tau, double stepsize, unsigned iter, bool debug);
RcppExport SEXP mcif_EB(SEXP ySEXP, SEXP bSEXP, SEXP sigmaSEXP, SEXP alphSEXP, SEXP dalphSEXP, SEXP tauSEXP, SEXP stepsizeSEXP, SEXP iterSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< mat >::type b(bSEXP);
    Rcpp::traits::input_parameter< mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< mat >::type alph(alphSEXP);
    Rcpp::traits::input_parameter< mat >::type dalph(dalphSEXP);
    Rcpp::traits::input_parameter< mat >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< unsigned >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    __result = Rcpp::wrap(EB(y, b, sigma, alph, dalph, tau, stepsize, iter, debug));
    return __result;
END_RCPP
}
