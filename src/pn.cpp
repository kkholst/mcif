#include "pn.h"
#include <mvtnormAPI.h>

static int _mvt_df = 0;
static double _mvt_abseps=0.00001;
static double _mvt_releps=0;
static int _mvt_maxpts=20000;
static int _mvt_inform;
static double _mvt_error[3];

/////////////////////////////////////////////////////////////////////////////////////////////////
/* Function that calculates normal cumulative distribution function */
double pn(double y, double mu, double sigma) {
  return(Rf_pnorm5(y,mu,sqrt(sigma),1,0));
}

// [[Rcpp::export]]
double pn(arma::mat y, arma::mat mu, arma::mat sigma) {
  /*int n = y.n_rows;*/
  int k = y.n_elem;
  double res;
  if (k==1) {
    res = Rf_pnorm5(y[0],mu[0],sqrt(sigma[0]),1,0);
    return(res);
  }

  arma::mat L = arma::mat(2,2); L.fill(0.0);
  L(0,0) = 1/sqrt(sigma(0,0));
  L(1,1) = 1/sqrt(sigma(1,1));
  y = L*(y-mu);
  double r = sigma(0,1)/(sqrt(sigma(0,0)*sigma(1,1)));

  int rand = 1;
  Rcpp::IntegerVector infin(k); // Infinity argument (all 0 since CDF)
  Rcpp::NumericVector _mvt_delta(k); // non-centrality parameter

  for (int i=0; i<k; i++) {
    infin[i]=0;
    _mvt_delta[i]=0;
  }

  double val;
  mvtnorm_C_mvtdst(
		   &k,            // dim
		   &_mvt_df,      // df
		   &y[0],        // lower
		   &y[0],        // upper
		   &infin[0],     // integration type
		   &r,            // correlation
		   &_mvt_delta[0],// non-centrality
		   &_mvt_maxpts,
		   &_mvt_abseps,
		   &_mvt_releps,
		   &_mvt_error[0],
		   &val,
		   &_mvt_inform,
		   &rand);
  res = val;
  return(res);
}
