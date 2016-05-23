// [[Rcpp::depends(RcppArmadillo,mvtnorm)]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <mvtnormAPI.h>

using namespace Rcpp;
using namespace arma;

int _mvt_df = 0;
double _mvt_abseps=0.0001;
double _mvt_releps=0;
int _mvt_maxpts=25000;
int _mvt_inform;
double _mvt_error[3];


/*
Funtion that calculates cumulative distribution function for 

1) zero-mean bivariate normal distribution (r:=vector of correlation
coefficients corresponding to the rows of y)

2) zero-mean univariate normal distribution (r:=vector of standard
deviations)
 */
// [[Rcpp::export]]
vec pn(mat y, vec r) {
  int n = y.n_rows;
  int k = y.n_cols;
  vec res(n);

  if (k==1) {    
    for (int i=0; i<n; i++) {
      res(i) = Rf_pnorm5(y[i],0.0,r[i],1,0);
    }
    return(res);
  }
  
  int rand = 1;    
  IntegerVector infin(k); // Infinity argument (all 0 since CDF)
  NumericVector _mvt_delta(k); // non-centrality parameter
  for (int i=0; i<k; i++) {
    infin[i]=0; 
    _mvt_delta[i]=0;
  }
  
  for (int i = 0; i < n; i++) {
    rowvec y0 = y.row(i);
    double val;
    mvtnorm_C_mvtdst(
    		     &k,            // dim
    		     &_mvt_df,      // df
    		     &y0[0],        // lower
    		     &y0[0],        // upper
    		     &infin[0],     // integration type
    		     &r[i],         // correlation
    		     &_mvt_delta[0],// non-centrality
    		     &_mvt_maxpts,
    		     &_mvt_abseps,
    		     &_mvt_releps,
    		     &_mvt_error[0],
    		     &val,
    		     &_mvt_inform,
    		     &rand);
    res[i] = val;
  }
  return res;
}


/*
Small example calling 'pn' 
*/
// [[Rcpp::export]]
NumericVector loglik_ex(mat y, vec theta) {
  int n = y.n_rows;
  /* int k = y.n_cols; */  
  Rcpp::Rcout << "theta=" << theta << std::endl << y;
  vec r(1);
 
  NumericVector res(n);
  for (int i=0; i<n; i++) {
    r(0) = theta[0];
    mat rr = y.row(i);
    res[i] = log(pn(y.row(i),r)[0]);
  }
  
  return res;  
}

