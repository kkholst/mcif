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

struct vecmat
{
  vec V;
  mat M1;
  mat M2;
};

const double twopi = 2*datum::pi;
const double h = 1e-8;

/*
Funtion that calculates cumulative distribution function for

1) zero-mean bivariate normal distribution (r:=vector of correlation
coefficients corresponding to the rows of y)

2) zero-mean univariate normal distribution (r:=vector of standard
deviations)
 */

// [[Rcpp::export]]
double pn(mat y, double r) {
  /*int n = y.n_rows;*/
  int k = y.n_elem;

  double res(1);

  if (k==1) {
    res = Rf_pnorm5(y[0],0.0,r,1,0);
    return(res);
  }

  int rand = 1;
  IntegerVector infin(k); // Infinity argument (all 0 since CDF)
  NumericVector _mvt_delta(k); // non-centrality parameter

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

/*
Funtion that calculates cumulative distribution function for

1) zero-mean bivariate normal distribution (r:=vector of correlation
coefficients corresponding to the rows of y)

2) zero-mean univariate normal distribution (r:=vector of standard
deviations)

[[Rcpp::export]]
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
*/


/*
Small example calling 'pn'
// [[Rcpp::export]]
NumericVector loglik_ex(mat y, vec theta) {
int n = y.n_rows;
// int k = y.n_cols;
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
*/

/*
Conditional mean and variance-covariance matrix
*/
vecmat conMuSig(mat sigma, vec mu, uvec rc1, uvec rc2) {

  mat sig11 = sigma.submat(rc1,rc1);
  mat sig12 = sigma.submat(rc1,rc2);
  mat sig21 = sigma.submat(rc2,rc1);
  mat sig22 = sigma.submat(rc2,rc2);

  mat isig22 = sig22.i();
  mat sigS = sig12*isig22;

  vec mu2 = mu.elem(rc2);

  vec c_mu = sig12*isig22*mu2;
  mat c_sig = sig11-sig12*isig22*sig21;

  vecmat out;
  out.V = c_mu;
  out.M1 = c_sig;
  out.M2 = sigS;

  return(out);
}

/* logfy */
// [[Rcpp::export]]
vec loglik(mat y, mat b, mat u, mat sigma, mat alph, mat dalph){
  /* y: nx2 matrix with event type (0, 1 or 2) of family member 1 and 2
     b: nx4 matrix with XB for event type 1 and 2 (b1 and b2) for family member 1 and 2
        the order is b1_1, b1_2, b2_1 and b2_2
     u: nx2 matrix with the random effects u1 and u2 affecting pi1 and pi2 (cluster-specific risk levels)
        the random effects are shared by family member 1 and 2
     sigma: 6x6 matrix. variance-covariance matrix, order: event1_1, event1_2, event2_1, event2_2, u1, u2
     alph: nx4 matrix, inside the Probit link, order a1_1, a1_2, a2_1, a2_2
     dalph: nx4 matrix, derivative of alph wrt. t, order da1_1, da1_2, da2_1, da2_2
  */

  /* No. of pairs */
  int n = y.n_rows;

  /* Estimatimg cluster-specific risk levels */
  vec pi1_1 = exp(b.col(0)+u.col(0))/(1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1)));
  vec pi1_2 = exp(b.col(1)+u.col(0))/(1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1)));
  vec pi2_1 = exp(b.col(2)+u.col(1))/(1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1)));
  vec pi2_2 = exp(b.col(3)+u.col(1))/(1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1)));

  /* Initialising loglik vector */
  vec res(n);

  for (int i=0; i<n; i++) {
    /* Both family members experience event 1, estimating ddF11 */
    if((y(i,0) == 1) & (y(i,1) == 1)){

      /* Specifying which parts of sigma apply */
      uvec rc1(2); rc1(0) = 0; rc1(1) = 1;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat ic_sig = c_sig.i(); // the inverse

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub(n,2);
      alph_sub.col(0) = alph.col(0); // alph1_1
      alph_sub.col(1) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 2xn */
      mat alph_f = alph_sub.t();

      /* Centering the alphas */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf11 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Loglikelihood contribution */
      double ddF11 = pi1_1(i)*pi1_2(i)*dalph(i,0)*dalph(i,1)*pdf11; // pi1_1, pi1_2, dalph1_1, dalph1_2
      res(i) = log(ddF11);
    }
    /* Family member 1 experience event 1, family member 2 experience event 2, estimating ddF12 */
    else if((y(i,0) == 1) & (y(i,1) == 2)){

      /* Specifying which parts of sigma apply */
      uvec rc1(2); rc1(0) = 0; rc1(1) = 3;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat ic_sig = c_sig.i(); // the inverse

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub(n,2);
      alph_sub.col(0) = alph.col(0); // alph1_1
      alph_sub.col(1) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 2xn */
      mat alph_f = alph_sub.t();

      /* Centering the alphas */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf12 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Loglikelihood contribution */
      double ddF12 = pi1_1(i)*pi2_2(i)*dalph(i,0)*dalph(i,3)*pdf12; // pi1_1, pi2_2, dalph1_1, dalph2_2
      res(i) = log(ddF12);
    }
    /* Family member 1 experience event 2, family member 2 experience event 1, estimating ddF21 */
    else if((y(i,0) == 2) & (y(i,1) == 1)){

      /* Specifying which parts of sigma apply */
      uvec rc1(2); rc1(0) = 2; rc1(1) = 1;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat ic_sig = c_sig.i(); // the inverse

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub(n,2);
      alph_sub.col(0) = alph.col(2); // alph2_1
      alph_sub.col(1) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 2xn */
      mat alph_f = alph_sub.t();

      /* Centering the alphas */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf21 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Loglikelihood contribution */
      double ddF21 = pi2_1(i)*pi1_2(i)*dalph(i,2)*dalph(i,1)*pdf21; // pi2_1, pi1_2, dalph2_1, dalph1_2
      res(i) = log(ddF21);
    }
    /* Both family members experience event 2, estimating ddF22 */
    else if((y(i,0) == 2) & (y(i,1) == 2)){

      /* Specifying which parts of sigma apply */
      uvec rc1(2); rc1(0) = 2; rc1(1) = 3;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat ic_sig = c_sig.i(); // the inverse

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub(n,2);
      alph_sub.col(0) = alph.col(2); // alph2_1
      alph_sub.col(1) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 2xn */
      mat alph_f = alph_sub.t();

      /* Centering the alphas */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf22 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Loglikelihood contribution */
      double ddF22 = pi2_1(i)*pi2_2(i)*dalph(i,2)*dalph(i,3)*pdf22; // pi2_1, pi2_2, dalph2_1, dalph2_2
      res(i) = log(ddF22);
    }
    /* Family member 1 experience event 0, family member 2 experience event 1, estimating dF01 */
    else if((y(i,0) == 0) & (y(i,1) == 1)){

      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 0;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 1 ; rc4(1) = 4; rc4(2) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal dF1_2 */
      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub(n,1);
      alph_sub.col(0) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f = alph_sub.t();

      /* Centering the alpha */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf1_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Estimation of dF1_2 */
      double dF1_2 = pi1_2(i)*dalph(i,1)*pdf1_2; // pi1_2, dalph1_2

      /* Conditional F1c1_1 */
      /* Conditional mean and variance-covariance matrix, conditional on alph1_2, u1 and u2 */
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub2(n,1);
      alph_sub2.col(0) = alph.col(0); // alph1_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f2 = alph_sub2.t();

      /* Centering the alpha */
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      /* Estimation of F1c1_1 */
      double F1c1_1 = pi1_1(i)*pn(alph_c2,sqrt(c_sig2[0]));

      /* Conditional F2c1_1 */
      /* Conditional mean and variance-covariance matrix, conditional on alph1_2, u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub3(n,1);
      alph_sub3.col(0) = alph.col(2); // alph2_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f3 = alph_sub3.t();

      /* Centering the alpha */
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      /* Estimation of F2c1_1 */
      double F2c1_1 = pi2_1(i)*pn(alph_c3,sqrt(c_sig3[0]));

      /* Loglikelihood contribution */
      double dF01 = dF1_2*(1-F1c1_1-F2c1_1);
      res(i) = log(dF01);
    }
    /* Family member 1 experience event 1, family member 2 experience event 0, estimating dF10 */
    else if((y(i,0) == 1) & (y(i,1) == 0)){

      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 1;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 0 ; rc4(1) = 4; rc4(2) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal dF1_1 */
      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub(n,1);
      alph_sub.col(0) = alph.col(0); // alph1_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f = alph_sub.t();

      /* Centering the alpha */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf1_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Estimation of dF1_1 */
      double dF1_1 = pi1_1(i)*dalph(i,0)*pdf1_1; // pi1_1, dalph1_1

      /* Conditional F1c1_2 */
      /* Conditional mean and variance-covariance matrix, conditional on alph1_1, u1 and u2 */
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub2(n,1);
      alph_sub2.col(0) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f2 = alph_sub2.t();

      /* Centering the alpha */
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      /* Estimation of F1c1_2 */
      double F1c1_2 = pi1_2(i)*pn(alph_c2,sqrt(c_sig2[0]));

      /* Conditional F2c1_2 */
      /* Conditional mean and variance-covariance matrix, conditional on alph1_1, u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub3(n,1);
      alph_sub3.col(0) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f3 = alph_sub3.t();

      /* Centering the alpha */
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      /* Estimation of F2c1_2 */
      double F2c1_2 = pi2_2(i)*pn(alph_c3,sqrt(c_sig3[0]));

      /* Loglikelihood contribution */
      double dF10 = dF1_1*(1-F1c1_2-F2c1_2);
      res(i) = log(dF10);
    }
    /* Family member 1 experience event 0, family member 2 experience event 2, estimating dF02 */
    else if((y(i,0) == 0) & (y(i,1) == 2)){

      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 0;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 3 ; rc4(1) = 4; rc4(2) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal dF2_2 */
      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc2, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub(n,1);
      alph_sub.col(0) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f = alph_sub.t();

      /* Centering the alpha */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf2_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Estimation of dF2_2 */
      double dF2_2 = pi2_2(i)*dalph(i,3)*pdf2_2; // pi2_2, dalph2_2

      /* Conditional F1c2_1 */
      /* Conditional mean and variance-covariance matrix, conditional on alph2_2, u1 and u2 */
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub2(n,1);
      alph_sub2.col(0) = alph.col(0); // alph1_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f2 = alph_sub2.t();

      /* Centering the alpha */
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      /* Estimation of F1c2_1 */
      double F1c2_1 = pi1_1(i)*pn(alph_c2,sqrt(c_sig2[0]));

      /* Conditional F2c2_1 */
      /* Conditional mean and variance-covariance matrix, conditional on alph2_2, u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub3(n,1);
      alph_sub3.col(0) = alph.col(2); // alph2_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f3 = alph_sub3.t();

      /* Centering the alpha */
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      /* Estimation of F2c2_1 */
      double F2c2_1 = pi2_1(i)*pn(alph_c3,sqrt(c_sig3[0]));

      /* Loglikelihood contribution */
      double dF02 = dF2_2*(1-F1c2_1-F2c2_1);
      res(i) = log(dF02);
    }
    /* Family member 1 experience event 2, family member 2 experience event 0, estimating dF20 */
    else if((y(i,0) == 2) & (y(i,1) == 0)){

      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 1;
      uvec rc2(1); rc2(0) = 3;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 2 ; rc4(1) = 4; rc4(2) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal dF2_1 */
      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc2, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub(n,1);
      alph_sub.col(0) = alph.col(2); // alph2_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f = alph_sub.t();

      /* Centering the alpha */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf2_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Estimation of dF2_1 */
      double dF2_1 = pi2_1(i)*dalph(i,2)*pdf2_1; // pi2_1, dalph2_1

      /* Conditional F1c2_2 */
      /* Conditional mean and variance-covariance matrix, conditional on alph2_1, u1 and u2 */
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub2(n,1);
      alph_sub2.col(0) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f2 = alph_sub2.t();

      /* Centering the alpha */
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      /* Estimation of F1c2_2 */
      double F1c2_2 = pi1_2(i)*pn(alph_c2,sqrt(c_sig2[0]));

      /* Conditional F2c2_2 */
      /* Conditional mean and variance-covariance matrix, conditional on alph2_1, u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub3(n,1);
      alph_sub3.col(0) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f3 = alph_sub3.t();

      /* Centering the alpha */
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      /* Estimation of F2c2_2 */
      double F2c2_2 = pi2_2(i)*pn(alph_c3,sqrt(c_sig3[0]));

      /* Loglikelihood contribution */
      double dF20 = dF2_1*(1-F1c2_2-F2c2_2);
      res(i) = log(dF20);
    }
    /* Family member 1 experience event 0, family member 2 experience event 0, estimating F00 */
    else{
      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 0;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;

      uvec rc4(2); rc4(0) = 0 ; rc4(1) = 1;
      uvec rc5(2); rc5(0) = 0 ; rc5(1) = 3;
      uvec rc6(2); rc6(0) = 2 ; rc6(1) = 1;
      uvec rc7(2); rc7(0) = 2 ; rc7(1) = 3;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal F1_1, F1_2, F2_1 and F2_2 */
      /* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
      vecmat out1 = conMuSig(sigma, mu, rc1, rc3);
      vec c_mu1 = out1.V;
      mat c_sig1 = out1.M1;

      vecmat out2 = conMuSig(sigma, mu, rc2, rc3);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub1_1(n,1);
      mat alph_sub1_2(n,1);
      mat alph_sub2_1(n,1);
      mat alph_sub2_2(n,1);
      alph_sub1_1.col(0) = alph.col(0); // alph1_1
      alph_sub1_2.col(0) = alph.col(1); // alph1_2
      alph_sub2_1.col(0) = alph.col(2); // alph2_1
      alph_sub2_2.col(0) = alph.col(3); // alph2_2

      /* Transposing matrices to be of the form 1xn */
      mat alph_f1_1 = alph_sub1_1.t();
      mat alph_f1_2 = alph_sub1_2.t();
      mat alph_f2_1 = alph_sub2_1.t();
      mat alph_f2_2 = alph_sub2_2.t();

      /* Centering the alphas */
      mat alph_c1_1 = alph_f1_1.col(i)-c_mu1;
      mat alph_c1_2 = alph_f1_2.col(i)-c_mu1;
      mat alph_c2_1 = alph_f2_1.col(i)-c_mu2;
      mat alph_c2_2 = alph_f2_2.col(i)-c_mu2;

      /* Estimation of F1_1, F1_2, F2_1 and F2_2 */
      double F1_1 = pi1_1(i)*pn(alph_c1_1,sqrt(c_sig1[0]));
      double F1_2 = pi1_2(i)*pn(alph_c1_2,sqrt(c_sig1[0]));
      double F2_1 = pi2_1(i)*pn(alph_c2_1,sqrt(c_sig2[0]));
      double F2_2 = pi2_2(i)*pn(alph_c2_2,sqrt(c_sig2[0]));

      /* Joint probabilities F11, F12, F21 and F22 */
      /* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc4, rc3);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;

      vecmat out4 = conMuSig(sigma, mu, rc5, rc3);
      vec c_mu4 = out4.V;
      mat c_sig4 = out4.M1;

      vecmat out5 = conMuSig(sigma, mu, rc6, rc3);
      vec c_mu5 = out5.V;
      mat c_sig5 = out5.M1;

      vecmat out6 = conMuSig(sigma, mu, rc7, rc3);
      vec c_mu6 = out6.V;
      mat c_sig6 = out6.M1;

      /* Correlation coefficient */
      double r1 = c_sig3(0,1)/(sqrt(c_sig3(0,0))*sqrt(c_sig3(1,1)));
      double r2 = c_sig4(0,1)/(sqrt(c_sig4(0,0))*sqrt(c_sig4(1,1)));
      double r3 = c_sig5(0,1)/(sqrt(c_sig5(0,0))*sqrt(c_sig5(1,1)));
      double r4 = c_sig6(0,1)/(sqrt(c_sig6(0,0))*sqrt(c_sig6(1,1)));

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub11(n,2);
      alph_sub11.col(0) = alph.col(0); // alph1_1
      alph_sub11.col(1) = alph.col(1); // alph1_2

      mat alph_sub12(n,2);
      alph_sub12.col(0) = alph.col(0); // alph1_1
      alph_sub12.col(1) = alph.col(3); // alph2_2

      mat alph_sub21(n,2);
      alph_sub21.col(0) = alph.col(2); // alph2_1
      alph_sub21.col(1) = alph.col(1); // alph1_2

      mat alph_sub22(n,2);
      alph_sub22.col(0) = alph.col(2); // alph2_1
      alph_sub22.col(1) = alph.col(3); // alph2_2

      /* Transposing matrices to be of the form 2xn */
      mat alph_f11 = alph_sub11.t();
      mat alph_f12 = alph_sub12.t();
      mat alph_f21 = alph_sub21.t();
      mat alph_f22 = alph_sub22.t();

      /* Centering the alphas */
      mat alph_c11 = alph_f11.col(i)-c_mu3;
      mat alph_c12 = alph_f12.col(i)-c_mu4;
      mat alph_c21 = alph_f21.col(i)-c_mu5;
      mat alph_c22 = alph_f22.col(i)-c_mu6;

      /* Estimating F11, F12, F21 and F22 */
      double F11 = pi1_1(i)*pi1_2(i)*pn(alph_c11,r1);
      double F12 = pi1_1(i)*pi2_2(i)*pn(alph_c12,r2);
      double F21 = pi2_1(i)*pi1_2(i)*pn(alph_c21,r3);
      double F22 = pi2_1(i)*pi2_2(i)*pn(alph_c22,r4);

      /* Calculating F01, F10, F02 and F20 */
      double F01 = F1_2-(F11+F21);
      double F10 = F1_1-(F11+F12);
      double F02 = F2_2-(F12+F22);
      double F20 = F2_1-(F21+F22);

      /* Loglikelihood contribution */
      double F00 = 1-(F11+F12+F21+F22+F01+F10+F02+F20);
      res(i) = log(F00);
    }
  }
  return(res);
}


//Rcpp::Rcout << "a12" << std::endl << alph_c12 ;
//Rcpp::Rcout << "everywhere";

/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/* Dlogfy */
// [[Rcpp::export]]
mat Dloglik(mat y, mat b, mat u, mat sigma, mat alph, mat dalph){
  /* y: nx2 matrix with event type (0, 1 or 2) of family member 1 and 2
     b: nx4 matrix with XB for event type 1 and 2 (b1 and b2) for family member 1 and 2
        the order is b1_1, b1_2, b2_1 and b2_2
     u: nx2 matrix with the random effects u1 and u2 affecting pi1 and pi2 (cluster-specific risk levels)
        the random effects are shared by family member 1 and 2
     sigma: 6x6 matrix. variance-covariance matrix, order: event1_1, event1_2, event2_1, event2_2, u1, u2
     alph: nx4 matrix, inside the Probit link, order a1_1, a1_2, a2_1, a2_2
     dalph: nx4 matrix, derivative of alph wrt. t, order da1_1, da1_2, da2_1, da2_2
  */

  /* No. of pairs */
  int n = y.n_rows;

  /* Estimating cluster-specific risk levels */
  vec pi1_1 = exp(b.col(0)+u.col(0))/(1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1)));
  vec pi1_2 = exp(b.col(1)+u.col(0))/(1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1)));
  vec pi2_1 = exp(b.col(2)+u.col(1))/(1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1)));
  vec pi2_2 = exp(b.col(3)+u.col(1))/(1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1)));

  /* Derivatives of the pis wrt. u1 and u2 */
  vec dpi1_1_u1 = (exp(b.col(0)+u.col(0)+b.col(2)+u.col(1))+exp(b.col(0)+u.col(0)))/pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2);
  vec dpi1_1_u2 = -exp(b.col(0)+u.col(0)+b.col(2)+u.col(1))/pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2);
  vec dpi2_1_u1 = -exp(b.col(0)+u.col(0)+b.col(2)+u.col(1))/pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2);
  vec dpi2_1_u2 = (exp(b.col(0)+u.col(0)+b.col(2)+u.col(1))+exp(b.col(2)+u.col(1)))/pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2);

  vec dpi1_2_u1 = (exp(b.col(1)+u.col(0)+b.col(3)+u.col(1))+exp(b.col(1)+u.col(0)))/pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2);
  vec dpi1_2_u2 = -exp(b.col(1)+u.col(0)+b.col(3)+u.col(1))/pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2);
  vec dpi2_2_u1 = -exp(b.col(1)+u.col(0)+b.col(3)+u.col(1))/pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2);
  vec dpi2_2_u2 = (exp(b.col(1)+u.col(0)+b.col(3)+u.col(1))+exp(b.col(3)+u.col(1)))/pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2);

  /* Derivatives of the product of two pis wrt. u1 and u2 */
  vec dpi11_u1 = (exp(b.col(0)+3*u.col(0)+2*b.col(1)+b.col(2)+u.col(1))+exp(2*b.col(0)+3*u.col(0)+b.col(1)+b.col(3)+u.col(1))+2*exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(3)+2*u.col(1)+b.col(2))+exp(b.col(0)+3*u.col(0)+2*b.col(1))+2*exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(3)+u.col(1))+exp(2*b.col(0)+3*u.col(0)+b.col(1))+2*exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(2)+u.col(1))+2*exp(b.col(0)+2*u.col(0)+b.col(1)))/(pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2)%pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2));

  vec dpi11_u2 = -(exp(b.col(0)+3*u.col(0)+2*b.col(1)+b.col(2)+u.col(1))+2*exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(3)+2*u.col(1)+b.col(2))+exp(2*b.col(0)+3*u.col(0)+b.col(1)+b.col(3)+u.col(1))+exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(2)+u.col(1))+exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(3)+u.col(1)))/(pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2)%pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2));

  vec dpi12_u1 = (-exp(2*b.col(0)+3*u.col(0)+b.col(1)+b.col(3)+u.col(1))+exp(b.col(0)+u.col(0)+2*b.col(3)+3*u.col(1)+b.col(2))+exp(b.col(0)+u.col(0)+2*b.col(3)+2*u.col(1))+exp(b.col(0)+u.col(0)+b.col(3)+2*u.col(1)+b.col(2))+exp(b.col(0)+u.col(0)+b.col(3)+u.col(1)))/(pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2)%pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2));

  vec dpi12_u2 = -(-exp(2*b.col(0)+3*u.col(0)+b.col(1)+b.col(3)+u.col(1))+exp(b.col(0)+u.col(0)+2*b.col(3)+3*u.col(1)+b.col(2))-exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(3)+u.col(1))-exp(2*b.col(0)+2*u.col(0)+b.col(3)+u.col(1))-exp(b.col(0)+u.col(0)+b.col(3)+u.col(1)))/(pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2)%pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2));

  vec dpi21_u1 = (-exp(b.col(0)+3*u.col(0)+2*b.col(1)+b.col(2)+u.col(1))+exp(2*b.col(2)+3*u.col(1)+b.col(3)+b.col(1)+u.col(0))+exp(2*b.col(2)+2*u.col(1)+b.col(1)+u.col(0))+exp(b.col(3)+2*u.col(1)+b.col(2)+b.col(1)+u.col(0))+exp(b.col(2)+u.col(1)+b.col(1)+u.col(0)))/(pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2)%pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2));

  vec dpi21_u2 = (exp(b.col(0)+3*u.col(0)+2*b.col(1)+b.col(2)+u.col(1))-exp(2*b.col(2)+3*u.col(1)+b.col(3)+b.col(1)+u.col(0))+exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(2)+u.col(1))+exp(b.col(2)+u.col(1)+2*b.col(1)+2*u.col(0))+exp(b.col(2)+u.col(1)+b.col(1)+u.col(0)))/(pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2)%pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2));

  vec dpi22_u1 = -(2*exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(3)+2*u.col(1)+b.col(2))+exp(b.col(0)+u.col(0)+2*b.col(3)+3*u.col(1)+b.col(2))+exp(2*b.col(2)+3*u.col(1)+b.col(3)+b.col(1)+u.col(0))+exp(b.col(0)+u.col(0)+b.col(3)+2*u.col(1)+b.col(2))+exp(b.col(3)+2*u.col(1)+b.col(2)+b.col(1)+u.col(0)))/(pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2)%pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2));

  vec dpi22_u2 = (2*exp(b.col(0)+2*u.col(0)+b.col(1)+b.col(3)+2*u.col(1)+b.col(2))+exp(b.col(0)+u.col(0)+2*b.col(3)+3*u.col(1)+b.col(2))+exp(2*b.col(2)+3*u.col(1)+b.col(3)+b.col(1)+u.col(0))+2*exp(b.col(0)+u.col(0)+b.col(3)+2*u.col(1)+b.col(2))+exp(2*b.col(2)+3*u.col(1)+b.col(3))+2*exp(b.col(3)+2*u.col(1)+b.col(2)+b.col(1)+u.col(0))+exp(b.col(2)+3*u.col(1)+2*b.col(3))+2*exp(b.col(2)+2*u.col(1)+b.col(3)))/(pow((1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1))),2)%pow((1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1))),2));

  /* Initialising Dloglik matrix */
  mat res(n,2);

  for (int i=0; i<n; i++) {
    /* Both family members experience event 1, estimating score contribution from ddF11 */
    if((y(i,0) == 1) & (y(i,1) == 1)){

      /* Specifying which parts of sigma apply */
      uvec rc1(2); rc1(0) = 0; rc1(1) = 1;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat c_sigX = out.M2;

      /* Correlation */
      double r = c_sig(0,1)/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1)));

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub(n,2);
      alph_sub.col(0) = alph.col(0); // alph1_1
      alph_sub.col(1) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 2xn */
      mat alph_f = alph_sub.t();

      /* Centering the alphas */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf11 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Derivative of the pdf wrt. u1 and u2 */
      double dpdf11_u1 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0))*(-1/(2*pow(1-r,2)))*((2*(pow(c_sigX(0,0),2)*sqrt(c_sig(1,1))*u(i,0)+c_sigX(0,0)*c_sigX(0,1)*sqrt(c_sig(1,1))*u(i,1)-2*c_sigX(0,0)*c_sigX(1,0)*r*u(i,0)-c_sigX(0,0)*c_sigX(1,1)*r*u(i,1)-c_sigX(0,1)*c_sigX(1,0)*r*u(i,1)+pow(c_sigX(1,0),2)*sqrt(c_sig(0,0))*u(i,0)+c_sigX(1,0)*c_sigX(1,1)*sqrt(c_sig(0,0))*u(i,1)+c_sigX(0,0)*r*alph(i,1)-c_sigX(0,0)*sqrt(c_sig(1,1))*alph(i,0)+c_sigX(1,0)*r*alph(i,0)-c_sigX(1,0)*sqrt(c_sig(0,0))*alph(i,1)))/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1))));

      double dpdf11_u2 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0))*(-1/(2*pow(1-r,2)))*((2*(c_sigX(0,0)*c_sigX(0,1)*sqrt(c_sig(1,1))*u(i,0)-c_sigX(0,0)*c_sigX(1,1)*r*u(i,0)+pow(c_sigX(0,1),2)*sqrt(c_sig(1,1))*u(i,1)-c_sigX(0,1)*c_sigX(1,0)*r*u(i,0)-2*c_sigX(0,1)*c_sigX(1,1)*r*u(i,1)+c_sigX(1,0)*c_sigX(1,1)*sqrt(c_sig(0,0))*u(i,0)+pow(c_sigX(1,1),2)*sqrt(c_sig(0,0))*u(i,1)+c_sigX(0,1)*r*alph(i,1)-c_sigX(0,1)*sqrt(c_sig(1,1))*alph(i,0)+c_sigX(1,1)*r*alph(i,0)-c_sigX(1,1)*sqrt(c_sig(0,0))*alph(i,1)))/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1))));

      /* Likelihood contribution */
      double ddF11 = pi1_1(i)*pi1_2(i)*dalph(i,0)*dalph(i,1)*pdf11; // pi1_1, pi1_2, dalph1_1, dalph1_2

      /* Score contributions from ddF11 wrt. u1 and u2 */
      double sc_u1 = (1/ddF11)*dalph(i,0)*dalph(i,1)*(dpi11_u1(i)*pdf11+pi1_1(i)*pi1_2(i)*dpdf11_u1);
      double sc_u2 = (1/ddF11)*dalph(i,0)*dalph(i,1)*(dpi11_u2(i)*pdf11+pi1_1(i)*pi1_2(i)*dpdf11_u2);

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;

    }
    /* Family member 1 experiences event 1, family member 2 experiences event 2, estimating score contribution from ddF12 */
    else if((y(i,0) == 1) & (y(i,1) == 2)){

      /* Specifying which parts of sigma apply */
      uvec rc1(2); rc1(0) = 0; rc1(1) = 3;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat c_sigX = out.M2;

      /* Correlation */
      double r = c_sig(0,1)/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1)));

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub(n,2);
      alph_sub.col(0) = alph.col(0); // alph1_1
      alph_sub.col(1) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 2xn */
      mat alph_f = alph_sub.t();

      /* Centering the alphas */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf12 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Derivative of the pdf wrt. u1 and u2 */
      double dpdf12_u1 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0))*(-1/(2*pow(1-r,2)))*((2*(pow(c_sigX(0,0),2)*sqrt(c_sig(1,1))*u(i,0)+c_sigX(0,0)*c_sigX(0,1)*sqrt(c_sig(1,1))*u(i,1)-2*c_sigX(0,0)*c_sigX(1,0)*r*u(i,0)-c_sigX(0,0)*c_sigX(1,1)*r*u(i,1)-c_sigX(0,1)*c_sigX(1,0)*r*u(i,1)+pow(c_sigX(1,0),2)*sqrt(c_sig(0,0))*u(i,0)+c_sigX(1,0)*c_sigX(1,1)*sqrt(c_sig(0,0))*u(i,1)+c_sigX(0,0)*r*alph(i,3)-c_sigX(0,0)*sqrt(c_sig(1,1))*alph(i,0)+c_sigX(1,0)*r*alph(i,0)-c_sigX(1,0)*sqrt(c_sig(0,0))*alph(i,3)))/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1))));

      double dpdf12_u2 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0))*(-1/(2*pow(1-r,2)))*((2*(c_sigX(0,0)*c_sigX(0,1)*sqrt(c_sig(1,1))*u(i,0)-c_sigX(0,0)*c_sigX(1,1)*r*u(i,0)+pow(c_sigX(0,1),2)*sqrt(c_sig(1,1))*u(i,1)-c_sigX(0,1)*c_sigX(1,0)*r*u(i,0)-2*c_sigX(0,1)*c_sigX(1,1)*r*u(i,1)+c_sigX(1,0)*c_sigX(1,1)*sqrt(c_sig(0,0))*u(i,0)+pow(c_sigX(1,1),2)*sqrt(c_sig(0,0))*u(i,1)+c_sigX(0,1)*r*alph(i,3)-c_sigX(0,1)*sqrt(c_sig(1,1))*alph(i,0)+c_sigX(1,1)*r*alph(i,0)-c_sigX(1,1)*sqrt(c_sig(0,0))*alph(i,3)))/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1))));

      /* Likelihood contribution */
      double ddF12 = pi1_1(i)*pi2_2(i)*dalph(i,0)*dalph(i,3)*pdf12; // pi1_1, pi2_2, dalph1_1, dalph2_2

      /* Score contributions from ddF12 wrt. u1 and u2 */
      double sc_u1 = (1/ddF12)*dalph(i,0)*dalph(i,3)*(dpi12_u1(i)*pdf12+pi1_1(i)*pi2_2(i)*dpdf12_u1);
      double sc_u2 = (1/ddF12)*dalph(i,0)*dalph(i,3)*(dpi12_u2(i)*pdf12+pi1_1(i)*pi2_2(i)*dpdf12_u2);

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;

    }
    /* Family member 1 experiences event 2, family member 2 experiences event 1, estimating score contribution from ddF21 */
    else if((y(i,0) == 2) & (y(i,1) == 1)){

      /* Specifying which parts of sigma apply */
      uvec rc1(2); rc1(0) = 2; rc1(1) = 1;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat c_sigX = out.M2;
      mat ic_sig = c_sig.i(); // the inverse

      /* Correlation */
      double r = c_sig(0,1)/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1)));

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub(n,2);
      alph_sub.col(0) = alph.col(2); // alph2_1
      alph_sub.col(1) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 2xn */
      mat alph_f = alph_sub.t();

      /* Centering the alphas */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf21 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Derivative of the pdf wrt. u1 and u2 */
      double dpdf21_u1 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0))*(-1/(2*pow(1-r,2)))*((2*(pow(c_sigX(0,0),2)*sqrt(c_sig(1,1))*u(i,0)+c_sigX(0,0)*c_sigX(0,1)*sqrt(c_sig(1,1))*u(i,1)-2*c_sigX(0,0)*c_sigX(1,0)*r*u(i,0)-c_sigX(0,0)*c_sigX(1,1)*r*u(i,1)-c_sigX(0,1)*c_sigX(1,0)*r*u(i,1)+pow(c_sigX(1,0),2)*sqrt(c_sig(0,0))*u(i,0)+c_sigX(1,0)*c_sigX(1,1)*sqrt(c_sig(0,0))*u(i,1)+c_sigX(0,0)*r*alph(i,1)-c_sigX(0,0)*sqrt(c_sig(1,1))*alph(i,2)+c_sigX(1,0)*r*alph(i,2)-c_sigX(1,0)*sqrt(c_sig(0,0))*alph(i,1)))/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1))));

      double dpdf21_u2 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0))*(-1/(2*pow(1-r,2)))*((2*(c_sigX(0,0)*c_sigX(0,1)*sqrt(c_sig(1,1))*u(i,0)-c_sigX(0,0)*c_sigX(1,1)*r*u(i,0)+pow(c_sigX(0,1),2)*sqrt(c_sig(1,1))*u(i,1)-c_sigX(0,1)*c_sigX(1,0)*r*u(i,0)-2*c_sigX(0,1)*c_sigX(1,1)*r*u(i,1)+c_sigX(1,0)*c_sigX(1,1)*sqrt(c_sig(0,0))*u(i,0)+pow(c_sigX(1,1),2)*sqrt(c_sig(0,0))*u(i,1)+c_sigX(0,1)*r*alph(i,1)-c_sigX(0,1)*sqrt(c_sig(1,1))*alph(i,2)+c_sigX(1,1)*r*alph(i,2)-c_sigX(1,1)*sqrt(c_sig(0,0))*alph(i,1)))/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1))));

      /* Likelihood contribution */
      double ddF21 = pi2_1(i)*pi1_2(i)*dalph(i,2)*dalph(i,1)*pdf21; // pi2_1, pi1_2, dalph2_1, dalph1_2

      /* Score contributions from ddF21 wrt. u1 and u2 */
      double sc_u1 = (1/ddF21)*dalph(i,2)*dalph(i,1)*(dpi21_u1(i)*pdf21+pi2_1(i)*pi1_2(i)*dpdf21_u1);
      double sc_u2 = (1/ddF21)*dalph(i,2)*dalph(i,1)*(dpi21_u2(i)*pdf21+pi2_1(i)*pi1_2(i)*dpdf21_u2);

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;

    }
/* Both family members experience event 2, estimating score contribution from ddF22 */
    else if((y(i,0) == 2) & (y(i,1) == 2)){

      /* Specifying which parts of sigma apply */
      uvec rc1(2); rc1(0) = 2; rc1(1) = 3;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat c_sigX = out.M2;
      mat ic_sig = c_sig.i(); // the inverse

      /* Correlation */
      double r = c_sig(0,1)/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1)));

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub(n,2);
      alph_sub.col(0) = alph.col(2); // alph2_1
      alph_sub.col(1) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 2xn */
      mat alph_f = alph_sub.t();

      /* Centering the alphas */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf22 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Derivative of the pdf wrt. u1 and u2 */
      double dpdf22_u1 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0))*(-1/(2*pow(1-r,2)))*((2*(pow(c_sigX(0,0),2)*sqrt(c_sig(1,1))*u(i,0)+c_sigX(0,0)*c_sigX(0,1)*sqrt(c_sig(1,1))*u(i,1)-2*c_sigX(0,0)*c_sigX(1,0)*r*u(i,0)-c_sigX(0,0)*c_sigX(1,1)*r*u(i,1)-c_sigX(0,1)*c_sigX(1,0)*r*u(i,1)+pow(c_sigX(1,0),2)*sqrt(c_sig(0,0))*u(i,0)+c_sigX(1,0)*c_sigX(1,1)*sqrt(c_sig(0,0))*u(i,1)+c_sigX(0,0)*r*alph(i,3)-c_sigX(0,0)*sqrt(c_sig(1,1))*alph(i,2)+c_sigX(1,0)*r*alph(i,2)-c_sigX(1,0)*sqrt(c_sig(0,0))*alph(i,3)))/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1))));

      double dpdf22_u2 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0))*(-1/(2*pow(1-r,2)))*((2*(c_sigX(0,0)*c_sigX(0,1)*sqrt(c_sig(1,1))*u(i,0)-c_sigX(0,0)*c_sigX(1,1)*r*u(i,0)+pow(c_sigX(0,1),2)*sqrt(c_sig(1,1))*u(i,1)-c_sigX(0,1)*c_sigX(1,0)*r*u(i,0)-2*c_sigX(0,1)*c_sigX(1,1)*r*u(i,1)+c_sigX(1,0)*c_sigX(1,1)*sqrt(c_sig(0,0))*u(i,0)+pow(c_sigX(1,1),2)*sqrt(c_sig(0,0))*u(i,1)+c_sigX(0,1)*r*alph(i,3)-c_sigX(0,1)*sqrt(c_sig(1,1))*alph(i,2)+c_sigX(1,1)*r*alph(i,2)-c_sigX(1,1)*sqrt(c_sig(0,0))*alph(i,3)))/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1))));

      /* Likelihood contribution */
      double ddF22 = pi2_1(i)*pi2_2(i)*dalph(i,2)*dalph(i,3)*pdf22; // pi2_1, pi2_2, dalph2_1, dalph2_2

      /* Score contributions from ddF22 wrt. u1 and u2 */
      double sc_u1 = (1/ddF22)*dalph(i,2)*dalph(i,3)*(dpi22_u1(i)*pdf22+pi2_1(i)*pi2_2(i)*dpdf22_u1);
      double sc_u2 = (1/ddF22)*dalph(i,2)*dalph(i,3)*(dpi22_u2(i)*pdf22+pi2_1(i)*pi2_2(i)*dpdf22_u2);

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;
    }
    // /* Family member 1 experience event 0, family member 2 experience event 1, estimating score contribution from dF01 */
    else if((y(i,0) == 0) & (y(i,1) == 1)){

      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 0;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 1 ; rc4(1) = 4; rc4(2) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal dF1_2 */
      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat c_sigX = out.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub(n,1);
      alph_sub.col(0) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f = alph_sub.t();

      /* Centering the alpha */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf1_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig[0]))*exp(-0.5*inner(0));

      /* Estimation of dF1_2 */
      double dF1_2 = pi1_2(i)*dalph(i,1)*pdf1_2; // pi1_2, dalph1_2, pdf1_2

      /* Score contribution from dF1_1 */
      double ddF1_2_u1 = dalph(i,1)*(dpi1_2_u1(i)*pdf1_2+pi1_2(i)*pdf1_2*(-1/(2*c_sig[0]))*2*(c_sigX(0)*u(i,0)+c_sigX(1)*u(i,1)-alph(i,1))*c_sigX(0));

      double ddF1_2_u2 = dalph(i,1)*(dpi1_2_u2(i)*pdf1_2+pi1_2(i)*pdf1_2*(-1/(2*c_sig[0]))*2*(c_sigX(0)*u(i,0)+c_sigX(1)*u(i,1)-alph(i,1))*c_sigX(1));

      /* Conditional F1c1_1 */
      /* Conditional mean and variance-covariance matrix, conditional on alph1_2, u1 and u2 */
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;
      mat c_sigX2 = out2.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub2(n,1);
      alph_sub2.col(0) = alph.col(0); // alph1_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f2 = alph_sub2.t();

      /* Centering the alpha */
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      /* Estimation of F1c1_1 */
      double F1c1_1 = pi1_1(i)*pn(alph_c2,sqrt(c_sig2[0]));

      /* Calculating the pdf */
      mat inner2 = alph_c2.t()*c_sig2.i()*alph_c2;
      double pdf1c1_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig2[0]))*exp(-0.5*inner2(0));

      /* Score contribution from F1c1_1 */
      double dF1c1_1_u1 = dpi1_1_u1(i)*pn(alph_c2,sqrt(c_sig2[0]))+pi1_1(i)*pdf1c1_1*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*alph(i,1)+c_sigX2(1)*u(i,0)+c_sigX2(2)*u(i,1)-alph(i,0))*c_sigX2(1);
      double dF1c1_1_u2 = dpi1_1_u2(i)*pn(alph_c2,sqrt(c_sig2[0]))+pi1_1(i)*pdf1c1_1*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*alph(i,1)+c_sigX2(1)*u(i,0)+c_sigX2(2)*u(i,1)-alph(i,0))*c_sigX2(2);

      /* Conditional F2c1_1 */
      /* Conditional mean and variance-covariance matrix, conditional on alph1_2, u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;
      mat c_sigX3 = out3.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub3(n,1);
      alph_sub3.col(0) = alph.col(2); // alph2_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f3 = alph_sub3.t();

      /* Centering the alpha */
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      /* Estimation of F2c1_1 */
      double F2c1_1 = pi2_1(i)*pn(alph_c3,sqrt(c_sig3[0]));

      /* Calculating the pdf */
      mat inner3 = alph_c3.t()*c_sig3.i()*alph_c3;
      double pdf2c1_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig3[0]))*exp(-0.5*inner3(0));

      /* Score contribution from F2c1_1 */
      double dF2c1_1_u1 = dpi2_1_u1(i)*pn(alph_c3,sqrt(c_sig3[0]))+pi2_1(i)*pdf2c1_1*(-1/(2*c_sig3[0]))*2*(c_sigX3(0)*alph(i,1)+c_sigX3(1)*u(i,0)+c_sigX3(2)*u(i,1)-alph(i,2))*c_sigX3(1);
      double dF2c1_1_u2 = dpi2_1_u2(i)*pn(alph_c3,sqrt(c_sig3[0]))+pi2_1(i)*pdf2c1_1*(-1/(2*c_sig3[0]))*2*(c_sigX3(0)*alph(i,1)+c_sigX3(1)*u(i,0)+c_sigX3(2)*u(i,1)-alph(i,2))*c_sigX3(2);

      /* Likelihood contribution */
      double dF01 = dF1_2*(1-F1c1_1-F2c1_1);

      /* Score contributions from dF01 wrt. u1 and u2 */
      double sc_u1 = (1/dF01)*(ddF1_2_u1*(1-F1c1_1-F2c1_1)+dF1_2*(-dF1c1_1_u1-dF2c1_1_u1));
      double sc_u2 = (1/dF01)*(ddF1_2_u2*(1-F1c1_1-F2c1_1)+dF1_2*(-dF1c1_1_u2-dF2c1_1_u2));

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;
    }
    /* Family member 1 experience event 1, family member 2 experience event 0, estimating score contribution from dF10 */
    else if((y(i,0) == 1) & (y(i,1) == 0)){

      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 1;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 0 ; rc4(1) = 4; rc4(2) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal dF1_1 */
      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc1, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat c_sigX = out.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub(n,1);
      alph_sub.col(0) = alph.col(0); // alph1_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f = alph_sub.t();

      /* Centering the alpha */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf1_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Estimation of dF1_1 */
      double dF1_1 = pi1_1(i)*dalph(i,0)*pdf1_1; // pi1_1, dalph1_1, pdf1_1

      /* Score contribution from dF1_1 */
      double ddF1_1_u1 = dalph(i,0)*(dpi1_1_u1(i)*pdf1_1+pi1_1(i)*pdf1_1*(-1/(2*c_sig[0]))*2*(c_sigX(0)*u(i,0)+c_sigX(1)*u(i,1)-alph(i,0))*c_sigX(0));

      double ddF1_1_u2 = dalph(i,0)*(dpi1_1_u2(i)*pdf1_1+pi1_1(i)*pdf1_1*(-1/(2*c_sig[0]))*2*(c_sigX(0)*u(i,0)+c_sigX(1)*u(i,1)-alph(i,0))*c_sigX(1));

      /* Conditional F1c1_2 */
      /* Conditional mean and variance-covariance matrix, conditional on alph1_1, u1 and u2 */
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;
      mat c_sigX2 = out2.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub2(n,1);
      alph_sub2.col(0) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f2 = alph_sub2.t();

      /* Centering the alpha */
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      /* Estimation of F1c1_2 */
      double F1c1_2 = pi1_2(i)*pn(alph_c2,sqrt(c_sig2[0]));

      /* Calculating the pdf */
      mat inner2 = alph_c2.t()*c_sig2.i()*alph_c2;
      double pdf1c1_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig2[0]))*exp(-0.5*inner2(0));

      /* Score contribution from F1c1_1 */
      double dF1c1_2_u1 = dpi1_2_u1(i)*pn(alph_c2,sqrt(c_sig2[0]))+pi1_2(i)*pdf1c1_2*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*alph(i,0)+c_sigX2(1)*u(i,0)+c_sigX2(2)*u(i,1)-alph(i,1))*c_sigX2(1);
      double dF1c1_2_u2 = dpi1_2_u2(i)*pn(alph_c2,sqrt(c_sig2[0]))+pi1_2(i)*pdf1c1_2*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*alph(i,0)+c_sigX2(1)*u(i,0)+c_sigX2(2)*u(i,1)-alph(i,1))*c_sigX2(2);

      /* Conditional F2c1_2 */
      /* Conditional mean and variance-covariance matrix, conditional on alph1_1, u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;
      mat c_sigX3 = out3.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub3(n,1);
      alph_sub3.col(0) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f3 = alph_sub3.t();

      /* Centering the alpha */
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      /* Estimation of F2c1_2 */
      double F2c1_2 = pi2_2(i)*pn(alph_c3,sqrt(c_sig3[0]));

      /* Calculating the pdf */
      mat inner3 = alph_c3.t()*c_sig3.i()*alph_c3;
      double pdf2c1_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig3[0]))*exp(-0.5*inner3(0));

      /* Score contribution from F2c1_1 */
      double dF2c1_2_u1 = dpi2_2_u1(i)*pn(alph_c3,sqrt(c_sig3[0]))+pi2_2(i)*pdf2c1_2*(-1/(2*c_sig3[0]))*2*(c_sigX3(0)*alph(i,0)+c_sigX3(1)*u(i,0)+c_sigX3(2)*u(i,1)-alph(i,3))*c_sigX3(1);
      double dF2c1_2_u2 = dpi2_2_u2(i)*pn(alph_c3,sqrt(c_sig3[0]))+pi2_2(i)*pdf2c1_2*(-1/(2*c_sig3[0]))*2*(c_sigX3(0)*alph(i,0)+c_sigX3(1)*u(i,0)+c_sigX3(2)*u(i,1)-alph(i,3))*c_sigX3(2);

      /* Likelihood contribution */
      double dF10 = dF1_1*(1-F1c1_2-F2c1_2);

      /* Score contributions from dF10 wrt. u1 and u2 */
      double sc_u1 = (1/dF10)*(ddF1_1_u1*(1-F1c1_2-F2c1_2)+dF1_1*(-dF1c1_2_u1-dF2c1_2_u1));
      double sc_u2 = (1/dF10)*(ddF1_1_u2*(1-F1c1_2-F2c1_2)+dF1_1*(-dF1c1_2_u2-dF2c1_2_u2));

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;
    }
    /* Family member 1 experience event 0, family member 2 experience event 2, estimating score contribution from dF02 */
    else if((y(i,0) == 0) & (y(i,1) == 2)){

      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 0;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 3 ; rc4(1) = 4; rc4(2) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal dF2_2 */
      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc2, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat c_sigX = out.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub(n,1);
      alph_sub.col(0) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f = alph_sub.t();

      /* Centering the alpha */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf2_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Estimation of dF2_2 */
      double dF2_2 = pi2_2(i)*dalph(i,3)*pdf2_2; // pi2_2, dalph2_2, pdf2_2

      /* Score contribution from dF2_2 */
      double ddF2_2_u1 = dalph(i,3)*(dpi2_2_u1(i)*pdf2_2+pi2_2(i)*pdf2_2*(-1/(2*c_sig[0]))*2*(c_sigX(0)*u(i,0)+c_sigX(1)*u(i,1)-alph(i,3))*c_sigX(0));
      double ddF2_2_u2 = dalph(i,3)*(dpi2_2_u2(i)*pdf2_2+pi2_2(i)*pdf2_2*(-1/(2*c_sig[0]))*2*(c_sigX(0)*u(i,0)+c_sigX(1)*u(i,1)-alph(i,3))*c_sigX(1));

      /* Conditional F1c2_1 */
      /* Conditional mean and variance-covariance matrix, conditional on alph2_2, u1 and u2 */
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;
      mat c_sigX2 = out2.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub2(n,1);
      alph_sub2.col(0) = alph.col(0); // alph1_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f2 = alph_sub2.t();

      /* Centering the alpha */
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      /* Estimation of F1c2_1 */
      double F1c2_1 = pi1_1(i)*pn(alph_c2,sqrt(c_sig2[0]));

      /* Calculating the pdf */
      mat inner2 = alph_c2.t()*c_sig2.i()*alph_c2;
      double pdf1c2_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig2[0]))*exp(-0.5*inner2(0));

      /* Score contribution from F1c2_1 */
      double dF1c2_1_u1 = dpi1_1_u1(i)*pn(alph_c2,sqrt(c_sig2[0]))+pi1_1(i)*pdf1c2_1*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*alph(i,3)+c_sigX2(1)*u(i,0)+c_sigX2(2)*u(i,1)-alph(i,0))*c_sigX2(1);
      double dF1c2_1_u2 = dpi1_1_u2(i)*pn(alph_c2,sqrt(c_sig2[0]))+pi1_1(i)*pdf1c2_1*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*alph(i,3)+c_sigX2(1)*u(i,0)+c_sigX2(2)*u(i,1)-alph(i,0))*c_sigX2(2);

      /* Conditional F2c2_1 */
      /* Conditional mean and variance-covariance matrix, conditional on alph2_2, u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;
      mat c_sigX3 = out3.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub3(n,1);
      alph_sub3.col(0) = alph.col(2); // alph2_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f3 = alph_sub3.t();

      /* Centering the alpha */
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      /* Estimation of F2c2_1 */
      double F2c2_1 = pi2_1(i)*pn(alph_c3,sqrt(c_sig3[0]));

      /* Calculating the pdf */
      mat inner3 = alph_c3.t()*c_sig3.i()*alph_c3;
      double pdf2c2_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig3[0]))*exp(-0.5*inner3(0));

      /* Score contribution from F2c2_1 */
      double dF2c2_1_u1 = dpi2_1_u1(i)*pn(alph_c3,sqrt(c_sig3[0]))+pi2_1(i)*pdf2c2_1*(-1/(2*c_sig3[0]))*2*(c_sigX3(0)*alph(i,3)+c_sigX3(1)*u(i,0)+c_sigX3(2)*u(i,1)-alph(i,2))*c_sigX3(1);
      double dF2c2_1_u2 = dpi2_1_u2(i)*pn(alph_c3,sqrt(c_sig3[0]))+pi2_1(i)*pdf2c2_1*(-1/(2*c_sig3[0]))*2*(c_sigX3(0)*alph(i,3)+c_sigX3(1)*u(i,0)+c_sigX3(2)*u(i,1)-alph(i,2))*c_sigX3(2);

      /* Likelihood contribution */
      double dF02 = dF2_2*(1-F1c2_1-F2c2_1);

      /* Score contributions from dF02 wrt. u1 and u2 */
      double sc_u1 = (1/dF02)*(ddF2_2_u1*(1-F1c2_1-F2c2_1)+dF2_2*(-dF1c2_1_u1-dF2c2_1_u1));
      double sc_u2 = (1/dF02)*(ddF2_2_u2*(1-F1c2_1-F2c2_1)+dF2_2*(-dF1c2_1_u2-dF2c2_1_u2));

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;
    }
    /* Family member 1 experience event 2, family member 2 experience event 0, estimating score contribution from dF20 */
    else if((y(i,0) == 2) & (y(i,1) == 0)){

      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 1;
      uvec rc2(1); rc2(0) = 3;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 2 ; rc4(1) = 4; rc4(2) = 5;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal dF2_1 */
      /* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
      vecmat out = conMuSig(sigma, mu, rc2, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat c_sigX = out.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub(n,1);
      alph_sub.col(0) = alph.col(2); // alph2_1

      /* Transposing matrix to be of the form 1xn */
      mat alph_f = alph_sub.t();

      /* Centering the alpha */
      mat alph_c = alph_f.col(i)-c_mu;

      /* Calculating the pdf */
      mat inner = alph_c.t()*c_sig.i()*alph_c;
      double pdf2_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      /* Estimation of dF2_1 */
      double dF2_1 = pi2_1(i)*dalph(i,2)*pdf2_1; // pi2_1, dalph2_1, pdf2_1

      /* Score contribution from dF2_1 */
      double ddF2_1_u1 = dalph(i,2)*(dpi2_1_u1(i)*pdf2_1+pi2_1(i)*pdf2_1*(-1/(2*c_sig[0]))*2*(c_sigX(0)*u(i,0)+c_sigX(1)*u(i,1)-alph(i,2))*c_sigX(0));
      double ddF2_1_u2 = dalph(i,2)*(dpi2_1_u2(i)*pdf2_1+pi2_1(i)*pdf2_1*(-1/(2*c_sig[0]))*2*(c_sigX(0)*u(i,0)+c_sigX(1)*u(i,1)-alph(i,2))*c_sigX(1));

      /* Conditional F1c2_2 */
      /* Conditional mean and variance-covariance matrix, conditional on alph2_1, u1 and u2 */
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;
      mat c_sigX2 = out2.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub2(n,1);
      alph_sub2.col(0) = alph.col(1); // alph1_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f2 = alph_sub2.t();

      /* Centering the alpha */
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      /* Estimation of F1c2_2 */
      double F1c2_2 = pi1_2(i)*pn(alph_c2,sqrt(c_sig2[0]));

      /* Calculating the pdf */
      mat inner2 = alph_c2.t()*c_sig2.i()*alph_c2;
      double pdf1c2_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig2[0]))*exp(-0.5*inner2(0));

      /* Score contribution from F1c2_2 */
      double dF1c2_2_u1 = dpi1_2_u1(i)*pn(alph_c2,sqrt(c_sig2[0]))+pi1_2(i)*pdf1c2_2*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*alph(i,2)+c_sigX2(1)*u(i,0)+c_sigX2(2)*u(i,1)-alph(i,1))*c_sigX2(1);
      double dF1c2_2_u2 = dpi1_2_u2(i)*pn(alph_c2,sqrt(c_sig2[0]))+pi1_2(i)*pdf1c2_2*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*alph(i,2)+c_sigX2(1)*u(i,0)+c_sigX2(2)*u(i,1)-alph(i,1))*c_sigX2(2);

      /* Conditional F2c2_2 */
      /* Conditional mean and variance-covariance matrix, conditional on alph2_1, u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;
      mat c_sigX3 = out3.M2;

      /* Pulling out the appropriate alpha from alph */
      mat alph_sub3(n,1);
      alph_sub3.col(0) = alph.col(3); // alph2_2

      /* Transposing matrix to be of the form 1xn */
      mat alph_f3 = alph_sub3.t();

      /* Centering the alpha */
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      /* Estimation of F2c2_2 */
      double F2c2_2 = pi2_2(i)*pn(alph_c3,sqrt(c_sig3[0]));

      /* Calculating the pdf */
      mat inner3 = alph_c3.t()*c_sig3.i()*alph_c3;
      double pdf2c2_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig3[0]))*exp(-0.5*inner3(0));

      /* Score contribution from F2c2_2 */
      double dF2c2_2_u1 = dpi2_2_u1(i)*pn(alph_c3,sqrt(c_sig3[0]))+pi2_2(i)*pdf2c2_2*(-1/(2*c_sig3[0]))*2*(c_sigX3(0)*alph(i,2)+c_sigX3(1)*u(i,0)+c_sigX3(2)*u(i,1)-alph(i,3))*c_sigX3(1);
      double dF2c2_2_u2 = dpi2_2_u2(i)*pn(alph_c3,sqrt(c_sig3[0]))+pi2_2(i)*pdf2c2_2*(-1/(2*c_sig3[0]))*2*(c_sigX3(0)*alph(i,2)+c_sigX3(1)*u(i,0)+c_sigX3(2)*u(i,1)-alph(i,3))*c_sigX3(2);

      /* Likelihood contribution */
      double dF20 = dF2_1*(1-F1c2_2-F2c2_2);

      /* Score contributions from dF20 wrt. u1 and u2 */
      double sc_u1 = (1/dF20)*(ddF2_1_u1*(1-F1c2_2-F2c2_2)+dF2_1*(-dF1c2_2_u1-dF2c2_2_u1));
      double sc_u2 = (1/dF20)*(ddF2_1_u2*(1-F1c2_2-F2c2_2)+dF2_1*(-dF1c2_2_u2-dF2c2_2_u2));

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;
    }
    /* Both family members experience event 0, estimating score contribution from F00 */
    else{
      /* Specifying which parts of sigma apply */
      uvec rc1(1); rc1(0) = 0;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;

      uvec rc4(2); rc4(0) = 0 ; rc4(1) = 1;
      uvec rc5(2); rc5(0) = 0 ; rc5(1) = 3;
      uvec rc6(2); rc6(0) = 2 ; rc6(1) = 1;
      uvec rc7(2); rc7(0) = 2 ; rc7(1) = 3;

      /* Mean vector (used for estimation of conditional mean) */
      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      /* Marginal F1_1, F1_2, F2_1 and F2_2 */
      /* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
      vecmat out1 = conMuSig(sigma, mu, rc1, rc3);
      vec c_mu1 = out1.V;
      mat c_sig1 = out1.M1;
      mat c_sigX1 = out1.M2;

      vecmat out2 = conMuSig(sigma, mu, rc2, rc3);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;
      mat c_sigX2 = out2.M2;

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub1_1(n,1);
      mat alph_sub1_2(n,1);
      mat alph_sub2_1(n,1);
      mat alph_sub2_2(n,1);
      alph_sub1_1.col(0) = alph.col(0); // alph1_1
      alph_sub1_2.col(0) = alph.col(1); // alph1_2
      alph_sub2_1.col(0) = alph.col(2); // alph2_1
      alph_sub2_2.col(0) = alph.col(3); // alph2_2

      /* Transposing matrices to be of the form 1xn */
      mat alph_f1_1 = alph_sub1_1.t();
      mat alph_f1_2 = alph_sub1_2.t();
      mat alph_f2_1 = alph_sub2_1.t();
      mat alph_f2_2 = alph_sub2_2.t();

      /* Centering the alphas */
      mat alph_c1_1 = alph_f1_1.col(i)-c_mu1;
      mat alph_c1_2 = alph_f1_2.col(i)-c_mu1;
      mat alph_c2_1 = alph_f2_1.col(i)-c_mu2;
      mat alph_c2_2 = alph_f2_2.col(i)-c_mu2;

      /* Estimation of F1_1, F1_2, F2_1 and F2_2 */
      double F1_1 = pi1_1(i)*pn(alph_c1_1,sqrt(c_sig1[0]));
      double F1_2 = pi1_2(i)*pn(alph_c1_2,sqrt(c_sig1[0]));
      double F2_1 = pi2_1(i)*pn(alph_c2_1,sqrt(c_sig2[0]));
      double F2_2 = pi2_2(i)*pn(alph_c2_2,sqrt(c_sig2[0]));

      /* Calculating the pdf of F1_1, F1_2, F2_1 and F2_2 */
      mat inner1 = alph_c1_1.t()*c_sig1.i()*alph_c1_1;
      mat inner2 = alph_c1_2.t()*c_sig1.i()*alph_c1_2;
      mat inner3 = alph_c2_1.t()*c_sig2.i()*alph_c2_1;
      mat inner4 = alph_c2_2.t()*c_sig2.i()*alph_c2_2;

      double pdf1_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig1))*exp(-0.5*inner1(0));
      double pdf1_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig1))*exp(-0.5*inner2(0));
      double pdf2_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig2))*exp(-0.5*inner3(0));
      double pdf2_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig2))*exp(-0.5*inner4(0));

      /* Score contributions from F1_1, F1_2, F2_1 and F2_2 */
      double dF1_1_u1 = dpi1_1_u1(i)*pn(alph_c1_1,sqrt(c_sig1[0]))+pi1_1(i)*pdf1_1*(-1/(2*c_sig1[0]))*2*(c_sigX1(0)*u(i,0)+c_sigX1(1)*u(i,1)-alph(i,0))*c_sigX1(0);
      double dF1_2_u1 = dpi1_2_u1(i)*pn(alph_c1_2,sqrt(c_sig1[0]))+pi1_2(i)*pdf1_2*(-1/(2*c_sig1[0]))*2*(c_sigX1(0)*u(i,0)+c_sigX1(1)*u(i,1)-alph(i,1))*c_sigX1(0);
      double dF2_1_u1 = dpi2_1_u1(i)*pn(alph_c2_1,sqrt(c_sig2[0]))+pi2_1(i)*pdf2_1*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*u(i,0)+c_sigX2(1)*u(i,1)-alph(i,2))*c_sigX2(0);
      double dF2_2_u1 = dpi2_2_u1(i)*pn(alph_c2_2,sqrt(c_sig2[0]))+pi2_2(i)*pdf2_2*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*u(i,0)+c_sigX2(1)*u(i,1)-alph(i,3))*c_sigX2(0);

      double dF1_1_u2 = dpi1_1_u2(i)*pn(alph_c1_1,sqrt(c_sig1[0]))+pi1_1(i)*pdf1_1*(-1/(2*c_sig1[0]))*2*(c_sigX1(0)*u(i,0)+c_sigX1(1)*u(i,1)-alph(i,0))*c_sigX1(1);
      double dF1_2_u2 = dpi1_2_u2(i)*pn(alph_c1_2,sqrt(c_sig1[0]))+pi1_2(i)*pdf1_2*(-1/(2*c_sig1[0]))*2*(c_sigX1(0)*u(i,0)+c_sigX1(1)*u(i,1)-alph(i,1))*c_sigX1(1);
      double dF2_1_u2 = dpi2_1_u2(i)*pn(alph_c2_1,sqrt(c_sig2[0]))+pi2_1(i)*pdf2_1*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*u(i,0)+c_sigX2(1)*u(i,1)-alph(i,2))*c_sigX2(1);
      double dF2_2_u2 = dpi2_2_u2(i)*pn(alph_c2_2,sqrt(c_sig2[0]))+pi2_2(i)*pdf2_2*(-1/(2*c_sig2[0]))*2*(c_sigX2(0)*u(i,0)+c_sigX2(1)*u(i,1)-alph(i,3))*c_sigX2(1);

      /* Joint probabilities F11, F12, F21 and F22 */
      /* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
      vecmat out3 = conMuSig(sigma, mu, rc4, rc3);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;
      mat c_sigX3 = out3.M2;

      vecmat out4 = conMuSig(sigma, mu, rc5, rc3);
      vec c_mu4 = out4.V;
      mat c_sig4 = out4.M1;
      mat c_sigX4 = out4.M2;

      vecmat out5 = conMuSig(sigma, mu, rc6, rc3);
      vec c_mu5 = out5.V;
      mat c_sig5 = out5.M1;
      mat c_sigX5 = out5.M2;

      vecmat out6 = conMuSig(sigma, mu, rc7, rc3);
      vec c_mu6 = out6.V;
      mat c_sig6 = out6.M1;
      mat c_sigX6 = out6.M2;

      /* Correlation coefficient */
      double r1 = c_sig3(0,1)/(sqrt(c_sig3(0,0))*sqrt(c_sig3(1,1)));
      double r2 = c_sig4(0,1)/(sqrt(c_sig4(0,0))*sqrt(c_sig4(1,1)));
      double r3 = c_sig5(0,1)/(sqrt(c_sig5(0,0))*sqrt(c_sig5(1,1)));
      double r4 = c_sig6(0,1)/(sqrt(c_sig6(0,0))*sqrt(c_sig6(1,1)));

      /* Pulling out the appropriate alphas from alph */
      mat alph_sub11(n,2);
      alph_sub11.col(0) = alph.col(0); // alph1_1
      alph_sub11.col(1) = alph.col(1); // alph1_2

      mat alph_sub12(n,2);
      alph_sub12.col(0) = alph.col(0); // alph1_1
      alph_sub12.col(1) = alph.col(3); // alph2_2

      mat alph_sub21(n,2);
      alph_sub21.col(0) = alph.col(2); // alph2_1
      alph_sub21.col(1) = alph.col(1); // alph1_2

      mat alph_sub22(n,2);
      alph_sub22.col(0) = alph.col(2); // alph2_1
      alph_sub22.col(1) = alph.col(3); // alph2_2

      /* Transposing matrices to be of the form 2xn */
      mat alph_f11 = alph_sub11.t();
      mat alph_f12 = alph_sub12.t();
      mat alph_f21 = alph_sub21.t();
      mat alph_f22 = alph_sub22.t();

      /* Centering the alphas */
      mat alph_c11 = alph_f11.col(i)-c_mu3;
      mat alph_c12 = alph_f12.col(i)-c_mu4;
      mat alph_c21 = alph_f21.col(i)-c_mu5;
      mat alph_c22 = alph_f22.col(i)-c_mu6;

      /* Estimating F11, F12, F21 and F22 */
      double F11 = pi1_1(i)*pi1_2(i)*pn(alph_c11,r1);
      double F12 = pi1_1(i)*pi2_2(i)*pn(alph_c12,r2);
      double F21 = pi2_1(i)*pi1_2(i)*pn(alph_c21,r3);
      double F22 = pi2_1(i)*pi2_2(i)*pn(alph_c22,r4);

      /* Calculating the pdf of F1_1, F1_2, F2_1 and F2_2 */
      mat inner11 = alph_c11.t()*c_sig3.i()*alph_c11;
      mat inner12 = alph_c12.t()*c_sig4.i()*alph_c12;
      mat inner21 = alph_c21.t()*c_sig5.i()*alph_c21;
      mat inner22 = alph_c22.t()*c_sig6.i()*alph_c22;

      double pdf11 = 1/(twopi)*1/sqrt(det(c_sig3))*exp(-0.5*inner11(0));
      double pdf12 = 1/(twopi)*1/sqrt(det(c_sig4))*exp(-0.5*inner12(0));
      double pdf21 = 1/(twopi)*1/sqrt(det(c_sig5))*exp(-0.5*inner21(0));
      double pdf22 = 1/(twopi)*1/sqrt(det(c_sig6))*exp(-0.5*inner22(0));

      /* Score contributions from pdf11, pdf12, pdf21 and pdf22 */
      double dpdf11_u1 = pdf11*(-1/(2*pow(1-r1,2)))*((2*(pow(c_sigX3(0,0),2)*sqrt(c_sig3(1,1))*u(i,0)+c_sigX3(0,0)*c_sigX3(0,1)*sqrt(c_sig3(1,1))*u(i,1)-2*c_sigX3(0,0)*c_sigX3(1,0)*r1*u(i,0)-c_sigX3(0,0)*c_sigX3(1,1)*r1*u(i,1)-c_sigX3(0,1)*c_sigX3(1,0)*r1*u(i,1)+pow(c_sigX3(1,0),2)*sqrt(c_sig3(0,0))*u(i,0)+c_sigX3(1,0)*c_sigX3(1,1)*sqrt(c_sig3(0,0))*u(i,1)+c_sigX3(0,0)*r1*alph(i,1)-c_sigX3(0,0)*sqrt(c_sig3(1,1))*alph(i,0)+c_sigX3(1,0)*r1*alph(i,0)-c_sigX3(1,0)*sqrt(c_sig3(0,0))*alph(i,1)))/(sqrt(c_sig3(0,0))*sqrt(c_sig3(1,1))));
      double dpdf12_u1 = pdf12*(-1/(2*pow(1-r2,2)))*((2*(pow(c_sigX4(0,0),2)*sqrt(c_sig4(1,1))*u(i,0)+c_sigX4(0,0)*c_sigX4(0,1)*sqrt(c_sig4(1,1))*u(i,1)-2*c_sigX4(0,0)*c_sigX4(1,0)*r2*u(i,0)-c_sigX4(0,0)*c_sigX4(1,1)*r2*u(i,1)-c_sigX4(0,1)*c_sigX4(1,0)*r2*u(i,1)+pow(c_sigX4(1,0),2)*sqrt(c_sig4(0,0))*u(i,0)+c_sigX4(1,0)*c_sigX4(1,1)*sqrt(c_sig4(0,0))*u(i,1)+c_sigX4(0,0)*r2*alph(i,3)-c_sigX4(0,0)*sqrt(c_sig4(1,1))*alph(i,0)+c_sigX4(1,0)*r2*alph(i,0)-c_sigX4(1,0)*sqrt(c_sig4(0,0))*alph(i,3)))/(sqrt(c_sig4(0,0))*sqrt(c_sig4(1,1))));
      double dpdf21_u1 = pdf21*(-1/(2*pow(1-r3,2)))*((2*(pow(c_sigX5(0,0),2)*sqrt(c_sig5(1,1))*u(i,0)+c_sigX5(0,0)*c_sigX5(0,1)*sqrt(c_sig5(1,1))*u(i,1)-2*c_sigX5(0,0)*c_sigX5(1,0)*r3*u(i,0)-c_sigX5(0,0)*c_sigX5(1,1)*r3*u(i,1)-c_sigX5(0,1)*c_sigX5(1,0)*r4*u(i,1)+pow(c_sigX5(1,0),2)*sqrt(c_sig5(0,0))*u(i,0)+c_sigX5(1,0)*c_sigX5(1,1)*sqrt(c_sig5(0,0))*u(i,1)+c_sigX5(0,0)*r3*alph(i,1)-c_sigX5(0,0)*sqrt(c_sig5(1,1))*alph(i,2)+c_sigX5(1,0)*r3*alph(i,2)-c_sigX5(1,0)*sqrt(c_sig5(0,0))*alph(i,1)))/(sqrt(c_sig5(0,0))*sqrt(c_sig5(1,1))));
      double dpdf22_u1 = pdf22*(-1/(2*pow(1-r4,2)))*((2*(pow(c_sigX6(0,0),2)*sqrt(c_sig6(1,1))*u(i,0)+c_sigX6(0,0)*c_sigX6(0,1)*sqrt(c_sig6(1,1))*u(i,1)-2*c_sigX6(0,0)*c_sigX6(1,0)*r4*u(i,0)-c_sigX6(0,0)*c_sigX6(1,1)*r4*u(i,1)-c_sigX6(0,1)*c_sigX6(1,0)*r4*u(i,1)+pow(c_sigX6(1,0),2)*sqrt(c_sig6(0,0))*u(i,0)+c_sigX6(1,0)*c_sigX6(1,1)*sqrt(c_sig6(0,0))*u(i,1)+c_sigX6(0,0)*r4*alph(i,3)-c_sigX6(0,0)*sqrt(c_sig6(1,1))*alph(i,2)+c_sigX6(1,0)*r4*alph(i,2)-c_sigX6(1,0)*sqrt(c_sig6(0,0))*alph(i,3)))/(sqrt(c_sig6(0,0))*sqrt(c_sig6(1,1))));

      double dpdf11_u2 = pdf11*(-1/2*pow(1-r1,2))*((2*(c_sigX3(0,0)*c_sigX3(0,1)*sqrt(c_sig3(1,1))*u(i,0)-c_sigX3(0,0)*c_sigX3(1,1)*r1*u(i,0)+pow(c_sigX3(0,1),2)*sqrt(c_sig3(1,1))*u(i,1)-c_sigX3(0,1)*c_sigX3(1,0)*r1*u(i,0)-2*c_sigX3(0,1)*c_sigX3(1,1)*r1*u(i,1)+c_sigX3(1,0)*c_sigX3(1,1)*sqrt(c_sig3(0,0))*u(i,0)+pow(c_sigX3(1,1),2)*sqrt(c_sig3(0,0))*u(i,1)+c_sigX3(0,1)*r1*alph(i,1)-c_sigX3(0,1)*sqrt(c_sig3(1,1))*alph(i,0)+c_sigX3(1,1)*r1*alph(i,0)-c_sigX3(1,1)*sqrt(c_sig3(0,0))*alph(i,1)))/(sqrt(c_sig3(0,0))*sqrt(c_sig3(1,1))));
      double dpdf12_u2 = pdf12*(-1/2*pow(1-r2,2))*((2*(c_sigX4(0,0)*c_sigX4(0,1)*sqrt(c_sig4(1,1))*u(i,0)-c_sigX4(0,0)*c_sigX4(1,1)*r2*u(i,0)+pow(c_sigX4(0,1),2)*sqrt(c_sig4(1,1))*u(i,1)-c_sigX4(0,1)*c_sigX4(1,0)*r2*u(i,0)-2*c_sigX4(0,1)*c_sigX4(1,1)*r2*u(i,1)+c_sigX4(1,0)*c_sigX4(1,1)*sqrt(c_sig4(0,0))*u(i,0)+pow(c_sigX4(1,1),2)*sqrt(c_sig4(0,0))*u(i,1)+c_sigX4(0,1)*r2*alph(i,3)-c_sigX4(0,1)*sqrt(c_sig4(1,1))*alph(i,0)+c_sigX4(1,1)*r2*alph(i,0)-c_sigX4(1,1)*sqrt(c_sig4(0,0))*alph(i,3)))/(sqrt(c_sig4(0,0))*sqrt(c_sig4(1,1))));
      double dpdf21_u2 = pdf21*(-1/2*pow(1-r3,2))*((2*(c_sigX5(0,0)*c_sigX5(0,1)*sqrt(c_sig5(1,1))*u(i,0)-c_sigX5(0,0)*c_sigX5(1,1)*r3*u(i,0)+pow(c_sigX5(0,1),2)*sqrt(c_sig5(1,1))*u(i,1)-c_sigX5(0,1)*c_sigX5(1,0)*r3*u(i,0)-2*c_sigX5(0,1)*c_sigX5(1,1)*r3*u(i,1)+c_sigX5(1,0)*c_sigX5(1,1)*sqrt(c_sig5(0,0))*u(i,0)+pow(c_sigX5(1,1),2)*sqrt(c_sig5(0,0))*u(i,1)+c_sigX5(0,1)*r3*alph(i,1)-c_sigX5(0,1)*sqrt(c_sig5(1,1))*alph(i,2)+c_sigX5(1,1)*r3*alph(i,2)-c_sigX5(1,1)*sqrt(c_sig5(0,0))*alph(i,1)))/(sqrt(c_sig5(0,0))*sqrt(c_sig5(1,1))));
      double dpdf22_u2 = pdf22*(-1/2*pow(1-r4,2))*((2*(c_sigX6(0,0)*c_sigX6(0,1)*sqrt(c_sig6(1,1))*u(i,0)-c_sigX6(0,0)*c_sigX6(1,1)*r4*u(i,0)+pow(c_sigX6(0,1),2)*sqrt(c_sig6(1,1))*u(i,1)-c_sigX6(0,1)*c_sigX6(1,0)*r4*u(i,0)-2*c_sigX6(0,1)*c_sigX6(1,1)*r4*u(i,1)+c_sigX6(1,0)*c_sigX6(1,1)*sqrt(c_sig6(0,0))*u(i,0)+pow(c_sigX6(1,1),2)*sqrt(c_sig6(0,0))*u(i,1)+c_sigX6(0,1)*r4*alph(i,3)-c_sigX6(0,1)*sqrt(c_sig6(1,1))*alph(i,2)+c_sigX6(1,1)*r4*alph(i,2)-c_sigX6(1,1)*sqrt(c_sig6(0,0))*alph(i,3)))/(sqrt(c_sig6(0,0))*sqrt(c_sig6(1,1))));

      /* Score contributions from F11, F12, F21 and F22 */
      double dF11_u1 = dpi11_u1(i)*pn(alph_c11,r1)+pi1_1(i)*pi1_2(i)*dpdf11_u1;
      double dF12_u1 = dpi12_u1(i)*pn(alph_c12,r2)+pi1_1(i)*pi2_2(i)*dpdf12_u1;
      double dF21_u1 = dpi21_u1(i)*pn(alph_c21,r3)+pi2_1(i)*pi1_2(i)*dpdf21_u1;
      double dF22_u1 = dpi22_u1(i)*pn(alph_c22,r4)+pi2_1(i)*pi2_2(i)*dpdf22_u1;

      double dF11_u2 = dpi11_u2(i)*pn(alph_c11,r1)+pi1_1(i)*pi1_2(i)*dpdf11_u2;
      double dF12_u2 = dpi12_u2(i)*pn(alph_c12,r2)+pi1_1(i)*pi2_2(i)*dpdf12_u2;
      double dF21_u2 = dpi21_u2(i)*pn(alph_c21,r3)+pi2_1(i)*pi1_2(i)*dpdf21_u2;
      double dF22_u2 = dpi22_u2(i)*pn(alph_c22,r4)+pi2_1(i)*pi2_2(i)*dpdf22_u2;

      /* Calculating F01, F10, F02 and F20 */
      double F01 = F1_2-(F11+F21);
      double F10 = F1_1-(F11+F12);
      double F02 = F2_2-(F12+F22);
      double F20 = F2_1-(F21+F22);

      /* Score contributions from F01, F10, F02 and F20 */
      double dF01_u1 = dF1_2_u1-(dF11_u1+dF21_u1);
      double dF10_u1 = dF1_1_u1-(dF11_u1+dF12_u1);
      double dF02_u1 = dF2_2_u1-(dF22_u1+dF12_u1);
      double dF20_u1 = dF2_1_u1-(dF22_u1+dF21_u1);

      double dF01_u2 = dF1_2_u2-(dF11_u2+dF21_u2);
      double dF10_u2 = dF1_1_u2-(dF11_u2+dF12_u2);
      double dF02_u2 = dF2_2_u2-(dF22_u2+dF12_u2);
      double dF20_u2 = dF2_1_u2-(dF22_u2+dF21_u2);

      /* Likelihood contribution */
      double F00 = 1-(F11+F12+F21+F22+F01+F10+F02+F20);

      /* Score contribution from F00 */
      double sc_u1 = (1/F00)*-(dF11_u1+dF12_u1+dF21_u1+dF22_u1+dF01_u1+dF10_u1+dF02_u1+dF20_u1);
      double sc_u2 = (1/F00)*-(dF11_u2+dF12_u2+dF21_u2+dF22_u2+dF01_u2+dF10_u2+dF02_u2+dF20_u2);

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;
    }
  }
  return(res);
}

/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------------*/
/* Hes */
// [[Rcpp::export]]
mat Hes(mat y, mat b, mat u, mat sigma, mat alph, mat dalph){
  /* y: 1x2 matrix with event type (0, 1 or 2) of family member 1 and 2
     b: 1x4 matrix with XB for event type 1 and 2 (b1 and b2) for family member 1 and 2
        the order is b1_1, b1_2, b2_1 and b2_2
     u: 1x2 matrix with the random effects u1 and u2 affecting pi1 and pi2 (cluster-specific risk levels)
        the random effects are shared by family member 1 and 2
     sigma: 6x6 matrix. variance-covariance matrix, order: event1_1, event1_2, event2_1, event2_2, u1, u2
     alph: 1x4 matrix, inside the Probit link, order a1_1, a1_2, a2_1, a2_2
     dalph: 1x4 matrix, derivative of alph wrt. t, order da1_1, da1_2, da2_1, da2_2
  */

  /* Initialising Hessian matrix */
  mat res(2,2);

  mat u1_p(1,2); u1_p(0,0) = u(0,0)+h; u1_p(0,1) = u(0,1);
  mat u1_m(1,2); u1_m(0,0) = u(0,0)-h; u1_m(0,1) = u(0,1);
  mat u2_p(1,2); u2_p(0,0) = u(0,0); u2_p(0,1) = u(0,1)+h;
  mat u2_m(1,2); u2_m(0,0) = u(0,0); u2_m(0,1) = u(0,1)-h;

  res.row(0) = (Dloglik(y, b, u1_p, sigma, alph, dalph)-Dloglik(y, b, u1_m, sigma, alph, dalph))/(2*h);
  res.row(1) = (Dloglik(y, b, u2_p, sigma, alph, dalph)-Dloglik(y, b, u2_m, sigma, alph, dalph))/(2*h);

  return(res);
}

//Rcpp::Rcout << "a12" << std::endl << alph_c12 ;
//Rcpp::Rcout << "everywhere";
