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


//   Rcpp::Rcout << "everywhere";
