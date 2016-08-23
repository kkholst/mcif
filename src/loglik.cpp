// [[Rcpp::depends(RcppArmadillo,mvtnorm)]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <mvtnormAPI.h>
#include "quadrule.h"

using namespace Rcpp;
using namespace arma;

static int _mvt_df = 0;
static double _mvt_abseps=0.00001;
static double _mvt_releps=0;
static int _mvt_maxpts=20000;
static int _mvt_inform;
static double _mvt_error[3];

const double twopi = 2*datum::pi;
const double sq_twopi = sqrt(twopi);
const double h = 1e-8;

const double _inner_NR_abseps=0.00001;

class matfour {
public:
  mat M1;
  mat M2;
  double M3;
  mat M4;
};

class ss {
public:
  mat e1_1s; /* Individual 1 event 1 conditional on u*/
  mat e2_1s; /* Individual 1 event 2 conditional on u*/
  mat e1_2s; /* Individual 2 event 1 conditional on u*/
  mat e2_2s; /* Individual 2 event 2 conditional on u*/

  mat e1_1i; /* Individual 1 event 1 conditional on u, the inverse*/
  mat e2_1i; /* Individual 1 event 2 conditional on u, the inverse*/
  mat e1_2i; /* Individual 2 event 1 conditional on u, the inverse*/
  mat e2_2i; /* Individual 2 event 2 conditional on u, the inverse*/

  double e1_1d; /* Individual 1 event 1 conditional on u, the determinant*/
  double e2_1d; /* Individual 1 event 2 conditional on u, the determinant*/
  double e1_2d; /* Individual 2 event 1 conditional on u, the determinant*/
  double e2_2d; /* Individual 2 event 2 conditional on u, the determinant*/

  mat e1_1x; /* Individual 1 event 1 conditional on u, the sigx*/
  mat e2_1x; /* Individual 1 event 2 conditional on u, the sigx*/
  mat e1_2x; /* Individual 2 event 1 conditional on u, the sigx*/
  mat e2_2x; /* Individual 2 event 2 conditional on u, the sigx*/

  mat e1c1_1s; /* Individual 1 event 1 conditional on u and individual 2 event 1*/
  mat e2c1_1s; /* Individual 1 event 2 conditional on u and individual 2 event 1*/
  mat e1c1_2s; /* Individual 2 event 1 conditional on u and individual 1 event 1*/
  mat e2c1_2s; /* Individual 2 event 2 conditional on u and individual 1 event 1*/

  mat e1c1_1i; /* Individual 1 event 1 conditional on u and individual 2 event 1, the inverse*/
  mat e2c1_1i; /* Individual 1 event 2 conditional on u and individual 2 event 1, the inverse*/
  mat e1c1_2i; /* Individual 2 event 1 conditional on u and individual 1 event 1, the inverse*/
  mat e2c1_2i; /* Individual 2 event 2 conditional on u and individual 1 event 1, the inverse*/

  double e1c1_1d; /* Individual 1 event 1 conditional on u and individual 2 event 1, the determinant*/
  double e2c1_1d; /* Individual 1 event 2 conditional on u and individual 2 event 1, the determinant*/
  double e1c1_2d; /* Individual 2 event 1 conditional on u and individual 1 event 1, the determinant*/
  double e2c1_2d; /* Individual 2 event 2 conditional on u and individual 1 event 1, the determinant*/

  mat e1c1_1x; /* Individual 1 event 1 conditional on u and individual 2 event 1, the sigx*/
  mat e2c1_1x; /* Individual 1 event 2 conditional on u and individual 2 event 1, the sigx*/
  mat e1c1_2x; /* Individual 2 event 1 conditional on u and individual 1 event 1, the sigx*/
  mat e2c1_2x; /* Individual 2 event 2 conditional on u and individual 1 event 1, the sigx*/

  mat e1c2_1s; /* Individual 1 event 1 conditional on u and individual 2 event 2*/
  mat e2c2_1s; /* Individual 1 event 2 conditional on u and individual 2 event 2*/
  mat e1c2_2s; /* Individual 2 event 1 conditional on u and individual 1 event 2*/
  mat e2c2_2s; /* Individual 2 event 2 conditional on u and individual 1 event 2*/

  mat e1c2_1i; /* Individual 1 event 1 conditional on u and individual 2 event 2, the inverse*/
  mat e2c2_1i; /* Individual 1 event 2 conditional on u and individual 2 event 2, the inverse*/
  mat e1c2_2i; /* Individual 2 event 1 conditional on u and individual 1 event 2, the inverse*/
  mat e2c2_2i; /* Individual 2 event 2 conditional on u and individual 1 event 2, the inverse*/

  double e1c2_1d; /* Individual 1 event 1 conditional on u and individual 2 event 2, the determinant*/
  double e2c2_1d; /* Individual 1 event 2 conditional on u and individual 2 event 2, the determinant*/
  double e1c2_2d; /* Individual 2 event 1 conditional on u and individual 1 event 2, the determinant*/
  double e2c2_2d; /* Individual 2 event 2 conditional on u and individual 1 event 2, the determinant*/

  mat e1c2_1x; /* Individual 1 event 1 conditional on u and individual 2 event 2, the sigx*/
  mat e2c2_1x; /* Individual 1 event 2 conditional on u and individual 2 event 2, the sigx*/
  mat e1c2_2x; /* Individual 2 event 1 conditional on u and individual 1 event 2, the sigx*/
  mat e2c2_2x; /* Individual 2 event 2 conditional on u and individual 1 event 2, the sigx*/

  mat e11s; /* Individual 1 event 1 and individual 2 event 1 conditional on u*/
  mat e12s; /* Individual 1 event 1 and individual 2 event 2 conditional on u*/
  mat e21s; /* Individual 1 event 2 and individual 2 event 1 conditional on u*/
  mat e22s; /* Individual 1 event 2 and individual 2 event 2 conditional on u*/

  mat e11i; /* Individual 1 event 1 and individual 2 event 1 conditional on u, the inverse*/
  mat e12i; /* Individual 1 event 1 and individual 2 event 2 conditional on u, the inverse*/
  mat e21i; /* Individual 1 event 2 and individual 2 event 1 conditional on u, the inverse*/
  mat e22i; /* Individual 1 event 2 and individual 2 event 2 conditional on u, the inverse*/

  double e11d; /* Individual 1 event 1 and individual 2 event 1 conditional on u, the determinant*/
  double e12d; /* Individual 1 event 1 and individual 2 event 2 conditional on u, the determinant*/
  double e21d; /* Individual 1 event 2 and individual 2 event 1 conditional on u, the determinant*/
  double e22d; /* Individual 1 event 2 and individual 2 event 2 conditional on u, the determinant*/

  mat e11x; /* Individual 1 event 1 and individual 2 event 1 conditional on u, the sigx*/
  mat e12x; /* Individual 1 event 1 and individual 2 event 2 conditional on u, the sigx*/
  mat e21x; /* Individual 1 event 2 and individual 2 event 1 conditional on u, the sigx*/
  mat e22x; /* Individual 1 event 2 and individual 2 event 2 conditional on u, the sigx*/

  mat us;
  mat ui;
  double squd;
};

/*
Funtion that calculates normal cumulative distribution function
 */
double pn(mat y, double mu, double sigma) {
  return(Rf_pnorm5(y[0],mu,sqrt(sigma),1,0));
}

// [[Rcpp::export]]
double pn(mat y, mat mu, mat sigma) {
  /*int n = y.n_rows;*/
  int k = y.n_elem;
  double res;
  if (k==1) {
    res = Rf_pnorm5(y[0],mu[0],sqrt(sigma[0]),1,0);
    return(res);
  }

  mat L = mat(2,2); L.fill(0.0);
  L(0,0) = 1/sqrt(sigma(0,0));
  L(1,1) = 1/sqrt(sigma(1,1));
  y = L*(y-mu);
  double r = sigma(0,1)/(sqrt(sigma(0,0)*sigma(1,1)));

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
Miscellaneous things related to the conditional variance-covariance matrix
*/
matfour condsig(mat sigma, uvec rc1, uvec rc2) {
  /* Decomposing sigma */
  mat sig11 = sigma.submat(rc1,rc1);
  mat sig12 = sigma.submat(rc1,rc2);
  mat sig21 = sigma.submat(rc2,rc1);
  mat sig22 = sigma.submat(rc2,rc2);

  mat isig22 = sig22.i();
  mat xc_sig = sig12*isig22;

  /* Conditional variance-covariance matrix */
  mat c_sig = sig11-xc_sig*sig21;

  /* Inverse */
  mat ic_sig = c_sig.i();

  /* The determinant */
  double dc_sig = det(c_sig);

  /* Output */
  matfour out;
  out.M1 = c_sig;
  out.M2 = ic_sig;
  out.M3 = dc_sig;
  out.M4 = xc_sig;
  return(out);
}

/*
Full loglikelihood
*/
/*// [[Rcpp::export]]*/
vec loglikfull(mat y, mat b, mat u, ss condsigma, mat alph, mat dalph, mat tau, bool full=1){
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

  /* Specifying which parts of mu to use */
  uvec rc1_1u(3); rc1_1u(0)=0; rc1_1u(1)=4; rc1_1u(2)=5;
  uvec rc1_2u(3); rc1_2u(0)=1; rc1_2u(1)=4; rc1_2u(2)=5;
  uvec rc2_1u(3); rc2_1u(0)=2; rc2_1u(1)=4; rc2_1u(2)=5;
  uvec rc2_2u(3); rc2_2u(0)=3; rc2_2u(1)=4; rc2_2u(2)=5;
  uvec rcu(2); rcu(0)=4; rcu(1)=5;

  /* Initialising loglik vector */
  vec res(n);

  for (int i=0; i<n; i++) {

    /* Mean vector (used for estimation of conditional mean) */
    vec mu(6) ;
    mu(0) = alph(i,0);
    mu(1) = alph(i,1);
    mu(2) = alph(i,2);
    mu(3) = alph(i,3);
    mu(4) = u(i,0);
    mu(5) = u(i,1);

    if ((tau(i,0) == 1) & (tau(i,1) == 1)){
      res(i) = log(1-pi1_1(i)-pi2_1(i))+log(1-pi1_2(i)-pi2_2(i));
    }
    else if ((tau(i,0) == 1) & (tau(i,1) == 0)){
      /* Family member 1 experience event 0, family member 2 experience event 1, estimating dF01 */
      if((y(i,0) == 0) & (y(i,1) == 1)){

	/* Marginal dF1_2 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e1_2s;
	mat ic_sig = condsigma.e1_2i; // the inverse
	double dc_sig = condsigma.e1_2d; // the determinant
	mat sigX = condsigma.e1_2x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub(n,1);
	alph_sub.col(0) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f = alph_sub.t();

	/* Centering the alpha */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf1_2 = log(1/sq_twopi)+log(1/sqrt(dc_sig))+log(dalph(i,1))-0.5*inner(0);

	/* Estimation of dF1_2 */
	double logdF1_2 = log(pi1_2(i))+logpdf1_2; // pi1_2, dalph1_2

	/* Loglikelihood contribution */
	double logdF01 = logdF1_2+log(1-pi1_1(i)-pi2_1(i));
	res(i) = logdF01;
      }
      /* Family member 1 experience event 0, family member 2 experience event 2, estimating dF02 */
      else if((y(i,0) == 0) & (y(i,1) == 2)){

	/* Marginal dF2_2 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e2_2s;
	mat ic_sig = condsigma.e2_2i; // the inverse
	double dc_sig = condsigma.e2_2d; // the determinant
	mat sigX = condsigma.e2_2x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub(n,1);
	alph_sub.col(0) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f = alph_sub.t();

	/* Centering the alpha */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf2_2 = log(1/sq_twopi)+log(1/sqrt(dc_sig))+log(dalph(i,3))-0.5*inner(0);

	/* Estimation of dF2_2 */
	double logdF2_2 = log(pi2_2(i))+logpdf2_2; // pi2_2, dalph2_2

	/* Loglikelihood contribution */
	double logdF02 = logdF2_2+log(1-pi1_1(i)-pi2_1(i));
	res(i) = logdF02;
      }
      else {
	/* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_2s;
	mat sigX1 = condsigma.e1_2x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	mat c_sig2 = condsigma.e2_2s;
	mat sigX2 = condsigma.e2_2x; // the sigx
	vec a2 = mu.elem(rcu);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub1_2(n,1);
	mat alph_sub2_2(n,1);
	alph_sub1_2.col(0) = alph.col(1); // alph1_2
	alph_sub2_2.col(0) = alph.col(3); // alph2_2

	/* Transposing matrices to be of the form 1xn */
	mat alph_f1_2 = alph_sub1_2.t();
	mat alph_f2_2 = alph_sub2_2.t();

	/* Centering the alphas */
	mat alph_c1_2 = alph_f1_2.col(i)-c_mu1;
	mat alph_c2_2 = alph_f2_2.col(i)-c_mu2;

	/* Estimation of F1_2 and F2_2 */
	double F1_2 = pi1_2(i)*pn(alph_c1_2,0.0,c_sig1[0]);
	double F2_2 = pi2_2(i)*pn(alph_c2_2,0.0,c_sig2[0]);

	/* Loglikelihood contribution */
	double logF00 = log(1-pi1_1(i)-pi2_1(i))+log(1-F1_2-F2_2);
	res(i) = logF00;
      }
    }
    else if ((tau(i,0) == 0) & (tau(i,1) == 1)){
      /* Family member 1 experience event 1, family member 2 experience event 0, estimating dF10 */
      if((y(i,0) == 1) & (y(i,1) == 0)){

	/* Marginal dF1_1 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e1_1s;
	mat ic_sig = condsigma.e1_1i; // the inverse
	double dc_sig = condsigma.e1_1d; // the determinant
	mat sigX = condsigma.e1_1x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub(n,1);
	alph_sub.col(0) = alph.col(0); // alph1_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f = alph_sub.t();

	/* Centering the alpha */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf1_1 = log(1/sq_twopi)+log(1/sqrt(dc_sig))+log(dalph(i,0))-0.5*inner(0);

	/* Estimation of dF1_1 */
	double logdF1_1 = log(pi1_1(i))+logpdf1_1; // pi1_1, dalph1_1

	/* Loglikelihood contribution */
	double logdF10 = logdF1_1+log(1-pi1_2(i)-pi2_2(i));
	res(i) = logdF10;
      }
      /* Family member 1 experience event 2, family member 2 experience event 0, estimating dF20 */
      else if((y(i,0) == 2) & (y(i,1) == 0)){

	/* Marginal dF2_1 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e2_1s;
	mat ic_sig = condsigma.e2_1i; // the inverse
	double dc_sig = condsigma.e2_1d; // the determinant
	mat sigX = condsigma.e2_1x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub(n,1);
	alph_sub.col(0) = alph.col(2); // alph2_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f = alph_sub.t();

	/* Centering the alpha */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf2_1 = log(1/sq_twopi)+log(1/sqrt(dc_sig))+log(dalph(i,2))-0.5*inner(0);

	/* Estimation of dF2_1 */
	double logdF2_1 = log(pi2_1(i))+logpdf2_1; // pi2_1, dalph2_1

	/* Loglikelihood contribution */
	double logdF20 = logdF2_1+log(1-pi1_2(i)-pi2_2(i));
	res(i) = logdF20;
      }
      else {
	/* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_1s;
	mat sigX1 = condsigma.e1_1x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	mat c_sig2 = condsigma.e2_1s;
	mat sigX2 = condsigma.e2_1x; // the sigx
	vec a2 = mu.elem(rcu);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub1_1(n,1);
	mat alph_sub2_1(n,1);
	alph_sub1_1.col(0) = alph.col(0); // alph1_1
	alph_sub2_1.col(0) = alph.col(2); // alph2_1

	/* Transposing matrices to be of the form 1xn */
	mat alph_f1_1 = alph_sub1_1.t();
	mat alph_f2_1 = alph_sub2_1.t();

	/* Centering the alphas */
	mat alph_c1_1 = alph_f1_1.col(i)-c_mu1;
	mat alph_c2_1 = alph_f2_1.col(i)-c_mu2;

	/* Estimation of F1_1, F1_2, F2_1 and F2_2 */
	double F1_1 = pi1_1(i)*pn(alph_c1_1,0.0,c_sig1[0]);
	double F2_1 = pi2_1(i)*pn(alph_c2_1,0.0,c_sig2[0]);

	/* Loglikelihood contribution */
	double logF00 = log(1-pi1_2(i)-pi2_2(i))+log(1-F1_1-F2_1);
	res(i) = logF00;
      }
    }
    else {
      /* Both family members experience event 1, estimating ddF11 */
      if((y(i,0) == 1) & (y(i,1) == 1)){

	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e11s;
	mat ic_sig = condsigma.e11i; // the inverse
	double dc_sig = condsigma.e11d; // the determinant
	mat sigX = condsigma.e11x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub(n,2);
	alph_sub.col(0) = alph.col(0); // alph1_1
	alph_sub.col(1) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 2xn */
	mat alph_f = alph_sub.t();

	/* Centering the alphas */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf of alpha */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf11 = log(1/twopi)+log(1/sqrt(dc_sig))+log(dalph(i,0))+log(dalph(i,1))-0.5*inner(0);

	/* Loglikelihood contribution */
	double logddF11 = log(pi1_1(i))+log(pi1_2(i))+logpdf11; // pi1_1, pi1_2, dalph1_1, dalph1_2

	res(i) = logddF11;
      }
      /* Family member 1 experience event 1, family member 2 experience event 2, estimating ddF12 */
      else if((y(i,0) == 1) & (y(i,1) == 2)){

	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e12s;
	mat ic_sig = condsigma.e12i; // the inverse
	double dc_sig = condsigma.e12d; // the determinant
	mat sigX = condsigma.e12x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub(n,2);
	alph_sub.col(0) = alph.col(0); // alph1_1
	alph_sub.col(1) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 2xn */
	mat alph_f = alph_sub.t();

	/* Centering the alphas */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf of alpha*/
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf12 = log(1/twopi)+log(1/sqrt(dc_sig))+log(dalph(i,0)*dalph(i,3))-0.5*inner(0);

	/* Loglikelihood contribution */
	double logddF12 = log(pi1_1(i))+log(pi2_2(i))+logpdf12; // pi1_1, pi2_2, dalph1_1, dalph2_2
	res(i) = logddF12;
      }
      /* Family member 1 experience event 2, family member 2 experience event 1, estimating ddF21 */
      else if((y(i,0) == 2) & (y(i,1) == 1)){

	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e21s;
	mat ic_sig = condsigma.e21i; // the inverse
	double dc_sig = condsigma.e21d; // the determinant
	mat sigX = condsigma.e21x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub(n,2);
	alph_sub.col(0) = alph.col(2); // alph2_1
	alph_sub.col(1) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 2xn */
	mat alph_f = alph_sub.t();
	/* Centering the alphas */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf21 = log(1/twopi)+log(1/sqrt(dc_sig))+log(dalph(i,2)*dalph(i,1))-0.5*inner(0);

	/* Loglikelihood contribution */
	double logddF21 = log(pi2_1(i))+log(pi1_2(i))+logpdf21; // pi2_1, pi1_2, dalph2_1, dalph1_2
	res(i) = logddF21;
      }
      /* Both family members experience event 2, estimating ddF22 */
      else if((y(i,0) == 2) & (y(i,1) == 2)){

	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e22s;
	mat ic_sig = condsigma.e22i; // the inverse
	double dc_sig = condsigma.e22d; // the determinant
	mat sigX = condsigma.e22x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub(n,2);
	alph_sub.col(0) = alph.col(2); // alph2_1
	alph_sub.col(1) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 2xn */
	mat alph_f = alph_sub.t();
	/* Centering the alphas */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf of alpha */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf22 = log(1/twopi)+log(1/sqrt(dc_sig))+log(dalph(i,2)*dalph(i,3))-0.5*inner(0);

	/* Loglikelihood contribution */
	double logddF22 = log(pi2_1(i))+log(pi2_2(i))+logpdf22; // pi2_1, pi2_2, dalph2_1, dalph2_2
	res(i) = logddF22;
      }
      /* Family member 1 experience event 0, family member 2 experience event 1, estimating dF01 */
      else if((y(i,0) == 0) & (y(i,1) == 1)){

	/* Marginal dF1_2 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e1_2s;
	mat ic_sig = condsigma.e1_2i; // the inverse
	double dc_sig = condsigma.e1_2d; // the determinant
	mat sigX = condsigma.e1_2x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub(n,1);
	alph_sub.col(0) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f = alph_sub.t();

	/* Centering the alpha */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf1_2 = log(1/sq_twopi)+log(1/sqrt(dc_sig))+log(dalph(i,1))-0.5*inner(0);

	/* Estimation of dF1_2 */
	double logdF1_2 = log(pi1_2(i))+logpdf1_2; // pi1_2, dalph1_2

	/* Conditional F1c1_1 */
	/* Conditional mean and variance-covariance matrix, conditional on alph1_2, u1 and u2 */
	mat c_sig2 = condsigma.e1c1_1s;
	mat sigX2 = condsigma.e1c1_1x; // the sigx
	vec a2 = mu.elem(rc1_2u);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub2(n,1);
	alph_sub2.col(0) = alph.col(0); // alph1_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f2 = alph_sub2.t();

	/* Centering the alpha */
	mat alph_c2 = alph_f2.col(i)-c_mu2;

	/* Estimation of F1c1_1 */
	double F1c1_1 = pi1_1(i)*pn(alph_c2,0.0,c_sig2[0]);

	/* Conditional F2c1_1 */
	/* Conditional mean and variance-covariance matrix, conditional on alph1_2, u1 and u2 */
	mat c_sig3 = condsigma.e2c1_1s;
	mat sigX3 = condsigma.e2c1_1x; // the sigx
	vec a3 = mu.elem(rc1_2u);
	vec c_mu3 = sigX3*a3;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub3(n,1);
	alph_sub3.col(0) = alph.col(2); // alph2_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f3 = alph_sub3.t();

	/* Centering the alpha */
	mat alph_c3 = alph_f3.col(i)-c_mu3;

	/* Estimation of F2c1_1 */
	double F2c1_1 = pi2_1(i)*pn(alph_c3,0.0,c_sig3[0]);

	/* Loglikelihood contribution */
	double logdF01 = logdF1_2+log(1-F1c1_1-F2c1_1);
	res(i) = logdF01;
      }
      /* Family member 1 experience event 1, family member 2 experience event 0, estimating dF10 */
      else if((y(i,0) == 1) & (y(i,1) == 0)){

	/* Marginal dF1_1 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e1_1s;
	mat ic_sig = condsigma.e1_1i; // the inverse
	double dc_sig = condsigma.e1_1d; // the determinant
	mat sigX = condsigma.e1_1x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub(n,1);
	alph_sub.col(0) = alph.col(0); // alph1_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f = alph_sub.t();

	/* Centering the alpha */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf1_1 = log(1/sq_twopi)+log(1/sqrt(dc_sig))+log(dalph(i,0))-0.5*inner(0);

	/* Estimation of dF1_1 */
	double logdF1_1 = log(pi1_1(i))+logpdf1_1; // pi1_1, dalph1_1

	/* Conditional F1c1_2 */
	/* Conditional mean and variance-covariance matrix, conditional on alph1_1, u1 and u2 */
	mat c_sig2 = condsigma.e1c1_2s;
	mat sigX2 = condsigma.e1c1_2x; // the sigx
	vec a2 = mu.elem(rc1_1u);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub2(n,1);
	alph_sub2.col(0) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f2 = alph_sub2.t();

	/* Centering the alpha */
	mat alph_c2 = alph_f2.col(i)-c_mu2;

	/* Estimation of F1c1_2 */
	double F1c1_2 = pi1_2(i)*pn(alph_c2,0.0,c_sig2[0]);

	/* Conditional F2c1_2 */
	/* Conditional mean and variance-covariance matrix, conditional on alph1_1, u1 and u2 */
	mat c_sig3 = condsigma.e2c1_2s;
	mat sigX3 = condsigma.e2c1_2x; // the sigx
	vec a3 = mu.elem(rc1_1u);
	vec c_mu3 = sigX3*a3;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub3(n,1);
	alph_sub3.col(0) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f3 = alph_sub3.t();

	/* Centering the alpha */
	mat alph_c3 = alph_f3.col(i)-c_mu3;

	/* Estimation of F2c1_2 */
	double F2c1_2 = pi2_2(i)*pn(alph_c3,0.0,c_sig3[0]);

	/* Loglikelihood contribution */
	double logdF10 = logdF1_1+log(1-F1c1_2-F2c1_2);
	res(i) = logdF10;
      }
      /* Family member 1 experience event 0, family member 2 experience event 2, estimating dF02 */
      else if((y(i,0) == 0) & (y(i,1) == 2)){

	/* Marginal dF2_2 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e2_2s;
	mat ic_sig = condsigma.e2_2i; // the inverse
	double dc_sig = condsigma.e2_2d; // the determinant
	mat sigX = condsigma.e2_2x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub(n,1);
	alph_sub.col(0) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f = alph_sub.t();

	/* Centering the alpha */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf2_2 = log(1/sq_twopi)+log(1/sqrt(dc_sig))+log(dalph(i,3))-0.5*inner(0);

	/* Estimation of dF2_2 */
	double logdF2_2 = log(pi2_2(i))+logpdf2_2; // pi2_2, dalph2_2

	/* Conditional F1c2_1 */
	/* Conditional mean and variance-covariance matrix, conditional on alph2_2, u1 and u2 */
	mat c_sig2 = condsigma.e1c2_1s;
	mat sigX2 = condsigma.e1c2_1x; // the sigx
	vec a2 = mu.elem(rc2_2u);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub2(n,1);
	alph_sub2.col(0) = alph.col(0); // alph1_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f2 = alph_sub2.t();

	/* Centering the alpha */
	mat alph_c2 = alph_f2.col(i)-c_mu2;

	/* Estimation of F1c2_1 */
	double F1c2_1 = pi1_1(i)*pn(alph_c2,0.0,c_sig2[0]);

	/* Conditional F2c2_1 */
	/* Conditional mean and variance-covariance matrix, conditional on alph2_2, u1 and u2 */
	mat c_sig3 = condsigma.e2c2_1s;
	mat sigX3 = condsigma.e2c2_1x; // the sigx
	vec a3 = mu.elem(rc2_2u);
	vec c_mu3 = sigX3*a3;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub3(n,1);
	alph_sub3.col(0) = alph.col(2); // alph2_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f3 = alph_sub3.t();

	/* Centering the alpha */
	mat alph_c3 = alph_f3.col(i)-c_mu3;

	/* Estimation of F2c2_1 */
	double F2c2_1 = pi2_1(i)*pn(alph_c3,0.0,c_sig3[0]);

	/* Loglikelihood contribution */
	double logdF02 = logdF2_2+log(1-F1c2_1-F2c2_1);
	res(i) = logdF02;
      }
      /* Family member 1 experience event 2, family member 2 experience event 0, estimating dF20 */
      else if((y(i,0) == 2) & (y(i,1) == 0)){

	/* Marginal dF2_1 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e2_1s;
	mat ic_sig = condsigma.e2_1i; // the inverse
	double dc_sig = condsigma.e2_1d; // the determinant
	mat sigX = condsigma.e2_1x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub(n,1);
	alph_sub.col(0) = alph.col(2); // alph2_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f = alph_sub.t();

	/* Centering the alpha */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Calculating the pdf */
	mat inner = alph_c.t()*ic_sig*alph_c;
	double logpdf2_1 = log(1/sq_twopi)+log(1/sqrt(dc_sig))+log(dalph(i,2))-0.5*inner(0);

	/* Estimation of dF2_1 */
	double logdF2_1 = log(pi2_1(i))+logpdf2_1; // pi2_1, dalph2_1

	/* Conditional F1c2_2 */
	/* Conditional mean and variance-covariance matrix, conditional on alph2_1, u1 and u2 */
	mat c_sig2 = condsigma.e1c2_2s;
	mat sigX2 = condsigma.e1c2_2x; // the sigx
	vec a2 = mu.elem(rc2_1u);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub2(n,1);
	alph_sub2.col(0) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f2 = alph_sub2.t();

	/* Centering the alpha */
	mat alph_c2 = alph_f2.col(i)-c_mu2;

	/* Estimation of F1c2_2 */
	double F1c2_2 = pi1_2(i)*pn(alph_c2,0.0,c_sig2[0]);

	/* Conditional F2c2_2 */
	/* Conditional mean and variance-covariance matrix, conditional on alph2_1, u1 and u2 */
	mat c_sig3 = condsigma.e2c2_2s;
	mat sigX3 = condsigma.e2c2_2x; // the sigx
	vec a3 = mu.elem(rc2_1u);
	vec c_mu3 = sigX3*a3;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub3(n,1);
	alph_sub3.col(0) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f3 = alph_sub3.t();

	/* Centering the alpha */
	mat alph_c3 = alph_f3.col(i)-c_mu3;

	/* Estimation of F2c2_2 */
	double F2c2_2 = pi2_2(i)*pn(alph_c3,0.0,c_sig3[0]);

	/* Loglikelihood contribution */
	double logdF20 = logdF2_1+log(1-F1c2_2-F2c2_2);
	res(i) = logdF20;
      }
      /* Family member 1 experience event 0, family member 2 experience event 0, estimating F00 */
      else{
	/* Marginal F1_1, F1_2, F2_1 and F2_2 */
	/* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_1s;
	mat sigX1 = condsigma.e1_1x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	mat c_sig2 = condsigma.e1_2s;
	mat sigX2 = condsigma.e1_2x; // the sigx
	vec a2 = mu.elem(rcu);
	vec c_mu2 = sigX2*a2;

	mat c_sig3 = condsigma.e2_1s;
	mat sigX3 = condsigma.e2_1x; // the sigx
	vec a3 = mu.elem(rcu);
	vec c_mu3 = sigX3*a3;

	mat c_sig4 = condsigma.e2_2s;
	mat sigX4 = condsigma.e2_2x; // the sigx
	vec a4 = mu.elem(rcu);
	vec c_mu4 = sigX4*a4;

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
	mat alph_c1_2 = alph_f1_2.col(i)-c_mu2;
	mat alph_c2_1 = alph_f2_1.col(i)-c_mu3;
	mat alph_c2_2 = alph_f2_2.col(i)-c_mu4;

	/* Estimation of F1_1, F1_2, F2_1 and F2_2 */
	double F1_1 = pi1_1(i)*pn(alph_c1_1,0.0,c_sig1[0]);
	double F1_2 = pi1_2(i)*pn(alph_c1_2,0.0,c_sig2[0]);
	double F2_1 = pi2_1(i)*pn(alph_c2_1,0.0,c_sig3[0]);
	double F2_2 = pi2_2(i)*pn(alph_c2_2,0.0,c_sig4[0]);

	/* Joint probabilities F11, F12, F21 and F22 */
	/* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
	mat c_sig5 = condsigma.e11s;
	mat sigX5 = condsigma.e11x; // the sigx
	vec a5 = mu.elem(rcu);
	vec c_mu5 = sigX5*a5;

	mat c_sig6 = condsigma.e12s;
	mat sigX6 = condsigma.e12x; // the sigx
	vec a6 = mu.elem(rcu);
	vec c_mu6 = sigX6*a6;

	mat c_sig7 = condsigma.e21s;
	mat sigX7 = condsigma.e21x; // the sigx
	vec a7 = mu.elem(rcu);
	vec c_mu7 = sigX7*a7;

	mat c_sig8 = condsigma.e22s;
	mat sigX8 = condsigma.e22x; // the sigx
	vec a8 = mu.elem(rcu);
	vec c_mu8 = sigX8*a8;

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

	/* Estimating F11, F12, F21 and F22 */
	double F11 = pi1_1(i)*pi1_2(i)*pn(alph_f11.col(i),c_mu5,c_sig5);
	double F12 = pi1_1(i)*pi2_2(i)*pn(alph_f12.col(i),c_mu6,c_sig6);
	double F21 = pi2_1(i)*pi1_2(i)*pn(alph_f21.col(i),c_mu7,c_sig7);
	double F22 = pi2_1(i)*pi2_2(i)*pn(alph_f22.col(i),c_mu8,c_sig8);

	/* Loglikelihood contribution */
	double F00 = (1-F1_1-F1_2-F2_1-F2_2+F11+F12+F21+F22);
	res(i) = log(F00);
      }
    }
    if (full) {
      /* For the pdf of u */
      mat sigu = condsigma.us;
      mat isigu = condsigma.ui;
      double sq_dsigu = condsigma.squd;
      mat pu = u.row(i);
      mat inu = pu*isigu*pu.t();

      /* The pdf of u */
      double logpdfu = log(1/twopi)+log(1/sq_dsigu)-0.5*inu(0);

      /* Adding to loglik */
      res(i) += logpdfu;
    }
  }
  return(res);
}

/*
Score function of full loglikelihood
*/
/*// [[Rcpp::export]]*/
mat Dloglikfull(mat y, mat b, mat u, ss condsigma, mat alph, mat dalph, mat tau, bool full=1){
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

  /* Defining denominator */
  vec denom1 = (1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1)));
  vec denom2 = (1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1)));

  /* Derivatives of the pis wrt. u1 and u2 */
  vec dpi1_1_u1 = (exp(b.col(0)+u.col(0)+b.col(2)+u.col(1))+exp(b.col(0)+u.col(0)))/pow(denom1,2);
  vec dpi1_2_u1 = (exp(b.col(1)+u.col(0)+b.col(3)+u.col(1))+exp(b.col(1)+u.col(0)))/pow(denom2,2);
  vec dpi2_1_u1 = -exp(b.col(0)+u.col(0)+b.col(2)+u.col(1))/pow(denom1,2);
  vec dpi2_2_u1 = -exp(b.col(1)+u.col(0)+b.col(3)+u.col(1))/pow(denom2,2);

  vec dpi1_1_u2 = -exp(b.col(0)+u.col(0)+b.col(2)+u.col(1))/pow(denom1,2);
  vec dpi1_2_u2 = -exp(b.col(1)+u.col(0)+b.col(3)+u.col(1))/pow(denom2,2);
  vec dpi2_1_u2 = (exp(b.col(0)+u.col(0)+b.col(2)+u.col(1))+exp(b.col(2)+u.col(1)))/pow(denom1,2);
  vec dpi2_2_u2 = (exp(b.col(1)+u.col(0)+b.col(3)+u.col(1))+exp(b.col(3)+u.col(1)))/pow(denom2,2);

  /* Derivatives of the product of two pis wrt. u1 and u2 */
  vec dpi11_u1 = dpi1_1_u1%pi1_2 + pi1_1%dpi1_2_u1;
  vec dpi12_u1 = dpi1_1_u1%pi2_2 + pi1_1%dpi2_2_u1;
  vec dpi21_u1 = dpi2_1_u1%pi1_2 + pi2_1%dpi1_2_u1;
  vec dpi22_u1 = dpi2_1_u1%pi2_2 + pi2_1%dpi2_2_u1;

  vec dpi11_u2 = dpi1_1_u2%pi1_2 + pi1_1%dpi1_2_u2;
  vec dpi12_u2 = dpi1_1_u2%pi2_2 + pi1_1%dpi2_2_u2;
  vec dpi21_u2 = dpi2_1_u2%pi1_2 + pi2_1%dpi1_2_u2;
  vec dpi22_u2 = dpi2_1_u2%pi2_2 + pi2_1%dpi2_2_u2;

  /* Specifying which parts of mu to use */
  uvec rc1_1u(3); rc1_1u(0)=0; rc1_1u(1)=4; rc1_1u(2)=5;
  uvec rc1_2u(3); rc1_2u(0)=1; rc1_2u(1)=4; rc1_2u(2)=5;
  uvec rc2_1u(3); rc2_1u(0)=2; rc2_1u(1)=4; rc2_1u(2)=5;
  uvec rc2_2u(3); rc2_2u(0)=3; rc2_2u(1)=4; rc2_2u(2)=5;
  uvec rcu(2); rcu(0)=4; rcu(1)=5;

  /* Initialising Dloglik matrix */
  mat res(n,2);

  for (int i=0; i<n; i++) {

    /* Mean vector (used for estimation of conditional mean) */
    vec mu(6) ;
    mu(0) = alph(i,0);
    mu(1) = alph(i,1);
    mu(2) = alph(i,2);
    mu(3) = alph(i,3);
    mu(4) = u(i,0);
    mu(5) = u(i,1);

    if ((tau(i,0) == 1) & (tau(i,1) == 1)){

      /* Likelihood contribution */
      double F00 = (1-pi1_1(i)-pi2_1(i))*(1-pi1_2(i)-pi2_2(i));

      /* Score contribution from F00 */
      double sc_u1 = (1/F00)*(-dpi1_1_u1(i)-dpi1_2_u1(i)-dpi2_1_u1(i)-dpi2_2_u1(i)+dpi11_u1(i)+dpi12_u1(i)+dpi21_u1(i)+dpi22_u1(i));
      double sc_u2 = (1/F00)*(-dpi1_1_u2(i)-dpi1_2_u2(i)-dpi2_1_u2(i)-dpi2_2_u2(i)+dpi11_u2(i)+dpi12_u2(i)+dpi21_u2(i)+dpi22_u2(i));

      /* Adding to return vector */
      res(i,0) = sc_u1;
      res(i,1) = sc_u2;
    }
    else if ((tau(i,0) == 1) & (tau(i,1) == 0)){
      /* Family member 1 experience event 0, family member 2 experience event 1, estimating score contribution from dF01 */
      if((y(i,0) == 0) & (y(i,1) == 1)){

	/* Marginal dF1_2 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_2s;
	mat ic_sig1 = condsigma.e1_2i; // the inverse
	mat sigX1 = condsigma.e1_2x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub1(n,1);
	alph_sub1.col(0) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f1 = alph_sub1.t();

	/* Centering the alpha */
	mat alph_c1 = alph_f1.col(i)-c_mu1;

	/* Likelihood contribution form F1c1_1 and F2c1_1 */
	double F1c1_1 = pi1_1(i);
	double F2c1_1 = pi2_1(i);

	/* Score contribution from F1c1_1 */
	double dF1c1_1_u1 = dpi1_1_u1(i);
	double dF1c1_1_u2 = dpi1_1_u2(i);

	/* Score contribution from F2c1_1 */
	double dF2c1_1_u1 = dpi2_1_u1(i);
	double dF2c1_1_u2 = dpi2_1_u2(i);

	/* Score contributions from dF01 wrt. u1 and u2 */
	double sc_u1 = (dpi1_2_u1(i)+pi1_2(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(0))/pi1_2(i)-(dF1c1_1_u1+dF2c1_1_u1)/(1-F1c1_1-F2c1_1);
	double sc_u2 = (dpi1_2_u2(i)+pi1_2(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(1))/pi1_2(i)-(dF1c1_1_u2+dF2c1_1_u2)/(1-F1c1_1-F2c1_1);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Family member 1 experience event 0, family member 2 experience event 2, estimating dF02 */
      else if((y(i,0) == 0) & (y(i,1) == 2)){

	/* Marginal dF2_2 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e2_2s;
	mat ic_sig1 = condsigma.e2_2i; // the inverse
	mat sigX1 = condsigma.e2_2x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub1(n,1);
	alph_sub1.col(0) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f1 = alph_sub1.t();

	/* Centering the alpha */
	mat alph_c1 = alph_f1.col(i)-c_mu1;

	/* Likelihood contributions from F1c2_1 and F2c2_1 */
	double F1c2_1 = pi1_1(i);
	double F2c2_1 = pi2_1(i);

	/* Score contribution from F1c2_1 */
	double dF1c2_1_u1 = dpi1_1_u1(i);
	double dF1c2_1_u2 = dpi1_1_u2(i);

	/* Score contribution from F2c2_1 */
	double dF2c2_1_u1 = dpi2_1_u1(i);
	double dF2c2_1_u2 = dpi2_1_u2(i);

	/* Score contributions from dF01 wrt. u1 and u2 */
	double sc_u1 = (dpi2_2_u1(i)+pi2_2(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(0))/pi2_2(i)-(dF1c2_1_u1+dF2c2_1_u1)/(1-F1c2_1-F2c2_1);
	double sc_u2 = (dpi2_2_u2(i)+pi2_2(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(1))/pi2_2(i)-(dF1c2_1_u2+dF2c2_1_u2)/(1-F1c2_1-F2c2_1);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      else {
	/* Marginal F1_2 and F2_2 */
	/* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_2s;
	mat ic_sig1 = condsigma.e1_2i; // the inverse
	mat sigX1 = condsigma.e1_2x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;
	double sd1 = sqrt(c_sig1[0]);

	mat c_sig2 = condsigma.e2_2s;
	mat ic_sig2 = condsigma.e2_2i; // the inverse
	mat sigX2 = condsigma.e2_2x; // the sigx
	vec a2 = mu.elem(rcu);
	vec c_mu2 = sigX2*a2;
	double sd2 = sqrt(c_sig2[0]);

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub1_2(n,1);
	mat alph_sub2_2(n,1);
	alph_sub1_2.col(0) = alph.col(1); // alph1_2
	alph_sub2_2.col(0) = alph.col(3); // alph2_2

	/* Transposing matrices to be of the form 1xn */
	mat alph_f1_2 = alph_sub1_2.t();
	mat alph_f2_2 = alph_sub2_2.t();

	/* Centering the alphas */
	mat alph_c1_2 = alph_f1_2.col(i)-c_mu1;
	mat alph_c2_2 = alph_f2_2.col(i)-c_mu2;

	/* Estimation of F1_2 and F2_2 */
	double p1 = pn(alph_c1_2,0.0,sd1*sd1);
	double p2 = pn(alph_c2_2,0.0,sd2*sd2);

	double F1_2 = pi1_2(i)*p1;
	double F2_2 = pi2_2(i)*p2;

	/* Calculating the pdf of F1_2 and F2_2 */
	double inner1 = pow(alph_c1_2[0],2)/c_sig1[0];
	double inner2 = pow(alph_c2_2[0],2)/c_sig2[0];

	double pdf1_2 = 1/sq_twopi*1/sd1*exp(-0.5*inner1);
	double pdf2_2 = 1/sq_twopi*1/sd2*exp(-0.5*inner2);

	/* Score contributions from F1_2 and F2_2 */
	double dF1_2_u1 = dpi1_2_u1(i)*p1+pi1_2(i)*pdf1_2*(-sigX1(0));
	double dF2_2_u1 = dpi2_2_u1(i)*p2+pi2_2(i)*pdf2_2*(-sigX2(0));

	double dF1_2_u2 = dpi1_2_u2(i)*p1+pi1_2(i)*pdf1_2*(-sigX1(1));
	double dF2_2_u2 = dpi2_2_u2(i)*p2+pi2_2(i)*pdf2_2*(-sigX2(1));

	/* Likelihood contributions from F1_1 and F2_1 */
	double F1_1 = pi1_1(i);
	double F2_1 = pi2_1(i);

	/* Score contribution from F1_1 */
	double dF1_1_u1 = dpi1_1_u1(i);
	double dF1_1_u2 = dpi1_1_u2(i);

	/* Score contribution from F2_1 */
	double dF2_1_u1 = dpi2_1_u1(i);
	double dF2_1_u2 = dpi2_1_u2(i);

	/* Score contributions from dF00 wrt. u1 and u2 */
	double sc_u1 = 1/(1-F1_1-F2_1)*(-dF1_1_u1-dF2_1_u1)+1/(1-F1_2-F2_2)*(-dF1_2_u1-dF2_2_u1);
	double sc_u2 = 1/(1-F1_1-F2_1)*(-dF1_1_u2-dF2_1_u2)+1/(1-F1_2-F2_2)*(-dF1_2_u2-dF2_2_u2);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
    }
    else if ((tau(i,0) == 0) & (tau(i,1) == 1)){
      /* Family member 1 experience event 1, family member 2 experience event 0, estimating dF10 */
      if((y(i,0) == 1) & (y(i,1) == 0)){

	/* Marginal dF1_1 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_1s;
	mat ic_sig1 = condsigma.e1_1i; // the inverse
	mat sigX1 = condsigma.e1_1x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub1(n,1);
	alph_sub1.col(0) = alph.col(0); // alph1_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f1 = alph_sub1.t();

	/* Centering the alpha */
	mat alph_c1 = alph_f1.col(i)-c_mu1;

	/* Likelihood contribution from F1c1_2 and F2c1_2 */
	double F1c1_2 = pi1_2(i);
	double F2c1_2 = pi2_2(i);

	/* Score contribution from F1c1_2 */
	double dF1c1_2_u1 = dpi1_2_u1(i);
	double dF1c1_2_u2 = dpi1_2_u2(i);

	/* Score contribution from F2c1_2 */
	double dF2c1_2_u1 = dpi2_2_u1(i);
	double dF2c1_2_u2 = dpi2_2_u2(i);

	/* Score contributions from dF01 wrt. u1 and u2 */
	double sc_u1 = (dpi1_1_u1(i)+pi1_1(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(0))/pi1_1(i)-(dF1c1_2_u1+dF2c1_2_u1)/(1-F1c1_2-F2c1_2);
	double sc_u2 = (dpi1_1_u2(i)+pi1_1(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(1))/pi1_1(i)-(dF1c1_2_u2+dF2c1_2_u2)/(1-F1c1_2-F2c1_2);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Family member 1 experience event 2, family member 2 experience event 0, estimating dF20 */
      else if((y(i,0) == 2) & (y(i,1) == 0)){

	/* Marginal dF2_1 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e2_1s;
	mat ic_sig1 = condsigma.e2_1i; // the inverse
	mat sigX1 = condsigma.e2_1x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub1(n,1);
	alph_sub1.col(0) = alph.col(2); // alph2_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f1 = alph_sub1.t();

	/* Centering the alpha */
	mat alph_c1 = alph_f1.col(i)-c_mu1;

	/* Likelihood contribution from F1c2_2 and F2c2_2 */
	double F1c2_2 = pi1_2(i);
	double F2c2_2 = pi2_2(i);

	/* Score contribution from F1c2_2 */
	double dF1c2_2_u1 = dpi1_2_u1(i);
	double dF1c2_2_u2 = dpi1_2_u2(i);

	/* Score contribution from F2c2_2 */
	double dF2c2_2_u1 = dpi2_2_u1(i);
	double dF2c2_2_u2 = dpi2_2_u2(i);

	/* Score contributions from dF01 wrt. u1 and u2 */
	double sc_u1 = (dpi2_1_u1(i)+pi2_1(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(0))/pi2_1(i)-(dF1c2_2_u1+dF2c2_2_u1)/(1-F1c2_2-F2c2_2);
	double sc_u2 = (dpi2_1_u2(i)+pi2_1(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(1))/pi2_1(i)-(dF1c2_2_u2+dF2c2_2_u2)/(1-F1c2_2-F2c2_2);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      else {
	/* Marginal F1_1 and F2_1 */
	/* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_1s;
	mat ic_sig1 = condsigma.e1_1i; // the inverse
	mat sigX1 = condsigma.e1_1x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;
	double sd1 = sqrt(c_sig1[0]);

	mat c_sig2 = condsigma.e2_1s;
	mat ic_sig2 = condsigma.e2_1i; // the inverse
	mat sigX2 = condsigma.e2_1x; // the sigx
	vec a2 = mu.elem(rcu);
	vec c_mu2 = sigX2*a2;
	double sd2 = sqrt(c_sig2[0]);

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub1_1(n,1);
	mat alph_sub2_1(n,1);
	alph_sub1_1.col(0) = alph.col(0); // alph1_1
	alph_sub2_1.col(0) = alph.col(2); // alph2_1

	/* Transposing matrices to be of the form 1xn */
	mat alph_f1_1 = alph_sub1_1.t();
	mat alph_f2_1 = alph_sub2_1.t();

	/* Centering the alphas */
	mat alph_c1_1 = alph_f1_1.col(i)-c_mu1;
	mat alph_c2_1 = alph_f2_1.col(i)-c_mu2;

	/* Estimation of F1_1 and F2_1 */
	double p1 = pn(alph_c1_1,0.0,sd1*sd1);
	double p2 = pn(alph_c2_1,0.0,sd2*sd2);

	double F1_1 = pi1_1(i)*p1;
	double F2_1 = pi2_1(i)*p2;

	/* Calculating the pdf of F1_1 and F2_1 */
	double inner1 = pow(alph_c1_1[0],2)/c_sig1[0];
	double inner2 = pow(alph_c2_1[0],2)/c_sig2[0];

	double pdf1_1 = 1/sq_twopi*1/sd1*exp(-0.5*inner1);
	double pdf2_1 = 1/sq_twopi*1/sd2*exp(-0.5*inner2);

	/* Score contributions from F1_1 and F2_1 */
	double dF1_1_u1 = dpi1_1_u1(i)*p1+pi1_1(i)*pdf1_1*(-sigX1(0));
	double dF2_1_u1 = dpi2_1_u1(i)*p2+pi2_1(i)*pdf2_1*(-sigX2(0));

	double dF1_1_u2 = dpi1_1_u2(i)*p1+pi1_1(i)*pdf1_1*(-sigX1(1));
	double dF2_1_u2 = dpi2_1_u2(i)*p2+pi2_1(i)*pdf2_1*(-sigX2(1));

	/* Likelihood contributions from F1_2 and F2_2 */
	double F1_2 = pi1_2(i);
	double F2_2 = pi2_2(i);

	/* Score contribution from F1_2 */
	double dF1_2_u1 = dpi1_2_u1(i);
	double dF1_2_u2 = dpi1_2_u2(i);

	/* Score contribution from F2_2 */
	double dF2_2_u1 = dpi2_2_u1(i);
	double dF2_2_u2 = dpi2_2_u2(i);

	/* Score contributions from F00 wrt. u1 and u2 */
	double sc_u1 = 1/(1-F1_2-F2_2)*(-dF1_2_u1-dF2_2_u1)+1/(1-F1_1-F2_1)*(-dF1_1_u1-dF2_1_u1);
	double sc_u2 = 1/(1-F1_2-F2_2)*(-dF1_2_u2-dF2_2_u2)+1/(1-F1_1-F2_1)*(-dF1_1_u2-dF2_1_u2);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
    }
    else {

      /* Both family members experience event 1, estimating score contribution from ddF11 */
      if((y(i,0) == 1) & (y(i,1) == 1)){

	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e11s;
	mat ic_sig = condsigma.e11i; // the inverse
	mat sigX = condsigma.e11x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Standard deviations */
	double sd1 = sqrt(c_sig(0,0));
	double sd2 = sqrt(c_sig(1,1));

	/* Correlation */
	double r = c_sig(0,1)/(sd1*sd2);

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub(n,2);
	alph_sub.col(0) = alph.col(0); // alph1_1
	alph_sub.col(1) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 2xn */
	mat alph_f = alph_sub.t();

	/* Centering the alphas */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Derivatives */
	double difflog = exp(-log(pi1_1(i))-log(pi1_2(i)));
	double a_u1 = (alph_c(0)*sigX(0,0)/pow(sd1,2)+alph_c(1)*sigX(1,0)/pow(sd2,2)-r*sigX(0,0)*alph_c(1)/(sd1*sd2)-r*alph_c(0)*sigX(1,0)/(sd1*sd2))/(1-pow(r,2));
	double a_u2 = (alph_c(0)*sigX(0,1)/pow(sd1,2)+alph_c(1)*sigX(1,1)/pow(sd2,2)-r*sigX(0,1)*alph_c(1)/(sd1*sd2)-r*alph_c(0)*sigX(1,1)/(sd1*sd2))/(1-pow(r,2));

	/* Score contributions from ddF11 wrt. u1 and u2 */
	double sc_u1 = dpi11_u1(i)*difflog+pi1_1(i)*pi1_2(i)*a_u1*difflog;
	double sc_u2 = dpi11_u2(i)*difflog+pi1_1(i)*pi1_2(i)*a_u2*difflog;

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Family member 1 experiences event 1, family member 2 experiences event 2, estimating score contribution from ddF12 */
      else if((y(i,0) == 1) & (y(i,1) == 2)){

	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e12s;
	mat ic_sig = condsigma.e12i; // the inverse
	mat sigX = condsigma.e12x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Standard deviations */
	double sd1 = sqrt(c_sig(0,0));
	double sd2 = sqrt(c_sig(1,1));

	/* Correlation */
	double r = c_sig(0,1)/(sd1*sd2);

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub(n,2);
	alph_sub.col(0) = alph.col(0); // alph1_1
	alph_sub.col(1) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 2xn */
	mat alph_f = alph_sub.t();

	/* Centering the alphas */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Derivatives */
	double difflog = exp(-log(pi1_1(i))-log(pi2_2(i)));
	double a_u1 = (alph_c(0)*sigX(0,0)/pow(sd1,2)+alph_c(1)*sigX(1,0)/pow(sd2,2)-r*sigX(0,0)*alph_c(1)/(sd1*sd2)-r*alph_c(0)*sigX(1,0)/(sd1*sd2))/(1-pow(r,2));
	double a_u2 = (alph_c(0)*sigX(0,1)/pow(sd1,2)+alph_c(1)*sigX(1,1)/pow(sd2,2)-r*sigX(0,1)*alph_c(1)/(sd1*sd2)-r*alph_c(0)*sigX(1,1)/(sd1*sd2))/(1-pow(r,2));

	/* Score contributions from ddF12 wrt. u1 and u2 */
	double sc_u1 = dpi12_u1(i)*difflog+pi1_1(i)*pi2_2(i)*a_u1*difflog;
	double sc_u2 = dpi12_u2(i)*difflog+pi1_1(i)*pi2_2(i)*a_u2*difflog;

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Family member 1 experiences event 2, family member 2 experiences event 1, estimating score contribution from ddF21 */
      else if((y(i,0) == 2) & (y(i,1) == 1)){

	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e21s;
	mat ic_sig = condsigma.e21i; // the inverse
	mat sigX = condsigma.e21x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Standard deviations */
	double sd1 = sqrt(c_sig(0,0));
	double sd2 = sqrt(c_sig(1,1));

	/* Correlation */
	double r = c_sig(0,1)/(sd1*sd2);

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub(n,2);
	alph_sub.col(0) = alph.col(2); // alph2_1
	alph_sub.col(1) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 2xn */
	mat alph_f = alph_sub.t();

	/* Centering the alphas */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Derivatives */
	double difflog = exp(-log(pi2_1(i))-log(pi1_2(i)));
	double a_u1 = (alph_c(0)*sigX(0,0)/pow(sd1,2)+alph_c(1)*sigX(1,0)/pow(sd2,2)-r*sigX(0,0)*alph_c(1)/(sd1*sd2)-r*alph_c(0)*sigX(1,0)/(sd1*sd2))/(1-pow(r,2));
	double a_u2 = (alph_c(0)*sigX(0,1)/pow(sd1,2)+alph_c(1)*sigX(1,1)/pow(sd2,2)-r*sigX(0,1)*alph_c(1)/(sd1*sd2)-r*alph_c(0)*sigX(1,1)/(sd1*sd2))/(1-pow(r,2));

	/* Score contributions from ddF21 wrt. u1 and u2 */
	double sc_u1 = dpi21_u1(i)*difflog+pi2_1(i)*pi1_2(i)*a_u1*difflog;
	double sc_u2 = dpi21_u2(i)*difflog+pi2_1(i)*pi1_2(i)*a_u2*difflog;

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Both family members experience event 2, estimating score contribution from ddF22 */
      else if((y(i,0) == 2) & (y(i,1) == 2)){

	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig = condsigma.e22s;
	mat ic_sig = condsigma.e22i; // the inverse
	mat sigX = condsigma.e22x; // the sigx
	vec a = mu.elem(rcu);
	vec c_mu = sigX*a;

	/* Standard deviations */
	double sd1 = sqrt(c_sig(0,0));
	double sd2 = sqrt(c_sig(1,1));

	/* Correlation */
	double r = c_sig(0,1)/(sd1*sd2);

	/* Pulling out the appropriate alphas from alph */
	mat alph_sub(n,2);
	alph_sub.col(0) = alph.col(2); // alph2_1
	alph_sub.col(1) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 2xn */
	mat alph_f = alph_sub.t();

	/* Centering the alphas */
	mat alph_c = alph_f.col(i)-c_mu;

	/* Derivatives */
	double difflog = exp(-log(pi2_1(i))-log(pi2_2(i)));
	double a_u1 = (alph_c(0)*sigX(0,0)/pow(sd1,2)+alph_c(1)*sigX(1,0)/pow(sd2,2)-r*sigX(0,0)*alph_c(1)/(sd1*sd2)-r*alph_c(0)*sigX(1,0)/(sd1*sd2))/(1-pow(r,2));
	double a_u2 = (alph_c(0)*sigX(0,1)/pow(sd1,2)+alph_c(1)*sigX(1,1)/pow(sd2,2)-r*sigX(0,1)*alph_c(1)/(sd1*sd2)-r*alph_c(0)*sigX(1,1)/(sd1*sd2))/(1-pow(r,2));

	/* Score contributions from ddF22 wrt. u1 and u2 */
	double sc_u1 = dpi22_u1(i)*difflog+pi2_1(i)*pi2_2(i)*a_u1*difflog;
	double sc_u2 = dpi22_u2(i)*difflog+pi2_1(i)*pi2_2(i)*a_u2*difflog;

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Family member 1 experience event 0, family member 2 experience event 1, estimating score contribution from dF01 */
      else if((y(i,0) == 0) & (y(i,1) == 1)){

	/* Marginal dF1_2 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_2s;
	mat ic_sig1 = condsigma.e1_2i; // the inverse
	mat sigX1 = condsigma.e1_2x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub1(n,1);
	alph_sub1.col(0) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f1 = alph_sub1.t();

	/* Centering the alpha */
	mat alph_c1 = alph_f1.col(i)-c_mu1;

	/* Conditional F1c1_1 */
	/* Conditional mean and variance-covariance matrix, conditional on alph1_2, u1 and u2 */
	mat c_sig2 = condsigma.e1c1_1s;
	mat ic_sig2 = condsigma.e1c1_1i; // the inverse
	mat sigX2 = condsigma.e1c1_1x; // the sigx
	vec a2 = mu.elem(rc1_2u);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub2(n,1);
	alph_sub2.col(0) = alph.col(0); // alph1_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f2 = alph_sub2.t();

	/* Centering the alpha */
	mat alph_c2 = alph_f2.col(i)-c_mu2;

	/* Estimation of F1c1_1 */
	double p1 = pn(alph_c2,0.0,c_sig2[0]);
	double F1c1_1 = pi1_1(i)*p1;

	/* Calculating the pdf */
	double sd2 = sqrt(c_sig2[0]);
	double inner2 = pow(alph_c2[0],2)/c_sig2[0];
	double pdf1c1_1 = 1/sq_twopi*1/sd2*exp(-0.5*inner2);

	/* Score contribution from F1c1_1 */
	double dF1c1_1_u1 = dpi1_1_u1(i)*p1+pi1_1(i)*pdf1c1_1*(-sigX2(1));
	double dF1c1_1_u2 = dpi1_1_u2(i)*p1+pi1_1(i)*pdf1c1_1*(-sigX2(2));

	/* Conditional F2c1_1 */
	/* Conditional mean and variance-covariance matrix, conditional on alph1_2, u1 and u2 */
	mat c_sig3 = condsigma.e2c1_1s;
	mat ic_sig3 = condsigma.e2c1_1i; // the inverse
	mat sigX3 = condsigma.e2c1_1x; // the sigx
	vec a3 = mu.elem(rc1_2u);
	vec c_mu3 = sigX3*a3;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub3(n,1);
	alph_sub3.col(0) = alph.col(2); // alph2_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f3 = alph_sub3.t();

	/* Centering the alpha */
	mat alph_c3 = alph_f3.col(i)-c_mu3;

	/* Estimation of F2c1_1 */
	double p3 = pn(alph_c3,0.0,c_sig3[0]);
	double F2c1_1 = pi2_1(i)*p3;

	/* Calculating the pdf */
	double sd3 = sqrt(c_sig3[0]);
	double inner3 = pow(alph_c3[0],2)/c_sig3[0];
	double pdf2c1_1 = 1/sq_twopi*1/sd3*exp(-0.5*inner3);

	/* Score contribution from F2c1_1 */
	double dF2c1_1_u1 = dpi2_1_u1(i)*p3+pi2_1(i)*pdf2c1_1*(-sigX3(1));
	double dF2c1_1_u2 = dpi2_1_u2(i)*p3+pi2_1(i)*pdf2c1_1*(-sigX3(2));

	/* Score contributions from dF01 wrt. u1 and u2 */
	double sc_u1 = (dpi1_2_u1(i)+pi1_2(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(0))/pi1_2(i)-(dF1c1_1_u1+dF2c1_1_u1)/(1-F1c1_1-F2c1_1);
	double sc_u2 = (dpi1_2_u2(i)+pi1_2(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(1))/pi1_2(i)-(dF1c1_1_u2+dF2c1_1_u2)/(1-F1c1_1-F2c1_1);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Family member 1 experience event 1, family member 2 experience event 0, estimating score contribution from dF10 */
      else if((y(i,0) == 1) & (y(i,1) == 0)){

	/* Marginal dF1_1 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_1s;
	mat ic_sig1 = condsigma.e1_1i; // the inverse
	mat sigX1 = condsigma.e1_1x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub1(n,1);
	alph_sub1.col(0) = alph.col(0); // alph1_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f1 = alph_sub1.t();

	/* Centering the alpha */
	mat alph_c1 = alph_f1.col(i)-c_mu1;

	/* Conditional F1c1_2 */
	/* Conditional mean and variance-covariance matrix, conditional on alph1_1, u1 and u2 */
	mat c_sig2 = condsigma.e1c1_2s;
	mat ic_sig2 = condsigma.e1c1_2i; // the inverse
	mat sigX2 = condsigma.e1c1_2x; // the sigx
	vec a2 = mu.elem(rc1_1u);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub2(n,1);
	alph_sub2.col(0) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f2 = alph_sub2.t();

	/* Centering the alpha */
	mat alph_c2 = alph_f2.col(i)-c_mu2;

	/* Estimation of F1c1_2 */
	double p2 = pn(alph_c2,0.0,c_sig2[0]);
	double F1c1_2 = pi1_2(i)*p2;

	/* Calculating the pdf */
	double sd2 = sqrt(c_sig2[0]);
	double inner2 = pow(alph_c2[0],2)/c_sig2[0];
	double pdf1c1_2 = 1/sq_twopi*1/sd2*exp(-0.5*inner2);

	/* Score contribution from F1c1_1 */
	double dF1c1_2_u1 = dpi1_2_u1(i)*p2+pi1_2(i)*pdf1c1_2*(-sigX2(1));
	double dF1c1_2_u2 = dpi1_2_u2(i)*p2+pi1_2(i)*pdf1c1_2*(-sigX2(2));

	/* Conditional F2c1_2 */
	/* Conditional mean and variance-covariance matrix, conditional on alph1_1, u1 and u2 */
	mat c_sig3 = condsigma.e2c1_2s;
	mat ic_sig3 = condsigma.e2c1_2i; // the inverse
	mat sigX3 = condsigma.e2c1_2x; // the sigx
	vec a3 = mu.elem(rc1_1u);
	vec c_mu3 = sigX3*a3;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub3(n,1);
	alph_sub3.col(0) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f3 = alph_sub3.t();

	/* Centering the alpha */
	mat alph_c3 = alph_f3.col(i)-c_mu3;

	/* Estimation of F2c1_2 */
	double p3 = pn(alph_c3,0.0,c_sig3[0]);
	double F2c1_2 = pi2_2(i)*p3;

	/* Calculating the pdf */
	double sd3 = sqrt(c_sig3[0]);
	double inner3 = pow(alph_c3[0],2)/c_sig3[0];
	double pdf2c1_2 = 1/sq_twopi*1/sd3*exp(-0.5*inner3);

	/* Score contribution from F2c1_1 */
	double dF2c1_2_u1 = dpi2_2_u1(i)*p3+pi2_2(i)*pdf2c1_2*(-sigX3(1));
	double dF2c1_2_u2 = dpi2_2_u2(i)*p3+pi2_2(i)*pdf2c1_2*(-sigX3(2));

	/* Score contributions from dF01 wrt. u1 and u2 */
	double sc_u1 = (dpi1_1_u1(i)+pi1_1(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(0))/pi1_1(i)-(dF1c1_2_u1+dF2c1_2_u1)/(1-F1c1_2-F2c1_2);
	double sc_u2 = (dpi1_1_u2(i)+pi1_1(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(1))/pi1_1(i)-(dF1c1_2_u2+dF2c1_2_u2)/(1-F1c1_2-F2c1_2);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Family member 1 experience event 0, family member 2 experience event 2, estimating score contribution from dF02 */
      else if((y(i,0) == 0) & (y(i,1) == 2)){

	/* Marginal dF2_2 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e2_2s;
	mat ic_sig1 = condsigma.e2_2i; // the inverse
	mat sigX1 = condsigma.e2_2x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub1(n,1);
	alph_sub1.col(0) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f1 = alph_sub1.t();

	/* Centering the alpha */
	mat alph_c1 = alph_f1.col(i)-c_mu1;

	/* Conditional F1c2_1 */
	/* Conditional mean and variance-covariance matrix, conditional on alph2_2, u1 and u2 */
	mat c_sig2 = condsigma.e1c2_1s;
	mat ic_sig2 = condsigma.e1c2_1i; // the inverse
	mat sigX2 = condsigma.e1c2_1x; // the sigx
	vec a2 = mu.elem(rc2_2u);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub2(n,1);
	alph_sub2.col(0) = alph.col(0); // alph1_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f2 = alph_sub2.t();

	/* Centering the alpha */
	mat alph_c2 = alph_f2.col(i)-c_mu2;

	/* Estimation of F1c2_1 */
	double p2 = pn(alph_c2,0.0,c_sig2[0]);
	double F1c2_1 = pi1_1(i)*p2;

	/* Calculating the pdf */
	double sd2 = sqrt(c_sig2[0]);
	double inner2 = pow(alph_c2[0],2)/c_sig2[0];
	double pdf1c2_1 = 1/sq_twopi*1/sd2*exp(-0.5*inner2);

	/* Score contribution from F1c2_1 */
	double dF1c2_1_u1 = dpi1_1_u1(i)*p2+pi1_1(i)*pdf1c2_1*(-sigX2(1));
	double dF1c2_1_u2 = dpi1_1_u2(i)*p2+pi1_1(i)*pdf1c2_1*(-sigX2(2));

	/* Conditional F2c2_1 */
	/* Conditional mean and variance-covariance matrix, conditional on alph2_2, u1 and u2 */
	mat c_sig3 = condsigma.e2c2_1s;
	mat ic_sig3 = condsigma.e2c2_1i; // the inverse
	mat sigX3 = condsigma.e2c2_1x; // the sigx
	vec a3 = mu.elem(rc2_2u);
	vec c_mu3 = sigX3*a3;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub3(n,1);
	alph_sub3.col(0) = alph.col(2); // alph2_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f3 = alph_sub3.t();

	/* Centering the alpha */
	mat alph_c3 = alph_f3.col(i)-c_mu3;

	/* Estimation of F2c2_1 */
	double p3 = pn(alph_c3,0.0,c_sig3[0]);
	double F2c2_1 = pi2_1(i)*p3;

	/* Calculating the pdf */
	double sd3 = sqrt(c_sig3[0]);
	double inner3 = pow(alph_c3[0],2)/c_sig3[0];
	double pdf2c2_1 = 1/sq_twopi*1/sd3*exp(-0.5*inner3);

	/* Score contribution from F2c2_1 */
	double dF2c2_1_u1 = dpi2_1_u1(i)*p3+pi2_1(i)*pdf2c2_1*(-sigX3(1));
	double dF2c2_1_u2 = dpi2_1_u2(i)*p3+pi2_1(i)*pdf2c2_1*(-sigX3(2));

	/* Score contributions from dF01 wrt. u1 and u2 */
	double sc_u1 = (dpi2_2_u1(i)+pi2_2(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(0))/pi2_2(i)-(dF1c2_1_u1+dF2c2_1_u1)/(1-F1c2_1-F2c2_1);
	double sc_u2 = (dpi2_2_u2(i)+pi2_2(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(1))/pi2_2(i)-(dF1c2_1_u2+dF2c2_1_u2)/(1-F1c2_1-F2c2_1);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Family member 1 experience event 2, family member 2 experience event 0, estimating score contribution from dF20 */
      else if((y(i,0) == 2) & (y(i,1) == 0)){

	/* Marginal dF2_1 */
	/* Conditional mean and variance-covariance matrix, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e2_1s;
	mat ic_sig1 = condsigma.e2_1i; // the inverse
	mat sigX1 = condsigma.e2_1x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub1(n,1);
	alph_sub1.col(0) = alph.col(2); // alph2_1

	/* Transposing matrix to be of the form 1xn */
	mat alph_f1 = alph_sub1.t();

	/* Centering the alpha */
	mat alph_c1 = alph_f1.col(i)-c_mu1;

	/* Conditional F1c2_2 */
	/* Conditional mean and variance-covariance matrix, conditional on alph2_1, u1 and u2 */
	mat c_sig2 = condsigma.e1c2_2s;
	mat ic_sig2 = condsigma.e1c2_2i; // the inverse
	mat sigX2 = condsigma.e1c2_2x; // the sigx
	vec a2 = mu.elem(rc2_1u);
	vec c_mu2 = sigX2*a2;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub2(n,1);
	alph_sub2.col(0) = alph.col(1); // alph1_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f2 = alph_sub2.t();

	/* Centering the alpha */
	mat alph_c2 = alph_f2.col(i)-c_mu2;

	/* Estimation of F1c2_2 */
	double sd2 = sqrt(c_sig2[0]);
	double p2 = pn(alph_c2,0.0,c_sig2[0]);
	double F1c2_2 = pi1_2(i)*p2;

	/* Calculating the pdf */
	double inner2 = pow(alph_c2[0],2)/c_sig2[0];
	double pdf1c2_2 = 1/sq_twopi*1/sd2*exp(-0.5*inner2);

	/* Score contribution from F1c2_2 */
	double dF1c2_2_u1 = dpi1_2_u1(i)*p2+pi1_2(i)*pdf1c2_2*(-sigX2(1));
	double dF1c2_2_u2 = dpi1_2_u2(i)*p2+pi1_2(i)*pdf1c2_2*(-sigX2(2));

	/* Conditional F2c2_2 */
	/* Conditional mean and variance-covariance matrix, conditional on alph2_1, u1 and u2 */
	mat c_sig3 = condsigma.e2c2_2s;
	mat ic_sig3 = condsigma.e2c2_2i; // the inverse
	mat sigX3 = condsigma.e2c2_2x; // the sigx
	vec a3 = mu.elem(rc2_1u);
	vec c_mu3 = sigX3*a3;

	/* Pulling out the appropriate alpha from alph */
	mat alph_sub3(n,1);
	alph_sub3.col(0) = alph.col(3); // alph2_2

	/* Transposing matrix to be of the form 1xn */
	mat alph_f3 = alph_sub3.t();

	/* Centering the alpha */
	mat alph_c3 = alph_f3.col(i)-c_mu3;

	/* Estimation of F2c2_2 */
	double sd3 = sqrt(c_sig3[0]);
	double p3 = pn(alph_c3,0.0,c_sig3[0]);
	double F2c2_2 = pi2_2(i)*p3;

	/* Calculating the pdf */
	double inner3 = pow(alph_c3[0],2)/c_sig3[0];
	double pdf2c2_2 = 1/sq_twopi*1/sd3*exp(-0.5*inner3);

	/* Score contribution from F2c2_2 */
	double dF2c2_2_u1 = dpi2_2_u1(i)*p3+pi2_2(i)*pdf2c2_2*(-sigX3(1));
	double dF2c2_2_u2 = dpi2_2_u2(i)*p3+pi2_2(i)*pdf2c2_2*(-sigX3(2));

	/* Score contributions from dF01 wrt. u1 and u2 */
	double sc_u1 = (dpi2_1_u1(i)+pi2_1(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(0))/pi2_1(i)-(dF1c2_2_u1+dF2c2_2_u1)/(1-F1c2_2-F2c2_2);
	double sc_u2 = (dpi2_1_u2(i)+pi2_1(i)*1/(2*c_sig1[0])*2*alph_c1[0]*sigX1(1))/pi2_1(i)-(dF1c2_2_u2+dF2c2_2_u2)/(1-F1c2_2-F2c2_2);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
      /* Both family members experience event 0, estimating score contribution from F00 */
      else{

	/* Marginal F1_1, F1_2, F2_1 and F2_2 */
	/* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
	mat c_sig1 = condsigma.e1_1s;
	mat ic_sig1 = condsigma.e1_1i; // the inverse
	mat sigX1 = condsigma.e1_1x; // the sigx
	vec a1 = mu.elem(rcu);
	vec c_mu1 = sigX1*a1;
	double sd1 = sqrt(c_sig1[0]);

	mat c_sig2 = condsigma.e1_2s;
	mat ic_sig2 = condsigma.e1_2i; // the inverse
	mat sigX2 = condsigma.e1_2x; // the sigx
	vec a2 = mu.elem(rcu);
	vec c_mu2 = sigX2*a2;
	double sd2 = sqrt(c_sig2[0]);

	mat c_sig3 = condsigma.e2_1s;
	mat ic_sig3 = condsigma.e2_1i; // the inverse
	mat sigX3 = condsigma.e2_1x; // the sigx
	vec a3 = mu.elem(rcu);
	vec c_mu3 = sigX3*a3;
	double sd3 = sqrt(c_sig3[0]);

	mat c_sig4 = condsigma.e2_2s;
	mat ic_sig4 = condsigma.e2_2i; // the inverse
	mat sigX4 = condsigma.e2_2x; // the sigx
	vec a4 = mu.elem(rcu);
	vec c_mu4 = sigX4*a4;
	double sd4 = sqrt(c_sig4[0]);

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
	mat alph_c1_2 = alph_f1_2.col(i)-c_mu2;
	mat alph_c2_1 = alph_f2_1.col(i)-c_mu3;
	mat alph_c2_2 = alph_f2_2.col(i)-c_mu4;

	/* Estimation of F1_1, F1_2, F2_1 and F2_2 */
	double p1 = pn(alph_c1_1,0.0,sd1*sd1);
	double p2 = pn(alph_c1_2,0.0,sd2*sd2);
	double p3 = pn(alph_c2_1,0.0,sd3*sd3);
	double p4 = pn(alph_c2_2,0.0,sd4*sd4);

	double F1_1 = pi1_1(i)*p1;
	double F1_2 = pi1_2(i)*p2;
	double F2_1 = pi2_1(i)*p3;
	double F2_2 = pi2_2(i)*p4;

	/* Calculating the pdf of F1_1, F1_2, F2_1 and F2_2 */
	double inner1 = pow(alph_c1_1[0],2)/c_sig1[0];
	double inner2 = pow(alph_c1_2[0],2)/c_sig2[0];
	double inner3 = pow(alph_c2_1[0],2)/c_sig3[0];
	double inner4 = pow(alph_c2_2[0],2)/c_sig4[0];

	double pdf1_1 = 1/sq_twopi*1/sd1*exp(-0.5*inner1);
	double pdf1_2 = 1/sq_twopi*1/sd2*exp(-0.5*inner2);
	double pdf2_1 = 1/sq_twopi*1/sd3*exp(-0.5*inner3);
	double pdf2_2 = 1/sq_twopi*1/sd4*exp(-0.5*inner4);

	/* Score contributions from F1_1, F1_2, F2_1 and F2_2 */
	double dF1_1_u1 = dpi1_1_u1(i)*p1+pi1_1(i)*pdf1_1*(-sigX1(0));
	double dF1_2_u1 = dpi1_2_u1(i)*p2+pi1_2(i)*pdf1_2*(-sigX2(0));
	double dF2_1_u1 = dpi2_1_u1(i)*p3+pi2_1(i)*pdf2_1*(-sigX3(0));
	double dF2_2_u1 = dpi2_2_u1(i)*p4+pi2_2(i)*pdf2_2*(-sigX4(0));

	double dF1_1_u2 = dpi1_1_u2(i)*p1+pi1_1(i)*pdf1_1*(-sigX1(1));
	double dF1_2_u2 = dpi1_2_u2(i)*p2+pi1_2(i)*pdf1_2*(-sigX2(1));
	double dF2_1_u2 = dpi2_1_u2(i)*p3+pi2_1(i)*pdf2_1*(-sigX3(1));
	double dF2_2_u2 = dpi2_2_u2(i)*p4+pi2_2(i)*pdf2_2*(-sigX4(1));

	/* Joint probabilities F11, F12, F21 and F22 */
	/* Conditional mean and variance-covariance matrices, conditional on u1 and u2 */
	mat c_sig5 = condsigma.e11s;
	mat ic_sig5 = condsigma.e11i; // the inverse
	mat sigX5 = condsigma.e11x; // the sigx
	vec a5 = mu.elem(rcu);
	vec c_mu5 = sigX5*a5;

	mat c_sig6 = condsigma.e12s;
	mat ic_sig6 = condsigma.e12i; // the inverse
	mat sigX6 = condsigma.e12x; // the sigx
	vec a6 = mu.elem(rcu);
	vec c_mu6 = sigX6*a6;

	mat c_sig7 = condsigma.e21s;
	mat ic_sig7 = condsigma.e21i; // the inverse
	mat sigX7 = condsigma.e21x; // the sigx
	vec a7 = mu.elem(rcu);
	vec c_mu7 = sigX7*a7;

	mat c_sig8 = condsigma.e22s;
	mat ic_sig8 = condsigma.e22i; // the inverse
	mat sigX8 = condsigma.e22x; // the sigx
	vec a8 = mu.elem(rcu);
	vec c_mu8 = sigX8*a8;

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

	/* Estimating F11, F12, F21 and F22 */
	double cdf11 = pn(alph_f11.col(i),c_mu5,c_sig5);
	double cdf12 = pn(alph_f12.col(i),c_mu6,c_sig6);
	double cdf21 = pn(alph_f21.col(i),c_mu7,c_sig7);
	double cdf22 = pn(alph_f22.col(i),c_mu8,c_sig8);

	double F11 = pi1_1(i)*pi1_2(i)*cdf11;
	double F12 = pi1_1(i)*pi2_2(i)*cdf12;
	double F21 = pi2_1(i)*pi1_2(i)*cdf21;
	double F22 = pi2_1(i)*pi2_2(i)*cdf22;

	/* Conditional mean and variance-covariance matrix for F1c1_1 and F1c1_2 */
	mat c_sig9 = condsigma.e1c1_1s;
	mat ic_sig9 = condsigma.e1c1_1i; // the inverse
	mat sigX9 = condsigma.e1c1_1x; // the sigx
	vec a9 = mu.elem(rc2_1u);
	vec c_mu9 = sigX9*a9;
	double sd9 = sqrt(c_sig9[0]);

	mat c_sig10 = condsigma.e1c1_2s;
	mat ic_sig10 = condsigma.e1c1_2i; // the inverse
	mat sigX10 = condsigma.e1c1_2x; // the sigx
	vec a10 = mu.elem(rc1_1u);
	vec c_mu10 = sigX10*a10;
	double sd10 = sqrt(c_sig10[0]);

	/* Centering the alphas */
	mat alph_c9 = alph_f1_1.col(i)-c_mu9;
	mat alph_c10 = alph_f1_2.col(i)-c_mu10;

	/* Conditional CDF */
	double cdf1c1_1 = pn(alph_c9,0.0,sd9*sd9);
	double cdf1c1_2 = pn(alph_c10,0.0,sd10*sd10);

	/* Score contributions */
	double dcdf11_u1 = cdf1c1_1*pdf1_2*(-sigX5(1,0)) + cdf1c1_2*pdf1_1*(-sigX5(0,0));
	double dcdf11_u2 = cdf1c1_1*pdf1_2*(-sigX5(1,1)) + cdf1c1_2*pdf1_1*(-sigX5(0,1));

	/* Conditional mean and variance-covariance matrix for F1c2_1 and F2c1_2 */
	mat c_sig11 = condsigma.e1c2_1s;
	mat ic_sig11 = condsigma.e1c2_1i; // the inverse
	mat sigX11 = condsigma.e1c2_1x; // the sigx
	vec a11 = mu.elem(rc2_2u);
	vec c_mu11 = sigX11*a11;
	double sd11 = sqrt(c_sig11[0]);

	mat c_sig12 = condsigma.e2c1_2s;
	mat ic_sig12 = condsigma.e2c1_2i; // the inverse
	mat sigX12 = condsigma.e2c1_2x; // the sigx
	vec a12 = mu.elem(rc1_1u);
	vec c_mu12 = sigX12*a12;
	double sd12 = sqrt(c_sig12[0]);

	/* Centering the alphas */
	mat alph_c11_n = alph_f1_1.col(i)-c_mu11;
	mat alph_c12_n = alph_f2_2.col(i)-c_mu12;

	/* Conditional CDF */
	double cdf1c2_1 = pn(alph_c11_n,0.0,sd11*sd11);
	double cdf2c1_2 = pn(alph_c12_n,0.0,sd12*sd12);

	/* Score contributions */
	double dcdf12_u1 = cdf1c2_1*pdf2_2*(-sigX6(1,0)) + cdf2c1_2*pdf1_1*(-sigX6(0,0));
	double dcdf12_u2 = cdf1c2_1*pdf2_2*(-sigX6(1,1)) + cdf2c1_2*pdf1_1*(-sigX6(0,1));

	/* Conditional mean and variance-covariance matrix for F2c1_1 and F1c2_2 */
	mat c_sig13 = condsigma.e2c1_1s;
	mat ic_sig13 = condsigma.e2c1_1i; // the inverse
	mat sigX13 = condsigma.e2c1_1x; // the sigx
	vec a13 = mu.elem(rc1_2u);
	vec c_mu13 = sigX13*a13;
	double sd13 = sqrt(c_sig13[0]);

	mat c_sig14 = condsigma.e1c2_2s;
	mat ic_sig14 = condsigma.e1c2_2i; // the inverse
	mat sigX14 = condsigma.e1c2_2x; // the sigx
	vec a14 = mu.elem(rc2_1u);
	vec c_mu14 = sigX14*a14;
	double sd14 = sqrt(c_sig14[0]);

	/* Centering the alphas */
	mat alph_c13 = alph_f2_1.col(i)-c_mu13;
	mat alph_c14 = alph_f1_2.col(i)-c_mu14;

	/* Conditional CDF */
	double cdf2c1_1 = pn(alph_c13,0.0,sd13*sd13);
	double cdf1c2_2 = pn(alph_c14,0.0,sd14*sd14);

	/* Score contributions */
	double dcdf21_u1 = cdf2c1_1*pdf1_2*(-sigX7(1,0)) + cdf1c2_2*pdf2_1*(-sigX7(0,0));
	double dcdf21_u2 = cdf2c1_1*pdf1_2*(-sigX7(1,1)) + cdf1c2_2*pdf2_1*(-sigX7(0,1));

	/* Conditional mean and variance-covariance matrix for F2c2_1 and F2c2_2 */
	mat c_sig15 = condsigma.e2c2_1s;
	mat ic_sig15 = condsigma.e2c2_1i; // the inverse
	mat sigX15 = condsigma.e2c2_1x; // the sigx
	vec a15 = mu.elem(rc2_2u);
	vec c_mu15 = sigX15*a15;
	double sd15 = sqrt(c_sig15[0]);

	mat c_sig16 = condsigma.e2c2_2s;
	mat ic_sig16 = condsigma.e2c2_2i; // the inverse
	mat sigX16 = condsigma.e2c2_2x; // the sigx
	vec a16 = mu.elem(rc2_1u);
	vec c_mu16 = sigX16*a16;
	double sd16 = sqrt(c_sig16[0]);

	/* Centering the alphas */
	mat alph_c15 = alph_f2_1.col(i)-c_mu15;
	mat alph_c16 = alph_f2_2.col(i)-c_mu16;

	/* Conditional CDF */
	double cdf2c2_1 = pn(alph_c15,0.0,sd15*sd15);
	double cdf2c2_2 = pn(alph_c16,0.0,sd16*sd16);

	/* Score contributions */
	double dcdf22_u1 = cdf2c2_1*pdf2_2*(-sigX8(1,0)) + cdf2c2_2*pdf2_1*(-sigX8(0,0));
	double dcdf22_u2 = cdf2c2_1*pdf2_2*(-sigX8(1,1)) + cdf2c2_2*pdf2_1*(-sigX8(0,1));

	/* Score contributions from F11, F12, F21 and F22 */
	double dF11_u1 = dpi11_u1(i)*cdf11+pi1_1(i)*pi1_2(i)*dcdf11_u1;
	double dF12_u1 = dpi12_u1(i)*cdf12+pi1_1(i)*pi2_2(i)*dcdf12_u1;
	double dF21_u1 = dpi21_u1(i)*cdf21+pi2_1(i)*pi1_2(i)*dcdf21_u1;
	double dF22_u1 = dpi22_u1(i)*cdf22+pi2_1(i)*pi2_2(i)*dcdf22_u1;

	double dF11_u2 = dpi11_u2(i)*cdf11+pi1_1(i)*pi1_2(i)*dcdf11_u2;
	double dF12_u2 = dpi12_u2(i)*cdf12+pi1_1(i)*pi2_2(i)*dcdf12_u2;
	double dF21_u2 = dpi21_u2(i)*cdf21+pi2_1(i)*pi1_2(i)*dcdf21_u2;
	double dF22_u2 = dpi22_u2(i)*cdf22+pi2_1(i)*pi2_2(i)*dcdf22_u2;

	/* Likelihood contribution */
	double F00 = 1-F1_1-F1_2-F2_1-F2_2+F11+F12+F21+F22;

	/* Score contribution from F00 */
	double sc_u1 = (1/F00)*(-dF1_1_u1-dF1_2_u1-dF2_1_u1-dF2_2_u1+dF11_u1+dF12_u1+dF21_u1+dF22_u1);
	double sc_u2 = (1/F00)*(-dF1_1_u2-dF1_2_u2-dF2_1_u2-dF2_2_u2+dF11_u2+dF12_u2+dF21_u2+dF22_u2);

	/* Adding to return vector */
	res(i,0) = sc_u1;
	res(i,1) = sc_u2;
      }
    }
    if (full) {
      /* u */
      mat sigu = condsigma.us;
      mat isigu = condsigma.ui;
      mat pu = u.row(i);

      /* Derivative of the pdf of the u's wrt u1 and u2 */
      double denom = sigu(0,0)*sigu(1,1)-sigu(0,1)*sigu(1,0);
      double dpdfu_u1 = (-0.5)*(2*sigu(1,1)*pu(0)-sigu(0,1)*pu(1)-sigu(1,0)*pu(1))/denom;
      double dpdfu_u2 = (-0.5)*(2*sigu(0,0)*pu(1)-sigu(0,1)*pu(0)-sigu(1,0)*pu(0))/denom;

      /* Adding to score */
      res(i,0) += dpdfu_u1;
      res(i,1) += dpdfu_u2;
    }
  }
  return(res);
}

/*
Hessian matrix of full loglikelihood
*/
/*// [[Rcpp::export]]*/
mat D2loglikfull(mat y, mat b, mat u, ss condsigma, mat alph, mat dalph, mat tau, bool full=1){
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

  /* Central difference */
  res.row(0) = (Dloglikfull(y, b, u1_p, condsigma, alph, dalph, tau, full)-Dloglikfull(y, b, u1_m, condsigma, alph, dalph, tau, full))/(2*h);
  res.row(1) = (Dloglikfull(y, b, u2_p, condsigma, alph, dalph, tau, full)-Dloglikfull(y, b, u2_m, condsigma, alph, dalph, tau, full))/(2*h);

  /* Return */
  return(res);
}

/*
Marginal likelihood via AGQ
*/
// [[Rcpp::export]]
vec loglik(mat y, mat b, mat sigma, mat alph, mat dalph, mat tau, mat eb0, int nq=1, double stepsize=0.7, unsigned iter=20, bool debug=false) {
  QuadRule gh(nq);
  double K = sqrt(2);
  vec z = gh.Abscissa();
  vec w = gh.Weight();
  int n = y.n_rows;
  vec warn(n); warn.fill(1);
  vec res(n); res.fill(0);

  /* Specifying components of sigma */
  uvec rc1(1); rc1(0) = 0;
  uvec rc2(1); rc2(0) = 1;
  uvec rc3(1); rc3(0) = 2;
  uvec rc4(1); rc4(0) = 3;

  uvec rc5(2); rc5(0) = 4; rc5(1) = 5;

  uvec rc6(2); rc6(0) = 0; rc6(1) = 1;
  uvec rc7(2); rc7(0) = 0; rc7(1) = 3;
  uvec rc8(2); rc8(0) = 2; rc8(1) = 1;
  uvec rc9(2); rc9(0) = 2; rc9(1) = 3;

  uvec rc10(3); rc10(0) = 0; rc10(1) = 4; rc10(2) = 5;
  uvec rc11(3); rc11(0) = 1; rc11(1) = 4; rc11(2) = 5;
  uvec rc12(3); rc12(0) = 2; rc12(1) = 4; rc12(2) = 5;
  uvec rc13(3); rc13(0) = 3; rc13(1) = 4; rc13(2) = 5;

  /* Estimation conditional sigmas etc */
  matfour e1_1 = condsig(sigma,rc1,rc5);
  matfour e1_2 = condsig(sigma,rc2,rc5);
  matfour e2_1 = condsig(sigma,rc3,rc5);
  matfour e2_2 = condsig(sigma,rc4,rc5);

  matfour e1c1_1 = condsig(sigma,rc1,rc11);
  matfour e2c1_1 = condsig(sigma,rc3,rc11);

  matfour e1c1_2 = condsig(sigma,rc2,rc10);
  matfour e2c1_2 = condsig(sigma,rc4,rc10);

  matfour e1c2_1 = condsig(sigma,rc1,rc13);
  matfour e2c2_1 = condsig(sigma,rc3,rc13);

  matfour e1c2_2 = condsig(sigma,rc2,rc12);
  matfour e2c2_2 = condsig(sigma,rc4,rc12);

  matfour e11 = condsig(sigma,rc6,rc5);
  matfour e12 = condsig(sigma,rc7,rc5);
  matfour e21 = condsig(sigma,rc8,rc5);
  matfour e22 = condsig(sigma,rc9,rc5);

  mat sigu = sigma.submat(rc5,rc5);
  mat isigu = sigu.i();
  double dsigu = det(sigu);
  double sq_dsigu = sqrt(dsigu);

  ss condsigma;
  condsigma.e1_1s = e1_1.M1;
  condsigma.e2_1s = e2_1.M1;
  condsigma.e1_2s = e1_2.M1;
  condsigma.e2_2s = e2_2.M1;

  condsigma.e1_1i = e1_1.M2;
  condsigma.e2_1i = e2_1.M2;
  condsigma.e1_2i = e1_2.M2;
  condsigma.e2_2i = e2_2.M2;

  condsigma.e1_1d = e1_1.M3;
  condsigma.e2_1d = e2_1.M3;
  condsigma.e1_2d = e1_2.M3;
  condsigma.e2_2d = e2_2.M3;

  condsigma.e1_1x = e1_1.M4;
  condsigma.e2_1x = e2_1.M4;
  condsigma.e1_2x = e1_2.M4;
  condsigma.e2_2x = e2_2.M4;

  condsigma.e1c1_1s = e1c1_1.M1;
  condsigma.e2c1_1s = e2c1_1.M1;
  condsigma.e1c1_2s = e1c1_2.M1;
  condsigma.e2c1_2s = e2c1_2.M1;

  condsigma.e1c1_1i = e1c1_1.M2;
  condsigma.e2c1_1i = e2c1_1.M2;
  condsigma.e1c1_2i = e1c1_2.M2;
  condsigma.e2c1_2i = e2c1_2.M2;

  condsigma.e1c1_1d = e1c1_1.M3;
  condsigma.e2c1_1d = e2c1_1.M3;
  condsigma.e1c1_2d = e1c1_2.M3;
  condsigma.e2c1_2d = e2c1_2.M3;

  condsigma.e1c1_1x = e1c1_1.M4;
  condsigma.e2c1_1x = e2c1_1.M4;
  condsigma.e1c1_2x = e1c1_2.M4;
  condsigma.e2c1_2x = e2c1_2.M4;

  condsigma.e1c2_1s = e1c2_1.M1;
  condsigma.e2c2_1s = e2c2_1.M1;
  condsigma.e1c2_2s = e1c2_2.M1;
  condsigma.e2c2_2s = e2c2_2.M1;

  condsigma.e1c2_1i = e1c2_1.M2;
  condsigma.e2c2_1i = e2c2_1.M2;
  condsigma.e1c2_2i = e1c2_2.M2;
  condsigma.e2c2_2i = e2c2_2.M2;

  condsigma.e1c2_1d = e1c2_1.M3;
  condsigma.e2c2_1d = e2c2_1.M3;
  condsigma.e1c2_2d = e1c2_2.M3;
  condsigma.e2c2_2d = e2c2_2.M3;

  condsigma.e1c2_1x = e1c2_1.M4;
  condsigma.e2c2_1x = e2c2_1.M4;
  condsigma.e1c2_2x = e1c2_2.M4;
  condsigma.e2c2_2x = e2c2_2.M4;

  condsigma.e11s = e11.M1;
  condsigma.e12s = e12.M1;
  condsigma.e21s = e21.M1;
  condsigma.e22s = e22.M1;

  condsigma.e11i = e11.M2;
  condsigma.e12i = e12.M2;
  condsigma.e21i = e21.M2;
  condsigma.e22i = e22.M2;

  condsigma.e11d = e11.M3;
  condsigma.e12d = e12.M3;
  condsigma.e21d = e21.M3;
  condsigma.e22d = e22.M3;

  condsigma.e11x = e11.M4;
  condsigma.e12x = e12.M4;
  condsigma.e21x = e21.M4;
  condsigma.e22x = e22.M4;

  condsigma.us = sigu;
  condsigma.ui = isigu;
  condsigma.squd = sq_dsigu;

  for (int i=0; i<n; i++) {
    mat y0 = y.row(i);
    mat b0 = b.row(i);
    mat alph0 = alph.row(i);
    mat dalph0 = dalph.row(i);
    mat tau0 = tau.row(i);
    mat u0(1,2); u0 = eb0.row(i);
    double conv = 1;
    mat H(2,2);
    mat U(1,2);
    /* Newton Raphson */
    unsigned j;
    for (j=0; j<iter; j++) {
      U = Dloglikfull(y0,b0,u0,condsigma,alph0,dalph0,tau0);
      H = D2loglikfull(y0,b0,u0,condsigma,alph0,dalph0,tau0);
      conv = (U(0)*U(0)+U(1)*U(1))/2;
      if (conv<_inner_NR_abseps) {
	warn(i) = 0;
	break;
      }
      u0 = u0-stepsize*U*H.i();
    }
    if (debug) {
      U = Dloglikfull(y0,b0,u0,condsigma,alph0,dalph0,tau0);
      vec L = loglikfull(y0,b0,u0,condsigma,alph0,dalph0,tau0);
      Rcpp::Rcout << "iter: " << j <<std::endl;
      Rcpp::Rcout << "conv: " << conv <<std::endl;
      Rcpp::Rcout << "L: " << L <<std::endl;
      Rcpp::Rcout << "U: " << U <<std::endl;
      Rcpp::Rcout << "i: " << i <<std::endl;
    }
    /* Laplace approximation */
    if (nq==0) {
      vec logf = loglikfull(y0,b0,u0,condsigma,alph0,dalph0,tau0);
      double lapl = log(twopi)-0.5*log(det(H))+logf(0);
      res(i) = lapl;
    } else {
      /* Adaptive Gaussian quadrature */
      bool useSVD = true;
      mat G = -H;
      double logdetG = 0;
      mat Bi;
      if (useSVD){
	mat W; vec s; mat V;
	svd(W,s,V,G);
	vec is=s;
	for (unsigned m=0; m<s.n_elem; m++){
	  if (s[m]<1e-9){
	    s[m] = 1e-9;
	    is[m] = 0;
	  } else {
	    is[m] = 1/sqrt(s[m]);
	  }
	}
	Bi = trans(W*diagmat(is)*V);
	logdetG = sum(log(s));
      } else {
	mat B = chol(G);
	Bi = B.i();
	logdetG = log(det(G));
      }
      double Sum = 0;
      for (unsigned k=0; k<z.n_elem; k++) {
      	for (unsigned l=0; l<z.n_elem; l++) {
      	  mat z0(2,1);
      	  z0(0) = z[k]; z0(1) = z[l];
      	  mat a0 = u0.t()+K*Bi*z0;
	  double w0 = w[k]*w[l]*exp(z0[0]*z0[0]+z0[1]*z0[1]);
	  double ll0 = loglikfull(y0,b0,a0.t(),condsigma,alph0,dalph0,tau0)[0];
	  Sum += exp(ll0)*w0;
      	}
      }
      res(i) = 2*log(K)-0.5*logdetG+log(Sum);
    }
  }
  return(res);
}

/*
EB - Empirical Bayes, intelligent starting values of u0 (for estimation of score and Hessian)
*/
// [[Rcpp::export]]
mat EB(mat y, mat b, mat sigma, mat alph, mat dalph, mat tau, double stepsize=0.7, unsigned iter=20, bool debug=false) {
  int n = y.n_rows;
  vec warn(n); warn.fill(1);
  mat eb0(n,2);
  for (int i=0; i<n; i++) {
    mat y0 = y.row(i);
    mat b0 = b.row(i);
    mat alph0 = alph.row(i);
    mat dalph0 = dalph.row(i);
    mat tau0 = tau.row(i);
    mat u0(1,2); u0.fill(0);
    double conv = 1;
    mat H(2,2);
    mat U(1,2);

    /* Specifying components of sigma */
    uvec rc1(1); rc1(0) = 0;
    uvec rc2(1); rc2(0) = 1;
    uvec rc3(1); rc3(0) = 2;
    uvec rc4(1); rc4(0) = 3;

    uvec rc5(2); rc5(0) = 4; rc5(1) = 5;

    uvec rc6(2); rc6(0) = 0; rc6(1) = 1;
    uvec rc7(2); rc7(0) = 0; rc7(1) = 3;
    uvec rc8(2); rc8(0) = 2; rc8(1) = 1;
    uvec rc9(2); rc9(0) = 2; rc9(1) = 3;

    uvec rc10(3); rc10(0) = 0; rc10(1) = 4; rc10(2) = 5;
    uvec rc11(3); rc11(0) = 1; rc11(1) = 4; rc11(2) = 5;
    uvec rc12(3); rc12(0) = 2; rc12(1) = 4; rc12(2) = 5;
    uvec rc13(3); rc13(0) = 3; rc13(1) = 4; rc13(2) = 5;

    /* Estimation conditional sigmas etc */
    matfour e1_1 = condsig(sigma,rc1,rc5);
    matfour e1_2 = condsig(sigma,rc2,rc5);
    matfour e2_1 = condsig(sigma,rc3,rc5);
    matfour e2_2 = condsig(sigma,rc4,rc5);

    matfour e1c1_1 = condsig(sigma,rc1,rc11);
    matfour e2c1_1 = condsig(sigma,rc3,rc11);

    matfour e1c1_2 = condsig(sigma,rc2,rc10);
    matfour e2c1_2 = condsig(sigma,rc4,rc10);

    matfour e1c2_1 = condsig(sigma,rc1,rc13);
    matfour e2c2_1 = condsig(sigma,rc3,rc13);

    matfour e1c2_2 = condsig(sigma,rc2,rc12);
    matfour e2c2_2 = condsig(sigma,rc4,rc12);

    matfour e11 = condsig(sigma,rc6,rc5);
    matfour e12 = condsig(sigma,rc7,rc5);
    matfour e21 = condsig(sigma,rc8,rc5);
    matfour e22 = condsig(sigma,rc9,rc5);

    mat sigu = sigma.submat(rc5,rc5);
    mat isigu = sigu.i();
    double dsigu = det(sigu);
    double sq_dsigu = sqrt(dsigu);

    ss condsigma;
    condsigma.e1_1s = e1_1.M1;
    condsigma.e2_1s = e2_1.M1;
    condsigma.e1_2s = e1_2.M1;
    condsigma.e2_2s = e2_2.M1;

    condsigma.e1_1i = e1_1.M2;
    condsigma.e2_1i = e2_1.M2;
    condsigma.e1_2i = e1_2.M2;
    condsigma.e2_2i = e2_2.M2;

    condsigma.e1_1d = e1_1.M3;
    condsigma.e2_1d = e2_1.M3;
    condsigma.e1_2d = e1_2.M3;
    condsigma.e2_2d = e2_2.M3;

    condsigma.e1_1x = e1_1.M4;
    condsigma.e2_1x = e2_1.M4;
    condsigma.e1_2x = e1_2.M4;
    condsigma.e2_2x = e2_2.M4;

    condsigma.e1c1_1s = e1c1_1.M1;
    condsigma.e2c1_1s = e2c1_1.M1;
    condsigma.e1c1_2s = e1c1_2.M1;
    condsigma.e2c1_2s = e2c1_2.M1;

    condsigma.e1c1_1i = e1c1_1.M2;
    condsigma.e2c1_1i = e2c1_1.M2;
    condsigma.e1c1_2i = e1c1_2.M2;
    condsigma.e2c1_2i = e2c1_2.M2;

    condsigma.e1c1_1d = e1c1_1.M3;
    condsigma.e2c1_1d = e2c1_1.M3;
    condsigma.e1c1_2d = e1c1_2.M3;
    condsigma.e2c1_2d = e2c1_2.M3;

    condsigma.e1c1_1x = e1c1_1.M4;
    condsigma.e2c1_1x = e2c1_1.M4;
    condsigma.e1c1_2x = e1c1_2.M4;
    condsigma.e2c1_2x = e2c1_2.M4;

    condsigma.e1c2_1s = e1c2_1.M1;
    condsigma.e2c2_1s = e2c2_1.M1;
    condsigma.e1c2_2s = e1c2_2.M1;
    condsigma.e2c2_2s = e2c2_2.M1;

    condsigma.e1c2_1i = e1c2_1.M2;
    condsigma.e2c2_1i = e2c2_1.M2;
    condsigma.e1c2_2i = e1c2_2.M2;
    condsigma.e2c2_2i = e2c2_2.M2;

    condsigma.e1c2_1d = e1c2_1.M3;
    condsigma.e2c2_1d = e2c2_1.M3;
    condsigma.e1c2_2d = e1c2_2.M3;
    condsigma.e2c2_2d = e2c2_2.M3;

    condsigma.e1c2_1x = e1c2_1.M4;
    condsigma.e2c2_1x = e2c2_1.M4;
    condsigma.e1c2_2x = e1c2_2.M4;
    condsigma.e2c2_2x = e2c2_2.M4;

    condsigma.e11s = e11.M1;
    condsigma.e12s = e12.M1;
    condsigma.e21s = e21.M1;
    condsigma.e22s = e22.M1;

    condsigma.e11i = e11.M2;
    condsigma.e12i = e12.M2;
    condsigma.e21i = e21.M2;
    condsigma.e22i = e22.M2;

    condsigma.e11d = e11.M3;
    condsigma.e12d = e12.M3;
    condsigma.e21d = e21.M3;
    condsigma.e22d = e22.M3;

    condsigma.e11x = e11.M4;
    condsigma.e12x = e12.M4;
    condsigma.e21x = e21.M4;
    condsigma.e22x = e22.M4;

    condsigma.us = sigu;
    condsigma.ui = isigu;
    condsigma.squd = sq_dsigu;

    /* Newton Raphson */
    unsigned j;
    for (j=0; j<iter; j++) {
      U = Dloglikfull(y0,b0,u0,condsigma,alph0,dalph0,tau0);
      H = D2loglikfull(y0,b0,u0,condsigma,alph0,dalph0,tau0);
      conv = (U(0)*U(0)+U(1)*U(1))/2;
      if (conv<_inner_NR_abseps) {
	warn(i) = 0;
	break;
      }
      u0 = u0-stepsize*U*H.i();
    }
    if (debug) {
      U = Dloglikfull(y0,b0,u0,condsigma,alph0,dalph0,tau0);
      vec L = loglikfull(y0,b0,u0,condsigma,alph0,dalph0,tau0);
      Rcpp::Rcout << "iter: " << j <<std::endl;
      Rcpp::Rcout << "conv: " << conv <<std::endl;
      Rcpp::Rcout << "L: " << L <<std::endl;
      Rcpp::Rcout << "U: " << U <<std::endl;
      Rcpp::Rcout << "i: " << i <<std::endl;
    }
    eb0.row(i) = u0;
  }
  return(eb0);
}
