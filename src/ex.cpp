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
  int n = y.n_rows;

  vec pi1_1 = exp(b.col(0)+u.col(0))/(1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1)));
  vec pi1_2 = exp(b.col(1)+u.col(0))/(1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1)));
  vec pi2_1 = exp(b.col(2)+u.col(1))/(1+exp(b.col(0)+u.col(0))+exp(b.col(2)+u.col(1)));
  vec pi2_2 = exp(b.col(3)+u.col(1))/(1+exp(b.col(1)+u.col(0))+exp(b.col(3)+u.col(1)));

  vec res(n);

  for (int i=0; i<n; i++) {
    // ddF11
    if((y(i,0) == 1) & (y(i,1) == 1)){
      uvec rc1(2); rc1(0) = 0; rc1(1) = 1;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat ic_sig = c_sig.i();

      mat alph_sub(n,2); // alph1_1 & alph1_2
      alph_sub.col(0) = alph.col(0);
      alph_sub.col(1) = alph.col(1);

      mat alph_f = alph_sub.t();
      mat alph_c = alph_f.col(i)-c_mu;

      mat inner = alph_c.t()*c_sig.i()*alph_c;

      double pdf11 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      double ddF11 = pi1_1(i)*pi1_2(i)*dalph(i,0)*dalph(i,1)*pdf11; // pi1_1, pi1_2, dalph1_1, dalph1_2
      res(i) = log(ddF11);
    }
    // ddF12
    else if((y(i,0) == 1) & (y(i,1) == 2)){
      uvec rc1(2); rc1(0) = 0; rc1(1) = 3;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat ic_sig = c_sig.i();

      mat alph_sub(n,2); // alph1_1 & alph2_2
      alph_sub.col(0) = alph.col(0);
      alph_sub.col(1) = alph.col(3);

      mat alph_f = alph_sub.t();
      mat alph_c = alph_f.col(i)-c_mu;

      mat inner = alph_c.t()*c_sig.i()*alph_c;

      double pdf12 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      double ddF12 = pi1_1(i)*pi2_2(i)*dalph(i,0)*dalph(i,3)*pdf12; // pi1_1, pi2_2, dalph1_1, dalph2_2
      res(i) = log(ddF12);
    }
    // ddF21
    else if((y(i,0) == 2) & (y(i,1) == 1)){
      uvec rc1(2); rc1(0) = 2; rc1(1) = 1;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat ic_sig = c_sig.i();

      mat alph_sub(n,2); // alph2_1 & alph1_2
      alph_sub.col(0) = alph.col(2);
      alph_sub.col(1) = alph.col(1);

      mat alph_f = alph_sub.t();
      mat alph_c = alph_f.col(i)-c_mu;

      mat inner = alph_c.t()*c_sig.i()*alph_c;

      double pdf21 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      double ddF21 = pi2_1(i)*pi1_2(i)*dalph(i,2)*dalph(i,1)*pdf21; // pi2_1, pi1_2, dalph2_1, dalph1_2
      res(i) = log(ddF21);
    }
    // ddF22
    else if((y(i,0) == 2) & (y(i,1) == 2)){
      uvec rc1(2); rc1(0) = 2; rc1(1) = 3;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      mat c_sig = out.M1;
      mat ic_sig = c_sig.i();

      mat alph_sub(n,2); // alph2_1 & alph2_2
      alph_sub.col(0) = alph.col(2);
      alph_sub.col(1) = alph.col(3);

      mat alph_f = alph_sub.t();
      mat alph_c = alph_f.col(i)-c_mu;

      mat inner = alph_c.t()*c_sig.i()*alph_c;

      double pdf22 = 1/(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      double ddF22 = pi2_1(i)*pi2_2(i)*dalph(i,2)*dalph(i,3)*pdf22; // pi2_1, pi2_2, dalph2_1, dalph2_2
      res(i) = log(ddF22);
    }
    // dF01
    else if((y(i,0) == 0) & (y(i,1) == 1)){
      uvec rc1(1); rc1(0) = 0;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 1 ; rc4(1) = 4; rc4(2) = 5;

      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      // Marginal dF1_2
      vecmat out = conMuSig(sigma, mu, rc1, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;

      mat alph_sub(n,1); // alph1_2
      alph_sub.col(0) = alph.col(1);

      mat alph_f = alph_sub.t();
      mat alph_c = alph_f.col(i)-c_mu;

      mat inner = alph_c.t()*c_sig.i()*alph_c;

      double pdf1_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      double dF1_2 = pi1_2(i)*dalph(i,1)*pdf1_2; // pi1_2, dalph1_2

      // Conditional F1c1_1
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;

      mat alph_sub2(n,1); // alph1_1
      alph_sub2.col(0) = alph.col(0);

      mat alph_f2 = alph_sub2.t();
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      double F1c1_1 = pi1_1(i)*pn(alph_c2,sqrt(c_sig2[0]));

      // Conditional F2c1_1
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;

      mat alph_sub3(n,1); // alph2_1
      alph_sub3.col(0) = alph.col(2);

      mat alph_f3 = alph_sub3.t();
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      double F2c1_1 = pi2_1(i)*pn(alph_c3,sqrt(c_sig3[0]));

      // dF01
      double dF01 = dF1_2*(1-F1c1_1-F2c1_1);
      res(i) = log(dF01);
    }
    // dF10
    else if((y(i,0) == 1) & (y(i,1) == 0)){
      uvec rc1(1); rc1(0) = 1;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 0 ; rc4(1) = 4; rc4(2) = 5;

      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      // Marginal dF1_1
      vecmat out = conMuSig(sigma, mu, rc1, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;

      mat alph_sub(n,1); // alph1_1
      alph_sub.col(0) = alph.col(0);

      mat alph_f = alph_sub.t();
      mat alph_c = alph_f.col(i)-c_mu;

      mat inner = alph_c.t()*c_sig.i()*alph_c;

      double pdf1_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      double dF1_1 = pi1_1(i)*dalph(i,0)*pdf1_1; // pi1_1, dalph1_1

      // Conditional F1c1_2
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;

      mat alph_sub2(n,1); // alph1_2
      alph_sub2.col(0) = alph.col(1);

      mat alph_f2 = alph_sub2.t();
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      double F1c1_2 = pi1_2(i)*pn(alph_c2,sqrt(c_sig2[0]));

      // Conditional F2c1_2
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;

      mat alph_sub3(n,1); // alph2_2
      alph_sub3.col(0) = alph.col(3);

      mat alph_f3 = alph_sub3.t();
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      double F2c1_2 = pi2_2(i)*pn(alph_c3,sqrt(c_sig3[0]));

      // dF10
      double dF10 = dF1_1*(1-F1c1_2-F2c1_2);
      res(i) = log(dF10);
    }
    // dF02
    else if((y(i,0) == 0) & (y(i,1) == 2)){
      uvec rc1(1); rc1(0) = 0;
      uvec rc2(1); rc2(0) = 2;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 3 ; rc4(1) = 4; rc4(2) = 5;

      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      // Marginal dF2_2
      vecmat out = conMuSig(sigma, mu, rc2, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;

      mat alph_sub(n,1); // alph2_2
      alph_sub.col(0) = alph.col(3);

      mat alph_f = alph_sub.t();
      mat alph_c = alph_f.col(i)-c_mu;

      mat inner = alph_c.t()*c_sig.i()*alph_c;

      double pdf2_2 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      double dF2_2 = pi2_2(i)*dalph(i,3)*pdf2_2; // pi2_2, dalph2_2

      // Conditional F1c2_1
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;

      mat alph_sub2(n,1); // alph1_1
      alph_sub2.col(0) = alph.col(0);

      mat alph_f2 = alph_sub2.t();
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      double F1c2_1 = pi1_1(i)*pn(alph_c2,sqrt(c_sig2[0]));

      // Conditional F2c2_1
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;

      mat alph_sub3(n,1); // alph2_1
      alph_sub3.col(0) = alph.col(2);

      mat alph_f3 = alph_sub3.t();
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      double F2c2_1 = pi2_1(i)*pn(alph_c3,sqrt(c_sig3[0]));

      // dF02
      double dF02 = dF2_2*(1-F1c2_1-F2c2_1);
      res(i) = log(dF02);
    }
    // dF20
    else if((y(i,0) == 2) & (y(i,1) == 0)){
      uvec rc1(1); rc1(0) = 1;
      uvec rc2(1); rc2(0) = 3;
      uvec rc3(2); rc3(0) = 4; rc3(1) = 5;
      uvec rc4(3); rc4(0) = 2 ; rc4(1) = 4; rc4(2) = 5;

      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(5) = u(i,1);

      // Marginal dF2_1
      vecmat out = conMuSig(sigma, mu, rc1, rc3);
      vec c_mu = out.V;
      mat c_sig = out.M1;

      mat alph_sub(n,1); // alph2_1
      alph_sub.col(0) = alph.col(2);

      mat alph_f = alph_sub.t();
      mat alph_c = alph_f.col(i)-c_mu;

      mat inner = alph_c.t()*c_sig.i()*alph_c;

      double pdf2_1 = 1/sqrt(twopi)*1/sqrt(det(c_sig))*exp(-0.5*inner(0));

      double dF2_1 = pi2_1(i)*dalph(i,2)*pdf2_1; // pi2_1, dalph2_1

      // Conditional F1c2_1
      vecmat out2 = conMuSig(sigma, mu, rc1, rc4);
      vec c_mu2 = out2.V;
      mat c_sig2 = out2.M1;

      mat alph_sub2(n,1); // alph1_2
      alph_sub2.col(0) = alph.col(1);

      mat alph_f2 = alph_sub2.t();
      mat alph_c2 = alph_f2.col(i)-c_mu2;

      double F1c2_2 = pi1_2(i)*pn(alph_c2,sqrt(c_sig2[0]));

      // Conditional F2c2_2
      vecmat out3 = conMuSig(sigma, mu, rc2, rc4);
      vec c_mu3 = out3.V;
      mat c_sig3 = out3.M1;

      mat alph_sub3(n,1); // alph2_2
      alph_sub3.col(0) = alph.col(3);

      mat alph_f3 = alph_sub3.t();
      mat alph_c3 = alph_f3.col(i)-c_mu3;

      double F2c2_2 = pi2_2(i)*pn(alph_c3,sqrt(c_sig3[0]));

      // dF20
      double dF20 = dF2_1;//*(1-F1c2_2-F2c2_2);
      res(i) = log(dF20);
    }
    else{
      res(i) = 0;
    }
  }
  return(res);
}

  /* return(c_mu);


  vec res(n);

  for (int i=0; i<n; i++) {
    if((y(i,0) == 1) & (y(i,1) == 1)){
      uvec rc1(2); rc1(0) = 0; rc1(1) = 1;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;
      uvec rc3(2); rc3(0) = 0; rc3(1) = 1;
      vecmat out = conMuSig(sigma, mu, rc1, rc2);
      vec c_mu = out.V;
      vec tc_mu = c_mu.t;
      mat c_sig = out.M1;
      double r = c_sig(0,0)/(sqrt(c_sig(0,1))*sqrt(c_sig(1,0)));
      vec alph_c = alph(i,rc3)-tc_mu;
      double ddF11 = pi1_1(i)*pi1_2(i)*dalph(i,0)*dalph(i,1)*pn(alph_c,r);
      res(i) = log(ddF11);
    } else {
      res(i) = 0;
    }
    return res;
  }
  */


/*{
      uvec rc1(2); rc1(0) = 0; rc1(1) = 1;
      uvec rc2(2); rc2(0) = 4; rc2(1) = 5;

      Rcpp::Rcout << "here";

      vec mu(6) ;
      mu(0) = alph(i,0);
      mu(1) = alph(i,1);
      mu(2) = alph(i,2);
      mu(3) = alph(i,3);
      mu(4) = u(i,0);
      mu(4) = u(i,1);

      Rcpp::Rcout << "there";

      vecmat out = conMuSig(sigma, mu, rc1, rc2);

      Rcpp::Rcout << "everywhere";

      vec c_mu = out.V;
      rowvec tc_mu = c_mu.t();
      mat c_sig = out.M1;

      mat alph_sub(n,2); // alph1_1 & alph1_2
      alph_sub.col(0) = alph.col(0);
      alph_sub.col(1) = alph.col(1);

      double r = c_sig(0,1)/(sqrt(c_sig(0,0))*sqrt(c_sig(1,1)));

      mat alph_c = alph_sub.row(i)-tc_mu;

      double ddF11 = pi1_1(i)*pi1_2(i)*dalph(i,0)*dalph(i,1)*pn(alph_c,r); //dalph1_1 & dalph1_2
      res(i) = log(ddF11);
    }*/
