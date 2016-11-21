// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rmath.h>

#include <vector>
#include <cassert>

#include "vmat.h"
#include "gmat.h"
#include "DataPairs.h"
#include "quadrule.h"
#include "pn.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

const double twopi = 2*datum::pi;
const double sq_twopi = sqrt(twopi);
const double loginvtwopi = log(1/twopi);
const double loginvsqtwopi = log(1/sq_twopi);

/////////////////////////////////////////////////////////////////////////////////////////////////

double logdF1(unsigned row, const unsigned &cause, const unsigned &indiv, const DataPairs &data, const gmat &sigma, vec u){

  vmat cond_sig = sigma(cause); // Attaining appropriate conditional sigma etc. (conditional on u)
  double cond_mean = as_scalar(cond_sig.proj*u); // Conditional mean (conditional on u)

  double alp = data.alphaMarg_get(row, cause, indiv);
  double gam = data.gammaMarg_get(row, cause, indiv);
  double alpgam = alp - gam;
  double inner = pow((alpgam-cond_mean),2)*as_scalar(cond_sig.inv);

  double logpdf = loginvsqtwopi + cond_sig.loginvsqdet + log(data.dalphaMarg_get(row, cause, indiv)) - 0.5*inner;

  double logdF1 = log(data.piMarg_get(row, cause, indiv)) + logpdf;
  return(logdF1);
};

/////////////////////////////////////////////////////////////////////////////////////////////////

double logdF2(unsigned row, const irowvec &causes, const DataPairs &data, const gmat &sigma, vec u){

  vmat cond_sig = sigma(causes); // Attaining appropriate conditional sigma etc. (conditional on u)
  vec cond_mean = cond_sig.proj*u; // Conditional mean (conditional on u)

  vec alp = data.alpha_get(row, causes);
  vec gam = data.gamma_get(row, causes);
  vec c_alpgam = (alp - gam) - cond_mean;
  double inner = as_scalar(c_alpgam.t()*cond_sig.inv*c_alpgam);

  double logpdf = loginvtwopi + cond_sig.loginvsqdet + log(data.dalphaMarg_get(row, causes(0), 1)) + log(data.dalphaMarg_get(row, causes(1), 2)) - 0.5*inner;

  double logdF2 = log(data.piMarg_get(row, causes(0), 1)) + log(data.piMarg_get(row, causes(1), 2)) + logpdf;
  return(logdF2);
};

/////////////////////////////////////////////////////////////////////////////////////////////////

double F1(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data, const gmat &sigma, vec u){

  vmat cond_sig = sigma(cause); // Attaining appropriate conditional sigma etc. (conditional on u)
  double cond_mean = as_scalar(cond_sig.proj*u); // Conditional mean (conditional on u and potentially alpgam of other individual)

  double alp = data.alphaMarg_get(row, cause, indiv);
  double gam = data.gammaMarg_get(row, cause, indiv);
  double alpgam = alp - gam;

  double F1 = data.piMarg_get(row, cause, indiv)*pn(alpgam, cond_mean, as_scalar(cond_sig.vcov));
  return(F1);
};

double F1(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data){
  double F1 = data.piMarg_get(row, cause, indiv);
  return(F1);
};

double F1(unsigned row, unsigned cause, unsigned indiv, unsigned cond_cause, const DataPairs &data, const gmat &sigma, vec u){

  // Other individual
  unsigned cond_indiv;
  if (indiv==1){
    cond_indiv=2;
  }
  else {
    cond_indiv=1;
  }

  // Alpgam of other individual
  double cond_alp = data.alphaMarg_get(row, cond_cause, cond_indiv);
  double cond_gam = data.gammaMarg_get(row, cond_cause, cond_indiv);
  double cond_alpgam = cond_alp - cond_gam;

  vec vcond_alpgam(1); vcond_alpgam(0) = cond_alpgam;

  // Joining u vector and alpgam from other individual
  vec alpgamu = join_cols(vcond_alpgam, u);

  vmat cond_sig = sigma(cause,cond_cause); // Attaining appropriate conditional sigma etc. (conditional on u and variance of other individual)
  double cond_mean = as_scalar(cond_sig.proj*alpgamu); // Conditional mean (conditional on alpgam of other individual and u)

  double alp = data.alphaMarg_get(row, cause, indiv);
  double gam = data.gammaMarg_get(row, cause, indiv);
  double alpgam = alp - gam;

  double F1 = data.piMarg_get(row, cause, indiv)*pn(alpgam, cond_mean, as_scalar(cond_sig.vcov));
  return(F1);
};

/////////////////////////////////////////////////////////////////////////////////////////////////

double F2(unsigned row, irowvec causes, const DataPairs &data, const gmat &sigma, vec u){
  vmat cond_sig = sigma(causes); // Attaining appropriate conditional sigma etc. (conditional on u)
  vec cond_mean = cond_sig.proj*u; // Conditional mean (conditional on u)

  vec alp = data.alpha_get(row, causes);
  vec gam = data.gamma_get(row, causes);
  vec alpgam = alp - gam;

  double F2 = data.piMarg_get(row, causes(0), 1)*data.piMarg_get(row, causes(1), 2)*pn(alpgam, cond_mean,cond_sig.vcov);
  return(F2);
};

/////////////////////////////////////////////////////////////////////////////////////////////////

double loglikfull(unsigned row, DataPairs &data, const gmat &sigmaMarg, const gmat &sigmaJoint, const gmat &sigmaMargCond, vmat sigmaU, vec u, bool full=1){

  data.pi_gen(row, u); // Estimation of pi based on u

  irowvec causes = data.causes_get(row); // Failure causes for pair in question

  double res = 0; // Initialising output (loglik contribution)

  if ((causes(0) > 0) & (causes(1) > 0)){
    /* Both individuals experience failure */
    res = logdF2(row, causes, data, sigmaJoint, u);
  }
  else if((causes(0) <= 0) & (causes(1) <= 0)){
    /* Neither individual experience failure */

    if ((causes(0) < 0) & (causes(1) < 0)){
      // Full follow-up for both individuals
      for (unsigned i=1; i<=2; i++){ // Over individuals
	double lik = 1;
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  double prob = F1(row, j, i, data);
	  lik -= prob;
	}
	res += log(lik);
      }
    }
    else if (((causes(0) < 0) & (causes(1) == 0)) | ((causes(0) == 0) & (causes(1) < 0))){
      // Full follow-up for only one individual
      for (unsigned i=1; i<=2; i++){ // Over individuals
	double lik = 1;
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  if (causes(i-1) < 0){
	    double prob = F1(row, j, i, data);
	    lik -= prob;
	  }
	  else {
	    double prob = F1(row, j, i, data, sigmaMarg, u);
	    lik -= prob;
	  }
	}
	res += log(lik);
      }
    }
    else {
      // Full follow-up for neither individual
      double lik = 1;
      // Marginal probabilities
      for (unsigned i=1; i<=2; i++){ // Over individuals
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  double prob = F1(row, j, i, data, sigmaMarg, u);
	  lik -= prob; // Subtracting
	}
      }
      // Bivariate probabilities
      for (unsigned k=1; k<=data.ncauses; k++){ // Over failure causes
	for (unsigned l=1; l<=data.ncauses; l++){
	  irowvec vcauses(2);
	  vcauses(0) = k; vcauses(1) = l;
	  double prob = F2(row, vcauses, data, sigmaJoint, u);
	  lik += prob; // Adding
	}
      }
      res = log(lik);
    }
  }
  else {
    /* One individual experiences failure the other does not */
    for (unsigned i=1; i<=2; i++){ // Over individuals
      unsigned cause = causes(i-1);
      if (cause > 0){
	// Marginal probability of failure
	res += logdF1(row, cause, i, data, sigmaMarg, u);
      }
      else {
	// Marginal probability of no failure
	unsigned cond_cause;
	if (i==1){
	  cond_cause = causes(1);
	}
	else {
	  cond_cause = causes(0);
	}
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  double lik = 1;
	  if (cause < 0){
	    // Unconditional
	    double prob = F1(row, j, i, data);
	    lik -= prob;
	  }
	  else {
	    // Conditional
	    double prob = F1(row, j, i, cond_cause, data, sigmaMargCond, u);
	    lik -= prob;
	  }
	  res += log(lik);
	}
      }
    }
  }
  /* Contribution from u */
  if (full){

    vmat sig = sigmaU; // Variance-covariance matrix etc. of u
    double inner = as_scalar(u*sig.inv*u.t());

    // PDF of u
    double logpdfu = log(pow(twopi,-(data.ncauses/2))) + sig.loginvsqdet - 0.5*inner;

    // Adding to the loglik
    res += logpdfu;
  }

  /* Return */
  return(res);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


vec Dloglikfull(unsigned row, DataPairs &data, const gmat &sigmaMarg, const gmat &sigmaJoint, const gmat &sigmaMargCond, vmat sigmaU, vec u, bool full=1){

  data.pi_gen(row, u); // Estimation of pi based on u

  // Estimation of the derivative of pi and pi*pi


  irowvec causes = data.causes_get(row); // Failure causes for pair in question

  vec res(data.ncauses); // Initialising output (score contribution)

  if ((causes(0) > 0) & (causes(1) > 0)){
    /* Both individuals experience failure */
    res = logdF2(row, causes, data, sigmaJoint, u);
  }
  else if((causes(0) <= 0) & (causes(1) <= 0)){
    /* Neither individual experience failure */

    if ((causes(0) < 0) & (causes(1) < 0)){
      // Full follow-up for both individuals
      for (unsigned i=1; i<=2; i++){ // Over individuals
	double lik = 1;
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  double prob = F1(row, j, i, data);
	  lik -= prob;
	}
	res += log(lik);
      }
    }
    else if (((causes(0) < 0) & (causes(1) == 0)) | ((causes(0) == 0) & (causes(1) < 0))){
      // Full follow-up for only one individual
      for (unsigned i=1; i<=2; i++){ // Over individuals
	double lik = 1;
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  if (causes(i-1) < 0){
	    double prob = F1(row, j, i, data);
	    lik -= prob;
	  }
	  else {
	    double prob = F1(row, j, i, data, sigmaMarg, u);
	    lik -= prob;
	  }
	}
	res += log(lik);
      }
    }
    else {
      // Full follow-up for neither individual
      double lik = 1;
      // Marginal probabilities
      for (unsigned i=1; i<=2; i++){ // Over individuals
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  double prob = F1(row, j, i, data, sigmaMarg, u);
	  lik -= prob; // Subtracting
	}
      }
      // Bivariate probabilities
      for (unsigned k=1; k<=data.ncauses; k++){ // Over failure causes
	for (unsigned l=1; l<=data.ncauses; l++){
	  irowvec vcauses(2);
	  vcauses(0) = k; vcauses(1) = l;
	  double prob = F2(row, vcauses, data, sigmaJoint, u);
	  lik += prob; // Adding
	}
      }
      res = log(lik);
    }
  }
  else {
    /* One individual experiences failure the other does not */
    for (unsigned i=1; i<=2; i++){ // Over individuals
      unsigned cause = causes(i-1);
      if (cause > 0){
	// Marginal probability of failure
	res += logdF1(row, cause, i, data, sigmaMarg, u);
      }
      else {
	// Marginal probability of no failure
	unsigned cond_cause;
	if (i==1){
	  cond_cause = causes(1);
	}
	else {
	  cond_cause = causes(0);
	}
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  double lik = 1;
	  if (cause < 0){ // Uncondtional
	    double prob = F1(row, j, i, data);
	    lik -= prob;
	  }
	  else { // Conditional
	    double prob = F1(row, j, i, cond_cause, data, sigmaMargCond, u);
	    lik -= prob;
	  }
	  res += log(lik);
	}
      }
    }
  }
  /* Contribution from u */
  if (full){

    vmat sig = sigmaU; // Variance-covariance matrix etc. of u
    double inner = as_scalar(u*sig.inv*u.t());

    // PDF of u
    double logpdfu = log(pow(twopi,-(data.ncauses/2))) + sig.loginvsqdet - 0.5*inner;

    // Adding to the loglik
    res += logpdfu;
  }

}

