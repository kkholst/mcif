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
#include "functions.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

const double twopi = 2*datum::pi;

/////////////////////////////////////////////////////////////////////////////////////////////////
/* Full loglikelihood */
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
    vmat sig = sigmaU; // Variance-covariance matrix of u
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
/* Score function of full loglikelihood */
rowvec Dloglikfull(unsigned row, DataPairs &data, const gmat &sigmaMarg, const gmat &sigmaJoint, const gmat &sigmaMargCond, vmat sigmaU, vec u, bool full=1){

  /* Estimation of pi, dpidu and dlogpidu */
  data.pi_gen(row, u);
  data.dpidu_gen(row, u);
  data.dlogpidu_gen(row, u);

  irowvec causes = data.causes_get(row); // Failure causes for pair in question

  rowvec res = zeros<rowvec>(data.ncauses); // Initialising output (score contribution)

  if ((causes(0) > 0) & (causes(1) > 0)){
    /* Both individuals experience failure */
    res = dlogdF2du(row, causes, data, sigmaJoint, u);
  }
  else if((causes(0) <= 0) & (causes(1) <= 0)){
    /* Neither individual experience failure */

    if ((causes(0) < 0) & (causes(1) < 0)){
      // Full follow-up for both individuals
      for (unsigned i=1; i<=2; i++){ // Over individuals
	double lik = 1;
	rowvec likdu = zeros<rowvec>(data.ncauses);
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  double prob = F1(row, j, i, data);
	  rowvec probdu = dF1du(row, j, i, data);
	  lik -= prob;
	  likdu -= probdu;
	}
	res += (1/lik)*likdu;
      }
    }
    else if (((causes(0) < 0) & (causes(1) == 0)) | ((causes(0) == 0) & (causes(1) < 0))){
      // Full follow-up for only one individual
      for (unsigned i=1; i<=2; i++){ // Over individuals
	double lik = 1;
	rowvec likdu = zeros<rowvec>(data.ncauses);
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  if (causes(i-1) < 0){
	    double prob = F1(row, j, i, data);
	    rowvec probdu = dF1du(row, j, i, data);
	    lik -= prob;
	    likdu -= probdu;
	  }
	  else {
	    double prob = F1(row, j, i, data, sigmaMarg, u);
	    rowvec probdu = dF1du(row, j, i, data, sigmaMarg, u);
	    lik -= prob;
	    likdu -= probdu;
	  }
	}
	res += (1/lik)*likdu;
      }
    }
    else {
      // Full follow-up for neither individual
      double lik = 1;
      rowvec likdu = zeros<rowvec>(data.ncauses);
      // Marginal probabilities
      for (unsigned i=1; i<=2; i++){ // Over individuals
	for (unsigned j=1; j<=data.ncauses; j++){ // Over failure causes
	  double prob = F1(row, j, i, data, sigmaMarg, u);
	  rowvec probdu = dF1du(row, j, i, data, sigmaMarg, u);
	  lik -= prob; // Subtracting
	  likdu -= probdu;
	}
      }
      // Bivariate probabilities
      for (unsigned k=1; k<=data.ncauses; k++){ // Over failure causes
	for (unsigned l=1; l<=data.ncauses; l++){
	  irowvec vcauses(2);
	  vcauses(0) = k; vcauses(1) = l;
	  double prob = F2(row, vcauses, data, sigmaJoint, u);
	  rowvec probdu = dF2du(row, vcauses, data, sigmaJoint, u);
	  lik += prob; // Adding
	  likdu += probdu;
	}
      }
      res = (1/lik)*likdu;
    }
  }
  else {
    /* One individual experiences failure the other does not */
    for (unsigned i=1; i<=2; i++){ // Over individuals
      unsigned cause = causes(i-1);
      if (cause > 0){
	// Marginal probability of failure
	res += dlogdF1du(row, cause, i, data, sigmaMarg, u);
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
	  rowvec likdu = zeros<rowvec>(data.ncauses);
	  if (cause < 0){ // Uncondtional
	    double prob = F1(row, j, i, data);
	    rowvec probdu = dF1du(row, j, i, data);
	    lik -= prob;
	    likdu -= probdu;
	  }
	  else { // Conditional
	    double prob = F1(row, j, i, cond_cause, data, sigmaMargCond, u);
	    rowvec probdu = dF1du(row, j, i, cond_cause, data, sigmaMargCond, u);
	    lik -= prob;
	    likdu -= probdu;
	  }
	  res += (1/lik)*likdu;
	}
      }
    }
  }
  /* Contribution from u */
  if (full){

    vmat sig = sigmaU; // Variance-covariance matrix etc. of u

    // Adding to the score
    res += -u.t()*sig.inv;
  };
  return(res);
}

/////////////////////////////////////////////////////////////////////////////
// FOR TESTING

double loglikout(unsigned row, mat sigma, vec u, int ncauses, imat causes, mat alpha, mat dalpha, mat beta, mat gamma){

  // Initialising gmats of sigma (Joint, Cond)
  gmat sigmaJoint = gmat(ncauses, ncauses);
  gmat sigmaCond = gmat(ncauses, ncauses);

  // Vectors for extracting rows and columns from sigma
  uvec rcJ(2); /* for joint */
  uvec rc1(1); /* for conditional */
  uvec rc2(ncauses+1); /* for conditional */
  uvec rcu(ncauses);
  for (int h=0; h<ncauses; h++){
    rcu(h) = ncauses + h;
  };

  // Estimating and setting vmats sigmaJoint
  for (int h=0; h<ncauses; h++){
    for (int i=0; i<ncauses; i++){
      rcJ(0)=h;
      rcJ(1)=ncauses+i;
      vmat x = vmat(sigma, rcJ, rcu);
      sigmaJoint.set(h,i,x);
    };
  };

  // Estimating and setting vmats of sigmaCond
  for (int h=0; h<ncauses; h++){
    for (int i=0; i<ncauses; i++){
      rc1(1) = h;
      rc2(0) = i;
      for (int j=0; j<ncauses; j++){
	rc2(j+1) = rcu(j);
      };
      vmat x = vmat(sigma, rc1, rc2);
      sigmaCond.set(h,i,x);
    };
  };

  // vmat of the us
  mat matU = sigma.submat(rcu,rcu);
  vmat sigmaU = vmat(matU);

  // Generating DataPairs
  DataPairs data = DataPairs(ncauses, causes, alpha, dalpha, beta, gamma);

  // Estimating likelihood contribution
  double loglik = loglikfull(unsigned row, data, sigmaJoint, sigmaCond, sigmaU, u, bool full=1);

  // Return
  return loglik;
};

//rowvec Dloglikout(unsigned row, mat sigma, mat data, vec u){

  // Generating gmats of sigma (Marg, Joint, MargCond, sigU)

  // Generating DataPairs

  // Estimating score contribution
//  rowvec score = Dloglikfull(unsigned row, DataPairs data, gmat sigmaJoint, gmat sigmaMargCond, vmat sigmaU, vec u, bool full=1);

  // Return
//  return score;
//};
