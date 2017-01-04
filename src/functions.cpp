#include "functions.h"

const double twopi = 2*datum::pi;
const double sq_twopi = sqrt(twopi);
const double loginvtwopi = log(1/twopi);
const double loginvsqtwopi = log(1/sq_twopi);

/////////////////////////////////////////////////////////////////////////////////////////////////
double logdF1(unsigned row, const unsigned &cause, const unsigned &indiv, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(cause);
  double cond_mean = as_scalar(cond_sig.proj*u);

  double alp = data.alphaMarg_get(row, cause, indiv);
  double gam = data.gammaMarg_get(row, cause, indiv);
  double alpgam = alp - gam;
  double inner = pow((alpgam-cond_mean),2)*as_scalar(cond_sig.inv);

  double logpdf = loginvsqtwopi + cond_sig.loginvsqdet + log(data.dalphaMarg_get(row, cause, indiv)) - 0.5*inner;

  double logdF1 = log(data.piMarg_get(row, cause, indiv)) + logpdf;
  return(logdF1);
};

double dlogdF1du(unsigned row, const unsigned &cause, const unsigned &indiv, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(cause);
  double cond_mean = as_scalar(cond_sig.proj*u);

  double alp = data.alphaMarg_get(row, cause, indiv);
  double gam = data.gammaMarg_get(row, cause, indiv);
  double alpgam = (alp - gam);
  double dinnerdu = (alpgam - cond_mean)*as_scalar(cond_sig.inv);

  double dlogdF1du = data.dlogpiduMarg_get(row, cause, indiv)) + dinnerdu;
  return(dlogdF1du);
};

/////////////////////////////////////////////////////////////////////////////////////////////////
double logdF2(unsigned row, const irowvec &causes, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(causes);
  vec cond_mean = cond_sig.proj*u;

  vec alp = data.alpha_get(row, causes);
  vec gam = data.gamma_get(row, causes);
  vec c_alpgam = (alp - gam) - cond_mean;
  double inner = as_scalar(c_alpgam.t()*cond_sig.inv*c_alpgam);

  double logpdf = loginvtwopi + cond_sig.loginvsqdet + log(data.dalphaMarg_get(row, causes(0), 1)) + log(data.dalphaMarg_get(row, causes(1), 2)) - 0.5*inner;

  double logdF2 = log(data.piMarg_get(row, causes(0), 1)) + log(data.piMarg_get(row, causes(1), 2)) + logpdf;
  return(logdF2);
};

double dlogdF2du(unsigned row, const irowvec &causes, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(causes);
  vec cond_mean = cond_sig.proj*u;

  vec alp = data.alphaMarg_get(row, causes);
  vec gam = data.gammaMarg_get(row, causes);
  vec c_alpgam = (alp - gam) - cond_mean;

  double dinnerdu = c_alpgam.t()*cond_sig.inv*cond_sig.proj;

  double dlogdF2du = data.dlogpiduMarg_get(row, causes(0), 1)) + data.dlogpiduMarg_get(row, causes(1), 2)) + innerdu;
  return(dlogdF2du);
};

/////////////////////////////////////////////////////////////////////////////////////////////////

/* Marginal */
double F1(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(cause);
  double cond_mean = as_scalar(cond_sig.proj*u);

  double alp = data.alphaMarg_get(row, cause, indiv);
  double gam = data.gammaMarg_get(row, cause, indiv);
  double alpgam = alp - gam;

  double F1 = data.piMarg_get(row, cause, indiv)*pn(alpgam, cond_mean, as_scalar(cond_sig.vcov));
  return(F1);
};

/* Full follow-up */
double F1(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data){
  double F1 = data.piMarg_get(row, cause, indiv);
  return(F1);
};

/* Conditional on other individual */
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

  // Attaining variance covariance matrix etc. (conditional on u and other individual)
  vmat cond_sig = sigma(cause,cond_cause);
  double cond_mean = as_scalar(cond_sig.proj*alpgamu);

  double alp = data.alphaMarg_get(row, cause, indiv);
  double gam = data.gammaMarg_get(row, cause, indiv);
  double alpgam = alp - gam;

  double F1 = data.piMarg_get(row, cause, indiv)*pn(alpgam, cond_mean, as_scalar(cond_sig.vcov));
  return(F1);
};

/////////////////////////////////////////////////////////////////////////////////////////////////
double F2(unsigned row, irowvec causes, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(causes);
  vec cond_mean = cond_sig.proj*u;

  vec alp = data.alpha_get(row, causes);
  vec gam = data.gamma_get(row, causes);
  vec alpgam = alp - gam;

  double F2 = data.piMarg_get(row, causes(0), 1)*data.piMarg_get(row, causes(1), 2)*pn(alpgam, cond_mean,cond_sig.vcov);
  return(F2);
};
