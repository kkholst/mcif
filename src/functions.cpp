#include "functions.h"

const double twopi = 2*datum::pi;
const double sq_twopi = sqrt(twopi);
const double invsqtwopi = 1/sq_twopi;
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

rowvec dlogdF1du(unsigned row, const unsigned &cause, const unsigned &indiv, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(cause);
  double cond_mean = as_scalar(cond_sig.proj*u);

  double alp = data.alphaMarg_get(row, cause, indiv);
  double gam = data.gammaMarg_get(row, cause, indiv);
  double alpgam = (alp - gam);
  rowvec dinnerdu = (alpgam - cond_mean)*as_scalar(cond_sig.inv)*cond_sig.proj;

  rowvec dlogdF1du = data.dlogpiduMarg_get(row, cause, indiv) + dinnerdu;
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

rowvec dlogdF2du(unsigned row, const irowvec &causes, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(causes);
  vec cond_mean = cond_sig.proj*u;

  vec alp = data.alpha_get(row, causes);
  vec gam = data.gamma_get(row, causes);
  vec c_alpgam = (alp - gam) - cond_mean;

  rowvec dinnerdu = c_alpgam.t()*cond_sig.inv*cond_sig.proj;

  rowvec dlogdF2du = data.dlogpiduMarg_get(row, causes(0), 1) + data.dlogpiduMarg_get(row, causes(1), 2) + dinnerdu;
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

rowvec dF1du(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(cause);
  double cond_mean = as_scalar(cond_sig.proj*u);

  double alp = data.alphaMarg_get(row, cause, indiv);
  double gam = data.gammaMarg_get(row, cause, indiv);
  double alpgam = alp - gam;
  double c_alpgam = alpgam - cond_mean;

  double inner = pow((c_alpgam),2)*as_scalar(cond_sig.inv);

  double cdf = pn(alpgam, cond_mean, as_scalar(cond_sig.vcov));
  double pdf = invsqtwopi*sqrt(as_scalar(cond_sig.vcov))*exp(-0.5*inner);

  rowvec dcdfdu = pdf*(-cond_sig.proj);

  rowvec dF1du = data.dpiduMarg_get(row, cause, indiv)*cdf + data.piMarg_get(row, cause, indiv)*dcdfdu;
  return(dF1du);
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

rowvec dF1du(unsigned row, unsigned cause, unsigned indiv, unsigned cond_cause, const DataPairs &data, const gmat &sigma, vec u){

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
  double c_alpgam = alpgam - cond_mean;

  double inner = pow((c_alpgam),2)*as_scalar(cond_sig.inv);

  double cdf = pn(alpgam, cond_mean, as_scalar(cond_sig.vcov));
  double pdf = invsqtwopi*sqrt(as_scalar(cond_sig.vcov))*exp(-0.5*inner);

  rowvec dcdfdu = pdf*(-cond_sig.proj);

  rowvec dF1du = data.dpiduMarg_get(row, cause, indiv)*cdf + data.piMarg_get(row, cause, indiv)*dcdfdu;
  return(dF1du);
};

/* Full follow-up */
double F1(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data){
  double F1 = data.piMarg_get(row, cause, indiv);
  return(F1);
};

rowvec dF1du(unsigned row, unsigned cause, unsigned indiv, const DataPairs &data){
  rowvec dF1du = data.dpiduMarg_get(row, cause, indiv);
  return(dF1du);
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

rowvec dF2du(unsigned row, irowvec causes, const DataPairs &data, const gmat &sigma, vec u){

  // Attaining variance covariance matrix etc. (conditional on u)
  vmat cond_sig = sigma(causes);
  vec cond_mean = cond_sig.proj*u;

  vec alp = data.alpha_get(row, causes);
  vec gam = data.gamma_get(row, causes);
  vec alpgam = alp - gam;

  double cdf = pn(alpgam, cond_mean, cond_sig.vcov);

  vec ll = sqrt(diagvec(cond_sig.vcov));
  mat Lambda = diagmat(ll);
  mat iLambda = diagmat(1/ll);
  mat R = iLambda*cond_sig.vcov*iLambda;
  mat LR = Lambda*R;

  double r = R(0,1);
  rowvec ytilde = iLambda*(alpgam - cond_mean);

  vecmat D = Dbvn(ytilde(0),ytilde(1),r);
  mat M = -LR*D.V;

  rowvec dcdfdu = cond_sig.proj*cond_sig.inv*M;

  rowvec dF2du_1 = data.dpiduMarg_get(row, causes(0), 1)*data.piMarg_get(row, causes(1), 2)*cdf ;
  rowvec dF2du_2 = data.dpiduMarg_get(row, causes(1), 2)*data.piMarg_get(row, causes(0), 1)*cdf;
  rowvec dF2du_3 = data.piMarg_get(row, causes(0), 1)*data.piMarg_get(row, causes(1), 2)*dcdfdu;

  rowvec dF2du = dF2du_1 + dF2du_2 + dF2du_3;
  return(dF2du);
};

