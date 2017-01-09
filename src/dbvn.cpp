#include "dbvn.h"
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

/////////////////////////////////////////////////////////////////////////////////////////////////
/* Function that calculates ... */
double dbvnorm(double y1, double y2, double R) {
  double detR = 1-R*R;
  // inv(R) = [1 -r; -r 1]/detR (prove by gauss elim.)
  double res = 1/(2*M_PI*sqrt(detR))*exp(-0.5/detR*(y1*y1+y2*y2-2*R*y1*y2));
  return(res);
}

vecmat Dbvn(double y1, double y2, double R) {
  vec DP(2);
  double R2 = R*R;
  DP(0) = Rf_dnorm4(y1,0.0,1.0,0)*Rf_pnorm5(y2-R*y1,0.0,sqrt(1-R2),1,0);
  DP(1) = Rf_dnorm4(y2,0.0,1.0,0)*Rf_pnorm5(y1-R*y2,0.0,sqrt(1-R2),1,0);
  mat HP(2,2);
  HP(1,0) = HP(0,1) = dbvnorm(y1,y2,R);
  HP(0,0) = -y1*DP(0) - R*HP(1,0);
  HP(1,1) = -y2*DP(1) - R*HP(1,0);
  vecmat res;
  res.V = DP;
  res.M= HP;
  return(res);
}
