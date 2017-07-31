#ifndef _DBVN_H_
#define _DBVN_H_

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

struct vecmat
{
arma::vec V;
arma::mat M;
};

double dbvnorm(double y1, double y2, double R);
vecmat Dbvn(double y1, double y2, double R);

#endif /* _DBVN_H_ */
