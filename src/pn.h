#ifndef _PN_H_
#define _PN_H_

#include <RcppArmadillo.h>
#include <Rmath.h>

double pn(double y, double mu, double sigma);
double pn(arma::mat y, arma::mat mu, arma::mat sigma);

#endif /* _PN_H_ */
