// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>

#include <vector>
#include <cassert>

//#include "vmat.h"
//#include "gmat.h"
#include "DataPairs.h"
//#include "quadrule.h"
//#include "pn.h"
//#include "functions.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//const double twopi = 2*datum::pi;

// [[Rcpp::export]]
double test(int ncauses, imat causes, mat alpha, mat dalpha, mat beta, mat gamma){

  // Generating DataPairs
  DataPairs data = DataPairs(ncauses, causes, alpha, dalpha, beta, gamma);

  double out = -(data.ncauses/2);

  // Return
  return out;
};
