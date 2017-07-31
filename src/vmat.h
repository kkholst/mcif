#ifndef VMAT_H
#define VMAT_H

#include <RcppArmadillo.h>
#include <vector>

using namespace arma;
using namespace std;

class vmat {
public:
  mat vcov;
  mat inv;
  double loginvsqdet;
  mat proj;

  arma::uvec dim() {
    arma::uvec res(2);
    res[0] = (unsigned)(this->vcov).n_cols;
    res[1] = (unsigned)(this->vcov).n_rows;
    return(res);
  }
  // Assignment operator:
  vmat& operator=(const vmat &rhs);

  // Constructor
  vmat(const mat& x, const uvec& rc1, const uvec& rc2);
  vmat(const mat &x);
  vmat(unsigned dim=2);
};

#endif /* VMAT_H */
