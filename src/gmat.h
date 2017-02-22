#ifndef GMAT_H
#define GMAT_H

#include <RcppArmadillo.h>
#include <vector>
#include <cassert>
#include "vmat.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

class gmat {
private:
  unsigned _ncol;
  unsigned _nrow;
  vector<vmat> vecvmat;
public:
  gmat(unsigned nrow, unsigned ncol);
  vmat operator()(unsigned row, unsigned col) const {
    assert((row>0) & (row<_nrow));
    assert((col>0) & (col<_ncol));
    unsigned pos = row+col*_nrow;
    return(vecvmat[pos]);
  }
  vmat operator()(irowvec idx) const {
    assert(idx.n_elem==2);
    unsigned row = idx[0] - 1;
    unsigned col = idx[1] - 1;
    unsigned pos = row+col*_nrow;
    return(vecvmat[pos]);
    // (this->)(idx[0],idx[1]);
  }
  vmat operator()(unsigned i) const {
    unsigned pos = i - 1;
    return(vecvmat[pos]);
  }
  void set(unsigned row, unsigned col, const vmat &x) {
    unsigned pos = row+col*_nrow;
    this->vecvmat[pos] = x;
  }
  unsigned ncol() const {
    return(this->_ncol);
  }
  unsigned nrow() const {
    return(this->_nrow);
  }
};

#endif /* GMAT_H */
