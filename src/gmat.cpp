// [[Rcpp::depends(RcppArmadillo)]]

#include "vmat.h"
#include "gmat.h"

gmat::gmat(unsigned nrow, unsigned ncol) {
  vmat vm = vmat(2);
  this->_ncol = ncol;
  this->_nrow = nrow;
  double pos = 0;
  for (unsigned i=0; i<ncol; i++) {
    for (unsigned j=0; j<nrow; j++) {
      std::cerr << "i,j" << i << "," << j << std::endl;
      pos++;
      vm.loginvsqdet = pos;
      vecvmat.push_back(vm);
    }
  }
}
