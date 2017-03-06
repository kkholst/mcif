#include "vmat.h"

// Member functions definitions including constructor
vmat& vmat::operator=(const vmat &rhs) {
    if (this != &rhs) {
      this->vcov = rhs.vcov;
      this->inv = rhs.inv;
      this->loginvsqdet = rhs.loginvsqdet;
      this->proj = rhs.proj;
    }
    return *this;
}

vmat::vmat(const arma::mat &x, const arma::uvec &rc1, const arma::uvec &rc2) {
  arma::mat x11 = x.submat(rc1,rc1);
  arma::mat x12 = x.submat(rc1,rc2);
  arma::mat x21 = x.submat(rc2,rc1);
  arma::mat x22 = x.submat(rc2,rc2);

  proj = x12*arma::inv(x22);
  vcov = x11-proj*x21;
  inv = arma::inv_sympd(vcov);
  loginvsqdet = log(1/sqrt(arma::det(vcov)));
}

vmat::vmat(const mat &x) {
  vcov = x;
  inv = arma::inv_sympd(x);
  loginvsqdet = log(1/sqrt(arma::det(x)));
}

vmat::vmat(unsigned dim) {
  mat x0 = arma::mat(dim,dim);
  vcov = x0;
  inv = x0;
  loginvsqdet = 0;
  proj = arma::mat(1,dim);
}
