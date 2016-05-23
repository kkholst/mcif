library(Rcpp)
sourceCpp("../src/ex.cpp")

y <- rbind(c(1,0),c(2,0));
r <- rbind(0.5)
loglik_ex(y,r)
log(pn(y,rbind(0.5,0.5)))
