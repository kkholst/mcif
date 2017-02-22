library(devtools); library(Rcpp)
ncauses <- 2
causes <- matrix(c(1,1),1,2)
alpha <- matrix(c(0,1,4,3),1,4)
dalpha <- matrix(c(4,1,2,1),1,4)
beta <- matrix(c(0.3,3,2,1),1,4)
gamma <- matrix(c(3,2,1,2),1,4)
sourceCpp("P:/PhD/Scripts/mcif/src/test.cpp")
test(ncauses, causes, alpha, dalpha, beta, gamma)

