library(devtools); library(Rcpp); library(numDeriv)
source("helpfunctions.R")

u <- c(0.1,0.3)
ncauses <- 2
causes <- matrix(c(0,0),1,2)
alpha <- matrix(c(0,1,4,3),1,4)
dalpha <- matrix(c(4,1,2,1),1,4)
beta <- matrix(c(0.3,3,2,1),1,4)
gamma <- matrix(c(3,2,1,2),1,4)

# Generating sigma
vcv <- Posdef(n=4, ev=1:4)

#u <- c(0.1,0.3,0.2)
#ncauses <- 3
#causes <- matrix(c(3,3),1,2)
#alpha <- matrix(c(0,1,4,3,2,2),1,6)
#dalpha <- matrix(c(4,1,2,1,2,3),1,6)
#beta <- matrix(c(0.3,3,2,1,0.2,1),1,6)
#gamma <- matrix(c(3,2,1,2,0.4,1),1,6)

# Generating sigma
#vcv <- Posdef(n=6, ev=1:6)
sigma <- SigmaGen(vcv, ncauses, old=FALSE)
#sigma

# For old function
tau <- matrix(c(0,0),nrow=1,ncol=2) # 1 full follow-up, 0 not full follow-up
y <- matrix(c(0,0),nrow=1,ncol=2)
sigma2 <- SigmaGen(vcv,ncauses,old=TRUE)
#sigma2
um <- matrix(u,nrow=1,ncol=2)

alpha2 <- matrix(alpha[c(1,3,2,4)],nrow=1,ncol=4)
dalpha2 <- matrix(dalpha[c(1,3,2,4)],nrow=1,ncol=4)
beta2 <- matrix(beta[c(1,3,2,4)],nrow=1,ncol=4)
gamma2 <- matrix(gamma[c(1,3,2,4)],nrow=1,ncol=4)

# Numerical derivative
test <- function(u){
    loglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)
}
# Numerical derivative
Dtest <- function(u){
    Dloglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)
}

# Testing
sourceCpp("../../src/newall.cpp")
sourceCpp("../../src/loglik.cpp")

loglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)
logliktest(y, beta2, um, sigma2, alpha2-gamma2, dalpha2, tau)

Dloglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)
Dlogliktest(y, beta2, um, sigma2, alpha2-gamma2, dalpha2, tau)
jacobian(test,u)

D2loglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)
hessian(test,u)
jacobian(Dtest,u)

library(microbenchmark)

microbenchmark(loglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma),logliktest(y, beta2, um, sigma2, alpha2-gamma2, dalpha2, tau),times=1000)
microbenchmark(Dloglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma),Dlogliktest(y, beta2, um, sigma2, alpha2-gamma2, dalpha2, tau),times=1000)
