# Loading libraries
library(devtools); library(Rcpp); library(numDeriv);

# Sourcing help functions etc.
source("helpfunctions_test.R")
source("dataprep.R")

# Sourcing loglik
sourceCpp("P:/PhD/Scripts/mcif/src/loglik.cpp")
#sourceCpp("../../loglik.cpp")

# Loading data
load("data.Rdata")

# Number of causes (must match number of causes in loaded data)
ncauses <- 2

# Generating positive definite sigma
vcv <- Posdef(n=4, ev=1:4) # For events
sigma <- SigmaGen(vcv, ncauses, old=FALSE) # Adding u and getting right format for loglik: (ncauses*3)*(ncauses*3)

# Prepping data
datprep <- data.prep(data, time, status="event", cova=NULL)
causes <- datprep$causes
eb0 <- t(datprep$eb0) # Must for some reason be transposed...

# Prepping alpha
a1 <- 3
a2 <- 2
alpha <- cbind(datprep$gt1*a1, datprep$gt1*a2, datprep$gt2*a1, datprep$gt2*a2)

# Prepping dalpha
dalpha <- cbind(datprep$dgt1*a1, datprep$dgt1*a2, datprep$dgt2*a1, datprep$dgt2*a2)

# Prepping gamma
g1 <- -1
g2 <- -1.2
gamma <- cbind(datprep$x.1%*%g1, datprep$x.1%*%g2, datprep$x.2%*%g1, datprep$x.2%*%g2)

# Prepping beta
b1 <- -1.9
b2 <- -0.2
beta <- cbind(datprep$x.1%*%b1, datprep$x.1%*%b2, datprep$x.2%*%b1, datprep$x.2%*%b2)

#-----------------------------------------------------------------------
# Estimation of loglikelihood
#-----------------------------------------------------------------------
loglik(sigma, ncauses, causes, alpha, dalpha, beta, gamma, eb0, nq=1)

#----------------------------------------------------------------------------------
# Miscellaneous checks - operate on rows, takes first row in data
u <- c(0.1,0.3)

# Full loglik
test <- function(u){
    loglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)
}
# Full Dloglik
Dtest <- function(u){
    Dloglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)
}

# Full D2loglik
D2test <- function(u){
    D2loglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)
}

# Checking full loglik
test(u)

# Checking full Dloglik
Dtest(u)
jacobian(test,u)

# Checking full D2loglik
D2test(u)
hessian(test,u)
jacobian(Dtest,u)


