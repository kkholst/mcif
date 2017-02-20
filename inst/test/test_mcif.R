library(devtools); library(Rcpp)
install_github("kkholst/mcif", force=TRUE)

library(mcif)
source("helpfunctions.R")

ncauses <- 2
ncauses <- 3

if (ncauses==2){
    u <- c(5,2)
    causes <- matrix(1,2,c(0,1))
    alpha <- matrix(1,2,c(1,2))
    dalpha <- matrix(1,2,c(4,1))
    beta <- matrix(1,2,c(0.3,3))
    gamma <- matrix(1,2,c(3,2))
    vcv <- Posdef(n=4, ev=1:4)
}

if (ncauses==3){
    u <- c(5,2,3)
    causes <- matrix(1,3,c(0,1,4))
    alpha <- matrix(1,3,c(0,1,2))
    dalpha <- matrix(1,3,c(4,1,1))
    beta <- matrix(1,3,c(0.3,3,2))
    gamma <- matrix(1,3,c(3,2,3))
    vcv <- Posdef(n=6, ev=1:6)
}

# Generating sigma
sigma <- SigmaGen(vcv, ncauses)

# Testing
loglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)


