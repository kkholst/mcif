#install_github("kkholst/mcif", force=TRUE)
#install("P:/PhD/Scripts/mcif", repos = NULL, force=TRUE)
#compileAttributes("mcif")
#sourceCpp("P:/PhD/Scripts/mcif/src/newall.cpp")
#detach("package:mcif", unload=TRUE)
#library(mcif)

library(devtools); library(Rcpp)
source("helpfunctions.R")

u <- c(5,2)
ncauses <- 2
causes <- matrix(c(1,1),1,2)
alpha <- matrix(c(0,1,4,3),1,4)
dalpha <- matrix(c(4,1,2,1),1,4)
beta <- matrix(c(0.3,3,2,1),1,4)
gamma <- matrix(c(3,2,1,2),1,4)

# Generating sigma
vcv <- Posdef(n=4, ev=1:4)
sigma <- SigmaGen(vcv, ncauses)

sigma

# Testing
sourceCpp("P:/PhD/Scripts/mcif/src/newall.cpp")
loglikout(sigma, u, ncauses, causes, alpha, dalpha, beta, gamma)

rcJ <- c(3,6)
rcu <- c(7,8,9)

sigma11 <- sigma[rcJ,rcJ]
sigma22 <- sigma[rcu,rcu]
sigma12 <- sigma[rcJ,rcu]
sigma21 <- sigma[rcu,rcJ]

vcov.new <- sigma11 - sigma12%*%solve(sigma22)%*%sigma21
inv.new <- solve(vcov.new)
loginvsqdet.new <- log(1/sqrt(det(vcov.new)))
proj.new <- sigma12%*%solve(sigma22)

vcov.new
inv.new
loginvsqdet.new
proj.new
