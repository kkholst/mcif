cif.loglik <- function(datprep, par, nq, stepsize, iter, grad=FALSE){

##-----------------------------------------------------------------------
## Parameters
##-----------------------------------------------------------------------
b1 <- par[1]
b2 <- par[2]

g1 <- par[3]
g2 <- par[4]

a1 <- exp(par[5])
a2 <- exp(par[6])

## Parameters for variance-covariance
s1 <- par[7]
s2 <- par[8]
s3 <- par[9]
s4 <- par[10]
s5 <- par[11]

s6 <- par[12]
s7 <- par[13]
s8 <- par[14]
s9 <- par[15]
s10 <- par[16]

##-----------------------------------------------------------------------
## Sigma
##-----------------------------------------------------------------------
## Cholesky (upper triangle)
L.chol <- matrix(c(exp(s1),0,0,0,
                     s5,exp(s2),0,0,
                     s6,s8,exp(s3),0,
                     s7,s9,s10,exp(s4)),
                   nrow=4, ncol=4)

## Variance-covariance matrix of eta and u
vcv <- t(L.chol)%*%L.chol

## Sigma for likelihood
sigma <- SigmaGen(vcv, 2)

##-----------------------------------------------------------------------
## The rest
##-----------------------------------------------------------------------
## Causes
causes <- datprep$causes
ncauses <- 2

## Weights
w <- datprep$weights

## Empical Bayes
eb0 <- t(datprep$eb0) # Must for some reason (aka my coding [LC]) be transposed...

## Prepping beta
beta <- cbind(datprep$x.1%*%b1, datprep$x.1%*%b2, datprep$x.2%*%b1, datprep$x.2%*%b2)

## Prepping alpha
alpha <- cbind(datprep$gt1*a1, datprep$gt1*a2, datprep$gt2*a1, datprep$gt2*a2)

## Prepping dalpha
dalpha <- cbind(datprep$dgt1*a1, datprep$dgt1*a2, datprep$dgt2*a1, datprep$dgt2*a2)

## Prepping gamma
gamma <- cbind(datprep$x.1%*%g1, datprep$x.1%*%g2, datprep$x.2%*%g1, datprep$x.2%*%g2)

##-----------------------------------------------------------------------
## Estimation of loglikelihood
##-----------------------------------------------------------------------
ll <- loglik(sigma, ncauses, causes, alpha, dalpha, beta, gamma, eb0, nq, stepsize, iter)

if (grad==FALSE){
    llw <- w*ll
    k <- -sum(llw)
    k
    return(k)
}
if(grad==TRUE){
    dg <- as.data.frame(cbind(as.numeric(ll),as.numeric(ID)))
    pc1 <- melt(tapply(dg$V1, dg$V2, sum))[,2]
    return(pc1)
}
}
