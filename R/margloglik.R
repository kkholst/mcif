##-----------------------------------------------------------------------
## Nested loglikelihood (submodel)
##-----------------------------------------------------------------------
cif.subloglik <- function(datprep, par){

##-----------------------------------------------------------------------
## Parameters
##-----------------------------------------------------------------------
b1 <- par[1]
b2 <- par[2]

g1 <- par[3]
g2 <- par[4]

a1 <- exp(par[5])
a2 <- exp(par[6])

##-----------------------------------------------------------------------
## The rest
##-----------------------------------------------------------------------
## Causes
causes <- datprep$causes

## Weights
w <- datprep$weights

## Prepping beta
beta <- cbind(datprep$x.1%*%b1, datprep$x.1%*%b2, datprep$x.2%*%b1, datprep$x.2%*%b2)

## Prepping alpha
alpha <- cbind(datprep$gt1*a1, datprep$gt1*a2, datprep$gt2*a1, datprep$gt2*a2)

## Prepping dalpha
dalpha <- cbind(datprep$dgt1*a1, datprep$dgt1*a2, datprep$dgt2*a1, datprep$dgt2*a2)

## Prepping gamma
gamma <- cbind(datprep$x.1%*%g1, datprep$x.1%*%g2, datprep$x.2%*%g1, datprep$x.2%*%g2)

## Estimating the pis
num1 <- exp(b1)
num2 <- exp(b2)
denom <- 1 + num1 + num2
pi1 <- num1/denom
pi2 <- num2/denom

## Alpha - gamma
alpgam <- alpha-gamma

## One individual format
causes_1 <- c(causes[,1],causes[,2])
alpgam_1 <- rbind(alpgam[,1:2],alpgam[,3:4])
dalp_1 <- rbind(dalpha[,1:2],dalpha[,3:4])
dat <- cbind(causes_1, alpgam_1, dalp_1)
ll <- apply(dat, 1, logliksub, pi1, pi2)
nll <- sum(-ll)
return(nll)
}

logliksub <- function(x,pi1,pi2){
    if (x[1]==1){
        l <- log(pi1*x[4]*dnorm(x[2]))
    }
    if (x[1]==2){
        l <- log(pi2*x[5]*dnorm(x[3]))
    }
    if (x[1]==0){
        l <- log(1-pi1*pnorm(x[2])-pi2*pnorm(x[3]))
    }
    if (x[1]<0){
        l <- log(1-pi1-pi2)
    }
    return(l)
}
