# Loading library
library(optimx); library(cmprsk);library(fda);library(numDeriv)

#-----------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------
g <- function(t){
    gt <- atanh((t-0.5*tau)/(0.5*tau))
    return(gt)
}

dg <- function(t){
    dgt <- (1/2)*tau/(t*(tau-t))
    return(dgt)
}

loglik <- function(x,pi1,pi2,tau){
    if (x[1]==1){
        l <- log(pi1*x[4]*dnorm(x[2]))
    }
    if (x[1]==2){
        l <- log(pi2*x[5]*dnorm(x[3]))
    }
    if (x[1]==0){
        if (x[6]==tau){
            l <- log(1-pi1-pi2)
        }
        else {
            l <- log(1-pi1*pnorm(x[2])-pi2*pnorm(x[3]))
        }
    }
    return(l)
}


a_lin <- function(gt,a){
    ax <- a[1]+exp(a[2])*gt
    return(ax)
}

da_lin <- function(dgt,a){
    dax <- exp(a[2])*dgt
    return(dax)
}

a_spl <- function(gt,a){
    spl <- bsplineS(gt,breaks=c(-10,0,10),norder=3,nderiv=0)#[,2:4]
    at <- vector()
    at[1] <- a[1]
    at[2] <- at[1]+exp(a[2])#a[2]#
    at[3] <- at[2]+exp(a[3])#a[3]#
    at[4] <- at[3]+exp(a[4])#a[4]#
    ax <- spl%*%at
    #ax <- at[1]+spl%*%at[2:4]
    return(ax)
}

da_spl <- function(gt,dgt,a){
    spl <- bsplineS(gt,breaks=c(-10,0,10),norder=3,nderiv=1)#[,2:4]
    at <- vector()
    at[1] <- a[1]
    at[2] <- at[1]+exp(a[2])#a[2]#
    at[3] <- at[2]+exp(a[3])#a[3]#
    at[4] <- at[3]+exp(a[4])#a[4]#
    dax <- spl%*%at*dgt
    #dax <- spl%*%at[2:4]*dgt
    return(dax)
}


#-----------------------------------------------------------------------
# Marginal loglikelihood
#-----------------------------------------------------------------------
cif.margloglik <- function(par, data, time, status, lin=TRUE){
    y <- c(data[,paste(status,1,sep="")],data[,paste(status,2,sep="")])
    t <- c(data[,paste(time,1,sep="")],data[,paste(time,2,sep="")])
    tau <- max(t)
    y.t <- y[order(t)]
    t.t <- t[order(t)]
    n <- nrow(data)*2
    b <- par[1:2]
    num1 <- exp(b[1])
    num2 <- exp(b[2])
    denom <- 1+exp(b[1])+exp(b[2])
    pi1 <- num1/denom
    pi2 <- num2/denom
    if (lin==TRUE){
        gt <- g(t)
        dgt <- dg(t)
        a1 <- par[3:4]
        a2 <- par[5:6]
        alp1 <- a_lin(gt,a1)
        alp2 <- a_lin(gt,a2)
        dalp1 <- da_lin(dgt,a1)
        dalp2 <- da_lin(dgt,a2)
    }
    else {
        a1 <- par[3:6]
        a2 <- par[7:10]
        y1 <- y.t[t.t<tau]
        y2 <- y.t[t.t>=tau]
        t1 <- t.t[t.t<tau]
        t2 <- t.t[t.t>=tau]
        n2 <- length(y2)
        gt <- g(t1)
        dgt <- dg(t1)
        alp1 <- c(a_spl(gt,a1),rep(1,n2))
        alp2 <- c(a_spl(gt,a2),rep(1,n2))
        dalp1 <- c(da_spl(gt,dgt,a1),rep(1,n2))
        dalp2 <- c(da_spl(gt,dgt,a2),rep(1,n2))
    }
    dat <- cbind(y,alp1,alp2,dalp1,dalp2,t)
    ll <- apply(dat, 1, loglik, pi1, pi2, tau)
    nll <- sum(-ll)
    return(nll)
}
