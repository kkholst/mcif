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

#-----------------------------------------------------------------------
# Data
#-----------------------------------------------------------------------
#data <- unique(refres3)
tau <- 55

# Checking
#head(data)

# Censoring at tau
data$out1 <- ifelse(data$time1 >= tau, tau, data$time1)
data$status1 <- ifelse(data$time1!=data$out1, 0, data$event1)
data$out2 <- ifelse(data$time2 >= tau, tau, data$time2)
data$status2 <- ifelse(data$time2!=data$out2, 0, data$event2)

#-----------------------------------------------------------------------
# Optimisation + miscellaneaous
#-----------------------------------------------------------------------
par6 <- c(2,3.5,-0.5,0.5,-0.5,0.5)
par10 <- c(2,3.5,-1.5,1,2,3,-1.5,1,2,3)

f1 <- function(p){cif.margloglik(p,data,time="out",status="status")}
f2 <- function(p){cif.margloglik(p,data,time="out",status="status",lin=FALSE)}

f1(par6)
f2(par10)

# Optimisation
op1 <- optimx(par6,f1,method=c("nlminb"),control=list(trace=1,kkt=FALSE),itnmax=500)
oppar1 <- coef(op1)
op2 <- optimx(par10,f2,method=c("nlminb"),control=list(trace=1,kkt=FALSE),itnmax=500)
oppar2 <- coef(op2)

# Checking
jacobian(f1,oppar1)
jacobian(f2,oppar2)

# Making plots
t.cif <- c(data[,"out1"],data[,"out2"])
y.cif <- c(data[,"status1"],data[,"status2"])
cif <- cuminc(t.cif,y.cif)
plot(cif,ylim=c(0,0.25))

# Time points
t <- seq(0.1,54.9,0.4)

# Parameters
b.1 <- oppar1[1:2]
a1.1 <- oppar1[3:4]
a2.1 <- oppar1[5:6]

b.2 <- oppar2[1:2]
a1.2 <- oppar2[3:6]
a2.2 <- oppar2[7:10]

# Pi1 and pi2
num1.1 <- exp(b.1[1])
num2.1 <- exp(b.1[2])
denom.1 <- 1+exp(b.1[1])+exp(b.1[2])

num1.2 <- exp(b.2[1])
num2.2 <- exp(b.2[2])
denom.2 <- 1+exp(b.2[1])+exp(b.2[2])

pi1.1 <- num1.1/denom.1
pi2.1 <- num2.1/denom.1

pi1.2 <- num1.2/denom.2
pi2.2 <- num2.2/denom.2

# Alpha and dalpha
alp1.1 <- a_lin(g(t),a1.1)
alp2.1 <- a_lin(g(t),a2.1)

alp1.2 <- a_spl(g(t),a1.2)
alp2.2 <- a_spl(g(t),a2.2)

# CIF
ymod1.1 <- pi1.1*pnorm(alp1.1)
ymod2.1 <- pi2.1*pnorm(alp2.1)

ymod1.2 <- pi1.2*pnorm(alp1.2)
ymod2.2 <- pi2.2*pnorm(alp2.2)

# Adding lines to plot
lines(t,ymod1.1,col=2)
lines(t,ymod2.1,col=2)

lines(t,ymod1.2,col=3)
lines(t,ymod2.2,col=3)

#-----------------------------------------------------------------------------
# Extra check
tt <- c(0.1,25,54)

lin <- function(x){a_lin(g(x),a1.1)}
dlin <-  function(x){da_lin(dg(x),a1.1)}
lin(tt)
dlin(tt)
colSums(jacobian(lin,tt))

spline <- function(x){a_spl(g(x),a1.2)}
dspline <- function(x){da_spl(g(x),dg(x),a1.2)}
spline(tt)
dspline(tt)
colSums(jacobian(spline,tt))
