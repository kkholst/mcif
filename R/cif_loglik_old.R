cif.loglik.old <- function(par, y, x.1, x.2, spl, dspl, dgt, eb0, ID, nq, stepsize, grad=FALSE,  ...){
#-----------------------------------------------------------------------
# Parameters
#-----------------------------------------------------------------------
# Number of pairs
n <- nrow(y)

# Parameters
b1 <- par[1]
b2 <- par[2]

#gam1 <- par[3]
#gam2 <- par[4]

a1.1 <- par[3]
a1.2 <- par[4]
a1.3 <- par[5]
a1.4 <- par[6]
a1.5 <- par[7]
a1.6 <- par[8]
#a1.7 <- par[11]
#a1.8 <- par[12]

a2.1 <- par[9]
a2.2 <- par[10]
a2.3 <- par[11]
a2.4 <- par[12]
a2.5 <- par[13]
a2.6 <- par[14]
#a2.7 <- par[19]
#a2.8 <- par[20]

s1 <- par[15]
s2 <- par[16]
s3 <- par[17]
s4 <- par[18]
s5 <- par[19]

s6 <- par[20]
s7 <- par[21]
s8 <- par[22]
s9 <- par[23]
s10 <- par[24]
s11 <- par[25]
s12 <- par[26]

#-----------------------------------------------------------------------
# Transforming variables for splines and variance-covariance matrix
#-----------------------------------------------------------------------
# Cholesky (upper triangle)
L.chol <- matrix(c(exp(s1),0,0,0,
                     s5,exp(s2),0,0,
                     s6,s8,exp(s3),0,
                     s7,s9,s10,exp(s4)),
                   nrow=4, ncol=4)

vcv <- t(L.chol)%*%L.chol

sigma88 <- matrix(c(exp(s11),0,0,0,0,0,0,0,
                    0,exp(s11),0,0,0,0,0,0,
                    0,0,exp(s12),0,0,0,0,0,
                    0,0,0,exp(s12),0,0,0,0,
                    0,0,0,0,vcv[1,1],vcv[2,1],vcv[3,1],vcv[4,1],
                    0,0,0,0,vcv[1,2],vcv[2,2],vcv[3,2],vcv[4,2],
                    0,0,0,0,vcv[1,3],vcv[2,3],vcv[3,3],vcv[4,3],
                    0,0,0,0,vcv[1,4],vcv[2,4],vcv[3,4],vcv[4,4]),
                  nrow=8,ncol=8)

A <- matrix(c(1,0,0,0,0,0,
              0,1,0,0,0,0,
              0,0,1,0,0,0,
              0,0,0,1,0,0,
              1,1,0,0,0,0,
              0,0,1,1,0,0,
              0,0,0,0,1,0,
              0,0,0,0,0,1),
            nrow=6, ncol=8)

# Sigma
sigma <- A%*%sigma88%*%t(A)

# Spline parameters - have to be monotonically increasing
a1 <- vector()
a1[1] <- a1.1
a1[2] <- a1[1]+exp(a1.2)
a1[3] <- a1[2]+exp(a1.3)
a1[4] <- a1[3]+exp(a1.4)
a1[5] <- a1[4]+exp(a1.5)
a1[6] <- a1[5]+exp(a1.6)
#a1[7] <- a1[6]+exp(a1.7)
#a1[8] <- a1[7]+exp(a1.8)

a2 <- vector()
a2[1] <- a2.1
a2[2] <- a2[1]+exp(a2.2)
a2[3] <- a2[2]+exp(a2.3)
a2[4] <- a2[3]+exp(a2.4)
a2[5] <- a2[4]+exp(a2.5)
a2[6] <- a2[5]+exp(a2.6)
#a2[7] <- a2[6]+exp(a2.7)
#a2[8] <- a2[7]+exp(a2.8)

#-----------------------------------------------------------------------
# X^T*beta1 and X^T*beta2 for individual 1 and 2
#-----------------------------------------------------------------------
b1_1 <- x.1%*%b1
b1_2 <- x.2%*%b1
b2_1 <- x.1%*%b2
b2_2 <- x.2%*%b2

b <- cbind(b1_1,b1_2,b2_1,b2_2)

#-----------------------------------------------------------------------
# Spline functions and their 1st derivatives for individual 1 and 2
#-----------------------------------------------------------------------
# Alphas
spl.1 <- spl%*%a1
spl.2 <- spl%*%a2

# Dalphas
dspl.1 <- (dspl%*%a1)*dgt
dspl.2 <- (dspl%*%a2)*dgt

#-----------------------------------------------------------------------
# X^T*gamma1 and X^T*gamma2 for individual 1 and 2
#-----------------------------------------------------------------------
#gam1_1 <- x.1%*%gam1
#gam1_2 <- x.2%*%gam1
#gam2_1 <- x.1%*%gam2
#gam2_2 <- x.2%*%gam2

#-----------------------------------------------------------------------
# Sum of spline+X^T*gamma1 and spline+X^T*gamma2 for individual 1 and 2
#-----------------------------------------------------------------------
# Sum
alpha1 <- spl.1 #+ c(gam1_1, gam1_2)
alpha2 <- spl.2 #+ c(gam2_1, gam2_2)

# Splitting
alpha <- cbind(alpha1[1:n], alpha1[(n+1):(n*2)], alpha2[1:n], alpha2[(n+1):(n*2)])
dalpha <- cbind(dspl.1[1:n], dspl.1[(n+1):(n*2)], dspl.2[1:n], dspl.2[(n+1):(n*2)])

#-----------------------------------------------------------------------
# Loglikelihood
#-----------------------------------------------------------------------
print(par)
ll <- loglik(y, b, sigma, alpha, dalpha, eb0, nq, stepsize)
if (grad==FALSE){
k <- -sum(ll)
k
return(k)
}
if(grad==TRUE){
    dg <- cbind(ll,ID)
    pc1 <- melt(tapply(dg$ll, dg$ID, sum))[,2]
    return(pc1)
}
}
