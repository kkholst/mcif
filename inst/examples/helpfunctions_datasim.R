#--------------------------------------------------------------------
# Pi
pi.1 <- function(x, b){
    num <- exp(b[1]+x[,1])
    den <- exp(b[1]+x[,1])+exp(b[2]+x[,2])+1
    pi <- num/den
    return(pi)
}

pi.2 <- function(x, b){
    num <- exp(b[2]+x[,2])
    den <- exp(b[1]+x[,1])+exp(b[2]+x[,2])+1
    pi <- num/den
    return(pi)
}

#--------------------------------------------------------------------
# Function for finding time points
timenew <- function(x, a1, a2, delta){
    u <- runif(dim(x)[1])
    t1 <- 0.5*delta*tanh((qnorm(u)+x[,1]+a1[1])/a1[2])+0.5*delta
    t2 <- 0.5*delta*tanh((qnorm(u)+x[,2]+a2[1])/a2[2])+0.5*delta
    t <- (x[,3]==1)*t1+(x[,3]==2)*t2+(x[,3]==3)*delta
    return(t)
}
