#-----------------------------------------------------------------------
# Time transforming function g(t)
#-----------------------------------------------------------------------
# Transformation g
g <- function(x,delta){
     atanh((x-delta/2)/(delta/2))
}

# Derivative of g
dg <- function(x,delta){
     (1/2)*delta/(x*(delta-x))
}

#-----------------------------------------------------------------------
# Data preparing function
#-----------------------------------------------------------------------
data.prep <- function(data, time, status, cova=NULL){

#-----------------------------------------------------------------------
# The Ys
#-----------------------------------------------------------------------
y <- cbind(data[, paste(status, 1, sep="")], data[, paste(status, 2, sep="")])

#-----------------------------------------------------------------------
# Xs
#-----------------------------------------------------------------------
n <- nrow(data)

if (!is.null(cova)){
    x1 <- paste(cova, "1", sep="")
    x2 <- paste(cova, "2", sep="")
}
if (is.null(cova)){
    x1 <- cova
    x2 <- cova
}
x.1 <- as.matrix(cbind(rep(1,n),data[,x1]))
x.2 <- as.matrix(cbind(rep(1,n),data[,x2]))

#-----------------------------------------------------------------------
# Basis spline functions and derivatives
#-----------------------------------------------------------------------
# Time points
t1 <- data[, paste(time, 1, sep="")]
t2 <- data[, paste(time, 2, sep="")]

# Max. time
delta <- max(t1,t2)+1e-10

# Transforming t1 and t2
gt1 <- g(t1,delta)
gt2 <- g(t2,delta)

# Knots
#knots <- c(-50, -1, -0.5, 0, 0.5, 1, 50)
knots <- c(-50, -1, 0, 1, 50)

# Spline functions
spl <- bsplineS(c(gt1, gt2), breaks=knots, norder=3, nderiv=0) # 2nd degree

# Detivative of spline functions wrt. t
dspl <- bsplineS(c(gt1, gt2), breaks=knots, norder=3, nderiv=1) # 2nd degree, 1st derivative

#-----------------------------------------------------------------------
# Derivative of g(t)
#-----------------------------------------------------------------------
dgt <- dg(c(t1,t2),delta)

#-----------------------------------------------------------------------
# Return
#-----------------------------------------------------------------------
res <- list("y"=y, "x.1"=x.1, "x.2"=x.2, "spl"=spl, "dspl"=dspl, "delta"=delta, "dgt"=dgt)
return(res)
}
