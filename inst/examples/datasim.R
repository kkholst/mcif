################################################################################
# FUNCTION FOR CREATING DATA SET WITH TWO CAUSES OF FAILURE
################################################################################
# Loading libraries
library(mets); library(MASS); library(miscF); library(cmprsk)

# Sourcing help function
source("helpfunctions_datasim.R")

# Simulated pairs of 3
n <- 30000

# Setting seed
set.seed(1449)

# Save data as
file.out <- 'data.Rdata'

# Censoring probability (larger than actual proportion censored)
cp <- 0.80

# Delta (max time)
delta <- 80

################################################################################
# Choosing parameters

# Betas
b1 <- -1.9
b2 <- -0.2
b <- c(b1,b2)

# Alphas
a10 <- 1 # in reality gamma1
a11 <- 3 # must be positive
a1 <- c(a10,a11)
a20 <- 1.2 # in reality gamma2
a21 <- 2 # must be positive
a2 <- c(a20,a21)

# Variance
var.e1 <- 1
var.e2 <- 1
var.u1 <- 1
var.u2 <- 1

# Covariance
cov.e1e2 <- 0.4
cov.u1u2 <- 0.4

cov.e1u1 <- 0.5
cov.e1u2 <- 0.4

cov.e2u1 <- 0.4
cov.e2u2 <- 0.3

################################################################################
# Mu
mu4 <- rep(0,4)

# Variance-covariance matrices
vcv <- matrix(c(var.e1, cov.e1e2, cov.e1u1, cov.e1u2,
                cov.e1e2, var.e2, cov.e2u1, cov.e2u2,
                cov.e1u1, cov.e2u1, var.u1, cov.u1u2,
                cov.e1u2, cov.e2u2, cov.u1u2, var.u2),
              nrow=4, ncol=4)

# Checking vcv matrix
eigen(vcv)

# Random effects
RE <- as.data.frame(mvrnorm(n, mu=mu4, Sigma=vcv))
RE <- cbind(RE,RE,RE)

# Assigning family IDs
RE$ID <- 1:n

# Column names
colnames(RE) <- c("r1.effectf1","r2.effectf1","r3.effectf1","r4.effectf1","r1.effectf2","r2.effectf2","r3.effectf2","r4.effectf2","r1.effectf3","r2.effectf3","r3.effectf3","r4.effectf3","ID")

# Rearranging data
ref <- as.data.frame(fast.reshape(RE, var=c("r1.effect","r2.effect","r3.effect","r4.effect"), idname="ID", numname="mem"))

################################################################################
# Pis
ref$pi1 <- pi.1(ref[,c("r3.effect","r4.effect")], b)
ref$pi2 <- pi.2(ref[,c("r3.effect","r4.effect")], b)
ref$pi3 <- 1-ref$pi1-ref$pi2
head(ref)

# Event type
ref$event <- rMultinom(p=cbind(ref$pi1,ref$pi2,ref$pi3))

# Event time
ref$time <- timenew(ref[,c("r1.effect","r2.effect","event")], a1, a2, delta)

# Removing unnecessary variables and renaming columns
refnew <- ref[,c("ID","time","event")]
colnames(refnew) <- c("ID","time","event")

# Renaming group 3
refnew$event <- ifelse(refnew$event==3, 0, refnew$event)

# Changing event at delta
refnew$event <- ifelse(refnew$time>=delta,0,refnew$event)
refnew$time <- ifelse(refnew$time>=delta,delta,refnew$time)

# Censoring
refnew$cens <- rbinom(nrow(refnew),1,prob=cp)

# Censoring times
refnew$censtime <- runif(nrow(refnew),1/365.25,delta+30)
refnew$newtime <- (refnew$cens==1 & (refnew$censtime < refnew$time))*refnew$censtime + (refnew$cens==0 | (refnew$censtime >= refnew$time))*refnew$time
refnew$newevent <- (refnew$cens==1 & (refnew$censtime < refnew$time))*0 + (refnew$cens==0 | (refnew$censtime >= refnew$time))*refnew$event

# Proportion censored
propcens <- sum((refnew$cens==1 & (refnew$censtime < refnew$time)))/dim(refnew)[1]*100
print(paste("Proportion censored:", round(propcens,2), "%"))

# Family pairs
pair <- familycluster.index(clusters=refnew$ID)

# Generating pairs
refres <- refnew[t(pair$pair),]

# Pair ID
refres$pairID <- rep(1:(nrow(refres)/2),each=2)

# Reshaping data
refres2 <- fast.reshape(refres,id=refres$pairID)

# Removing unnecessary columns
refres3 <- refres2[,c("ID1","pairID1","newtime1","newevent1","newtime2","newevent2")]

# Renaming
colnames(refres3) <- c("ID","pairID","time1","event1","time2","event2")

# Regrouping if time1 or time2 is equal to delta
refres3$event1 <- ifelse(refres3$time1==delta, -1, refres3$event1)
refres3$event2 <- ifelse(refres3$time2==delta, -1, refres3$event2)

# Saving
data <- refres3
save(data,file=file.out)
