library(Rcpp)
sourceCpp("../src/ex.cpp")

source("P:/PhD/Scripts/M2/M2/Working/Working/step1_prep.R")

#data <- subset(data, kombo==11|kombo==12|kombo==21|kombo==22)

dim(data)

t <- cbind(data$event1,data$event2)
b <- cbind(data$b1_1,data$b1_2,data$b2_1,data$b2_2)
u <- cbind(rep(1,dim(data)[1]),rep(1,dim(data)[1]))
a <- cbind(data$alph1_1,data$alph1_2,data$alph2_1,data$alph2_2)
da <- cbind(data$dalph1_1,data$dalph1_2,data$dalph2_1,data$dalph2_2)

test <- loglik(y=t,b,u,Sigma,alph=a,dalph=da)
test1 <- test[which(data$kombo==11|data$kombo==12|data$kombo==21|data$kombo==22),]

source("P:/PhD/Scripts/M2/Working/Working/logfy.R")

test2 <- logfy(u)[1:26]

sum(test1)
sum(test2)


head(b)
test <- exp(data$b1_1+u[,1])/(1+exp(data$b1_1+u[,1])+exp(data$b2_1+u[,2]))



o1 <- loglik(t,b,u)
o2 <- exp(b[,1]+u[,1])/(1+exp(b[,1]+u[,1])+exp(b[,3]+u[,2]))
o1==o2


loglik(t,b,u)



y <- rbind(c(1,0),c(2,0));
r <- rbind(0.5)
loglik_ex(y,r)
log(pn(y,rbind(0.5,0.5)))

y <- Sigma

conSig(y,rc1=c(2,0),rc2=c(4,5))
conMu(y,rc1=c(2,0),rc2=c(4,5),mu=x)

rc1 <- c(3,1)
rc2 <- c(5,6)

y11 <- y[rc1,rc1]
y12 <- y[rc1,rc2]
y21 <- y[rc2,rc1]
y22 <- y[rc2,rc2]

y11-y12%*%solve(y22)%*%y21
y12%*%solve(y22)%*%x[rc2]

