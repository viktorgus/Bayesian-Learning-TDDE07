library(LaplacesDemon)
library(mvtnorm)
library(bayestestR)
library(ggplot2)

data = read.table("TempLinkoping.txt",header=TRUE)

X = cbind(rep(1,nrow(data)), data$time, data$time^2)
Y = data$temp

colnames(X) = c("1","time","time2")


## Uppgift 1 
#omega0 = diag(x=0.01,3)
omega0 = diag(x=0.5,3)

v0 = 4
mu0 = c(-10,140,-140)
#mu0 = c(-10,100,-100)
sigma02 = 1

  prior.sigma2 = rinvchisq(1,v0,sigma02)
  prior.b = rnorm(3,mu0,prior.sigma2*solve(omega0))
  yPrior = rnorm(nrow(X),X%*%prior.b,prior.sigma2 * diag(3))
  

  plot(round(365*X[,2]),yPrior,xlim =c(0,365),ylim=c(min(yPrior),max(yPrior)),xlab="Day",ylab="Temperature")
  





## Uppgift 2 

bhat = (solve(t(X) %*% X)) %*% t(X) %*% Y

s2 = 1/(nrow(X)-ncol(X))*t(Y-X %*% bhat)%*%(Y-X %*% bhat)
vn = v0 + nrow(X)

muN = solve(t(X)%*%X + omega0)%*%(t(X)%*%X%*%bhat + omega0%*%mu0)
omegaN = t(X)%*%X + omega0

sigma2N = (v0*sigma02 + (t(Y)%*%Y + t(mu0)%*%omega0%*%mu0 - t(muN)%*%omegaN%*%muN))/vn


## Compute list for histograms of parameters
draws = 500
list.sig = rep(0,draws)
list.b = matrix(nrow=draws,ncol=3)
list.y = matrix(nrow=nrow(X),ncol=draws)
list.xMax = rep(0,draws)


colnames(list.b) = c("b0","b1","b2")
for (i in 1:draws){
  post.sigma2 = rinvchisq(1,vn,sigma2N)
  post.b = rmvn(n=1,muN[,1],round(post.sigma2*solve(omegaN),3))
  list.sig[i]=post.sigma2
  list.b[i,]=post.b
  list.y[,i] = X%*%t(post.b)
  
  ##Upgift c
  list.xMax[i] = X[which.max(list.y[,i]),2]
}

#uppgift c
hist(round(365*list.xMax),xlab="day", ylab="density",main="day with max temperature")


hist(list.b[,1],xlab = "b0",main=sprintf("Histogram of b0 from %i draws",draws))
hist(list.b[,2],xlab = "b1",main=sprintf("Histogram of b1 from %i draws",draws))
hist(list.b[,3],xlab = "b2",main=sprintf("Histogram of b2 from %i draws",draws))
hist(list.sig,xlab = "sigma-squared",main=sprintf("Histogram of sigma from %i draws",draws))



## Compute CI interval for all X values
upper = rep(0,nrow(X))
lower = rep(0,nrow(X)) 
for (i in 1:nrow(list.y)){
interval <- ci(list.y[i,], method = "ETI",ci=0.95)
upper[i] = interval$CI_low
lower[i] = interval$CI_high
}

median.b =c(median(list.b[,1]),median(list.b[,2]),median(list.b[,3]))

median.y = X%*%(median.b)
df <- data.frame(
  median = median.y,
  x = round(365*X[,2]),
  upper=upper,
  lower=lower,
  y=Y
)
#sum(apply(list.y,1,median) - median.y)/nrow(median.y)

ggplot(df,aes(x=x,y=median)) + geom_line() + ylim(-15,25) + 
  labs(x="days",y="prediction") +
  geom_line(aes(x=x,y=upper),color="red") + 
  geom_line(aes(x=x,y=lower),color="red") +
  geom_point(aes(x=x,y=y))

