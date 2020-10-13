library("mvtnorm")
library("bayestestR")
library("ggplot2")
library("combinat")


#Assignment 1

# Loading data from file
Data<-read.table("WomenWork.dat",header=TRUE)  # Spam data from Hastie et al.
y <- as.vector(Data[,1]); # Data from the read.table function is a data frame. Let's convert y and X to vector and matrix.
X <- as.matrix(Data[,2:9]);
names(X) = names(data)[2:9]

nPara = dim(X)[2]

#Setup logistic prior
tau <- 10;
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara);

logpost = function(beta,y,X,mu,Sigma){
  nPara = length(beta)
  pred = X%*%beta
  
  likelihood = sum(y*pred - log(1+exp(pred)))
  if (abs(likelihood) == Inf) likelihood = -20000
  priorlike = dmvnorm(beta, matrix(0,nPara,1), Sigma, log=TRUE)
  return(likelihood+priorlike)
}


initVal <- as.vector(rep(0,nPara)); 
OptimResults<-optim(initVal,logpost,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

betaMode = OptimResults$par
Jinv = -solve(OptimResults$hessian)
names(betaMode) = names(X)

samples = 10000
samples = rmvnorm(samples,betaMode, Jinv)
betaList = samples[,7] 

betaList = sort(betaList)
dens = density(betaList)
plot(dens,main="density of b7 for NsmallChild" ,xlab="b7")
ci_eti <- ci(betaList, method = "ETI",ci=0.95)


ggplot(data.frame(x=dens$x,y=dens$y),aes(x=x,y=y)) + geom_line() + geom_area(fill="blue") 
+ geom_vline(xintercept=ci_eti$CI_low, color="red", size=1) +
  geom_vline(xintercept=ci_eti$CI_high, color="red", size=1)


## Assignment b
x = c(1,10,8,10,1,40,1,1)
preds = t(x%*%t(samples))[,1]
preds = 1/(1+exp(preds))
plot(density(preds),main="",sub="p(y=0|x)")
##Chance of working
length(which(preds<0.5))/length(preds)


## assignemnt c
preds = t(x%*%t(samples))[,1]
p = exp(preds)/(1+exp(preds))

sampleProbs = matrix(nrow=length(p),ncol=10)
for (i in 1:length(p)){
  
  for(j in 1:10){
    sampleProbs[i,j] = dbinom(j, 10, p[i])
  }
}
sumrow = colSums(sampleProbs)/nrow(sampleProbs)
sumrow = sumrow/sum(sumrow)
qplot(x=c(1:10),y=sumrow,ylab="probability", main="PDF for number of working women",geom=c("line","point"),xlab="Number of working out of ten identical women")
