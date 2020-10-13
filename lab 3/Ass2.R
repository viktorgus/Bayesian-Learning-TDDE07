library("mvtnorm")
library("MASS")
library("plotly")

##  a)
data.ebay = read.table("eBayNumberOfBidderData.dat", header=TRUE)
glm.fit = glm(nBids ~ .-Const, data=data.ebay ,family="poisson")
summary(glm.fit)

X = as.matrix(data.ebay[,2:10])
y = data.ebay$nBids
 
# b )
logPosterior = function(beta, y, X){
  
  #Given from the log of the likelihood function for poisson with lambda = exp(x*B)
  likelihood = sum(y * X %*% beta - exp(X %*% beta))
  
  if (abs(likelihood) == Inf) likelihood = -20000
  priorlike = dmvnorm(x = beta, mean = rep(0,9), sigma = 100*solve(t(X) %*% X), log=TRUE)
  return(likelihood+priorlike)
}

initVal <- as.vector(rep(0,9)); 

# optimize the logposterior to find the posterior mean and the associated hessian
OptimResults<-optim(initVal,logPosterior,gr=NULL,y,X,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)


## c)
metropolis = function(theta, c, iterations, logPosterior, ...){
  
  # Run the optim-algo to find the hessian
  initVal <- as.vector(rep(0,9)); 
  OptimResults<-optim(initVal,logPosterior,gr=NULL,y,X,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)
  sigma = solve(-OptimResults$hessian)
  
  #thetaList = rep(0,iterations)
  thetaM = matrix(nrow = iterations, ncol = length(theta))
  for(i in 1:iterations){
    thetaM[i,] = theta
    proposal = mvrnorm(1,theta,sigma*c)
    
    acceptenceProb = min(1, exp(logPosterior(proposal, ...) - logPosterior(theta, ...)))
    #print(sprintf("acceptence prob: %f iteration: %i",acceptenceProb, i))
    if( runif(1)>1-acceptenceProb ){
      theta = proposal
    }
  }
  return(thetaM)
}


plotBeta = function(beta,matrix){
  
  betas = matrix[-(1:350),beta]
  plot_ly(x = betas,
          type = "histogram",
          histnorm = "probability density") %>% 
    layout(title = sprintf("Sampled Beta%i Distrubution", beta),
           yaxis = list(title = "Density",
                        zeroline = FALSE),
           xaxis = list(title = sprintf("Beta%i",beta),
                        zeroline = FALSE,range=c(min(betas),max(betas))))

}

#function to get one poisson draw for each draw of beta for a given x value
getPossValues = function(x,betas){
  lambdas = exp(x %*% t(betas))
  valList = rep(0, length(lambdas))
  for (i in 1:length(lambdas)){
    valList[i] = rpois(1,lambdas[i])
  }
  return(valList)
}

# run the metropolis hastings algorithm with the zero-vector as initial start
initVal <- as.vector(rep(0,9)); 
matrix = metropolis(initVal,1,4000,logPosterior, y, X)

# plot one distrubution and one convergence
plotBeta(1,matrix)
plot(matrix[,9],xlab = "iteration", ylab = "beta9")

## d)

# x*B, remove burn in for distrubution of betas
x = c(1,1,1,1,0,0,0,1,0.5)


# remove 1:300 values of betadraws from M.H to remove burn-in
# For each betavalue given from M.H: Get a draw from poisson of exp(x*B) as lambda (from the given
# poisson regression expression in a )
values = getPossValues(x,matrix[-c(1:300)])

hist(values,xlab = "nBids", main ="Distrubution of nBids")
length(which(values==0))/length(values)
