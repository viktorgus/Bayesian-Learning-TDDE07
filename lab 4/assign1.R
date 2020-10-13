library("rstan")


# 1a) 
Ar = function(phi,mu, sigma2,T){
  # x is vector of x values to be returned from AR process
  X = rep(0,T)
  err = rnorm(T,0,sigma2)
  X[1] = mu + err[1]
  for(i in 2:T){
    X[i] = mu +phi*(X[i-1]-mu)+err[i]
  }
  
  return(X)
  
}


xT = Ar(0.95,10,2,200)
plot( xT , xlab="Iteration", main="Ar-process for phi = 0.8")


# 1b) 

st = function(x){
  data = list(T=200, xT=x) 
  burnin = 1000 
  niter = 2000 
  
  StanModel = '
data {
  int<lower=0> T;
  vector[T] xT;
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
}
model {
  phi ~ uniform(-1,1);
  for (n in 2:T)
    xT[n] ~ normal(mu + phi * (xT[n-1]-mu), sigma);
}'

fit = stan(model_code=StanModel,data=data, warmup=burnin,iter=niter,chains=4)
  return(fit)
  
}




fit1 = st(Ar(0.95,10,2,200))
fit2 = st(Ar(0.3,10,2,200))

print(fit1,digits_summary=3)
print(fit2,digits_summary=3)

# Extract posterior samples 
 postDraws <- extract(fit2)
plot(x=postDraws$mu,y=postDraws$phi,xlab = "Mu", ylab="phi")
 
traceplot(fit2)
# Bivariate posterior plots 
pairs(fit2)
