library("rstan")

dat = read.csv("campy.dat")

data = list(T=length(dat$c), ct=dat$c) 
burnin = 1000 
niter = 2000 

StanModel = '
data {
  int<lower=0> T;
  int ct[T];
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
  vector[T] xt;
}
model {
  sigma ~ scaled_inv_chi_square(T,0.03);
  phi ~ uniform(-1,1);
  for (n in 2:T)
    xt[n] ~ normal(mu + phi * (xt[n-1] - mu), sqrt(sigma));
  for (n in 1:T)
    ct[n] ~ poisson(exp(xt[n]));
}'



fit = stan(model_code=StanModel,data=data, warmup=burnin,iter=niter)
print(fit,digits_summary=3)

sum = summary(fit)$summary
theta =  exp(sum[5:nrow(sum)-1,])

df = data.frame(
  y = dat$c,
  x = 1:nrow(dat),
  mean = theta[,1],
  lower = theta[,4],
  upper = theta[,8]
  
)

ggplot(df, aes(x=x,y=y),) + geom_point() + geom_line(aes(y=mean), color="red") + geom_line(aes(y=upper), color ="orange") + geom_line(aes(y=lower),color="orange") + ggtitle("Intensity",subtitle="Red: mean, Orange: 95% CDI") + xlab("t")
  
                                                                                          