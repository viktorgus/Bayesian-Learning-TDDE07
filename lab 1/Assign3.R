# uppgift a

y = c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23, 2.07,2.02)

vonMise = function(mu,y,k){
  return ( prod(exp(k*cos(y-mu))/(2*pi*besselI(k,0)))*dexp(k) )
}


n = 1000
mu = 2.39
k = seq(0,10,0.001)
posterior = k

for ( i in 1:length(k)){
posterior[i] = vonMise(mu,y,k[i])
}

posterior = posterior/sum(posterior)
plot(x=k,y=posterior)


# uppgift b
k[which.max(posterior)]
