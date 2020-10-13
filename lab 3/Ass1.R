library(LaplacesDemon)
library(ggplot2)


dat = read.csv2("rainfall.dat")

X = dat$X136

mu0 = 0
sigma0 = 5
tau0 = 5
v0 = 20
n = length(X)[1]
xmean = sum(X)/n

mu = rnorm(1,mu0,tau0^2)
sigma = sqrt(rinvchisq(1,v0,sigma0^2)) 


### Part i)
draws = 700

sigmas = rep(0,draws)
mus = rep(0,draws)
meanMus = rep(0,draws)
meanSigmas = rep(0,draws)

for (i in 1:draws){
  tauN = sqrt((n/sigma^2+1/tau0^2)^-1)
  
  w = (n/sigma^2)/(n/sigma^2+1/tau0^2)
  muN = w*xmean +(1-w)*mu0
  
  mu = rnorm(1,muN,tauN^2)
  sigma = sqrt(rinvchisq(  1, 
                           v0+n,(v0*sigma0^2+sum((X-mu)^2))/(n+v0) ))
  mus[i] = mu
  sigmas[i] = sigma
  meanMus[i] = sum(mus[1:i])/i
  meanSigmas[i] = sum(sigmas[1:i])/i
}

plot(meanMus, xlab = "Iterations", ylab = "Mean of sampled Mus")
  plot(meanSigmas, xlab ="Iterations", ylab ="Mean of sampled Sigmas")

### Part ii)
plot_dat = data.frame(
  sigma = sigmas,
  mu = mus
)


ggplot(plot_dat) + geom_path( aes(x=sigma,y=mus))

#1 c) , need to run NormalMixtureGibbs first therefor, commented out

#hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Mixed Dens, vs Single Dens", ylim = ylim)
#lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
#dens = dnorm(xGrid,meanMus[draws],sd = meanSigmas[draws])
#lines(xGrid, dens, type = "l", lty = 2, lwd = 3, col = 'blue')
#legend("topright", box.lty = 1, legend = c("Data histogram",'Normal Density a)', 'Mixture b)'), 
#       col = c("black",'blue', 'red'), lwd = 2)

