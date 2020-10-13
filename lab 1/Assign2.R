library(LaplacesDemon)
library(bayestestR)
library(dplyr)
library(plotly)
library(ggplot2)


ChiInvPlot <- function(n,f){
  xGrid <- seq(0.001, 0.999, by=0.001)
  prior = dinvchisq(xGrid, n, f)
  maxDensity <- max(prior) # Use to make the y-axis high enough
  plot(xGrid, prior, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = "sigma2", 
       ylab = 'Density', main = 'inv-chitwo(n,f) density')
}



## Uppgift a
y =  c(25, 45, 52, 30, 63, 19, 50, 34, 67)
mu = 3.7
n = length(y)

t2 = sum((log(y)-mu)^2)/n
samples = rinvchisq(10000,n,t2)
#hist(samples,main="Sample Distrubution chitwo-inv",xlim=c(0,1),breaks=100)

plot_ly(x = samples,
        type = "histogram",
        histnorm = "probability density") %>% 
  layout(title = "Sample Distrubution chitwo-inv",
         yaxis = list(title = "Density",
                      zeroline = FALSE),
         xaxis = list(title = "Theta",
                      zeroline = FALSE,range=c(0,1)))

ChiInvPlot(n,t2)

## Uppgift b

G = 2*pnorm(sqrt(samples)/sqrt(2))-1
G = sort(G)

#hist(G,main="Sample Distrubution of G",breaks=100)
plot_ly(x = G,
        type = "histogram",
        histnorm = "probability density") %>% 
  layout(title = "Sample Distrubution of G",
         yaxis = list(title = "Density",
                      zeroline = FALSE),
         xaxis = list(title = "G",
                      zeroline = FALSE,range=c(0,1)))


# Generate a normal distribution
gDens = density(G)
#hist(density$y,main="Sample Distrubution of G",breaks=1000,xlim=c(0,1))

# Compute HDI and ETI
ci_eti <- ci(G, method = "ETI",ci=0.9)

# Plot the distribution and add the limits of the two CIs
G %>% 
  estimate_density(extend=TRUE) %>% 
  ggplot(aes(x=x, y=y)) +
  ggtitle("G equal tail credible interval")+
  geom_area(fill="orange") +
  theme_classic() +
  # Quantile in red
  geom_vline(xintercept=ci_eti$CI_low, color="red", size=1) +
  geom_vline(xintercept=ci_eti$CI_high, color="red", size=1)

tail = function(dens,a,b){
  
  return( (sum(dens$y[0:a])+sum(dens$y[b:length(dens$y)]))/sum(dens$y) )
}

leftTail = function(dens,a){
  return ( (sum(dens$y[0:a]))/sum(dens$y) )
}


rightTail = function(dens,b){
  return ( (sum(dens$y[b:length(dens$y)]))/sum(dens$y) )
}

getHDI = function(dens,conf){
  a = which.max(dens$y)
  b = which.max(dens$y)
  
  while( tail(dens,a,b)>conf ){
    if( dens$y[a-1] >= dens$y[b+1] ){
      a = a - 1
    }else{
      b = b + 1
    }
    
  }
  HDI = c(dens$x[a],dens$x[b],a,b)
  return(HDI)
  
}

getCDI = function(dens,conf){
  a = which.max(dens$y)
  b = which.max(dens$y)
  
  while(leftTail(dens,a) > conf/2){
    a = a - 1
  }
  while(rightTail(dens,b) > conf/2){
    b=b+1
  }
  CDI = c(dens$x[a],dens$x[b],a,b)
  return(CDI)
  
}

plotG = function(a,b,title){
  G %>% 
    estimate_density(extend=TRUE) %>% 
    ggplot(aes(x=x, y=y)) +
    ggtitle(title)+
    geom_area(fill="orange") +
    theme_classic() +
    # Quantile in red
    geom_vline(xintercept=a, color="red", size=1) +
    geom_vline(xintercept=b, color="red", size=1)
  
}

HDI = getHDI(gDens,0.1)
CDI = getCDI(gDens,0.1)

plotG(HDI[1],HDI[2],"HDI")
plotG(CDI[1],CDI[2],"CDI")
