library(plotly)

BetaPlot <- function(a,b){
  xGrid <- seq(0.001, 0.999, by=0.001)
  prior = dbeta(xGrid, a, b)
  maxDensity <- max(prior) # Use to make the y-axis high enough
  plot(xGrid, prior, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = "theta", 
       ylab = 'Density', main = 'Beta(7,17) density')
}

s=5;
f=15;
priorA = priorB = 2;
a = s+priorA;
b = f + priorB;
realMean = (a)/(a+b)
realVariance = a*b/((a+b)^2*(a+b+1))

## Uppgift a

samples=10

posterior = rbeta(samples,a,b)
sampleVar = var(posterior)
sampleMean = mean(posterior)
hist(posterior)

smeans = c()
svars = c()

i=2;
while (abs(sampleMean - realMean) > 0.00001 || abs(sampleVar-realVariance) > 0.00001 ){
  i=i+1;
  samples = i*10
  sampleVar = var(posterior)
  sampleMean = mean(posterior)
  posterior = rbeta(samples,a,b)
  smeans = c(smeans,sampleMean)
  svars = c(svars,sampleVar)
}

plot(smeans,col="red",main="Sample Mean vs Real",sub="Red: Sampled, Blue: Real",xlab="samplesammount/10",ylab="Mean")
abline(h=realMean,col="blue")

plot(svars,col="red",main="Sample Sdev vs Real",sub="Red: Sampled, Blue: Real",xlab="samplesammount/10",ylab="sdev")
abline(h=realVariance,col="blue")

samples
#hist(posterior,breaks=100,xlim=c(0,1),main="Sample Distrubution")

plot_ly(x = posterior,
              type = "histogram",
              histnorm = "probability density") %>% 
layout(title = "Sample Beta Distrubution",
       yaxis = list(title = "Density",
                    zeroline = FALSE),
       xaxis = list(title = "Theta",
                    zeroline = FALSE,range=c(0,1)))


BetaPlot(a,b)

## Uppgift b
nDraws = 10000
posterior = rbeta(nDraws,a,b)

sampleProb = length(which(posterior>0.3))/length(posterior)
q = (0.3-realMean)/realVariance

realProb = pbeta(0.3,a,b,lower.tail = FALSE)

realProb
sampleProb

## Uppgift c
posterior = rbeta(nDraws,a,b)
logOdds = log(posterior/(1-posterior))
plot_ly(x = logOdds,
        type = "histogram",
        histnorm = "probability density") %>% 
  layout(title = "Sample Logodds Distrubution",
         yaxis = list(title = "Density",
                      zeroline = FALSE),
         xaxis = list(title = "Theta",
                      zeroline = FALSE))

density(logOdds)
