
rm(list=ls(all=TRUE))
set.seed(27)

# metropolis-hasting to Uniform distribution.

x = runif(1000,-2,5)

posterior.uniform = function(x,a,b){
  n=length(x) 
  resultado=n*(-log(b-a)) 
  resultado=resultado+(((a^2)/1000)*((b^2)/1000))
  resultado 
}

nmc = 10500 
amc=array(0,c(nmc)) 
bmc=array(0,c(nmc)) 

# initial values
amc[1]=min(x)- 1 # the initial value for amc will always have to be less than the minimum of the generated vector 'x', so put 'min (x) - 2'.
bmc[1]=max(x)+ 1 # the initial value for bmc will always have to be greater than the maximum of the generated vector 'x', so put 'max (x) +2'.

# this condition for the initial kick occurs because of the indicator of the distribution function
# uniform, which is: a <x (1) <x (n) <b.

for (i in 2:nmc){
  Sigmaa=0.0001 # declared value for the variance of amc
  Sigmab=0.0001 # Declared value for the variance of bmc

  aest=rnorm(1,amc[i-1], sqrt(Sigmaa))
  best=rnorm(1,bmc[i-1], sqrt(Sigmab))
  # if command, used to only accept values that are between the indicator of the uniform
  # in the case, the aest always has to be smaller than the best, aest always smaller than the
  # minimum of the vector x, and best always greater than the maximum of x. a <x (1) <x (n) <b.
  if ((aest<best) && (aest<min(x)) && (best>max(x))){
    num = dnorm(amc[i-1], aest, sqrt(Sigmaa))*dnorm(bmc[i-1], best, sqrt(Sigmab))
    denom = dnorm(aest,amc[i-1], sqrt(Sigmaa))*dnorm(best, bmc[i-1], sqrt(Sigmab))
    ratio = exp(posterior.uniform(x,aest,best)-posterior.uniform(x,amc[i-1],bmc[i-1]))
    probac = min(1,ratio*(num/denom))
    
    # generates a point on the uniform (0,1)
    u = runif(1)
    
    if (u<probac){
      amc[i] = aest
      bmc[i] = best
    } else{
      amc[i] = amc[i-1]
      bmc[i] = bmc[i-1]
    }
  } else {
    # if the aliasing of the if of the beginning is not false, if a <x (1) <x (n) <b, it
    # return to the beginning, and generate new values of 'aest' and 'best' plus value
    # of the vector will continue the same in the observation in which the condition was false.
    amc[i]=amc[i-1]
    bmc[i]=bmc[i-1]
  }
}

par(mfrow=c(1,2)) 
plot.ts(amc,col=2, main="Plot of 'amc'",ylab="amc",xlab="time") 
plot.ts(bmc,col=4, main="Plot of 'bmc'",ylab="bmc",xlab="time")

# burn-in: exclude the first 500 observations (initail values effect).
amc = amc[c(501:10500)]
bmc = bmc[c(501:10500)]

# thinning by 10 
thin=10 
nmcb=length(amc)/thin

amcb=array(0,c(nmcb))  
bmcb=array(0,c(nmcb))

for (i in 1:nmcb) {
  amcb[i]=amc[thin*i]  
  # to always get the value of 10 in 10 of the 'amc'
  bmcb[i]=bmc[thin*i]  
  # to always get the value of 10 in 10 of the 'bmc'
}

par(mfrow=c(1,2)) 
plot.ts(amcb,col=2, main="Plot of 'amcb'",ylab="amcb",xlab="time")  
plot.ts(bmcb,col=4, main="Plot of 'bmcb'",ylab="bmcb",xlab="time")  

# histogram of the posterior distribution
hist(amcb, col=2, main="Histogram of 'amcb'",ylab="Frequency",xlab="amcb")  # histograma para 'amcb'
hist(bmcb, col=4, main="Histogram of 'bmcb'",ylab="Frequency",xlab="bmcb")  # histograma para 'bmcb'

# mean of 'amcb' e 'bmcb'
mean(amcb)
mean(bmcb)

# Credibility intervals of the posterior distributions
LIa=quantile(amcb,0.025)
LSa=quantile(amcb,0.975)

LIa
LSa

LIb=quantile(bmcb,0.025)
LSb=quantile(bmcb,0.975)

LIb
LSb
