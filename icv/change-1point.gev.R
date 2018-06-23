
rm(list=ls(all=TRUE))
set.seed(27)
library(evir)

# metropolis-hastings for posterior of gev distribution
posterior.gev.cp=function(data,mumca,sigmamca,ximca){
  n=length(data)
  lpost=function(mu,sigma,xi){
    logpost=-n*log(sigma)-sum((1+xi*(data-mu)/sigma)^(-1/xi))
    logpost=logpost-(1+1/xi)*sum(log((1+xi*(data-mu)/sigma)))
    logpost=logpost+(0.1-1)*log(sigma)-0.1*sigma-mu^2/10000-xi^2/1
    logpost
  }
  Vu=5
  Vsigma=5
  Vxi=0.5
  muest=rnorm(1,mumca,Vu)
  xiest=rnorm(1,ximca,Vxi)
  sigmaest=max(rgamma(1,sigmamca^2/Vsigma,sigmamca/Vsigma),0.01)
  while (min(1+xiest*(data-muest)/sigmaest)<0){
    muest=rnorm(1,mumca,Vu)
    xiest=rnorm(1,ximca,Vxi)
    sigmaest=max(rgamma(1,sigmamca^2/Vsigma,sigmamca/Vsigma),0.01)
  }
  if (min(1+ximca*(data-mumca)/sigmamca)<0){
    alpha=0
  } else{
    alpha=exp(lpost(muest,sigmaest,xiest)-lpost(mumca,sigmamca,ximca))
    alpha=alpha*dgamma(sigmamca,sigmaest^2/Vsigma,sigmaest/Vsigma)/(dgamma(sigmaest,sigmamca^2/Vsigma,sigmamca/Vsigma))
  }
  u=runif(1)
  if (u<alpha){
    mumc=muest
    sigmamc=sigmaest
    ximc=xiest
  } else{
    mumc=mumca
    sigmamc=sigmamca
    ximc=ximca
  }
  cbind(mumc,sigmamc,ximc)
}

# function to sample the change-point in the time series
change.1point.gev=function(int,data,kf,n){
  mumc1=array(0,c(int))
  sigmamc1=array(0,c(int))
  ximc1=array(0,c(int))
  mumc2=array(0,c(int))
  sigmamc2=array(0,c(int))
  ximc2=array(0,c(int))
  k=array(0,c(int))
  l=array(0,c(n-1))
  ln=array(0,c(n-1))
  lp=array(0,c(n-1))
  
  # initial values for the parameters
  k[1]=kf
  
  ajuste1=gev(data[1:k[1]],1)
  mumc1[1]=ajuste1$par.ests[3]
  sigmamc1[1]=ajuste1$par.ests[2]
  ximc1[1]=ajuste1$par.ests[1]
  
  ajuste2=gev(data[(k[1]+1):n],1)
  mumc2[1]=ajuste2$par.ests[3]
  sigmamc2[1]=ajuste2$par.ests[2]
  ximc2[1]=ajuste2$par.ests[1]
  
    for (i in 2:int){   
      # updating the parameters of GEV
      at1=posterior.gev.cp(data[1:k[i-1]],mumc1[i-1],sigmamc1[i-1],ximc1[i-1])
      at2=posterior.gev.cp(data[(k[i-1]+1):n],mumc2[i-1],sigmamc2[i-1],ximc2[i-1])
      mumc1[i]=at1[1]
      sigmamc1[i]=at1[2]
      ximc1[i]=at1[3]
      mumc2[i]=at2[1]
      sigmamc2[i]=at2[2]
      ximc2[i]=at2[3]
      
      # update k
      for(j in 1:(n-1)){
        ln[j] = sum((-log(sigmamc1[i])) - (((1/ximc1[i])+1)*log(1+ximc1[i]*((data[1:j] - mumc1[i])/sigmamc1[i]))) - ((1 + (ximc1[i]*((data[1:j] - mumc1[i])/sigmamc1[i])))^(-1/ximc1[i]))) + sum((-log(sigmamc2[i])) - (((1/ximc2[i])+1)*log(1+ximc2[i]*((data[(j+1):n] - mumc2[i])/sigmamc2[i]))) - ((1 + (ximc2[i]*((data[(j+1):n] - mumc2[i])/sigmamc2[i])))^(-1/ximc2[i])))
        if (is.nan(ln[j])){l[j]=min(l)}  else {l[j]=ln[j]} 
      }
      
      lp=(exp(l-max(l))/sum(exp(min(l)-l)))
      k[i]=sample(1:(n-1), size=1, prob=lp)
      cat(i, "/", int, "\r")
    }
    options(warn=-1)
    list(mumc1=mumc1,sigmamc1=sigmamc1,ximc1=ximc1,mumc2=mumc2,sigmamc2=sigmamc2,ximc2=ximc2,k=k)
}

# simulation
data.1 = rgev(100,-0.3,40,100)
data.2 = rgev(150,-0.3,80,100)
data=c(data.1,data.2)
int=1000 
n=length(data)

#initial kick for the change point
kf=90

post.values = change.1point.gev(int,data,kf,n)

sapply(post.values,mean)

sapply(post.values,quantile,probs=c(0.05, 0.95))

nn=length(post.values$k)
tm = table(post.values$k)/(nn)
plot(tm, xlab="k", ylab="Probabilities", col=2)

# aplication example for parnaiba river
dataP=read.table("C:\\Users\\WYARAVMS\\Google Drive\\graduação-ufpi-2011.2015\\ufpi-2015\\icv-valores-extremos(GEV)\\Análise II\\Parndados.csv", header=T, sep=";")
data=dataP[,2]
int=1000 
n=length(data)
kf=100

post.values = change.1point.gev(int,data,kf,n)

sapply(post.values,mean)

sapply(post.values,quantile,probs=c(0.05, 0.95))

nn=length(post.values$k)
tm = table(post.values$k)/(nn)
plot(tm, xlab="k", ylab="Probabilities", col=2)


