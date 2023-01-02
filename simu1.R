#setwd("C:/Users/Yiwei/Documents/Fanyiwei/Essays/Changepoint/Code")
source("myfunction.R")
library(gtools)
library(nloptr)

#simu2
#parameter setting

simu1<-function(K,nn){
  n=1000 #sample size
  J=3 #number of components
  #K: number of categories, K=10,30,50
  
  mset=c(1:8)#set of number of trials
  for (j in 1:length(mset)){
    m=mset[j] #number of trials
    Tc=c()
    while(length(Tc)==0){
      #generate data
      mydata=GenerateData(n,K,J,lambda2=nn)
      #true changepoints
      Tc=mydata$Tc
    }
    
    X=t(mydata$X)

    #burn-in period, estimation
    BF_fit=Estimate(X[1:100,],100,m,K,1)
    Jhat=which.min(BF_fit$BIC)
    Phat=unlist(BF_fit$Phat[Jhat])
    Alphahat=matrix(unlist(BF_fit$Alphahat[Jhat]),nrow=K)
    result=Bayes_detec(X,n,m,Phat,Alphahat,K,Jhat,2)
    
    #evaluation
    result=evaluation(n,Tc,m,result$bf_value,eta=2)
    print(result)
  }
}

