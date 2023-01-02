library(logOfGamma)
#compute Pr(y(s,m)|gamma_sj,alpha)
basic<-function(y,gamma,alpha,K){
  #y: the multinomial matrix of m*K
  #gamma: binary indicator vector
  #alpha: a matrix of K*J
  #K: the number of categories
  #J: the number of components
  
  j=which(gamma==1) #the component Pr(y(s,m)|gamma_sj,alpha) belongs to
  # the first multiplier
  temp0=sum(apply(y,1,function(x){gammaln(sum(x)+1)-sum(gammaln(x+1))}))
  # the second multiplier
  temp1=gammaln(sum(alpha[,j]))-gammaln(sum(colSums(y))+sum(alpha[,j]))
  #the third multiplier
  temp2=sum(gammaln(colSums(y)+alpha[,j])-gammaln(alpha[,j]))
  
  result=list(basic_comp=exp(temp0+temp1+temp2),
              bf_comp=exp(temp1+temp2))
  return (result)
}

#compute E_gamma_s(gamma_s|y(s,m),p,alpha) as a vector of length J
expectation<-function(y,p,alpha,K,J){
  #y: the multinomial matrix of m*K
  #p: the current estimator of p
  #alpha: the current estimator alpha
  
  temp=rep(0,J)
  for (j in 1:J){
    gamma=rep(0,J)
    gamma[j]=1
    temp[j]=basic(y,gamma,alpha,K)$basic_comp
  }
  temp=p*temp
  
  return(temp/sum(temp))
}

#compute loglikelihood
loglihood<-function(X,n,m,p,alpha,K,J){
  temp=matrix(0,nrow=n-m+1,ncol=J)
  for (i in 1:(n-m+1)){
    for (j in 1:J){
      gamma=rep(0,J)
      gamma[j]=1
      temp[i,j]=basic(matrix(X[i:(i+m-1),],nrow=m),gamma,alpha,K)$basic_comp
    }
  }
  return(sum(log(p%*%t(temp))))
}


Estimate<-function(X,n,m,K,J,epsilon=0.1){
  #J is the candidate set
  P_list=list()
  Alpha_list=list()
  BIC=rep(0,J)
  #choose Jhat
  for (j in 1:J){
    #Bayes detection
    #initialization
    OldP=rep(1/j,j)
    OldAlpha=matrix(1,nrow=K,ncol=j)
    errP=100
    errAl=100
    i=1
    Estep=matrix(0,nrow=n-m+1,ncol=j)
    #iteration
    while ((errP>epsilon|errAl>epsilon)&i<=(n-m+1)){
      Estep[i,]=expectation(matrix(X[i:(i+m-1),],nrow=m),OldP,OldAlpha,K,j)
      NewP=colMeans(matrix(Estep[1:i,],nrow=i))
      SolutionAl=lbfgs(rep(1/(K*j),K*j),function(alpha){
        tempalpha=matrix(alpha,nrow=K,ncol=j)
        tempmatrix=matrix(0,nrow=i,ncol=j)
        for (ii in 1:i){
          for (jj in 1:j){
            tempmatrix[ii,jj]=gammaln(sum(tempalpha[,jj]))-
              gammaln(sum(X[ii:(ii+m-1),])+sum(tempalpha[,jj]))+
              sum(gammaln(colSums(matrix(X[ii:(ii+m-1),],nrow=m))+
                            tempalpha[,jj])-
                    gammaln(tempalpha[,jj]))
          }
        }
        return(sum(Estep[1:i,]*tempmatrix))
      },
                       lower=rep(0,K*j))$par
      NewAlpha=matrix(SolutionAl,nrow=K,ncol=j)
      errP=max(abs(NewP-OldP)/OldP)
      errAl=max(abs(NewAlpha-OldAlpha)/OldAlpha)
      OldP=NewP
      OldAlpha=NewAlpha
      i=i+1
    }
    P_list[[j]]=NewP
    Alpha_list[[j]]=NewAlpha
    BIC[j]=-loglihood(X[1:(i+m-2),],(i+m-2),m,NewP,NewAlpha,K,j)+0.5*j*(2*K+1)*log(i)
  }
  result=list(Phat=P_list,Alphahat=Alpha_list,BIC=BIC)
  return(result)
}

#compute b_tm
bayes<-function(y,alpha,p,K,J){
  #y: the multinomial matrix of m*K
  #p: the result of p
  #alpha: the result of alpha
  #K: the number of categories
  temp=rep(0,J)
  for (j in 1:J){
    gamma=rep(0,J)
    gamma[j]=1
    temp[j]=basic(y,gamma,alpha,K)$basic_comp
  }
  return(sum(p*temp))
}

#detected change-points, non-change-points
Bayes_detec<-function(X,n,m,p,alpha,K,J,eta){
  #eta: threshold
  bf_value=rep(0,(n-m+1))
  for(i in 1:(n-2*m+1)){
    y1=matrix(X[i:(i+m-1),],nrow=m)
    y2=matrix(X[(i+m):(i+2*m-1),],nrow=m)
    y3=matrix(X[i:(i+2*m-1),],nrow=2*m)
    b1=bayes(y1,alpha,p,K,J)
    b2=bayes(y2,alpha,p,K,J)
    b3=bayes(y3,alpha,p,K,J)
    bf_value[i]=2*log(b1*b2/b3)
  }
  bf_value[is.na(bf_value)]=Inf
  That=which(bf_value>=eta)
  tempT=c()
  if(length(That)>0){
    temp=split(That,cumsum(seq_along(That) %in% (
      which(diff(That)>1)+1)))
    for(i in 1:length(temp)){
      tempT=c(tempT,
              temp[[i]][which.max(bf_value[temp[[i]]])])
    } 
  }
  That=tempT+m
  return(list(bf_value=bf_value,That=That))
}

#power analysis
evaluation<-function(N,Tc,m,bf_value,eta){
  #N: number of sample size
  #Tc: true change-points
  #bf_value: value of bayes factor
  #eta: threshold
  TP=0
  FP=0
  TN=0
  FN=0
  for (i in 1:length(Tc)){
    if (Tc[i]>m){
      if (bf_value[Tc[i]-m]>=eta){
        TP=TP+1}
      if (bf_value[Tc[i]-m]<=((-1)*eta)){
        FN=FN+1}
    }
  }
  
  NTc=setdiff(c(1:N),Tc)
  for (i in 1:length(NTc)){
    if (NTc[i]>m){
      if (bf_value[NTc[i]-m]>=eta){
        FP=FP+1}
      if (bf_value[NTc[i]-m]<=((-1)*eta)){
        TN=TN+1}
    }
  }
  
  Power0=min(1,TN/(length(NTc)))
  Power1=min(1,TP/(length(Tc)))
  
  result=list(P0=Power0, P1=Power1)
  return (result)
}


#evaluation metrics
evaluation2<-function(N,m=11,Tc,That){
  #Tc: true change-points
  #That: detected change-points

  TP=0
  for (i in 1:length(That)){
    temp=intersect(Tc,c((That[i]-m+1):(That[i]+m-1)))
    if (length(temp)!=0){
      TP=TP+1}
  }
  
  PT=0
  for (i in 1:length(Tc)){
    temp=intersect(That,c((Tc[i]-m+1):(Tc[i]+m-1)))
    if (length(temp)!=0){
      PT=PT+1}
  }
  
  Precision=TP/length(That)
  Recall=PT/length(Tc)
  Fvalue=(2*Precision*Recall)/(Precision+Recall)
  result=list(Pre=Precision, Rec=Recall, Fv=Fvalue)
  
  return (result)
}


#Generate data
GenerateData<-function(n,K,J,lambda1=20,lambda2=5){
  #n: length of the data stream
  #K: number of categories
  #J: number of components
  #lambda1: interval of change-points
  #lambda2: size of X_i
  
  p=rep(0,J)
  for (j in 1:(J-1)){
    p[j]=runif(1,0,1/(J-1))
  }
  p[J]=1-sum(p)
  
  Alpha=matrix(runif(K*J),nrow=K,ncol=J)
  for (j in 1:J){
    Alpha[j,j]=5
  }
  Beta=apply(Alpha,2,function(x){rdirichlet(1,x)})
  X=c()
  Tc=c()
  Gamma=c()
  d0=0
  while(max(nrow(X),1)<n){
    d=rpois(1,lambda1)
    d0=d0+d
    gamma=t(rmultinom(1,1,p))
    Gamma=c(Gamma,which(gamma==1))
    #data stream
    for(ii in 1:d){
      X=rbind(X,t(rmultinom(n=1,rpois(1,lambda2),gamma%*%t(Beta))))
    }
    Tc=c(Tc,d0+1)
  }
  Tc=Tc[1:(length(Tc)-1)]
  Tc=Tc[Gamma[1:(length(Gamma)-1)]!=Gamma[2:length(Gamma)]]
  result=list(Tc=Tc,X=t(X))
  return(result)
}