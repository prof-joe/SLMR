### Statistical Learning with Math and R, Springer ###

# Chapter 4 Resampling

## 4.1 Cross Validation

cv.linear=function(X,y,k){
  n=length(y); m=n/k; S=0 # Assume k divides n
  for(j in 1:k){
      test=((j-1)*m+1):(j*m)
      ## Specify test data out of n
      beta=solve(t(X[-test,])%*%X[-test,])%*%t(X[-test,])%*%y[-test]
      ## Estimate beta using data except test
      e=y[test]-X[test,]%*%beta; S=S+drop(t(e)%*%e)
      ## Evaluate the estimated beta by applying it to the test
  }
  return(S/n)
}
## Data Generation ##
n=100; p=5; X=matrix(rnorm(n*p),ncol=p); X=cbind(1,X)
beta=rnorm(p+1); beta[c(2,3)]=0; eps=rnorm(n); y=X%*%beta+eps
## Evaluation via CV ##
cv.linear(X[,c(1,4,5,6)], y, 10)
cv.linear(X[,c(1,2,3,4)], y, 10)
cv.linear(X,y,10)

n=100; p=5; X=matrix(rnorm(n*p),ncol=p); X=cbind(1,X); beta=rnorm(p+1); beta[c(2,3)]=0;
U=NULL; V=NULL
for(j in 1:100){
  eps=rnorm(n); y=X%*%beta+eps
  U=c(U,cv.linear(X[,c(1,4,5,6)], y, 10)); V=c(V,cv.linear(X,y,10))
}
plot(U,V,xlab="cv.linear(X[,c(1,4,5)], y, 10)", ylab="cv.linear(X, y, 10)",
     main="OverEstimation due to selecting too many")
abline(a=0,b=1,col="red")

## Data Generation ##
plot(0,0,xlab="k",ylab="CV values", xlim=c(2,n),ylim=c(0.3,1.5),type="n")
for(j in 2:11){
  n=100; p=5; X=matrix(rnorm(n*p),ncol=p);
  X=cbind(1,X); beta=rnorm(p+1); eps=rnorm(n); y=X%*%beta+eps
  U=NULL; V=NULL;
  for(k in 2:n)if(n%%k==0){U=c(U,k); V=c(V,cv.linear(X,y,k))}; lines(U,V, col=j)
}

knn.1=function(x,y,z,k){
    x=as.matrix(x); n=nrow(x); p=ncol(x); dis=array(dim=n)
    for(i in 1:n)dis[i]=norm(z-x[i,],"2")
    S=order(dis)[1:k]  ## The set of indices i that are the closest k
    u=sort(table(y[S]),decreasing=TRUE) ## The most often y[i] among the k and its Number
 ## Tie Breaking start
    while(length(u)>1 && u[1]==u[2]){ k=k-1; S=order(dis)[1:k]; u=sort(table(y[S]),decreasing=TRUE)}
 ## Tie Breaking end
    return(names(u)[1])
}
knn=function(x,y,z,k){
  n=nrow(z); w=array(dim=n)
  for(i in 1:n)w[i]=knn.1(x,y,z[i,],k)
  return(w)
}
df=iris; df=df[sample(1:150,150,replace=FALSE),]
n=nrow(df); U=NULL; V=NULL
for(k in 1:10){
  top.seq=1+seq(0,135,15); S=0
  for(top in top.seq){
      index= top:(top+14); knn.ans=knn(df[-index,1:4],df[-index,5],df[index,1:4],k)
      ans= df[index,5]; S=S+sum(knn.ans!=ans)
   }
  S=S/n; U=c(U,k);V=c(V,S)
}
plot(0,0,type="n", xlab="k", ylab="Error Rate", xlim=c(1,10),ylim=c(0,0.1),
    main="Error Probability evaluated by CV")
lines(U,V,col="red")

## 4.2 The formula for linear regression

cv.fast=function(X,y,k){
  n=length(y); m=n/k; H=X%*%solve(t(X)%*%X)%*%t(X);
  I=diag(rep(1,n)); e=(I-H)%*%y; I=diag(rep(1,m))
  S=0
  for(j in 1:k){ test=((j-1)*m+1):(j*m); S=S+norm(solve(I-H[test,test])%*%e[test],"2")^2 }
  return(S/n)
}
n=1000; p=5; beta=rnorm(p+1); x=matrix(rnorm(n*p),ncol=p)
X=cbind(rep(1,n),x); y=X%*%beta+rnorm(n) ## Data Generation
plot(0,0,xlab="k",ylab="Execution", xlim=c(2,n),ylim=c(0,0.5),type="n")
U=NULL; V=NULL
for(k in 10:n)if(n%%k==0){
  t=proc.time()[3]; cv.fast(X,y,k); U=c(U,k); V=c(V, (proc.time()[3]-t)) 
}
lines(U,V, col="blue")
U=NULL; V=NULL
for(k in 10:n)if(n%%k==0){
  t=proc.time()[3]; cv.linear(X,y,k); U=c(U,k); V=c(V,(proc.time()[3]-t))
}
lines(U,V, col="red")
legend("topleft",legend=c("cv.linear","cv.fast"),col=c("red","blue"), lty=1)

## 4.3 BootStrap

bt=function(df,f,r){
  m=nrow(df); org=f(df,1:m); u=array(dim=r)
  for(j in 1:r){index=sample(1:m,m,replace=TRUE); u[j]=f(df,index)}
  return(list(original=org, bias=mean(u)-org, stderr=sd(u)))
}
func.1=function(data,index){
  X=data$X[index]; Y=data$Y[index]
  return((var(Y)-var(X))/(var(X)+var(Y)-2*cov(X,Y)))
}
library(ISLR)
bt(Portfolio,func.1,1000)
df=read.table("crime.txt")
for(j in 1:3){
  func.2=function(data,index)coef(lm(V1~V3+V4,data=data,subset=index))[j]
  print(bt(df,func.2,1000))
}
