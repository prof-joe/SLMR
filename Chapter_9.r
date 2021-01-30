### Statistical Learning with Math and R, Springer ###

# Chapter 9 Support Vector Machine

## 9.3 The solution of Support Vector Machine

library(quadprog)
svm.1=function(X,y,C){
  eps=0.0001; n=nrow(X); meq=1; Dmat=matrix(nrow=n,ncol=n);
  for(i in 1:n)for(j in 1:n)Dmat[i,j]=sum(X[i,]*X[j,])*y[i]*y[j];
  Dmat=Dmat+eps*diag(n) ## Adding small values to the diagonal elements to make the matrix nonsingular
  dvec= rep(1,n)
  Amat=matrix(nrow=(2*n+1),ncol=n); Amat[1,]=y;
  Amat[2:(n+1),1:n]=-diag(n); Amat[(n+2):(2*n+1),1:n]= diag(n);
  Amat=t(Amat) ## This package requires to specify the transpose for Amat.
  bvec= c(0,rep(-C,n),rep(0,n))
  alpha=solve.QP(Dmat,dvec,Amat,bvec=bvec,meq=1)$solution
  beta=drop((alpha*y)%*%X);
  index=(1:n)[eps<alpha&alpha<C-eps]
  beta.0=mean(y[index]-X[index,]%*%beta)
  return(list(beta=beta,beta.0=beta.0))
}

a=rnorm(1); b=rnorm(1); n=100;
X=matrix(rnorm(n*2),ncol=2,nrow=n); 
y=sign(a*X[,1]+b*X[,2]+0.1*rnorm(n))
plot(-3:3,-3:3,xlab="第1成分",ylab="第2成分", type="n")
for(i in 1:n){
  if(y[i]==1)points(X[i,1],X[i,2],col="red") 
  else points(X[i,1],X[i,2],col="blue")
}
qq=svm.1(X,y,10); 
abline(-qq$beta.0/qq$beta[2],-qq$beta[1]/qq$beta[2])


# 9.4 Extension of Support Vector Machine using Kernels
K.linear <-function(x,y) {return(t(x)%*%y)}
K.poly <-function(x,y) {return((1+t(x)%*%y)^2)}
svm.2=function(X,y,C, K){ ## svm.2 is another function name
  eps=0.0001; n=nrow(X); Dmat=matrix(nrow=n,ncol=n)
  for(i in 1:n)for(j in 1:n) Dmat[i,j]=K(X[i,],X[j,])*y[i]*y[j]
  Dmat=Dmat+eps*diag(n); dvec=rep(1,n)
  Amat=matrix(nrow=(2*n+1),ncol=n); Amat[1,]=y; Amat[2:(n+1),1:n]=-diag(n);
  Amat[(n+2):(2*n+1),1:n]=diag(n) ; Amat=t(Amat)
  bvec=c(0,-C*rep(1,n),rep(0,n)); meq=1
  alpha=solve.QP(Dmat,dvec,Amat,bvec=bvec,meq=1)$solution
  index=(1:n)[eps<alpha&alpha<C-eps]
  beta=drop((alpha*y)%*%X); beta.0=mean(y[index]-X[index,]%*%beta)
  return(list(alpha=alpha,beta.0=beta.0))
}
# Function Definition 
plot.kernel=function(K, lty){ ## Specify the line property via the parameter lty
  qq=svm.2(X,y,1,K); alpha=qq$alpha; beta.0=qq$beta.0
  f=function(u,v){x=c(u,v); S=beta.0; for(i in 1:n)S=S+alpha[i]*y[i]*K(X[i,], x); return(S)}
  ## f gives the height z at (x,y), and the contour can be obtained.
  u=seq(-2,2,.1);v=seq(-2,2,.1);w=array(dim=c(41,41))
  for(i in 1:41)for(j in 1:41)w[i,j]=f(u[i],v[j])
  contour(u,v,w,level=0,add=TRUE,lty=lty)
}
# Execution
a=3; b=-1
n=200; X=matrix(rnorm(n*2),ncol=2,nrow=n); y=sign(a*X[,1]+b*X[,2]^2+0.3*rnorm(n))
plot(-3:3,-3:3,xlab="X[,1]",ylab="X[,2]", type="n")
for(i in 1:n){
  if(y[i]==1)points(X[i,1],X[i,2],col="red") else points(X[i,1],X[i,2],col="blue")
}
plot.kernel(K.linear,1); 
plot.kernel(K.poly,2)

library(e1071)
x=matrix(rnorm(200*2), ncol=2); x[1:100,]=x[1:100,]+2; x[101:150,]=x[101:150,]-2
y=c(rep(1,150), rep(2,50)); dat=data.frame(x=x,y=as.factor(y))
train=sample(200,100)
svmfit=svm(y~., data=dat[train,], kernel ="radial", gamma=1, cost=100)
plot(svmfit, dat[train,])

tune.out=tune(svm, y~., data=dat[train ,], kernel="radial",
ranges=list(cost=c(0.1,1,10,100,1000), gamma=c(0.5,1,2,3,4)))
summary (tune.out)
