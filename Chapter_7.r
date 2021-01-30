### Statistical Learning with Math and R, Springer ###

# Chapter 7 NonLinear Regression

## 7.1 Polynomial Regression

n=100; x=rnorm(n); y=sin(x)+rnorm(n)　　　　　　　　　　　## data generation
m=3;p.set=c(3,5,7); col.set=c("red","blue","green")
g=function(beta,u){S=beta[1]; for(j in 1:p)S=S+beta[j+1]*u^j; return(S)}
for(i in 1:m){
  p=p.set[i]; X=rep(1,n);for(j in 1:p)X=cbind(X,x^j)
  beta=drop(solve(t(X)%*%X)%*%t(X)%*%y); f=function(u)g(beta,u)
  curve(f(x),-3,3, col=col.set[i],yaxt="n"); par(new=TRUE)
}
legend("topleft",lty=1,paste0("p=",p.set),col=col.set); points(x,y)

## Data Generation Closed to Even Function
n=100; x=rnorm(n)*pi; y=ceiling(x)%%2*2-1+rnorm(n)*0.2
plot(x,y,xaxt="n",yaxt="n",ann=FALSE, main="Follow Random Numbers via sin and cos")

## The function f below chooses 1, cos x, cos 2x, cos 3x as the basis
X=cbind(1,cos(x),cos(2*x),cos(3*x)); beta=solve(t(X)%*%X)%*%t(X)%*%y
f=function(x)beta[1]+beta[2]*cos(x)+beta[3]*cos(2*x)+beta[4]*cos(3*x)
par(new=TRUE); 
curve(f(x),-5,5, col="red",yaxt="n",ann=FALSE)

## The function g below chooses 1, sin x, sin 2x, sin 3x as the basis
X=cbind(1,sin(x),sin(2*x),sin(3*x)); beta=solve(t(X)%*%X)%*%t(X)%*%y
g=function(x)beta[1]+beta[2]*sin(x)+beta[3]*sin(2*x)+beta[4]*sin(3*x)
par(new=TRUE); 
curve(g(x),-5,5,col="blue",yaxt="n",ann=FALSE)

## 7.2 Spline Regression

n=100; x=rnorm(n)*2*pi; y=sin(x)+0.2*rnorm(n) ## Data Generation via Random Numbers
col.set=c("red","green","blue"); K.set=c(5,7,9) ## Knots
for(k in 1:3){
  K=K.set[k]; knots=seq(-2*pi,2*pi,length=K)
  X=matrix(nrow=n,ncol=K+4)
  for(i in 1:n){
    X[i,1]= 1; X[i,2]= x[i]; X[i,3]= x[i]^2; X[i,4]= x[i]^3
    for(j in 1:K)X[i,j+4]=max((x[i]-knots[j])^3,0)
  }
  beta=solve(t(X)%*%X)%*%t(X)%*%y ## beta estimation
  f=function(x){S=beta[1]+beta[2]*x+beta[3]*x^2+beta[4]*x^3;
    for(j in 1:K)S=S+beta[j+4]*max((x-knots[j])^3,0)
    return(S)
  }　　　　　　　　　　　　　　　　　　　　　## obtaining Function f
  u.seq=seq(-5,5,0.02); v.seq=NULL; for(u in u.seq)v.seq=c(v.seq,f(u))
  plot(u.seq,v.seq,type="l",col=col.set[k], yaxt="n", xlab="x", ylab="f(x)")
  par(new=TRUE)
}
legend(-2.2,1,paste0("K=",K.set), lty=1, col=col.set); 
points(x,y)

## 7.3 Regression to Natural Spline Curves 

d=function(j,x,knots){
  K=length(knots);
  (max((x-knots[j])^3,0)-max((x-knots[K])^3,0))/(knots[K]-knots[j])
}
h=function(j,x,knots){
  K=length(knots);
  if(j==1) return(1) 
  else if(j==2)return(x) 
  else return(d(j-2,x,knots)-d(K-1,x,knots))
}
n=100; x=rnorm(n)*2*pi; y=sin(x)+0.2*rnorm(n); ## Data Generation
K=11; knots=seq(-5,5,length=K); X=matrix(nrow=n,ncol=K+4)
for(i in 1:n){
  X[i,1]= 1; X[i,2]= x[i]; X[i,3]= x[i]^2; X[i,4]= x[i]^3
  for(j in 1:K)X[i,j+4]=max((x[i]-knots[j])^3,0)
}
beta=solve(t(X)%*%X)%*%t(X)%*%y
f=function(x){ ## Spline Function
  S=beta[1]+beta[2]*x+beta[3]*x^2+beta[4]*x^3;
  for(j in 1:K)S=S+beta[j+4]*max((x-knots[j])^3,0)
  return(S)
}
X=matrix(nrow=n,ncol=K); X[,1]=1; for(j in 2:K)for(i in 1:n)X[i,j]=h(j,x[i],knots)
gamma=solve(t(X)%*%X)%*%t(X)%*%y
g=function(x){　　　　　　　　　　　　　　　## Natural Spline Function
  S=gamma[1]; for(j in 2:K)S=S+gamma[j]*h(j,x,knots); return(S)
}
u.seq=seq(-6,6,0.02); ## Draw the Functions as a Graph
v.seq=NULL; for(u in u.seq)v.seq=c(v.seq,f(u))
plot(u.seq,v.seq,type="l",col="blue", yaxt="n", xlab="x",ylab="f(x),g(x)"); par(new=TRUE);
w.seq=NULL; for(u in u.seq)w.seq=c(w.seq,g(u))
plot(u.seq,w.seq,type="l",col="red", yaxt="n", xlab="",ylab=""); par(new=TRUE)
legend(-3.7,1.1,c("スプライン","自然なスプライン"), lty=1, col=c("blue","red"))
points(x,y); abline(v=knots,lty=3); abline(v=c(-5,5),lwd=2); title("K=11")

## 7.4 Smoothing Spline

G=function(x){ ## The values of x is assumed to be in the ascending ordere
  n=length(x); g=matrix(0, nrow=n,ncol=n)
  for(i in 3:n)for(j in i:n){
    g[i,j]=12*(x[n]-x[n-1])*(x[n-1]-x[j-2])*(x[n-1]-x[i-2])/(x[n]-x[i-2])/(x[n]-x[j-2])+
      (12*x[n-1]+6*x[j-2]-18*x[i-2])*(x[n-1]-x[j-2])^2/(x[n]-x[i-2])/(x[n]-x[j-2])
    g[j,i]=g[i,j]
  }
  return(g)
}
n=100; x=runif(n,-5,5); y=x+sin(x)*2+rnorm(n) ## Data Generation
index=order(x); x=x[index];y=y[index]
X=matrix(nrow=n,ncol=n); X[,1]=1
for(j in 2:n)for(i in 1:n)X[i,j]=h(j,x[i],x) ## Generation of Matrix X
GG=G(x); ## Generation of Matrix G
lambda.set=c(1,30,80); col.set=c("red","blue","green")
for(i in 1:3){
  lambda=lambda.set[i]
  gamma=solve(t(X)%*%X+lambda*GG)%*%t(X)%*%y
  g=function(u){S=gamma[1]; for(j in 2:n)S=S+gamma[j]*h(j,u,x); return(S)}
  u.seq=seq(-8,8,0.02); v.seq=NULL; for(u in u.seq)v.seq=c(v.seq,g(u))
  plot(u.seq,v.seq,type="l",yaxt="n", xlab="x",ylab="g(x)",ylim=c(-8,8), col=col.set[i])
  par(new=TRUE)
}
points(x,y); legend("topleft", paste0("lambda=",lambda.set), col=col.set, lty=1)
title("Smoothing Spline (n=100)")

cv.ss.fast=function(X,y,lambda, G, k){
  n=length(y); m=n/k;
  H=X%*%solve(t(X)%*%X+lambda*G)%*%t(X); df=sum(diag(H))
  I=diag(rep(1,n)); e=(I-H)%*%y; I=diag(rep(1,m))
  S=0
  for(j in 1:k){
    test=((j-1)*m+1):(j*m);
    S=S+norm(solve(I-H[test,test])%*%e[test],"2")^2
  }
  return(list(score=S/n,df=df))
}
## Data Generation
n=100; x=runif(n,-5,5); y=x-0.02*sin(x)-0.1*rnorm(n)
index=order(x); x=x[index];y=y[index]
X=matrix(nrow=n,ncol=n); X[,1]=1; for(j in 2:n)for(i in 1:n)X[i,j]=h(j,x[i], x)
GG=G(x)
## Plot the Efficient Degree and Prediction Error for each lambda
u=seq(1,50); v=NULL; w=NULL
for(lambda in u){
  result=cv.ss.fast(X,y,lambda,GG,n); v=c(v,result$df); w=c(w,result$score)
}
plot(v,w,type="l",col="red",xlab="Efficient Degree",ylab="Prediction Error via CV")
title("The efficient degree and Predition Error via CV")

## 7.5 Local Regression

n=250; x=2*rnorm(n); y=sin(2*pi*x)+rnorm(n)/4 ## Data Generation
D=function(t) max(0.75*(1-t^2),0) ## Function Definition D
K=function(x,y,lambda) D(abs(x-y)/lambda) ## Function Definition K
f=function(z,lambda){ ## Function Definition f
  S=0; T=0;
  for(i in 1:n){S=S+K(x[i],z,lambda)*y[i]; T=T+K(x[i],z,lambda)}
  return(S/T)
}
plot(seq(-3,3,length=10),seq(-2,3,length=10),type="n",xlab="x", ylab="y"); points(x,y)
xx=seq(-3,3,0.1)
yy=NULL;for(zz in xx)yy=c(yy,f(zz,0.05)); lines(xx,yy,col="green")
yy=NULL;for(zz in xx)yy=c(yy,f(zz,0.25)); lines(xx,yy,col="blue")
## The curves for lambda=0.05, 0.25 have been drawn
m=n/10
lambda.seq=seq(0.05,1,0.01); SS.min=Inf
for(lambda in lambda.seq){
  SS=0
  for(k in 1:10){
    test=((k-1)*m+1):(k*m); train=setdiff(1:n,test)
      for(j in test){u=0; v=0;
        for(i in train){
          kk=K(x[i],x[j],lambda); u=u+kk*y[i]; v=v+kk
        }
        if(v==0){
          d.min=Inf;
          for(i in train){d=abs(x[j]-x[i]); if(d<d.min){d.min=d; index=i}};
          z=y[index]
        }
        else z=u/v
        SS=SS+(y[j]-z)^2
      }
  }
  if(SS<SS.min){SS.min=SS;lambda.best=lambda}
}
## The optimum lambda has been computed
yy=NULL;for(zz in xx)yy=c(yy,f(zz,lambda.best)); lines(xx,yy,col="red")
title("Nadaraya-Watson Estimator")
legend("topleft",legend=paste0("lambda=",c(0.05, 0.25, "lambda.best")),
      lwd=1,col=c("green","blue","red"))
local=function(x,y,z=x){
  X=cbind(rep(1,n),x); yy=NULL; beta.hat=array(dim=2)
  for(u in z){
    w=array(dim=n); for(i in 1:n)w[i]=K(x[i],u,lambda=1); W=diag(w)
    beta.hat= solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y; yy=c(yy,beta.hat[1]+beta.hat[2]*u)
  }
  return(yy)
}
n=30; x=runif(n)*2*pi-pi; 
y=sin(x)+rnorm(n); 
plot(x,y) ## Data Generaion
m=200; U=seq(-pi,pi,pi/m); V=local(x,y,U)
lines(U,V,col="red",type="l"); 
title("Local Linear Regression (p=1, N=30)")

## 7.6 Generalized Additive Model

poly=function(x,y,z=x){
  n=length(x);m=length(z); X=cbind(rep(1,n),x,x^2,x^3); yy=array(dim=n);
  beta.hat=array(dim=4); beta.hat= solve(t(X)%*%X)%*%t(X)%*%y;
  X=cbind(rep(1,m),z,z^2,z^3); yy= X%*% beta.hat
  return(yy)
}
## The same function local is used as in the previous section.
n=30; x=runif(n)*2*pi-pi; y=sin(x)+rnorm(n); plot(x,y) ## Data Generation
y.1=0; y.2=0; for(k in 1:10){y.1=poly(x,y-y.2); y.2= local(x,y-y.1)}
z=seq(-2,2,0.1); par(mfrow=c(1,2))
plot(z,poly(x,y.1,z),type="l", xlab="x", ylab="f(x)", main="Polynomial Regression （Order 3）", col="green")
plot(z,local(x,y.2,z),type="l", xlab="x", ylab="f(x)", main="Local Linear Regression",col="green")