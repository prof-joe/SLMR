# Chapter 3 Classification
## 3.1 Logistic Curve

f=function(x)exp(beta.0+beta*x)/(1+exp(beta.0+beta*x))
beta.0=0; beta.seq=c(0,0.2,0.5,1,2,10); m=length(beta.seq); beta=beta.seq[1]
plot(f,xlim=c(-10,10),ylim=c(0,1),xlab="x",ylab="P(Y=1|x)", col=1,
  main="Logistic Curve")
for(i in 2:m){
  beta=beta.seq[i]; par(new=TRUE);
  plot(f,xlim=c(-10,10),ylim=c(0,1),xlab="", ylab="", axes=FALSE, col=i)
 }
legend("topleft", legend=beta.seq, col=1:length(beta.seq), lwd=2, cex=.8)

## 3.2 Newton-Raphson
f=function(x) x^2-1; f.=function(x) 2*x
curve(f(x),-1,5); abline(h=0,col="blue")
x=4
for(i in 1:10){
  X=x; Y=f(x); x=x-f(x)/f.(x); y=f(x)
  segments(X,Y,x,0); segments(X,Y,X,0, lty=3)
  points(x,0,col="red",pch=16)
}

f=function(z)z[1]^2+z[2]^2-1; f.x=function(z) 2*z[1]; f.y=function(z)2*z[2];
g=function(z)z[1]+z[2]; g.x=function(z) 1; g.y=function(z) 1;
z=c(3,4)
for(i in 1:10){
  z=z-solve(matrix(c(f.x(z),f.y(z),g.x(z),g.y(z)),ncol=2,byrow=TRUE))%*%c(f(z),g(z))
}
z

## Data Generation　##
N=1000; p=2; X=matrix(rnorm(N*p),ncol=p); X=cbind(rep(1,N),X)
beta=rnorm(p+1); y=array(N); s=as.vector(X%*%beta); prob=1/(1+exp(s));
for(i in 1:N)if(runif(1)>prob[i])y[i]=1 else y[i]=-1
beta

## Maximum Likelihood　##
beta=Inf; gamma=rnorm(p+1)
while(sum((beta-gamma)^2)>0.001){
  beta=gamma
  s=as.vector(X%*%beta); 
  v=exp(-s*y); 
  u= y*v/(1+v)
  w= v/(1+v)^2; W=diag(w); 
  z= s+u/w
  gamma=as.vector(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z)
  print(gamma)
}

n=100; x=c(rnorm(n)+1, rnorm(n)-1); y=c(rep(1,n),rep(-1,n))
train=sample(1:(2*n),n,replace=FALSE); df=data.frame(x,y)
x=as.matrix(df[train,1]); y=as.vector(df[train,2])
p=1; X=cbind(1,x); beta=0; gamma=rnorm(p+1)
while(sum((beta-gamma)^2)>0.001){
  beta=gamma
  s=as.vector(X%*%beta); v=exp(-s*y); u= y*v/(1+v)
  w= v/(1+v)^2; W=diag(w); z= s+u/w
  gamma=as.vector(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z)
  print(gamma)
}

x=as.matrix(df[-train,1]); y=as.vector(df[-train,2]) ## y is the correct answer
z=2*as.integer(beta[1]+x*beta[2]>0)-1 ## z is used for discrimination
table(y,z) ## the number of diagonal elements out of 100

## 3.3 Linear and Quadratic Descrimination
mu.1=c(2,2); sigma.1=2; sigma.2=2; rho.1=0
mu.2=c(-3,-3); sigma.3=1; sigma.4=1; rho.2=-0.8
n=100
u=rnorm(n); v=rnorm(n); x.1=sigma.1*u+mu.1[1];
y.1=(rho.1*u+sqrt(1-rho.1^2)*v)*sigma.2+mu.1[2]
u=rnorm(n); v=rnorm(n); x.2=sigma.3*u+mu.2[1];
y.2=(rho.2*u+sqrt(1-rho.2^2)*v)*sigma.4+mu.2[2]
f=function(x,mu,inv,de)drop(-0.5*t(x-mu)%*%inv%*%(x-mu)-0.5*log(de))
mu.1=mean(c(x.1,y.1)); mu.2=mean(c(x.2,y.2));
df=data.frame(x.1,y.1); mat=cov(df); inv.1=solve(mat); de.1=det(mat) #
df=data.frame(x.2,y.2); mat=cov(df); inv.2=solve(mat); de.2=det(mat) #
f.1=function(u,v)f(c(u,v),mu.1,inv.1,de.1);
f.2=function(u,v)f(c(u,v),mu.2,inv.2,de.2)
pi.1=0.5; pi.2=0.5
u = v = seq(-6, 6, length=50); m=length(u); w=array(dim=c(m,m))
for(i in 1:m)for(j in 1:m)w[i,j]=log(pi.1)+f.1(u[i],v[j])-log(pi.2)-f.2(u[i],v[j])
# plot
contour(u,v,w,level=0)   
points(x.1,y.1,col="red"); points(x.2,y.2,col="blue")

df=data.frame(c(x.1,y.1)-mu.1, c(x.2,y.2)-mu.2); inv.1=solve(mat); de.1=det(mat)
inv.2=inv.1; de.2=de.1
#iris Data
f=function(w,mu,inv,de)-0.5*(w-mu)%*%inv%*%t(w-mu)-0.5*log(de)
df=iris; df[[5]]=c(rep(1,50),rep(2,50),rep(3,50))
n=nrow(df); train=sample(1:n,n/2,replace=FALSE); test=setdiff(1:n,train)
mat=as.matrix(df[train,])
mu=list(); covv=list()
for(j in 1:3){
    x=mat[mat[,5]==j,1:4];
    mu[[j]]=c(mean(x[,1]),mean(x[,2]),mean(x[,3]),mean(x[,4]))
    covv[[j]]=cov(x)
}
g=function(v,j)f(v,mu[[j]],solve(covv[[j]]),det(covv[[j]]))
z=array(dim=length(test))
for(i in test){
    u=as.matrix(df[i,1:4]); a=g(u,1);b=g(u,2); c=g(u,3)
    if(a<b){if(b<c)z[i]=3 else z[i]=2} 
    else {if(a<c)z[i]=3 else z[i]=1}
}
table(z[test],df[test,5])

## 3.4 K Nearest neighbour
knn.1=function(x,y,z,k){
    x=as.matrix(x); n=nrow(x); p=ncol(x); dis=array(dim=n)
    for(i in 1:n)dis[i]=norm(z-x[i,],"2")
    S=order(dis)[1:k];　　
    u=sort(table(y[S]),decreasing=TRUE) 
 ## Tie Breasking
    while(length(u)>1 && u[1]==u[2]){ k=k-1; S=order(dis)[1:k]; u=sort(table(y[S]),decreasing=TRUE)}
 ##
    return(names(u)[1])
}
knn=function(x,y,z,k){
  n=nrow(z); w=array(dim=n); for(i in 1:n)w[i]=knn.1(x,y,z[i,],k)
  return(w)
}

df=iris;
    n=150; train=sample(1:n,n/2,replace=FALSE); test=setdiff(1:n,train)
	x=as.matrix(df[train,1:4]); y=as.vector(df[train,5])
	z=as.matrix(df[test,1:4]); ans=as.vector(df[test,5])
	w=knn(x, y, z, k=3)
	table(w,ans)

## 3S.5 ROC Curves
N.0=10000;N.1=1000; mu.1=1; mu.0=-1; var.1=1; var.0=1
x=rnorm(N.0,mu.0,var.0); y=rnorm(N.1,mu.1,var.1) # x: healthy，y: sick
plot(1:1,1:1,xlim=c(0,1),ylim=c(0,1), xlab="False Positive", ylab="True Positive",
main="ROC curve", type="n")
theta.seq=exp(seq(-10,100,0.1))
U=NULL; V=NULL
for(theta in theta.seq){
    u=sum(dnorm(x,mu.1,var.1)/dnorm(x,mu.0,var.0)>theta)/N.0 #
	v=sum(dnorm(y,mu.1,var.1)/dnorm(y,mu.0,var.0)>theta)/N.1 # 
	U=c(U,u); V=c(V,v)
}
lines(U,V,col="blue")
M=length(theta.seq)-1; AUC=0; for(i in 1:M)AUC=AUC+abs(U[i+1]-U[i])*V[i]
text(0.5,0.5,paste("AUC=",AUC),col="red")
