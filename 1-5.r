# 第５章 正則化
## 5.1 Ridge

ridge=function(X, y, lambda=0){
  X=as.matrix(X); p=ncol(X); n=length(y); X.bar=array(dim=p); s=array(dim=p)
  for(j in 1:p){X.bar[j]=mean(X[,j]);X[,j]=X[,j]-X.bar[j];};
  for(j in 1:p){s[j]=sd(X[,j]);X[,j]=X[,j]/s[j]};
  y.bar=mean(y); y=y-y.bar
  beta=drop(solve(t(X)%*%X+n*lambda*diag(p))%*%t(X)%*%y)
  for(j in 1:p)beta[j]=beta[j]/s[j]
  beta.0= y.bar-sum(X.bar*beta)
  return(list(beta=beta, beta.0=beta.0))
}
df=read.table("crime.txt"); x=df[,3:7]; y=df[,1]; p=ncol(x);
lambda.seq=seq(0,50); 
plot(lambda.seq, xlim=c(0,50), ylim=c(-7.5,15),
      xlab="lambda",ylab="beta",main="各lambdaについての各係数の値", type="n", col="red")
for(j in 1:p){
  coef.seq=NULL; for(lambda in lambda.seq)coef.seq=c(coef.seq,ridge(x,y,lambda)$beta[j])
  par(new=TRUE); lines(lambda.seq,coef.seq, col=j)
}
legend("topright",legend=c("警察への年間資金","25歳以上で高校を卒業した人の割合",
    "16--19歳で高校に通っていない人の割合","18--24歳で大学生の割合",
    "25歳以上で4年制大学を卒業した人の割合"), col=1:p, lwd=2, cex =.8)

## 5.2 劣微分
curve(x^2-3*x+abs(x), -2,2, main="y=x^2-3x+|x|"); points(1,-1, col="red", pch=16)
curve(x^2+x+2*abs(x), -2,2, main="x^2+x+2|x|"); points(0,0, col="red", pch=16)


## 5.3 Lasso
soft.th=function(lambda,x)sign(x)*pmax(abs(x)-lambda,0)
curve(soft.th(5,x),-10,10, main="soft.th(lambda,x)")
segments(-5,-4,-5,4, lty=5, col="blue"); segments(5,-4,5,4, lty=5, col="blue")
text(-0.2,1,"lambda=5",cex=1.5)

lasso=function(X, y, lambda=0){
  X=as.matrix(X); p=ncol(X); n=length(y); X.bar=array(dim=p); s=array(dim=p)
  for(j in 1:p){X.bar[j]=mean(X[,j]);X[,j]=X[,j]-X.bar[j];}; 
  for(j in 1:p){s[j]=sd(X[,j]);X[,j]=X[,j]/s[j]} 
  y.bar=mean(y); y=y-y.bar
  eps=1; beta=array(0, dim=p); beta.old=array(0, dim=p)
  while(eps>0.001){
      for(j in 1:p){
        r= y-X[,-j]%*%beta[-j]
        beta[j]= soft.th(lambda,sum(r*X[,j])/n)
      }
      eps=max(abs(beta-beta.old)); beta.old=beta
  }
  for(j in 1:p)beta[j]=beta[j]/s[j] 
  beta.0= y.bar-sum(X.bar*beta)
  return(list(beta=beta, beta.0=beta.0))
}
df=read.table("crime.txt"); x=df[,3:7]; y=df[,1]; p=ncol(x);
lambda.seq=seq(0,200);
plot(lambda.seq, xlim=c(0,200), ylim=c(-7.5,15),
xlab="lambda",ylab="beta",main="各lambdaについての各係数の値", type="n", col="red")
for(j in 1:p){
   coef.seq=NULL; for(lambda in lambda.seq)coef.seq=c(coef.seq,lasso(x,y,lambda)$beta[j])
   par(new=TRUE); lines(lambda.seq,coef.seq, col=j)
}
legend("topright",legend=
    c("警察への年間資金","25歳以上で高校を卒業した人の割合",
    "16--19歳で高校に通っていない人の割合","18--24歳で大学生の割合",
    "25歳以上で4年制大学を卒業した人の割合"), col=1:p, lwd=2, cex =.8)

## 5.5 λ の値の設定
library(glmnet)
df=read.table("crime.txt"); X=as.matrix(df[,3:7])
cv.fit=cv.glmnet(X,y); 
plot(cv.fit)
lambda.min=cv.fit$lambda.min; 
lambda.min
fit=glmnet(X,y,lambda=lambda.min); 
fit$beta
