# Chapter 2 Linear Regression
## 1.1 Least Square Method

min.sq=function(x,y){ # 最小二乗法の切片と傾きを求める関数 min.sq
  x.bar=mean(x); y.bar=mean(y)
  beta.1=sum((x-x.bar)*(y-y.bar))/sum((x-x.bar)^2); beta.0=y.bar-beta.1*x.bar
  return(list(a=beta.0, b=beta.1))
  }
a=rnorm(1); b=rnorm(1); # 直線の係数をランダムに生成
N=100; x=rnorm(N); y=a*x+b+rnorm(N) # 直線の周りの点をランダムに生成
plot(x,y); abline(h=0); abline(v=0) # 点のプロット
abline(min.sq(x,y)$a, min.sq(x,y)$b,col="red") # 中心化前の直線
x=x-mean(x); y=y-mean(y) # 中心化
abline(min.sq(x,y)$a,min.sq(x,y)$b,col="blue") # 中心化後の直線
legend("topleft",c("中心化前","中心化後"),lty=1, col=c("red","blue")) # 凡例

## 1.2 Multiple Regression

n=100; p=2; beta=c(1,2,3)
x=matrix(rnorm(n*p),nrow=n,ncol=p)
y=beta[1]+beta[2]*x[,1]+beta[3]*x[,2]+rnorm(n)　　　# 標準正規分布にしたがう雑音
X=cbind(1,x) 　　　　　　　 　# 左側にすべて1の列をおく

solve(t(X)%*%X)%*%t(X)%*%y # betaを推定する

## 1.4 RSSの分布

i=1; curve(dchisq(i,x), 0, 20, col=i)
for(i in 2:10)curve(dchisq(x,i), 0, 20,col=i,add=TRUE,ann=FALSE)
legend("topright",legend=1:10,lty=1, col=1:10)

## 1.5 β^j≠0の仮説検定

curve(dnorm(x), -10,10, ann=FALSE, ylim=c(0,0.5), lwd=5)
for(i in 1:10)curve(dt(x,df=i), -10, 10, col=i, add=TRUE, ann=FALSE)
legend("topright",legend=1:10,lty=1, col=1:10)

n=100; x=rnorm(n)+2;
plot(1,1,xlim=c(0.5,1.5),ylim=c(0.5,1.5),xlab="beta.0",ylab="beta.1")
for(i in 1:100){
  y=1+x+rnorm(n); z=cbind(1,x); beta.est=solve(t(z)%*%z)%*%t(z)%*%y
  points(beta.est[1],beta.est[2],col=i)
  }
abline(v=1); abline(h=1)

sum(x)/n; sum(x^2)/n

N=100; x=rnorm(N); y=rnorm(N)
x.bar=mean(x); y.bar=mean(y)
beta.0=sum(y.bar*sum(x^2)-x.bar*sum(x*y))/sum((x-x.bar)^2)
beta.1=sum((x-x.bar)*(y-y.bar))/sum((x-x.bar)^2)
RSS=sum((y-beta.0-beta.1*x)^2); RSE=sqrt(RSS/(N-1-1))
B.0=sum(x^2)/N/sum((x-x.bar)^2); B.1=1/sum((x-x.bar)^2)
se.0=RSE*sqrt(B.0); se.1=RSE*sqrt(B.1)
t.0=beta.0/se.0; t.1=beta.1/se.1
p.0=2*(1-pt(abs(t.0),N-2)) # p値（その値より外側にある確率）
p.1=2*(1-pt(abs(t.1),N-2)) # p値（その値より外側にある確率）
beta.0;se.0;t.0;p.0;
beta.1;se.1;t.1;p.1
lm(y~x)

summary(lm(y~x))

N=100; r=1000
T=NULL
for(i in 1:r){
  x=rnorm(N); y=rnorm(N); x.bar=mean(x); y.bar=mean(y)
  fit=lm(y~x);beta=fit$coefficients
  RSS=sum((y-fit$fitted.values)^2); RSE=sqrt(RSS/(N-1-1))
  B.1=1/sum((x-x.bar)^2); se.1=RSE*sqrt(B.1)
  T=c(T,beta[2]/se.1)
 }
hist(T,breaks=sqrt(r),probability=TRUE, xlab="tの値",ylab="確率密度",
main="tの値のヒストグラムと理論値（赤）")
curve(dt(x, N-2),-3,3,type="l", col="red",add=TRUE)

## 1.6 決定係数と共線形性の検出

R2=function(x,y){
  y.hat=lm(y~x)$fitted.values; y.bar=mean(y)
  RSS=sum((y-y.hat)^2); TSS=sum((y-y.bar)^2)
  return(1-RSS/TSS)
}
N=100; m=2; x=matrix(rnorm(m*N),ncol=m); y=rnorm(N); R2(x,y)
N=100; m=1; x=matrix(rnorm(m*N),ncol=m); y=rnorm(N)
R2(x,y)
cor(x,y)^2

vif=function(x){
  p=ncol(x); values=array(dim=p); for(j in 1:p)values[j]=1/(1-R2(x[,-j],x[,j]))
  return(values)
  }
library(MASS) 
x=as.matrix(Boston) 
vif(x)

## 1.7 信頼区間と予測区間

#データを生成
N=100; p=1; X=matrix(rnorm(N*p),ncol=p); X=cbind(rep(1,N),X)
beta=c(1,1); epsilon=rnorm(N); y=X%*%beta+epsilon
 #関数f(x), g(x)を定義。 Uはt(X)%*%Xの逆行列
U=solve(t(X)%*%X); beta.hat=U%*%t(X)%*%y;
RSS=sum((y-X%*%beta.hat)^2); RSE=sqrt(RSS/(N-p-1)); alpha=0.05
f=function(x, a){ #a=0なら信頼区間，a=1なら予測区間
x=cbind(1,x); range=qt(df=N-p-1,1-alpha/2)*RSE*sqrt(a+x%*%U%*%t(x));
return(list(lower=x%*%beta.hat-range,upper=x%*%beta.hat+range))
}
x.seq=seq(-10,10,0.1)
# グラフで信頼区間を表示
lower.seq=NULL; for(x in x.seq)lower.seq=c(lower.seq, f(x,0)$lower)
upper.seq=NULL; for(x in x.seq)upper.seq=c(upper.seq, f(x,0)$upper)
x.lim=c(min(x.seq),max(x.seq)); y.lim=c(min(lower.seq),max(upper.seq))
plot(x.seq, lower.seq, col="blue",xlim=x.lim, ylim=y.lim, xlab="x",
     ylab="y", type="l")
par(new=TRUE); 
plot(x.seq, upper.seq,col="red", xlim=x.lim, ylim=y.lim, 
     xlab="",ylab="", type="l", axes=FALSE)
par(new=TRUE); 
# グラフで予測区間を表示
lower.seq=NULL; for(x in x.seq)lower.seq=c(lower.seq, f(x,1)$lower)
upper.seq=NULL; for(x in x.seq)upper.seq=c(upper.seq, f(x,1)$upper)
x.lim=c(min(x.seq),max(x.seq)); y.lim=c(min(lower.seq),max(upper.seq))
plot(x.seq, lower.seq, col="blue",xlim=x.lim, ylim=y.lim, 
     xlab="",ylab="", type="l", lty=4, axes=FALSE)
par(new=TRUE); 
plot(x.seq, upper.seq, col="red", xlim=x.lim, ylim=y.lim, 
     xlab="",ylab="", type="l", lty=4,
axes=FALSE)
abline(beta.hat[1],beta.hat[2])

