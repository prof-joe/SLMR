# 第３章 リサンプリング
## 3.1 クロスバリデーション

cv.linear=function(X,y,k){
  n=length(y); m=n/k; S=0 # kはnを割り切ることを仮定
  for(j in 1:k){
      test=((j-1)*m+1):(j*m)
      ## n組のデータのうち，どのデータをテスト用に使うかを指定
      beta=solve(t(X[-test,])%*%X[-test,])%*%t(X[-test,])%*%y[-test]
      ##テストで用いるデータ以外で，係数betaを推定する
      e=y[test]-X[test,]%*%beta; S=S+drop(t(e)%*%e)
      ## 係数betaに，テスト用のデータをあてはめて評価する
  }
  return(S/n)
}
## データ生成 ##
n=100; p=5; X=matrix(rnorm(n*p),ncol=p); X=cbind(1,X)
beta=rnorm(p+1); beta[c(2,3)]=0; eps=rnorm(n); y=X%*%beta+eps
## クロスバリデーションによる評価 ##
cv.linear(X[,c(1,4,5,6)], y, 10)
cv.linear(X[,c(1,2,3,4)], y, 10)
cv.linear(X,y,10)
1.05391062681677
5.30939219302796
1.08605237073673
n=100; p=5; X=matrix(rnorm(n*p),ncol=p); X=cbind(1,X); beta=rnorm(p+1); beta[c(2,3)]=0;
U=NULL; V=NULL
for(j in 1:100){
  eps=rnorm(n); y=X%*%beta+eps
  U=c(U,cv.linear(X[,c(1,4,5,6)], y, 10)); V=c(V,cv.linear(X,y,10))
}
plot(U,V,xlab="cv.linear(X[,c(1,4,5)], y, 10)", ylab="cv.linear(X, y, 10)",
     main="変数を多く選びすぎて過学習")
abline(a=0,b=1,col="red")

## データ生成 ##
plot(0,0,xlab="k",ylab="CVの値", xlim=c(2,n),ylim=c(0.3,1.5),type="n")
for(j in 2:11){
  n=100; p=5; X=matrix(rnorm(n*p),ncol=p);
  X=cbind(1,X); beta=rnorm(p+1); eps=rnorm(n); y=X%*%beta+eps
  U=NULL; V=NULL;
  for(k in 2:n)if(n%%k==0){U=c(U,k); V=c(V,cv.linear(X,y,k))}; lines(U,V, col=j)
}

knn.1=function(x,y,z,k){
    x=as.matrix(x); n=nrow(x); p=ncol(x); dis=array(dim=n)
    for(i in 1:n)dis[i]=norm(z-x[i,],"2")
    S=order(dis)[1:k];　　## 距離の小さいk個の添え字iの集合
    u=sort(table(y[S]),decreasing=TRUE) ## k個の中の最頻のy[i]と頻度
 ## タイブレーキングの処理
    while(length(u)>1 && u[1]==u[2]){ k=k-1; S=order(dis)[1:k]; u=sort(table(y[S]),decreasing=TRUE)}
 ## ここまで
    return(names(u)[1])
}
knn=function(x,y,z,k){
  n=nrow(z); w=array(dim=n); for(i in 1:n)w[i]=knn.1(x,y,z[i,],k)
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
plot(0,0,type="n", xlab="k", ylab="誤り率", xlim=c(1,10),ylim=c(0,0.1),
    main="cvによる誤り率の評価")
lines(U,V,col="red")

## 3.2 線形回帰の場合の公式

cv.fast=function(X,y,k){
  n=length(y); m=n/k; H=X%*%solve(t(X)%*%X)%*%t(X);
  I=diag(rep(1,n)); e=(I-H)%*%y; I=diag(rep(1,m))
  S=0
  for(j in 1:k){ test=((j-1)*m+1):(j*m); S=S+norm(solve(I-H[test,test])%*%e[test],"2")^2 }
  return(S/n)
}
n=1000; p=5; beta=rnorm(p+1); x=matrix(rnorm(n*p),ncol=p)
X=cbind(rep(1,n),x); y=X%*%beta+rnorm(n) ## データ生成
plot(0,0,xlab="k",ylab="実行時間", xlim=c(2,n),ylim=c(0,0.5),type="n")
U=NULL; V=NULL
for(k in 10:n)if(n%%k==0){
  t=proc.time()[3]; cv.fast(X,y,k); U=c(U,k); V=c(V, (proc.time()[3]-t)) 
# 2020-5-6 AF氏が指摘、2020-5-10 JSが修正
}
lines(U,V, col="blue")
U=NULL; V=NULL
for(k in 10:n)if(n%%k==0){
  t=proc.time()[3]; cv.linear(X,y,k); U=c(U,k); V=c(V,(proc.time()[3]-t))
# 2020-5-6 AF氏が指摘、2020-5-12 JSが修正
}
lines(U,V, col="red")
legend("topleft",legend=c("cv.linear","cv.fast"),col=c("red","blue"), lty=1)

## 3.3 ブートストラップ

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
