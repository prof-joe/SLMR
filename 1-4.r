# Chapter 5 Information Criteria
## 5.1 Information Criteria 

RSS.min=function(X,y,T){
  m=ncol(T); S.min=Inf
  for(j in 1:m){
    q=T[,j]; S=sum((lm(y~X[,q])$fitted.values-y)^2)
  if(S<S.min){S.min=S; set.q=q}
  }
return(list(value=S.min,set=set.q))
}
library(MASS)
df=Boston; X=as.matrix(df[,c(1,3,5,6,7,8,10,11,12,13)]); y=df[[14]];
p=ncol(X); n=length(y)
AIC.min=Inf
for(k in 1:p){
  T=combn(1:p,k); res=RSS.min(X,y,T)
  AIC= n*log(res$value/n)+2*k ##
  if(AIC<AIC.min){AIC.min=AIC; set.min= res$set}
}
AIC.min
set.min

y.bar=mean(y); TSS=sum((y-y.bar)^2); D.max=-Inf
for(k in 1:p){
  T=combn(1:p,k); res=RSS.min(X,y,T)
  D= 1-res$value/(n-k-1)/TSS/(n-1); if(D>D.max){D.max=D ; set.max= res$set}
}
D.max
set.max

library(MASS)
df=Boston; X=as.matrix(df[,c(1,3,5,6,7,8,10,11,12,13)]); y=df[[14]];
n=nrow(X); p=ncol(X)
IC=function(k){
  T=combn(1:p,k); 
  res=RSS.min(X,y,T)
  AIC= n*log(res$value/n)+2*k; 
  BIC= n*log(res$value/n)+k*log(n)
  return(list(AIC=AIC,BIC=BIC))
}
AIC.seq=NULL; BIC.seq=NULL;
for(k in 1:p){AIC.seq=c(AIC.seq, IC(k)$AIC); BIC.seq=c(BIC.seq, IC(k)$BIC)}
plot(1:p, ylim=c(min(AIC.seq),max(BIC.seq)), type="n",xlab="# of variables", ylab="IC values")
lines(AIC.seq,col="red"); lines(BIC.seq,col="blue")
legend("topright",legend=c("AIC","BIC"), col=c("red","blue"), lwd=1, cex =.8)

# 4.2 Efficient Estimator and Fisher Information Matrix
# 4.3 Kullback-Leibler divergence
# 4.4 Derivation of Akaike's Information Criterion
