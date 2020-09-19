# 第７章 決定木
## 7.1 回帰の決定木

sq.loss=function(y){y.bar=mean(y); return(sum((y-y.bar)^2))}
branch=function(x,y,f,S,m=ncol(x)){ ## Sはサンプル集合の添え字
  n=length(S); p=ncol(x); if(n==0)return(NULL); best.score=Inf
  for(j in 1:p)for(i in S){
    left=NULL; right=NULL
    for(k in S)if(x[k,j]<x[i,j])left=c(left,k) else right=c(right,k)
    L=f(y[left]); R=f(y[right]); score=L+R
    if(score<best.score){
      best.score=score;
      info=list(i=i, j=j, left=left, right=right, score=best.score,left.score=L, right.score=R)
    }
  }
  return(info)
}
dt=function(x,y,f="sq.loss",alpha=0, n.min=1, m=ncol(x)){
  if(f=="sq.loss")g=sq.loss 
  else if(f=="mis.match") g=mis.match 
  else if(f=="gini") g=gini 
  else g=entropy
  n=length(y); stack=list();
  stack[[1]]=list(parent=0, set=1:n, score=g(y)); vertex=list(); k=0
  while(length(stack)>0){
    r=length(stack); node=stack[[r]]; stack=stack[-r] ## POP
    k=k+1 ## PUSHされた頂点のid=kは，POPされてはじめて付与
    res=branch(x, y, g, node$set, m) ## 最適な分割に関する情報
    if(node$score-res$score<alpha||length(node$set)<n.min||length(res$left)==0||length(res$right)==0){
      vertex[[k]]=list(parent=node$parent, j=0, set=node$set)
    } 
    else{
      vertex[[k]]=list(parent=node$parent,set=node$set, th=x[res$i,res$j],j=res$j)
      stack[[r]]=list(parent=k, set=res$right, score=res$right.score) ## PUSH
      stack[[r+1]]=list(parent=k, set=res$left, score=res$left.score) ## PUSH
    }
  }
  ## 端点にその値(center)，内点にその左右の子のID(left,right)
  mode=function(y) names(sort(table(y),decreasing=TRUE))[1] ## 最頻値
  r=length(vertex)
  for(h in 1:r){vertex[[h]]$left=0; vertex[[h]]$right=0}
  for(h in r:2){
    pa=vertex[[h]]$parent
    if(vertex[[pa]]$right==0)vertex[[pa]]$right=h 
    else vertex[[pa]]$left=h
  }
  if(f=="sq.loss") g=mean else g=mode
  for(h in 1:r)if(vertex[[h]]$j==0)vertex[[h]]$center=g(y[vertex[[h]]$set])
  return(vertex)
}
library(MASS); library(igraph)
x=as.matrix(Boston[,1:13]); y=as.vector(Boston[,14])
vertex=dt(x,y,n.min=50)
## 以下，グラフの出力
r=length(vertex); col=array(dim=r); edge.list=matrix(nrow=r,ncol=2)
for(h in 1:r)col[h]=vertex[[h]]$j; for(h in 1:r)edge.list[h,]=c(vertex[[h]]$parent,h)
edge.list=edge.list[-1,]; g=graph_from_edgelist(edge.list); V(g)$color=col
plot(g, layout = layout.reingold.tilford(g, root=1))
## 以下，表の出力
VAR=NULL; TH=NULL
for(h in 1:r)if(vertex[[h]]$j!=0){
  i=vertex[[h]]$i; j=vertex[[h]]$j; VAR=c(VAR,j); TH=c(TH,x[i,j])
}
cbind(VAR,TH)
VAR

value=function(u,vertex){
  r=1
  while(vertex[[r]]$j!=0)if(u[vertex[[r]]$j]<vertex[[r]]$th)
  r=vertex[[r]]$left else r=vertex[[r]]$right
  return(vertex[[r]]$center)
}
df=Boston[1:100,]; n=nrow(df); p=ncol(df)
x=as.matrix(df[,1:13]); y=as.vector(df[,14])
alpha.seq=seq(0,1.5,0.1); s=floor(n/10); m=length(vertex); out=NULL
for(alpha in alpha.seq){
  SS=0
  for(h in 1:10){ ## 10-fold CV
    test=(h*s-s+1):(h*s); train=setdiff(1:n,test)
    vertex=dt(x[train,],y[train],alpha=alpha)
    for(t in test)SS=SS+(y[t]-value(x[t,],vertex))^2
  }
  out=c(out,SS/100)
}
plot(alpha.seq,out,type="l",ylim=c(10.5,12), xlab="alpha",
     ylab="二乗誤差", main="CVで最適なalpha (Bostonデータセット, N=100)")

## 7.2 分類の決定木
mode=function(y)names(sort(table(y),decreasing=TRUE))[1] ## 最頻値
## 誤り率
mis.match=function(y){y.hat=mode(y); return(sum(y!=y.hat))}
## Gini
gini=function(y){
  n=length(y); if(n==0)return(0); z=as.vector(table(y)); m=length(z);
  T=0;for(j in 1:m)T=T+z[j]*(n-z[j])/n;
  return(T/n)
}
## エントロピー
entropy=function(y){
  n=length(y); if(n==0)return(0); z=as.vector(table(y))
  m=length(z); T=0;for(j in 1:m)if(z[j]!=0)T=T+z[j]*log(n/z[j]);
  return(T)
}
df=iris; x=as.matrix(df[,1:4]); y=as.matrix(df[,5])
vertex=dt(x,y,"mis.match",n.min=4); m=length(vertex);
u=NULL; v=NULL
for(h in 1:m)if(vertex[[h]]$j==0){
  w=y[vertex[[h]]$set]; u=c(u,rep(mode(w),length(w))); v=c(v,w)
}
table(u,v)
col=array(dim=m); edge.list=matrix(nrow=m,ncol=2)
for(h in 1:m)col[h]=vertex[[h]]$j
for(h in 1:m)edge.list[h,]=c(vertex[[h]]$parent,h); 
edge.list=edge.list[-1,]
library(igraph)
g=graph_from_edgelist(edge.list); 
V(g)$color=col; 
V(g)$name=""
plot(g, pin=c(4,4), layout = layout.reingold.tilford(g, root=1)); title("誤り率")
## Giniとエントロピーも同様に実行した（mis.matchをgini, entropyに）

df=iris; n=nrow(df); p=ncol(df); index=sample(1:n,n,replace=FALSE)
x=as.matrix(df[index,1:4]); y=as.matrix(df[index,5])
n.min.seq=1:10; s=15; out=NULL
for(n.min in n.min.seq){
  SS=0
  for(h in 10:1){
    test=(15*h-14):(15*h); train=setdiff(1:n,test)
    vertex=dt(x[train,],y[train],"mis.match",n.min,p);
    for(i in test)SS=SS+as.integer(value(x[i,],vertex)!=y[i])
  }
  print(SS)
}

## 7.3 バギング

par(mfrow=c(2,4))
n=200; p=5; x=matrix(rnorm(n*p),nrow=n,ncol=p)
beta=rnorm(p); y=abs(round(x%*%beta+rnorm(n)))+1 ## データ生成
for(h in 1:8){
  index=sample(1:n,n,replace=TRUE); x=x[index,];y=y[index]
  vertex=dt(x,y,"mis.match",n.min=6)
  r=length(vertex); col=array(dim=m); edge.list=matrix(nrow=r,ncol=2)
  for(h in 1:r)col[h]=vertex[[h]]$j
  for(h in 1:r)edge.list[h,]=c(vertex[[h]]$parent,h); 
  edge.list=edge.list[-1,]
  g=graph_from_edgelist(edge.list); 
  V(g)$color=col; 
  V(g)$name=""
  plot(g, layout = layout.reingold.tilford(g, root=1))
}
par(mfrow=c(1,1))

## 7.4 ランダムフォレスト

branch=function(x,y,f,S, m=ncol(x)){ ## mの値の設定, デフォルトはp
  n=length(S); p=ncol(x); if(n==0)return(NULL); best.score=Inf
  if(m<p) T=sample(1:p, m, replace=FALSE) 
  else T=1:p ## ここが違う
  for(j in T)for(i in S){ ## Tの中で最適な変数を選ぶ
    left=NULL; right=NULL
    for(k in S)if(x[k,j]<x[i,j])left=c(left,k) else right=c(right,k)
    L=f(y[left]); R=f(y[right]); score=L+R
    if(score<best.score){
      best.score=score;
      info=list(i=i, j=j, left=left, right=right, score=best.score,
      left.score=L, right.score=R)
    }
  }
  return(info)
}
rf=function(z){
  zz=array(dim=c(B,50)); zzz=NULL
  for(b in 1:B){
    for(i in 1:50)zz[b,i]=(mode(z[1:b,i])==y[i+100])
    zzz=c(zzz,sum(zz[b,]))
  }
  return(zzz)
}
set.seed(101)
df=iris; n=nrow(df); p=ncol(df); 
index=sample(1:n,n,replace=FALSE)
x=as.matrix(df[index,1:4]); y=as.vector(df[index,5])
train=1:100; test=101:150; B=100; z=array(dim=c(B,50)); m=4
for(b in 1:B){
  index=sample(train, 100, replace=TRUE)
  vertex=dt(x[index,], y[index], "mis.match", n.min=2, m=m)
  for(i in test)z[b,i-100]=value(x[i,],vertex)
}
z4=z
## m=4をm=3,m=2,m=1に変えて, z4,z3,z2,z1に格納
m=3
for(b in 1:B){
  index=sample(train, 100, replace=TRUE)
  vertex=dt(x[index,], y[index], "mis.match", n.min=2, m=m)
  for(i in test)z[b,i-100]=value(x[i,],vertex)
}
z3=z
m=2
for(b in 1:B){
  index=sample(train, 100, replace=TRUE)
  vertex=dt(x[index,], y[index], "mis.match", n.min=2, m=m)
  for(i in test)z[b,i-100]=value(x[i,],vertex)
}
z2=z
m=1
for(b in 1:B){
  index=sample(train, 100, replace=TRUE)
  vertex=dt(x[index,], y[index], "mis.match", n.min=2, m=m)
  for(i in test)z[b,i-100]=value(x[i,],vertex)
}
z1=z
plot(1:B,rf(z4)-0.1,type="l",ylim=c(0,50),col=2,xlab="木の生成回数",ylab="正答回数/50回",
      main="ランダムフォレスト")
lines(1:B,rf(z3),col=3); 
lines(1:B,rf(z2)+0.1,col=4); 
lines(1:B,rf(z1)-0.1,col=5)
legend("bottomright",legend=c("m=4","m=3","m=2","m=1"),col=c(2,3,4,5),lty=1)

## 7.5 ブースティング

b.dt=function(x,y,d,f="sq.loss"){
  n=nrow(x)
  if(f=="sq.loss")g=sq.loss 
  else if(f=="mis.match") g=mis.match 
  else if(f=="gini")g=gini 
  else g=entropy
  vertex=list(); vertex[[1]]=list(parent=0, set=1:n, score=g(y), j=0)
  while(length(vertex)<=2*d-1){ #
    r=length(vertex); gain.max=-Inf #
    for(h in 1:r)if(vertex[[h]]$j==0){ #
      res=branch(x,y,g,vertex[[h]]$set) #
      gain=vertex[[h]]$score-res$score #
      if(gain>gain.max){gain.max=gain; h.max=h; res.max=res} #
    } #
    vertex[[h.max]]$th=x[res.max$i,res.max$j]; vertex[[h.max]]$j=res.max$j
    vertex[[r+1]]=list(parent=h.max, set=res.max$left,
    score=res.max$left.score, j=0)
    vertex[[r+2]]=list(parent=h.max, set=res.max$right,
    score=res.max$right.score, j=0)
  }
  r=2*d+1 #
  for(h in 1:r){vertex[[h]]$left=0; vertex[[h]]$right=0}
  for(h in r:2){
    pa=vertex[[h]]$parent
    if(vertex[[pa]]$right==0)vertex[[pa]]$right=h else vertex[[pa]]$left=h
    if(vertex[[h]]$right==0&vertex[[h]]$left==0)vertex[[h]]$j=0
  }
  mode=function(y) names(sort(table(y),decreasing=TRUE))[1]
  if(f=="sq.loss") g=mean else g=mode
  for(h in 1:r)if(vertex[[h]]$j==0)vertex[[h]]$center=g(y[vertex[[h]]$set])
  return(vertex)
}
library(MASS)
df=Boston; x=as.matrix(df[,1:13]); y=as.matrix(df[,14])
train=1:200; test=201:300; B=200; lambda=0.1; d=1
trees=list(); r=y[train]
for(b in 1:B){
  trees[[b]]=b.dt(x[train,],r,d)
  for(i in train)r[i]=r[i]-lambda*value(x[i,],trees[[b]])
}
z=array(0,dim=c(B,600))
for(i in test)z[1,i]=lambda*value(x[i,],trees[[1]])
for(b in 2:B)for(i in test)z[b,i]=z[b-1,i]+lambda*value(x[i,],trees[[b]])
out=NULL; 
for(b in 1:B)out=c(out,sum((y[test]-z[b,test])^2)/length(test))
out1=out

#d=2,d=3でも実行する
d=2;trees=list(); r=y[train]
for(b in 1:B){
  trees[[b]]=b.dt(x[train,],r,d)
  for(i in train)r[i]=r[i]-lambda*value(x[i,],trees[[b]])
}
z=array(0,dim=c(B,600))
for(i in test)z[1,i]=lambda*value(x[i,],trees[[1]])
for(b in 2:B)for(i in test)z[b,i]=z[b-1,i]+lambda*value(x[i,],trees[[b]])
out=NULL; 
for(b in 1:B)out=c(out,sum((y[test]-z[b,test])^2)/length(test))
out2=out

# d=3
d=3;trees=list(); r=y[train]
for(b in 1:B){
  trees[[b]]=b.dt(x[train,],r,d)
  for(i in train)r[i]=r[i]-lambda*value(x[i,],trees[[b]])
}
z=array(0,dim=c(B,600))
for(i in test)z[1,i]=lambda*value(x[i,],trees[[1]])
for(b in 2:B)for(i in test)z[b,i]=z[b-1,i]+lambda*value(x[i,],trees[[b]])
out=NULL; 
for(b in 1:B)out=c(out,sum((y[test]-z[b,test])^2)/length(test))
out3=out
plot(21:100, xlab="生成した木の個数",ylab="テストデータでの二乗誤差", type="n",
    xlim=c(20,100),ylim=c(0,35),main="ブースティング")
lines(21:100, out1[21:100],col="red")
lines(21:100, out2[21:100],col="blue")
lines(21:100, out3[21:100],col="green")
legend("topright",legend=c("d=1","d=2","d=3"),col=c("red","blue","green"),lty=1)

