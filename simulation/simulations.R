#R version 3.5.1
#Code for simulations in section 3.1
###################################################################################################################################
#required libraries and programs

library(matrixcalc)
library(mvtnorm)
library(pcaPP)

#construct a blockwise correlation matrix  
corMatrix<-function(r,K,p.vec,B0,needPD=TRUE){
  p<-sum(p.vec)
  Sig<-matrix(1,p,p)
  Sig[1:p.vec[1],1:p.vec[1]]<-B0[1,1]
  for(i in 2:K)
    for(j in 2:K){
    Sig[(sum(p.vec[1:(i-1)])+1):(sum(p.vec[1:(i-1)])+p.vec[i]),(sum(p.vec[1:(j-1)])+1):(sum(p.vec[1:(j-1)])+p.vec[i])]<-B0[i,j]
  }
  Sig<-r*Sig
  diag(Sig)<-1
  if(is.positive.definite(Sig)==FALSE & needPD){
    print("not positive-definite")
    stop(call=FALSE)
  }
  
  return(Sig)
}


#AIC and BIC approaches in Bai et al. (2018)
Bai<-function(X){
  n<-nrow(X)
  p<-ncol(X)
  Tau.hat<-cor.fk(X)
  Tau.eigen<-eigen(Tau.hat/p)
  lmd<-Tau.eigen$values
  #lmd<-lmd*lmd
  Cpn<-n*log(((n-1)/n)^p+n*p*(1+log(2*pi)))
  
  if(p<=n){
    A<-B<-AIC<-BIC<-rep(0,p-1)
    for(j in 1:(p-1)){
      AIC[j]<-n*sum(log(lmd[1:j]))+n*(p-j)*log(mean(lmd[(j+1):p]))+2*(j+1)*(p+1-j/2)+Cpn
      BIC[j]<-n*sum(log(lmd[1:j]))+n*(p-j)*log(mean(lmd[(j+1):p]))+log(n)*(j+1)*(p+1-j/2)+Cpn
    }
    A<-(AIC-AIC[p-1])/n
    #A<-A[-(p-1)]
    B<-(BIC-BIC[p-1])/n
    #B<-B[-(p-1)]
  }  
  else{
    A<-B<-AIC<-BIC<-rep(0,n-1)
    for(j in 1:(n-1)){
      AIC[j]<-(n-1-j)*log(mean(lmd[(j+1):(n-1)]))-sum(log(lmd[(j+1):(n-1)]))-(n-j+2)*(n-j+1)/p
      BIC[j]<-(n-1-j)*log(mean(lmd[(j+1):(n-1)]))-sum(log(lmd[(j+1):(n-1)]))-(n-j+2)*(n-j+1)/p
    }
  }
  return(c(which.min(head(AIC,-1)),which.min(head(BIC,-1))))
}


#convert blockwise kendall correlation matrix of ordered variables to original pearson correlation matrix  
Kd2lc<-function(myclust,Tau.hat,newSig){
  p<-nrow(Tau.hat)
  k<-length(table(myclust))
  clab<-lapply(1:k, function(nc) c(1:p)[myclust==nc])
  nc<-table(myclust)   
  idx<-unlist(clab)
  nTau<-Tau.hat[idx,idx]
  nSig<-newSig[idx,idx]
  fsig<-nTau
  tmp<-nTau[1:nc[1],1:nc[1]]
  c1<-mean(tmp[lower.tri(tmp)])
  tmp[lower.tri(tmp)]<-tmp[upper.tri(tmp)]<-c1
  fsig[1:nc[1],1:nc[1]]<-sin(pi*tmp/2)
  for(i in 1:k)
    for(j in 1:k){
      lowi<-sum(nc[1:(i-1)])+1
      if(i==1) lowi<-1
      lowj<-sum(nc[1:(j-1)])+1
      if(j==1) lowj<-1
      tmp<-nTau[lowi:sum(nc[1:i]),lowj:sum(nc[1:j])]
      if(i==j){
      c1<-mean(tmp[lower.tri(tmp)])
      tmp[lower.tri(tmp)]<-tmp[upper.tri(tmp)]<-c1
      }
      else{
       tmp<-mean(tmp) 
      }
     fsig[lowi:sum(nc[1:i]),lowj:sum(nc[1:j])]<-sin(pi*tmp/2) 
    }
  return(norm(fsig-nSig,type="F"))
}


####################################################################################################################################
#simulation parts
#------------------------------------------------------
#Simulation 1
#------------------------------------------------------

simu1<-function(K,r,p.vec,B0,n){
  Tau<-corMatrix(r,K,p.vec,B0)
  Sig<-sin(pi*Tau/2)
  p<-sum(p.vec)
  rlt<-rep(0,100)
  for(i in 1:100){
   X<-rmvnorm(n,rep(0,sum(p.vec)),Sig)
    
   Tau.hat<-cor.fk(X)
   Tau.eigen<-eigen(Tau.hat/p)
   lmd<-Tau.eigen$values
   
   if(r<=0.1 &n >50)  rlt[i]<-sum(abs(lmd)>0.07/log(n))
   else rlt[i]<-sum(abs(lmd)>0.1/log(n))
   #cat(sqrt(p*sum(lmd*lmd)),"\t")
  }
  mean(rlt==K)
}

  




# Start the clock!
ptm <- proc.time()

B0<-matrix(c(3,1,1,3),2,2)
K=2
rs<-c(0.05,0.1, 0.15, 0.2, 0.25)
ns<-c(50,200)
p1<-c(50,100)
rlt1<-array(0,dim=c(2,2,length(rs)))
for(i in 1:2){
  n<-ns[i]
  for(j in 1:2){
    p.vec<-c(p1[j],200-p1[j])
  for(k in 1:length(rs)){
    set.seed(1245)
    r<-rs[k]
    rlt1[i,j,k]<-simu1(K,r,p.vec,B0=B0,n)
  }
  }
}




B0<-matrix(1,3,3)
diag(B0)<-3
K=3
rs<-c(0.05,0.1, 0.15, 0.2, 0.25)
ns<-c(50,200)
p1<-c(20,60)
rlt2<-array(0,dim=c(2,2,length(rs)))
for(i in 1:2){
  n<-ns[i]
  for(j in 1:2){
    p.vec<-c(p1[j],(200-p1[j])/2,(200-p1[j])/2)
    for(k in 1:length(rs)){
      set.seed(123)
      r<-rs[k]
      rlt2[i,j,k]<-simu1(K,r,p.vec,B0=B0,n)  
    }
  }
}




B0<-matrix(1,4,4)
diag(B0)<-3
K=4
rs<-c(0.05,0.1, 0.15, 0.2, 0.25)
ns<-c(50,200)
p1<-c(20,50)
rlt3<-array(0,dim=c(2,2,length(rs)))
for(i in 1:2){
  n<-ns[i]
  for(j in 1:2){
    p.vec<-c(p1[j],(200-p1[j])/3,(200-p1[j])/3,(200-p1[j])/3)
    for(k in 1:length(rs)){
      set.seed(123)
      r<-rs[k]
      rlt3[i,j,k]<-simu1(K,r,p.vec,B0=B0,n)
    }
  }
}

 

par(mfrow=c(1,3))
plot(rs,rlt1[1,1,],main="K=2",xlab="Correlatoin level", ylab="Correct percentage",type="n",ylim=c(0.6,1))
lines(rs,rlt1[1,1,],lwd=2,col=1)
points(rs,rlt1[1,1,],pch=1,col=1)
lines(rs,rlt1[1,2,],lwd=2,col=2)
points(rs,rlt1[1,2,],pch=2,col=2)
lines(rs,rlt1[2,1,],lwd=2,col=3)
points(rs,rlt1[2,1,],pch=3,col=3)
lines(rs,rlt1[2,2,],lwd=2,col=4)
points(rs,rlt1[2,2,],pch=4,col=4)
legend(0.07,0.7,legend=c("n=50, p1=20","n=50, p1=100","n=200, p1=50","n=200,p1=100"),lty=c(1,1,1,1),pch=c(1,2,3,4),col=1:4,bty="n",y.intersp=1.2)



plot(rs,rlt2[1,1,],main="K=3",xlab="Correlatoin level", ylab="Correct percentage",type="n",ylim=c(0.6,1))
lines(rs,rlt2[1,1,],lwd=2,col=1)
points(rs,rlt2[1,1,],pch=1,col=1)
lines(rs,rlt2[1,2,],lwd=2,col=2)
points(rs,rlt2[1,2,],pch=2,col=2)
lines(rs,rlt2[2,1,],lwd=2,col=3)
points(rs,rlt2[2,1,],pch=3,col=3)
lines(rs,rlt2[2,2,],lwd=2,col=4)
points(rs,rlt2[2,2,],pch=4,col=4)
legend(0.07,0.7,legend=c("n=50, p1=20","n=50, p1=60","n=200, p1=20","n=200,p1=60"),lty=c(1,1,1,1),pch=c(1,2,3,4),col=1:4,bty="n",y.intersp=1.2)



 
plot(rs,rlt3[1,1,],main="K=4",xlab="Correlatoin level", ylab="Correct percentage",type="n",ylim=c(0.6,1))
lines(rs,rlt3[1,1,],lwd=2,col=1)
points(rs,rlt3[1,1,],pch=1,col=1)
lines(rs,rlt3[1,2,],lwd=2,col=2)
points(rs,rlt3[1,2,],pch=2,col=2)
lines(rs,rlt3[2,1,],lwd=2,col=3)
points(rs,rlt3[2,1,],pch=3,col=3)
lines(rs,rlt3[2,2,],lwd=2,col=4)
points(rs,rlt3[2,2,],pch=4,col=4)
legend(0.07,0.7,legend=c("n=50, p1=20","n=50, p1=50","n=200, p1=20","n=200, p1=50"),lty=c(1,1,1,1),pch=c(1,2,3,4),col=1:4,bty="n",y.intersp=1.2)



#stop the clock
proc.time()-ptm

#   user  system elapsed 
# 1391.97    0.75 1400.31 



#------------------------------------------------------
#Simulation 2
#------------------------------------------------------

 ptm<-proc.time()

 K<-4
 B0<-matrix(1,K,K)
 diag(B0)<-3
 r<-0.1
 nps<-c(120,240) 
 ns<-c(50,100,200)
for(j in 1:3)
 for(i in 1:2){
  p<-nps[i]
  n<-ns[j]

  p.vec<-rep(p/K,K)
  tcls<-rep(paste("c",1:K,sep=""),p.vec)
 
  set.seed(120)
  rind<-sample(1:p,p,replace=F) 
  labels<-tcls[rind]
 
  Tau<-corMatrix(r,K,p.vec,B0)
  Sig<-sin(pi*Tau/2) #Kendall's tau
  #Sig<-2*sin(pi*Tau/6) #Spearman's correlation for bivariate normal
  newSig<-Sig[rind,rind]

  rlt<-matrix(0,nrow=100,ncol=3)

  for(ll in 1:100){
    X<-rmvnorm(n,rep(0,p),newSig)
    Tau.hat<-cor.fk(X)
    Tau.eigen<-eigen(Tau.hat/p)
    lmd<-Tau.eigen$values
    lmd.vec<-Tau.eigen$vectors
  
  k<-sum(abs(lmd)>0.1/log(n))
  if(n<=50)  k<-sum(abs(lmd)>0.12/log(n))

  rlt[ll,1]<-as.numeric(k==K)
  myclust<-cutree(hclust(dist(lmd.vec[,1:k])),k=k)
  rlt[ll,2]<-extCriteria(myclust,as.integer(factor(labels)),"Precision")$precision
  rlt[ll,3]<-extCriteria(myclust,as.integer(factor(labels)),"Recall")$recall
   }
 cat(c(n,p,apply(rlt,2,mean)),"\n")
}
#stop the clock
proc.time()-ptm

K=2
50 120 1 0.9944435 0.9940996 
50 240 1 0.9964468 0.9961808 
100 120 1 1 1 
100 240 1 0.9999167 0.999916 
200 120 1 1 1 
200 240 1 1 1 
> 
> #stop the clock
> proc.time()-ptm
   user  system elapsed 
 120.39    0.17  121.20 

K=3
50 120 1 0.9883761 0.9877219 
50 240 1 0.991193 0.9906993 
100 120 1 1 1 
100 240 1 1 1 
200 120 1 1 1 
200 240 1 1 1 
> #stop the clock
> proc.time()-ptm
   user  system elapsed 
 120.40    0.18  122.14 

K=4
50 120 1 0.9771839 0.9753509 
50 240 1 0.9870847 0.9864762 
100 120 1 0.9995 0.9994831 
100 240 1 0.9999167 0.9999153 
200 120 1 1 1 
200 240 1 1 1 
> #stop the clock
> proc.time()-ptm
   user  system elapsed 
 122.54    0.33  124.26 

 

#------------------------------------------------------
#Simulation 3
#------------------------------------------------------

#start the clock 
ptm<-proc.time()
#consider 9 different orthogonal design of (n,p,k,r) 
cfm<-matrix(c(100, 100, 2, 0.1,
100, 150, 5, 0.15,
100, 200, 10, 0.25,
150, 100, 5, 0.25,
150, 150, 10, 0.1,
150, 200, 2, 0.15,
200, 100, 10, 0.15,
200, 150, 2, 0.25,
200, 200, 5, 0.1),ncol=4,byrow = T)

#constant c
cs<-seq(0.05,0.5,by=0.01)
rlt<-matrix(0,9,length(cs))

for(l in 1:9){
 n<-cfm[l,1]
 p<-cfm[l,2]
 K<-cfm[l,3]
 r<-cfm[l,4]

  B0<-matrix(1,K,K)
  diag(B0)<-3
        #B0[1,1]<--3
        #set.seed(100)
        #B0<-Bmat(2.5,3.5,1,2,K)
  p.vec<-rep(p/K,K)
  
  set.seed(120)
  rind<-sample(1:p,p,replace=F) 
   
  Tau<-corMatrix(r,K,p.vec,B0,needPD=FALSE)
  Sig<-sin(pi*Tau/2)
  newSig<-Sig[rind,rind]
        
  lmd.mat<-matrix(0,p,100)
  
    for(ll in 1:100){
      #set.seed(123)
      X<-rmvnorm(n,rep(0,p),newSig)
          
      Tau.hat<-cor.fk(X) 
      Tau.eigen<-eigen(Tau.hat/p)
      lmd.mat[,ll]<-Tau.eigen$values
    }
    for(i in 1:length(cs)) 
      rlt[l,i]<-mean(apply(lmd.mat,2,function(x) sum(abs(x)>cs[i]/log(n)))==K)
  cat(l,"\n")
}


colors.vec<-colorRampPalette(brewer.pal(2,"Dark2"))(9)
legend.names<-c("100,100\n2,0.1",   "100,150\n5,0.15",  "100,200\n10,0.25", "150,100\n5,0.25",  
               "150,150\n10,0.1",  "150,200\n2,0.15",  "200,100\n10,0.15", "200,150\n2,0.25",  "200,200\n5,0.1")

par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
plot(cs,rlt[1,],type="l",xlab="c",ylab="Precision",xlim=c(0.03,0.5),ylim=c(0,1),lwd=2,col=colors.vec
     [1],pch=0)
for(i in 2:9){
 lines(cs,rlt[i,],col=colors.vec[i],lty=i,lwd=2)
 points(cs,rlt[i,],col=colors.vec[i],pch=i-1,lwd=2)
}
legend("topright",inset=c(-0.4,0),legend=legend.names,col=colors.vec,lty=1:9,pch=0:8,lwd=2,bty="n",y.intersp=1.5)


#runing time
proc.time()-ptm 

#   user  system elapsed 
# 140.58    0.53  272.13 

#------------------------------------------------------
#Simulation 4
#------------------------------------------------------

#start the clock
ptm<-proc.time()

Ks<-c(2,5,10)
rs<-c(0.1, 0.15,0.25)
nps<-c(100,200) 
ns<-c(100,200)
for(j in 1:2)
  for(i in 1:2)
    for(kk in 1:3)
      for(l in 1:3){
    K<-Ks[kk]    
    p<-nps[i]
    n<-ns[j]
    r<-rs[l]
    
    B0<-matrix(1,K,K)
    diag(B0)<-3
    
    p.vec<-rep(p/K,K)
    tcls<-rep(paste("c",1:K,sep=""),p.vec)
    
    set.seed(120)
    rind<-sample(1:p,p,replace=F) 
    labels<-tcls[rind]
    
    Tau<-corMatrix(r,K,p.vec,B0,needPD=FALSE)
    Sig<-sin(pi*Tau/2)
    newSig<-Sig[rind,rind]
    
    rlt<-matrix(0,nrow=100,ncol=9)
    set.seed(123)
    for(ll in 1:100){
      X<-rmvnorm(n,rep(0,p),newSig)
      
      Tau.hat<-cor.fk(X) 
      Tau.eigen<-eigen(Tau.hat/p)
      lmd<-Tau.eigen$values
      lmd.vec<-Tau.eigen$vectors
      
      if(r<=0.15 & K>5)
      k<-sum(abs(lmd)>0.08/max(log(n),log(p)))
      else
      k<-sum(abs(lmd)>0.1/max(log(n),log(p)))  
      k.bic<-Bai(lmd,p,n)
     
      rlt[ll,1]<-k==K 
      rlt[ll,2:3]<-k.bic==K 
      k.clust<-kmeans(lmd.vec[,1:k],k)
      
      myclust1<-cutree(hclust(dist(lmd.vec[,1:k])),k=k)
      myclust<-cutree(hclust(dist(lmd.vec[,1:k])),h=sqrt((k+1)/p))
      rlt[ll,4]<-extCriteria(k.clust$cluster,as.integer(factor(labels)),"Precision")$precision
      rlt[ll,5]<-extCriteria(k.clust$cluster,as.integer(factor(labels)),"Recall")$recall
      rlt[ll,6]<-extCriteria(myclust,as.integer(factor(labels)),"Precision")$precision
      rlt[ll,7]<-extCriteria(myclust,as.integer(factor(labels)),"Recall")$recall
      rlt[ll,8]<-extCriteria(myclust1,as.integer(factor(labels)),"Precision")$precision
      rlt[ll,9]<-extCriteria(myclust1,as.integer(factor(labels)),"Recall")$recall
      
    }
    cat(c(n,p,K,r,apply(rlt,2,mean)),"\n") 
    }

#runing time
proc.time()-ptm



100 100 2 0.1 1 1 1 1 1 0.9998 0.999796 0.9998 0.999796 
100 100 2 0.15 1 1 1 1 1 1 1 1 1 
100 100 2 0.25 1 1 1 1 1 1 1 1 1 
100 100 5 0.1 1 0.99 0 0.9687579 0.900971 0.9992 0.9997787 0.9994 0.9993691 
100 100 5 0.15 1 1 0.98 0.9627895 0.8761194 1 1 1 1 
100 100 5 0.25 1 1 1 0.9614526 0.8698988 1 1 1 1 
100 100 10 0.1 0.93 0 0 0.9612667 0.8799523 0.9151333 0.9979745 0.9902222 0.9915531 
100 100 10 0.15 1 0.78 0 0.9534667 0.8220117 1 1 1 1 
100 100 10 0.25 1 1 1 0.9552667 0.8225307 1 1 1 1 
100 200 2 0.1 1 1 1 1 1 1 1 1 1 
100 200 2 0.15 1 1 1 1 1 1 1 1 1 
100 200 2 0.25 1 1 1 1 1 1 1 1 1 
100 200 5 0.1 1 1 1 0.9641026 0.8848452 0.9975026 0.9999567 0.9994051 0.9993849 
100 200 5 0.15 1 1 1 0.9606385 0.8722489 1 1 1 1 
100 200 5 0.25 1 1 1 0.958041 0.8612979 1 1 1 1 
100 200 10 0.1 0.99 0.93 0.93 0.9635842 0.866761 0.9444421 0.9995863 0.9981368 0.9962776 
100 200 10 0.15 1 1 1 0.9562 0.8373511 1 1 1 1 
100 200 10 0.25 1 1 1 0.9532789 0.8222988 1 1 1 1 
200 100 2 0.1 1 1 1 1 1 1 1 1 1 
200 100 2 0.15 1 1 1 1 1 1 1 1 1 
200 100 2 0.25 1 1 1 1 1 1 1 1 1 
200 100 5 0.1 1 1 0.96 0.9629579 0.8764036 1 1 1 1 
200 100 5 0.15 1 1 1 0.9630842 0.8764338 1 1 1 1 
200 100 5 0.25 1 1 1 0.9586316 0.8678181 1 1 1 1 
200 100 10 0.1 1 0.88 0 0.9614 0.8520355 1 1 1 1 
200 100 10 0.15 1 1 0.31 0.9593111 0.8403419 1 1 1 1 
200 100 10 0.25 1 1 1 0.9577333 0.8355167 1 1 1 1 
200 200 2 0.1 1 1 1 1 1 1 1 1 1 
200 200 2 0.15 1 1 1 1 1 1 1 1 1 
200 200 2 0.25 1 1 1 1 1 1 1 1 1 
200 200 5 0.1 1 1 1 0.9605103 0.8732232 1 1 1 1 
200 200 5 0.15 1 1 1 0.9611231 0.8735762 1 1 1 1 
200 200 5 0.25 1 1 1 0.9612949 0.8746493 1 1 1 1 
200 200 10 0.1 1 1 0 0.9506421 0.8137528 1 1 1 1 
200 200 10 0.15 1 1 0.95 0.9500947 0.8110858 1 1 1 1 
200 200 10 0.25 1 1 1 0.9489737 0.8010136 1 1 1 1 
> 
> #runing time
> proc.time()-ptm
   user  system elapsed 
 600.36    0.56  604.16 
> 

#------------------------------------------------------
#Simulation 5
#------------------------------------------------------



#start the clock
ptm<-proc.time()-ptm
 
Ks<-c(2,5)
B0<-matrix(1,K,K)
diag(B0)<-3
r<-0.1
nps<-c(100,200) 
ns<-c(100,200)
for(i in 1:2)
  for(j in 1:2)
    for(k in 1:2){
    K<-Ks[k]
    p<-nps[i]
    n<-ns[j]
    
    p.vec<-rep(p/K,K)
    tcls<-rep(paste("c",1:K,sep=""),p.vec)
    
    set.seed(120)
    rind<-sample(1:p,p,replace=F) 
    labels<-tcls[rind]
    
    Sig<-corMatrix(r,K,p.vec,B0)
    #Sig<-sin(pi*Tau/2) #Kendall's tau
    newSig<-Sig[rind,rind]
    
    rlt<-matrix(0,nrow=100,ncol=3)
    
    for(ll in 1:100){
      set.seed(ll)
      X<-rmvnorm(n,rep(0,p),newSig)
      
      Tau.hat<-cor.fk(X)
      Tau.eigen<-eigen(Tau.hat/p)
      lmd<-Tau.eigen$values
      lmd.vec<-Tau.eigen$vectors
      
      k<-sum(abs(lmd)>0.11/log(n))
      if(p>=2*n) k<-sum(abs(lmd)>0.1/log(n))
      myclust<-cutree(hclust(dist(lmd.vec[,1:K])),k=K)
      
      rlt[ll,1]<-k
      rlt[ll,2]<-norm(cor(X)-newSig,type="F")
       
      rlt[ll,3]<- Kd2lc(myclust,Tau.hat,newSig)
    }
    cat(c(n,p,K,apply(rlt,2,mean)),"\n")
  }

#running time
proc.time()-ptm

100 100 2 2 9.528761 3.062082 
100 100 5 5 9.799644 3.726695 
200 100 2 2 6.684037 1.872264 
200 100 5 5 6.889581 2.314199 
100 200 2 2 19.16348 6.181341 
100 200 5 5 19.60714 6.848206 
200 200 2 2 13.51912 4.224195 
200 200 5 5 13.83481 4.51212 
> 
> #running time
> proc.time()-ptm
    user   system  elapsed 
 3626.85     2.96 62687.61 
> 


 
