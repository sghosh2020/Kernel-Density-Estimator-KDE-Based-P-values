
skew<- function(u){
  return(mean((u-mean(u))^3))
         }


kurto<- function(u){
return(mean((x-mean(x))^4)-3*var(u))
}


KDEPV<- function(x,y){
  d<-dim(x)[[2]]
  n1<-dim(x)[1]
  n2<-dim(y)[1]
  n<-n1+n2
  la1<-n1/(n1+n2)
  la2<-n2/(n1+n2)
  sx<-((n1-1)/n1)*apply(x,2,var)
  sy<-((n2-2)/n1)*apply(y,2,var)
  gx.hat<-apply(x,2,skew)
  gy.hat<-apply(y,2,skew)
  kx<-apply(x,2,kurto)
  ky<-apply(y,2,kurto)

  a1<-sx/la1 +sy/la2
  a2<-gx.hat/la1^2
  a3<-gy.hat/la2^2
  a4<-kx/la1^3
  a5<-ky/la2^3
  a6<-sx^{2}/la1^3+sy^{2}/la2^3
  a7<-(sx*sy)/(la1^{2}*la2^{2})

  a6<-sx^{2}/la1^3+sy^{2}/la2^3
  a7<-(sx*sy)/(la1^{2}*la2^{2})

  t1<-sqrt(n)*(mean(x)-mean(y))/sqrt(a1)
  p1<-2*pnorm(abs(t1),0,1,lower.tail=FALSE)
  pr<-min(pnorm(abs(t1))+6^{-1}*n^{-1/2}*(a1^{-3/2}*(a2-a3)*(2*t1^{2}+1))*dnorm(t1),0.99999)
  w<-qnorm(pr,0,1)
  k1<- -(12^{-1}*a1^{-2}*(a4+a5)*(w^{2}-3)-18^{-1}*a1^{-3}*(a2-a3)^{2}*(w^{4}+2*w^{2}-3)-4^{-1}*a1^{-2}*(w^{2}*a6+3*a6+2*a7))
  k2<-  (12^{-1}*a1^{-2}*(a4+a5)*(w^{2}-3)-18^{-1}*a1^{-3}*(a2-a3)^{2}*(w^{4}+2*w^{2}-3)-4^{-1}*a1^{-2}*(w^{2}*a6+3*a6+2*a7))

h.hat<- apply(cbind(n^{-1/2}*(2*a1*k1*(1/la1+1/la2)^{-1})^{1/2},n^{-1/2}*(2*a1*k2*(1/la1+1/la2)^{-1})^{1/2}),1, na.omit)
t2<- sqrt(n)*(apply(x,2,mean)-apply(y,2,mean))/sqrt(a1+h.hat*(1/la1+1/la2))
p<-2*pnorm(abs(t2),0,1,F)
adj.pvalue<-cbind(p.adjust(cbind(p, seq(1,d,1))[order(cbind(p, seq(1,d,1))[,1]),][,1],"BH", d),cbind(p, seq(1,d,1))[order(cbind(p, seq(1,d,1))[,1]),][,2])
colnames(adj.pvalue) <- c("BH adjusted pvalues", "Hypothesis")
return(adj.pvalue)
}






