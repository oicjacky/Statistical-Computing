#est. E(X^4) by monte carlo 
f<-function(x){ 1/sqrt(2*pi)*(x^4)*exp(-x^2/2) }
I<-integrate(f,-Inf,Inf)$value
plot(x=seq(0.01,5,0.01),y=f(seq(0.01,5,0.01)))

#with U(0,1)
n<-1e04;b<-10 #integrate(f,10,Inf)
I_Muni<-2*sum(f(runif(n,0,b)))*b/n
confi.of.I_Muni<-function(N,alpha){
  x<-rep(0,N)
  for(i in 1:N){  x[i]<-x[i]+2*sum(f(runif(n,0,b)))*b/n }
  I.bar<-sum(x)/N
  I.sd<-sd(x)
  cat("The 95% confidence interval for estimating I is [",
      I.bar-I.sd/sqrt(N)*qnorm(1-alpha),"",
      I.bar+I.sd/sqrt(N)*qnorm(1-alpha),"]","\n",
      I.bar)
}
confi.of.I_Muni(100,0.05)

#with exp(1)
n<-1e04;lambda<-1
Exp<-function(n,lambda){
  X<--log(runif(n,0,1))/lambda
  return(X)
}
fex<-function(x){f(x)/exp(-x) }
I_Mexp<-2*sum(fex(Exp(n,1)))/n
confi.of.I_Mexp<-function(N,alpha){
  x<-rep(0,N)
  for(i in 1:N){  x[i]<-x[i]+2*sum(fex(rexp(n,1)))/n }
  I.bar<-sum(x)/N
  I.sd<-sd(x)
  cat("The 95% confidence interval for estimating I is [",
      I.bar-I.sd/sqrt(N)*qnorm(1-alpha),"",
      I.bar+I.sd/sqrt(N)*qnorm(1-alpha),"]","\n",I.bar)
}
confi.of.I_Mexp(100,0.05)


#acceptance rejection of N(0,1) with g(x)~exp(1)
plot(x=c(seq(0,3,0.001),seq(0,3,0.001)),y=c(dnorm(seq(0,3,0.001),0,1),dunif(seq(0,3,0.001),0,2)),
     col =c(rep(1,3001),rep(2,3001)) )
windows();plot(x=c(seq(0,3,0.001),seq(0,3,0.001)),y=c(dnorm(seq(0,3,0.001),0,1),2*dexp(seq(0,3,0.001),1)),
     col =c(rep(1,3001),rep(2,3001)) )
n<-1e04
envo<-function(x,c){  (sqrt(2/pi)*exp(-x^2/2)) / (c*exp(-x)) }
acc.rej.exp<-function(n,c){
  u1<-runif(n,0,1)
  Y<--log(u1)
  u2<-runif(n,0,1)
  X<-rep(0,n)
  N<-length(which(u2 <= envo(Y,c) ))
  for(i in 1:N){   X[i]<-Y[which(u2 <=envo(Y,c))][i]    }    #exp(-(Y-1)^2/2)
  while(N<n){
    uu1<-runif(n-N,0,1)
    Y<--log(uu1)
    uu2<-runif(n-N,0,1)
    if(length(which(uu2<=envo(Y,c)))>0){
      for(i in 1:length(which(uu2<=envo(Y,c)))){ X[N+i]<-Y[which(uu2<=envo(Y,c))][i] }}
    N<-N+length(which(uu2 <= envo(Y,c)))
  }
  u3<-runif(n,0,1)
  X[which(u3<=0.5)]<--X[which(u3<=0.5)]
  return(X)
}
A<-acc.rej.exp(10000,2)
plot(A,dnorm(A,0,1)) ; lines(x=c(0,0),y=c(0,0.4))
windows();qqnorm(A);qqline(A)
#with N(0,1)
fnor<-function(x){  f(x) / (1/sqrt(2*pi)*exp(-x^2/2)) }
I_Mnorm<-sum(fnor(acc.rej.exp(n,2)))/n
confi.of.I_Mnorm<-function(N,alpha){
  x<-rep(0,N)
  for(i in 1:N){  x[i]<-x[i]+sum(fnor(acc.rej.exp(n,2)))/n }
  I.bar<-sum(x)/N
  I.sd<-sd(x)
  cat("The 95% confidence interval for estimating I is [",
      I.bar-I.sd/sqrt(N)*qnorm(1-alpha),"",
      I.bar+I.sd/sqrt(N)*qnorm(1-alpha),"]","\n",I.bar)
}
confi.of.I_Mnorm(100,0.05)

