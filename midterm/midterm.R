rm(list=ls())
f_z<-function(z){ (1+( (z-0.25)/0.5 )^2 )^-1.5}
C.true<-integrate(f_z,0,1)$value
I_true<-integrate(f_z,0,0.7)$value /C.true
plot(seq(0,1,0.001), f_z(seq(0,1,0.001))/C )
par(mfrow=c(1,1))     
#est. C
C<-( max(f_z(seq(0,1,0.001)) )+min(f_z(seq(0,1,0.001)) ) )/2 
plot(c(seq(0,1,0.001),seq(0,1,0.001) ) ,c( f_z(seq(0,1,0.001))/C , f_z(seq(0,1,0.001))/C.true ),
     col= c(rep(3,1001),rep(1,1001)) )
#conv. rate
secdiff.f_z<-function(z){  60/C*( 1+(2*z-0.5)^2 )^(-3.5) *(2*z-0.5)^2 -12/C*( 1+(2*z-0.5)^2 )^(-2.5) }
plot(seq(0,1,0.001),secdiff.f_z(seq(0,1,0.001)));lines(c(0.6,0.6),c(-22,5))
n=sqrt( (0.7)^3*secdiff.f_z(0.6)/12* (10^5) ) 
########################
#Quadrature integration
a<-0 ;b<-0.7 ;n<-100 
h<-seq(a,b,(b-a)/n)
#Rectangle rule
I_upper<-0;I_lower<-0
for(i in 1:n){
  I_upper<-I_upper + f_z(h[i+1])*(b-a)/n/C
  I_lower<-I_lower + f_z(h[i])*(b-a)/n/C 
} 

a<-0 ;b<-0.7 ;n<-0
I_upper<-0
while(abs(I_upper-I_true)>=10^-4){
  cat("\n","Now is",n,"-th iteration ,bias is",abs(I_upper-I_true),"\n")
  n<-n+1 
  h<-seq(a,b,(b-a)/n)
  I_upper<-0
  for(i in 1:n){  I_upper<-I_upper + f_z(h[i+1])*(b-a)/n/C }
}
cat("\n","The sample size is",n," ,bias is",abs(I_upper-I_true),"\n")

a<-0 ;b<-0.7 ;n<-0
I_lower<-0
while(abs(I_lower-I_true)>=10^-4){
  cat("\n","Now is",n,"-th iteration ,bias is",abs(I_lower-I_true),"\n")
  n<-n+1 
  h<-seq(a,b,(b-a)/n)
  I_lower<-0
  for(i in 1:n){  I_lower<-I_lower + f_z(h[i])*(b-a)/n/C }
}
cat("\n","The sample size is",n," ,bias is",abs(I_lower-I_true),"\n")

#Trapezoidal rule
I_T<-0
for(i in 1:n ){  I_T<-I_T+( f_z(h[i])+f_z(h[i+1]) )/2  }
I_T<-I_T*(b-a)/n/C

a<-0 ;b<-0.7 ;n<-0
I_T<-0
while(abs(I_T-I_true)>=10^-4){
  cat("\n","Now is",n,"-th iteration ,bias is",abs(I_T-I_true),"\n")
  n<-n+1 
  h<-seq(a,b,(b-a)/n)
  I_T<-0
  for(i in 1:n){  I_T<-I_T+( f_z(h[i])+f_z(h[i+1]) )/2 }
  I_T<-I_T*(b-a)/n/C
}
cat("\n","The sample size is",n," ,bias is",abs(I_T-I_true),"\n")

#Simpsons rule
if(n%%2==0){
  z1<-0;z2<-0
  for(i in 2:(n/2)){   z1<-z1+2*f_z(h[2*i-1])  }
  for(i in 1:(n/2)){   z2<-z2+4*f_z(h[2*i])  }
  I_S<-( f_z(h[1])+f_z(h[n+1])+z1+z2 )*(b-a)/n/3/C
  rm(z1,z2)
} else { print("n need to be even") }

a<-0 ;b<-0.7 ;n<-0
I_S<-0
while(abs(I_S-I_true)>=10^-4){
  cat("\n","Now is",n,"-th iteration ,bias is",abs(I_S-I_true),"\n")
  n<-n+2 
  h<-seq(a,b,(b-a)/n)
  I_S<-0
  z1<-0;z2<-0
  for(i in 2:(n/2)){   z1<-z1+2*f_z(h[2*i-1])  }
  for(i in 1:(n/2)){   z2<-z2+4*f_z(h[2*i])  }
  I_S<-( f_z(h[1])+f_z(h[n+1])+z1+z2 )*(b-a)/n/3/C
};rm(z1,z2)
cat("\n","The sample size is",n," ,bias is",abs(I_S-I_true),"\n")

#Monte Carlo method with U(0,0.7)
G<-function(z,n){
  z<-c(z,rep(0,n))
  for(i in 1:n){z[i+1]<-(16807*z[i])%%(2^31-1)}
  u<-z[-1]/(2^31-1)
  return(u)
}
confi.of.I_M<-function(n,N,alpha){
  #n=100 ;seed<-1;N=50 ;alpha=0.05
  x<-rep(0,N) ;seed<-1
  for(i in 1:N){  
    u<-rep(0,n) 
    for(j in 1:n){
      u[j]<-f_z(0.7*G(seed,1))/(10/7)
      seed<-G(seed,1)*(2^31-1) }
    x[i]<-x[i]+sum(u)/n/C              #I_M<-sum(u)/n/C
  }
  
  I.bar<-sum(x)/N
  I.sd<-sd(x)
  cat("The 95% confidence interval for estimating I is [",
      I.bar-I.sd/sqrt(N)*qnorm(1-alpha),"",
      I.bar+I.sd/sqrt(N)*qnorm(1-alpha),"]","\n",I.bar)
}
confi.of.I_M(6,50,0.05)
seed<-2
I_M<-0 ;n<-0 
while(abs(I_M-I_true)>=10^-4){
  cat("\n","Now is",n,"-th iteration ,bias is",abs(I_M-I_true),"\n")
  n<-n+1
  u<-rep(0,n)
  for(j in 1:n){
    u[j]<-f_z(0.7*G(seed,1))/(10/7)
    seed<-G(seed,1)*(2^31-1) }
  I_M<-sum(u)/n/C
}
cat("\n","The sample size is",n," ,bias is",abs(I_M-I_true),"\n")

#Acce. Rejection
windows();plot(x=c(seq(0,1,0.001),seq(0,1,0.001)),y=c(f_z(seq(0,1,0.001))/C,1.8*dunif(seq(0,1,0.001),0,1)),
               col =c(rep(3,1001),rep(2,1001)) )
c<-1.8 ;C*c
envo<-function(x,c){  f_z(x)/C / c }
acc.rej.exp<-function(n,c=1.8){
  u1<-runif(n,0,1)
  Y<-u1
  u2<-runif(n,0,1)
  X<-rep(0,n)
  N<-length(which(u2 <= envo(Y,c) ))
  for(i in 1:N){   X[i]<-Y[which(u2 <=envo(Y,c))][i]    }    #exp(-(Y-1)^2/2)
  while(N<n){
    uu1<-runif(n-N,0,1)
    Y<-uu1
    uu2<-runif(n-N,0,1)
    if(length(which(uu2<=envo(Y,c)))>0){
      for(i in 1:length(which(uu2<=envo(Y,c)))){ X[N+i]<-Y[which(uu2<=envo(Y,c))][i] }}
    N<-N+length(which(uu2 <= envo(Y,c)))
  }
  return(X)
}
PDF<-function(a,b,n,band){
  #a=0;b=1;n=1000;band=0.2
  x<-a
  pdf<-rep(0,length(seq(x,b,0.005)))
  j<-1
  while(a<=b){
    I<-0
    for(i in 1:n){ if( abs(a-acc.rej.exp(1))<=band){ I<-I+1  }  }
    pdf[j]<-I/(2*n*band)
    a<-a+0.005
    j<-j+1
  }
  windows();plot(seq(x,b,0.005),pdf) #;lines(c(0,0),c(0,1)) ;lines(c(-x,b),c(0.5,0.5))
  return(pdf)
}
PDF(0,1,5000,0.2)
# F(0.7)
cdf<-0 ;n<-1
while(abs(cdf-I_true)>=10^-4){
  cat("\n","bias is",abs(cdf-I_true),"\n")
  if(n<10^5){n<-n*10}
  I<-0
  for(i in 1:n){ if(acc.rej.exp(1)<=0.7){ I<-I+1  }  }
  cdf<-I/n
}
cat("\n","The sample size is",n," ,bias is",abs(cdf-I_true),"\n") 

