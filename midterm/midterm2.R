rm(list=ls())
px.theta<-function(x,theta){ 1/sqrt(2*pi)*exp(-(x-theta)^2/2) }
p.theta<-function(theta){ 1/(pi*(1+theta^2)) }
px<-function(theta){ px.theta(x,theta)*p.theta(theta) }         #plot(seq(-10,10,0.01),px(seq(-10,10,0.01)))  #plot p(x) 
post<-function(theta){ px.theta(x,theta)*p.theta(theta)/I_Spx }  #p(theta|x)      (posterior)  #x<-rnorm(1,5,1);plot(seq(-10,10,0.01),post(seq(-10,10,0.01)))#plot p(theta|x)
f<-function(theta){  theta*post(theta) }                         #theta*p(theta|x)
#plot p(theta|x)
theta<-5 ;x<-rnorm(1,theta,1)
a<-theta-10 ;b<-theta+10 ;n<-100 
h<-seq(a,b,(b-a)/n)
if(n%%2==0){
  z1<-0;z2<-0
  for(i in 2:(n/2)){   z1<-z1+2*px(h[2*i-1])  }
  for(i in 1:(n/2)){   z2<-z2+4*px(h[2*i])  }
  I_Spx<-( px(h[1])+px(h[n+1])+z1+z2 )*(b-a)/n/3
  rm(z1,z2)
} else { print("n need to be even") }
plot(seq(-5,10,0.01),post(seq(-5,10,0.01))) ;lines(c(theta,theta),c(-0.1,0.6));lines(c(seq(-5,10,0.01)[which(post(seq(-5,10,0.01))==max(post(seq(-5,10,0.01))))],seq(-5,10,0.01)[which(post(seq(-5,10,0.01))==max(post(seq(-5,10,0.01))))]),c(-0.1,max(post(seq(-5,10,0.01)))) )
#estmator with theta=5
N=10000 ;theta=5
y<-rep(0,N)
for(j in 1:N){
  x<-rnorm(1,theta,1)
  a<-theta-10 ;b<-theta+10 ;n<-100 
  h<-seq(a,b,(b-a)/n)
  if(n%%2==0){
    z1<-0;z2<-0
    for(i in 2:(n/2)){   z1<-z1+2*px(h[2*i-1])  }
    for(i in 1:(n/2)){   z2<-z2+4*px(h[2*i])  }
    I_Spx<-( px(h[1])+px(h[n+1])+z1+z2 )*(b-a)/n/3
    rm(z1,z2)
  } else { print("n need to be even") }
  if(n%%2==0){
    z1<-0;z2<-0
    for(i in 2:(n/2)){   z1<-z1+2*f(h[2*i-1])  }
    for(i in 1:(n/2)){   z2<-z2+4*f(h[2*i])  }
    I_Sposterior<-( f(h[1])+f(h[n+1])+z1+z2 )*(b-a)/n/3
    rm(z1,z2)
  } else { print("n need to be even") }
  y[j]<-I_Sposterior
};rm(i,I_Sposterior)
I_bar<-sum(y)/N
I_var<-(sum(y^2)-N*I_bar^2)/(N-1)
list(theta=theta,I_bar=I_bar,I_var=I_var)

#use simpsons to est px's integral
a<--5 ;b<-10 ;n<-100 
h<-seq(a,b,(b-a)/n)
if(n%%2==0){
  z1<-0;z2<-0
  for(i in 2:(n/2)){   z1<-z1+2*px(h[2*i-1])  }
  for(i in 1:(n/2)){   z2<-z2+4*px(h[2*i])  }
  I_Spx<-( px(h[1])+px(h[n+1])+z1+z2 )*(b-a)/n/3
  rm(z1,z2)
} else { print("n need to be even") }          #integrate(px,-Inf,Inf)
#use simpsons to est [theta*p(theta|x)]'s integral
a<--5 ;b<-10 ;n<-1000 
h<-seq(a,b,(b-a)/n)
if(n%%2==0){
  z1<-0;z2<-0
  for(i in 2:(n/2)){   z1<-z1+2*f(h[2*i-1])  }
  for(i in 1:(n/2)){   z2<-z2+4*f(h[2*i])  }
  I_Sposterior<-( f(h[1])+f(h[n+1])+z1+z2 )*(b-a)/n/3
  rm(z1,z2)
} else { print("n need to be even") }        #integrate(f,-Inf,Inf)  

#######################
f.1<-function(theta){ 
  a<-1#; n=50;theta=1
  for(i in 1:20){
    x<-rnorm(1,theta,1)
    a<-a*px.theta(x,theta)
  }
  a*p.theta(theta)
}
integrate(f.1,-Inf,Inf)
f.2<-function(theta){
  a<-1#; n=50;theta=1
  for(i in 1:20){
    x<-rnorm(1,theta,1)
    a<-a*px.theta(x,theta)
  }
  theta*a*p.theta(theta)
}
integrate(f.2,-Inf,Inf)$value/integrate(f.1,-Inf,Inf)$value

