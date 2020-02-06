rm(list=ls())
px.theta<-function(x,theta){ 1/sqrt(2*pi)*exp(-(x-theta)^2/2) }
p.theta<-function(theta){ 1/(pi*(1+theta^2)) }
px<-function(theta){ px.theta(x,theta)*p.theta(theta) }         #plot(seq(-10,10,0.01),px(seq(-10,10,0.01)))  #plot p(x) 
post<-function(theta){ px.theta(x,theta)*p.theta(theta)/I_Spx }  #p(theta|x)      (posterior)  #x<-rnorm(1,5,1);plot(seq(-10,10,0.01),post(seq(-10,10,0.01)))#plot p(theta|x)
#f<-function(theta){  theta*post(theta) }                         #theta*p(theta|x)
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
#estmator 
E<-rep(0,100)
for(k in 1:100){
  y<-c()
  for(j in -1000:1000){
    theta<-j
    x<-rnorm(1,theta,1)
    a<-theta-10 ;b<-theta+10 ;n<-100 
    h<-seq(a,b,(b-a)/n)
    z1<-0;z2<-0
    for(i in 2:(n/2)){   z1<-z1+2*px(h[2*i-1])  }
    for(i in 1:(n/2)){   z2<-z2+4*px(h[2*i])  }
    I_Spx<-( px(h[1])+px(h[n+1])+z1+z2 )*(b-a)/n/3
    y<-c(y,theta*post(theta))
    }
  E[k]<-sum(y)
}
mean(E)
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




