rm(list = ls())
f<-function(x){exp(-x^2)}
a<-0;b<-1;n<-100
I<-integrate(f,0,1)$value
#Monte Carlo method
F_u<-function(n){  (b-a)/n*sum(f(runif(n,0,1))) }         #f(runif(n,0,1))%*%rep(1,n)
I_U<-F_u(n)

N<-50 ;x<-c()
for(i in 1:N){  x<-c(x,F_u(n)) }
I_U.bar<-mean(x) 
I_U.var<-sd(x)^2
bias<-abs(I_U.bar-I)
rm(x)