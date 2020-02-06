rm(list = ls())
f<-function(x){x^(-1/3)+x/10}
a<-0;b<-1;n<-1000;N<-100
I<-integrate(f,0,1)$value
#Monte Carlo method with U(0,1)
F_u<-function(n){  (b-a)/n*sum(f(runif(n,0,1))) }         #f(runif(n,0,1))%*%rep(1,n)
I_U<-F_u(n)
confi.of.I_U<-function(N,alpha){
  x<-rep(0,N)
  for(i in 1:N){  x[i]<-x[i]+F_u(n) }
  I.bar<-sum(x)/N
  I.sd<-sd(x)
  cat("The 95% confidence interval for estimating I is [",
      I.bar-I.sd/sqrt(N)*qnorm(1-alpha),"",
      I.bar+I.sd/sqrt(N)*qnorm(1-alpha),"]","\n")
}
confi.of.I_U(100,0.05)

####################
# confidence of estimator
# confi.est<-function(N,alpha){
#   x<-rep(0,N)
#   for(i in 1:N){  x[i]<-x[i]+ }
#   I.bar<-sum(x)/N
#   I.sd<-sd(x)
#   cat("The 95% confidence interval for estimating I is [",
#       I.bar-I.sd/sqrt(N)*qnorm(1-alpha),"",
#       I.bar+I.sd/sqrt(N)*qnorm(1-alpha),"]","\n")
# }
####################
#p_1(x) ~ U(0,1)
var1<-(3+3/25+1/300-(31/20)^2)/n
x1<-rep(0,N)
for(i in 1:N){
  x1[i]<-x1[i]+F_u(n)
}
var.hat1<-sd(x1)^2;rm(x1)
#p_2(x)=2/3*x^(-1/3)
var2<-(9/4+3/20+9/2000-(31/20)^2)/n
F_2<-function(n){ sum(3/2+3/20*(runif(n,0,1))^(4/3))/n }
I_2<-F_2(n)
x2<-rep(0,N)
for(i in 1:N){
  x2[i]<-x2[i]+F_2(n)
}
var.hat2<-sd(x2)^2;rm(x2)




