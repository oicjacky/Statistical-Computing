theta_0<- 0.5
theta<- 0
theta_2<- theta_0
n<- 0   # by increasing sample size to control monte carlo estimation

Theta<-c() # storge theta

g<-function(x){ x*factorial(125)/factorial(x)/factorial(125-x) *( theta/(2+theta) )^x *( 2/(2+theta) )^(125-x) }
Q<-function(x){ 38*log(1/2-x/2) +34*log(x/4) +log(x/4)*I_y3 }
while(abs(theta-theta_2) >= 10^-7){
  n<-n+1000
  cat("now theta is",theta_2,"\n")
  theta <-theta_2
  # estimate E(Y3|X) 
  I_y3<-sum( g(runif(n,0,125)) ) /n *125
  # max Q(theta)
  h<-10^-6
  x1<-0 ;x2<-0.5  
  
  while(abs(x2-x1)>= 10^-5 ){
    x1<-x2
    firdif<-(Q(x1+h)-Q(x1))/h
    secdif<-(Q(x1+2*h)-2*Q(x1+h)+Q(x1) ) / h^2
    x2<-x1 -firdif/secdif
  }
  theta_2<-x2
  Theta<-c(Theta,theta_2)
}

# observed data (Y1,Y2) ~ multnomial(72 ,p1' ,p2')
log.obser<-function(x){ 38*log( (2-2*x)/(2-x) ) +34*log( x/(2-x) ) }
plot(log.obser(Theta), type = "l" ,main = "MC-step" ,xlab = "i-th iteration")
