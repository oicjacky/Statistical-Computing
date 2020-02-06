#(1)  f(x1,x2)=(x1-2)^4+(x1-2*x2)^2
rm(list = ls())
f<-function(x){ (x[1]-2)^4+(x[1]-2*x[2])^2 } 
#ff<-function(x1,x2){ (x1-2)^4+(x1-2*x2)^2 }  
#library(rgl);plot3d(ff)
#gradient & Hessian 
h<-10^-6
gradient.f<-function(x){
  y<-c( ( f(c(x[1]+h,x[2]))-f(x) ) / h ,  ( f(c(x[1],x[2]+h))-f(x) ) / h )
  return(y)
}
Hessian.f<-function(x){
  y<-matrix(c( ( f(c(x[1]+2*h,x[2])) - 2*f(c(x[1]+h,x[2])) + f(x) ) / h^2,
               ( f(c(x[1]+h,x[2]+h)) - f(c(x[1]+h,x[2])) - f(c(x[1],x[2]+h)) + f(x) ) / h^2,
               ( f(c(x[1]+h,x[2]+h)) - f(c(x[1]+h,x[2])) - f(c(x[1],x[2]+h)) + f(x) ) / h^2,
               ( f(c(x[1],x[2]+2*h)) - 2*f(c(x[1],x[2]+h)) + f(x) ) / h^2) ,2,2 ) 
  return(y)
}

#Newton method
#initial value 
x0<-c(0,3)
epslon<-10^-5
loop<-0
#step
x2<-x0 ;x1<-c(1,1)
while( sqrt(sum((x2-x1)^2)) >= epslon){
  loop<-loop+1
  x1<-x2
  x2<-x1 - solve(Hessian.f(x1)) %*% gradient.f(x1)
  cat("The ",loop," step ","f(x) is ",f(x2)," ,and x is (",x2,") \n")
};cat("The minimum of f(x) is ",f(x2),", at (",x2,")")

#Steepest Descent method
#initial value 
x0<-c(0,3)
epslon<-10^-5
loop1<-0
#step
x2<-x0 
x1<-c(1,1)
A<-c()
g<-function(alpha){ f(x1-alpha*gradient.f(x1)) }
while(sqrt(sum((x2-x1)^2)) >= epslon &&loop1<40){
  loop1<-loop1+1
  x1<-x2
  
  #alpha.hat = argmin g(alpha)
  alpha0<-0
  loop2<-0
  h<-10^-6
  
  alpha2<-alpha0 ;alpha1<-1
  while( abs(alpha2-alpha1) >= epslon ){
    loop2<-loop2+1
    alpha1<-alpha2
    firdif<-(g(alpha1+h)-g(alpha1))/h
    secdif<-(g(alpha1+2*h)-2*g(alpha1+h)+g(alpha1) ) / h^2
    alpha2<-alpha1 - firdif/secdif             
    #cat("The ",loop2," step ","g(alpha) is ",g(alpha2),"\n")
  };#cat("The minimum of g(alpha) is ",g(alpha2),", the alpha.hat is ",alpha2,"\n")
  alpha.hat<-alpha2
  #step
  x2<-x1 - alpha.hat*grad(f,x1)
  
  cat("The",loop1,"step f(x) is",f(x2),",and x is (",x2,"),alpha is",alpha.hat,"\n")
};cat("The minimum of f(x) is ",f(x2),", at (",x2,")")
