rm(list = ls())
library(numDeriv)
#(1)  f(x1,x2)=(x1-2)^4+(x1-2*x2)^2
f<-function(x){ (x[1]-2)^4+(x[1]-2*x[2])^2 } 
ff<-function(x1,x2){ (x1-2)^4+(x1-2*x2)^2 }  
plot3d(ff)
#gradient & Hessian 
h<-10^-6
gradient.f<-function(x){
  y<-c( ( f(c(x[1]+h,x[2]))-f(x) ) / h ,  ( f(c(x[1],x[2]+h))-f(x) ) / h )
  return(y)
}
gradient.f(x0)
#dxdx 
( f(c(x[1]+2*h,x[2])) - 2*f(c(x[1]+h,x[2])) + f(x) ) / h^2
#dxdy=dydx
( f(c(x[1]+h,x[2]+h)) - f(c(x[1]+h,x[2])) - f(c(x[1],x[2]+h)) + f(x) ) / h^2
#dydy
( f(c(x[1],x[2]+2*h)) - 2*f(c(x[1],x[2]+h)) + f(x) ) / h^2
Hessian<-function(x){
  y<-matrix(c( ( f(c(x[1]+2*h,x[2])) - 2*f(c(x[1]+h,x[2])) + f(x) ) / h^2,
               ( f(c(x[1]+h,x[2]+h)) - f(c(x[1]+h,x[2])) - f(c(x[1],x[2]+h)) + f(x) ) / h^2,
               ( f(c(x[1]+h,x[2]+h)) - f(c(x[1]+h,x[2])) - f(c(x[1],x[2]+h)) + f(x) ) / h^2,
               ( f(c(x[1],x[2]+2*h)) - 2*f(c(x[1],x[2]+h)) + f(x) ) / h^2) ,2,2 ) 
  return(y)
}
Hessian(x0)
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
  x2<-x1 - solve(hessian(f,x1)) %*% grad(f,x1)
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
g<-function(alpha){ f(x1-alpha*grad(f,x1)) }
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
    alpha2<-alpha1 - firdif/secdif             #solve(hessian(g,alpha1)) %*% grad(g,alpha1)
    #cat("The ",loop2," step ","g(alpha) is ",g(alpha2),"\n")
  };#cat("The minimum of g(alpha) is ",g(alpha2),", the alpha.hat is ",alpha2,"\n")
  alpha.hat<-alpha2
  #step
  x2<-x1 - alpha.hat*grad(f,x1)
  #A<-rbind(A,c(loop1,f(x2),x2,alpha.hat))
  cat("The",loop1,"step f(x) is",f(x2),",and x is (",x2,"),alpha is",alpha.hat,"\n")
};cat("The minimum of f(x) is ",f(x2),", at (",x2,")")
#colnames(A)<-c("i-th step","f(x)","x1","x2","alpha.hat")

#(2)  f(x1,x2)=(x1-2)^2+(x1-2*x2)^2
rm(list = ls())
ff<-function(x){ (x[1]-2)^2+(x[1]-2*x[2])^2 }

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
  x2<-x1 - solve(hessian(ff,x1)) %*% grad(ff,x1)
  cat("The ",loop," step ","f(x) is ",ff(x2)," ,and x is (",x2,") \n")
};cat("The minimum of f(x) is ",ff(x2),", at (",x2,")")

#Steepest Descent method
#initial value 
x0<-c(0,3)
epslon<-10^-5
loop1<-0
#step
x2<-x0 
x1<-c(1,1)
#A<-c()
g<-function(alpha){ ff(x1-alpha*grad(ff,x1)) }
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
    firdiff<-(g(alpha1+h)-g(alpha1))/h
    secdiff<-(g(alpha1+2*h)-2*g(alpha1+h)+g(alpha1) ) / h^2
    alpha2<-alpha1 - firdiff/secdiff             #solve(hessian(g,alpha1)) %*% grad(g,alpha1)
    #cat("The ",loop2," step ","g(alpha) is ",g(alpha2),"\n")
  };#cat("The minimum of g(alpha) is ",g(alpha2),", the alpha.hat is ",alpha2,"\n")
  alpha.hat<-alpha2
  #step
  x2<-x1 - alpha.hat*grad(ff,x1)
  #A<-rbind(A,c(loop1,ff(x2),x2,alpha.hat))
  cat("The",loop1,"step f(x) is",ff(x2),",and x is (",x2,"),alpha is",alpha.hat,"\n")
};cat("The minimum of f(x) is ",ff(x2),", at (",x2,")")
#colnames(A)<-c("i-th step","f(x)","x1","x2","alpha.hat")

#Conjugate gradient method
rm(list = ls())
ff<-function(x){ (x[1]-2)^2+(x[1]-2*x[2])^2 }

Q<-matrix(c(4,-4,-4,8),2,2)
y1<-c(0,3)
d1<--grad(ff,y1) ;d3<-d1

loop<-0
y3<-y1 
while( loop<dim(Q)[1] ){
  loop<-loop+1
  d2<-d3      
  y2<-y3         #"2" is now ,"3" is next
  
  alpha.hat<-- grad(ff,y2) %*% d2 / d2 %*% Q %*% d2
  y3<-y2 + alpha.hat*d2
  if(loop<dim(Q)[1]){
    lambda2<- - grad(ff,y3) %*% Q %*% d2 / d2 %*% Q %*% d2
    d3<- - grad(ff,y3) + lambda2*d2
  }
}
