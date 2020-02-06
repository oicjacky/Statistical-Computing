#(2)  f(x1,x2)=(x1-2)^2+(x1-2*x2)^2
rm(list = ls())
ff<-function(x){ (x[1]-2)^2+(x[1]-2*x[2])^2 }

#gradient & Hessian 
h<-10^-6
gradient.ff<-function(x){
  y<-c( ( ff(c(x[1]+h,x[2]))-ff(x) ) / h ,  ( ff(c(x[1],x[2]+h))-ff(x) ) / h )
  return(y)
}

Hessian.ff<-function(x){
  y<-matrix(c( ( ff(c(x[1]+2*h,x[2])) - 2*ff(c(x[1]+h,x[2])) + ff(x) ) / h^2,
               ( ff(c(x[1]+h,x[2]+h)) - ff(c(x[1]+h,x[2])) - ff(c(x[1],x[2]+h)) + ff(x) ) / h^2,
               ( ff(c(x[1]+h,x[2]+h)) - ff(c(x[1]+h,x[2])) - ff(c(x[1],x[2]+h)) + ff(x) ) / h^2,
               ( ff(c(x[1],x[2]+2*h)) - 2*ff(c(x[1],x[2]+h)) + ff(x) ) / h^2) ,2,2 ) 
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
  x2<-x1 - solve(Hessian.ff(x1)) %*% gradient.ff(x1)
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
g<-function(alpha){ ff(x1-alpha*gradient.ff(x1)) }
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
    alpha2<-alpha1 - firdiff/secdiff             
    #cat("The ",loop2," step ","g(alpha) is ",g(alpha2),"\n")
  };#cat("The minimum of g(alpha) is ",g(alpha2),", the alpha.hat is ",alpha2,"\n")
  alpha.hat<-alpha2
  #step
  x2<-x1 - alpha.hat*gradient.ff(x1)
  
  cat("The",loop1,"step f(x) is",ff(x2),",and x is (",x2,"),alpha is",alpha.hat,"\n")
};cat("The minimum of f(x) is ",ff(x2),", at (",x2,")")


#Conjugate gradient method
Q<-matrix(c(4,-4,-4,8),2,2)
y1<-c(0,3)
d1<--gradient.ff(y1) ;d3<-d1
D<-c()

loop<-0
y3<-y1 
while( loop<dim(Q)[1] ){
  loop<-loop+1
  d2<-d3      
  y2<-y3         #"2" is now ,"3" is next
  
  alpha.hat<-- gradient.ff(y2) %*% d2 / d2 %*% Q %*% d2
  y3<-y2 + alpha.hat*d2
  D<-rbind(D,d2)
  #if(loop<dim(Q)[1]){
    lambda2<- - gradient.ff(y3) %*% Q %*% d2 / d2 %*% Q %*% d2
    d3<- - gradient.ff(y3) + lambda2*d2
  #}
  cat("The",loop,"step f(x) is",ff(y3),",and x is (",y3,") \n")
}
cat("The minimum of f(x) is ",ff(y3),", at (",y3,")")
row.names(D)<-c("d1","d2")
