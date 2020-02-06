#newton method
f<-function(x){  -x^4+3*x^2+2 }
plot(f ,xlim = c(-3,3) ,ylim = c(-4,4) ,lwd= 2) ;lines(c(-3,3),c(0,0))
epslon<-10^-9 ;h<-10^-8 ;loop<-0
seed<--0.5 ;x2<-seed
while(abs(f(x2))>=epslon && loop<500 ){
  cat("the ",loop,"-th step, f(x) is ",abs(f(x2)),"\n")
  loop<-loop+1
  x1<-x2
  dif<-(f(x1+h)-f(x1)) /h
  x2<-x1 - f(x1)/dif
  lines(c(x1,x1),c(0,f(x1)),col=2) ;lines(c(x1,x2),c(f(x1),0) ,col=3)
};cat("the root of f(x) is",x2)

g<-function(x){ x-(-x^4+3*x^2+2)/(-4*x^3+6*x)}
windows();plot(g,xlim = c(-3,3), ylim=c(-4,4),lwd=2,type="p") ;lines(c(-3,3),c(1,1)) ;lines(c(-3,3),c(-1,-1)) ;lines(c(-3,3),c(-3,3),lwd=2)
x<-seq(0,2,0.001)
plot(x,abs((g(x+h)-g(x))/h), ylim = c(0,2))
x[abs((g(x+h)-g(x))/h) < 1]

#bisection
f<-function(x){ -x^4+3*x^2+2 }
plot(f ,xlim = c(-3,3) ,ylim = c(-4,4) ,lwd= 2) ;lines(c(-3,3),c(0,0))
epslon<-10^-9 ;loop<-0
xl<--51 ;xr<--1 ;xm<-(xl+xr)/2
while(abs(xl-xr)>=epslon){
  cat("the ",loop,"-th step, f(x) is ",abs(f(xm)),"\n")
  loop<-loop+1
  xm<-(xl+xr)/2
  if(f(xl)*f(xr)>0){
    break
    }else if(f(xl)*f(xm)<0){
      xr<-xm
      }else if(f(xr)*f(xm)<0){
        xl<-xm }
};cat("the root of f(x) is",xm)
