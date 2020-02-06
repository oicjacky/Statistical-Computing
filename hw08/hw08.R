#newton method
f<-function(x){ x^4-14*x^3+60*x^2-70*x }
seed<-1.9
for(i in 1:20){
epslon<-10^-5 ;h<-10^-6
x1<-0;x2<-seed ;loop<-1

while(abs(x2-x1)>=epslon &&loop<500){
  cat("the ",loop,"-th step,","x is ",x2,", f(x) is ",f(x2),"\n")
  loop<-loop+1
  x1<-x2
  firdif<-(f(x1+h)-f(x1))/h
  secdif<-(f(x1+2*h)-2*f(x1+h)+f(x1) ) / h^2
  x2<-x1-firdif/secdif
}
seed<-seed+0.1
}
#golden-section
g<-function(x){ -(x^4-14*x^3+60*x^2-70*x) }
g(seq(0,100,1))
xl<-0 ;g(xl)
xm<-1 ;g(xm)
xr<-2 ;g(xr)
loop<-1
while((xr-xl)>=epslon){#f(xm)-f(xl)>=epslon && f(xm)-f(xr)>=epslon
  if(abs(xl-xm)>=abs(xr-xm)){                #choose large interval
    x4<-(xl+xm)/2
    if(g(x4)>g(xm)){ xr<-xm ;xm<-x4 }else{ xl<-x4 }
    }else{
      x4<-(xr+xm)/2
      if(g(x4)>g(xm)){ xl<-xm ;xm<-x4 }else{ xr<-x4 }
    }
  loop<-loop+1
  cat("left g(x):",g(xl),", middle g(x):",g(xm),", right g(x):",g(xr),"\n")
};cat("left point:",xl,", middle point:",xm,", right point:",xr,"\n")
