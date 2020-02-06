g1<-function(x){ x^(3/2) } ;plot(g1,xlim = c(0,3),ylim =c(0,5)) ;lines(c(0,3),c(0,3))
g2<-function(x){ x^(1/2) } ;plot(g2,xlim = c(0,3),ylim =c(0,5)) ;lines(c(0,3),c(0,3))

#for g1
epslon<-10^-9
seed<-2.2
x<-seed ;x2<-1 ;x1<-0
plot(g1,xlim = c(0,3),ylim =c(0,5),lwd=3);lines(c(0,3),c(0,3),lwd =3)
while(abs(x2-x1)>=epslon && abs(x2-x1)!=Inf){
  cat("length of step is",abs(x2-x1),"\n")
  x1<-x
  x2<-g1(x1);print(x2);lines(c(x1,x1),c(0,g1(x1)),lwd =1,col=2);lines(c(x1,x2),c(g1(x1),x2),col=3)
  x<-x2
};rm(x)
#for g2
epslon<-10^-9
seed<-2.2
x<-seed ;x2<-1 ;x1<-0
plot(g2,xlim = c(0,3),ylim =c(0,5),lwd=3);lines(c(0,3),c(0,3),lwd =3)
while(abs(x2-x1)>=epslon){
  cat("length of step is",abs(x2-x1),"\n")
  x1<-x
  x2<-g2(x1);print(x2);lines(c(x1,x1),c(0,g2(x1)),lwd =1,col=2);lines(c(x1,x2),c(g2(x1),x2),col=3)
  x<-x2
};rm(x)

#g3 do not exist fixed point
g3<-function(x){ x+1/x }
epslon<-10^-9
seed<-1
x<-seed ;x2<-1 ;x1<-0
plot(g3,xlim = c(0,3),ylim =c(0,5),lwd=3);lines(c(0,3),c(0,3),lwd =3)
while(abs(x2-x1)>=epslon){
  cat("length of step is",abs(x2-x1),"\n")
  x1<-x
  x2<-g3(x1);print(x2);lines(c(x1,x1),c(0,g3(x1)),lwd =1,col=2);lines(c(x1,x2),c(g3(x1),x2),col=3)
  x<-x2
};rm(x)
