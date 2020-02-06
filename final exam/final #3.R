#############################################
#C1: x_0 = 0 ; n = 1000 ; different epslon  #
#############################################
#case: x_0 = 0 ; epslon = 1
x_0<- 0
epslon<- 1
x_2<- x_0
X<- c()

n<-1000
for(i in 1:n){
x_1<- x_2
# generate y ~ U( x_1-epslon ,x_1+epslon )
y<- runif(1, x_1-epslon ,x_1+epslon)
# acceptance prob.
A<-function(x,y){ min( exp(-y^2/2 +x^2/2), 1 ) }
u<- runif(1, 0, 1)
if( u <= A(x_1,y) ){ x_2<- y }else{ x_2<- x_1 }
X<- c(X,x_2)
}
windows();par(mfrow=c(3,2))
plot(X ,type = "l" ,xlab = "n=1000" ,ylab = "X_n" ,col = "1",main ="x_0=0 ,epslon=1" )
qqnorm(X);qqline(X)
#case: x_0 = 0 ; epslon = 0.2
x_0<- 0
epslon<- 0.2
x_2<- x_0
X<- c()

n<-1000
for(i in 1:n){
  x_1<- x_2
  # generate y ~ U( x_1-epslon ,x_1+epslon )
  y<- runif(1, x_1-epslon ,x_1+epslon)
  # acceptance prob.
  A<-function(x,y){ min( exp(-y^2/2 +x^2/2), 1 ) }
  u<- runif(1, 0, 1)
  if( u <= A(x_1,y) ){ x_2<- y }else{ x_2<- x_1 }
  X<- c(X,x_2)
}
plot(X ,type = "l" ,xlab = "n=1000" ,ylab = "X_n" ,col = "1",main ="x_0=0 ,epslon=0.2" )
qqnorm(X);qqline(X)
#case: x_0 = 0 ; epslon = 20
x_0<- 0
epslon<- 20
x_2<- x_0
X<- c()

n<-1000
for(i in 1:n){
  x_1<- x_2
  # generate y ~ U( x_1-epslon ,x_1+epslon )
  y<- runif(1, x_1-epslon ,x_1+epslon)
  # acceptance prob.
  A<-function(x,y){ min( exp(-y^2/2 +x^2/2), 1 ) }
  u<- runif(1, 0, 1)
  if( u <= A(x_1,y) ){ x_2<- y }else{ x_2<- x_1 }
  X<- c(X,x_2)
}
plot(X ,type = "l" ,xlab = "n=1000" ,ylab = "X_n" ,col = "1",main ="x_0=0 ,epslon=20" )
qqnorm(X);qqline(X)
###########################################
#C2: x_0 = -50 ; epslon = 1 ; different n # 
###########################################
#case: x_0 = -50 ; epslon = 1 ; n = 1000 
x_0<- -50
epslon<- 1
x_2<- x_0
X<- c()

n<-1000
for(i in 1:n){
  x_1<- x_2
  # generate y ~ U( x_1-epslon ,x_1+epslon )
  y<- runif(1, x_1-epslon ,x_1+epslon)
  # acceptance prob.
  A<-function(x,y){ min( exp(-y^2/2 +x^2/2), 1 ) }
  u<- runif(1, 0, 1)
  if( u <= A(x_1,y) ){ x_2<- y }else{ x_2<- x_1 }
  X<- c(X,x_2)
}
windows();par(mfrow=c(4,2))
plot(X ,type = "l" ,xlab = "n=1000" ,ylab = "X_n" ,col = "1",main = "x_0=-50, epslon=1, without burn-in")
qqnorm(X);qqline(X)
plot(X[(n/2):n] ,type = "l" ,xlab = "n=1000" ,ylab = "X_n" ,main = "x_0=-50, epslon=1, burn-in with 500 samples" ,cex.main =1)
qqnorm(X[(n/2):n]);qqline(X[(n/2):n])
#case: x_0 = -50 ; epslon = 10 ; n = 1000 
x_0<- -50
epslon<- 10
x_2<- x_0
X<- c()

n<-1000
for(i in 1:n){
  x_1<- x_2
  # generate y ~ U( x_1-epslon ,x_1+epslon )
  y<- runif(1, x_1-epslon ,x_1+epslon)
  # acceptance prob.
  A<-function(x,y){ min( exp(-y^2/2 +x^2/2), 1 ) }
  u<- runif(1, 0, 1)
  if( u <= A(x_1,y) ){ x_2<- y }else{ x_2<- x_1 }
  X<- c(X,x_2)
}
plot(X ,type = "l" ,xlab = "n=1000" ,ylab = "X_n" ,col = "1",main = "x_0=-50, epslon=10, without burn-in")
qqnorm(X);qqline(X)
plot(X[(n/4):n] ,type = "l" ,xlab = "n=1000" ,ylab = "X_n" ,main = "x_0=-50, epslon=10, burn-in with 250 samples" ,cex.main =1)
qqnorm(X[(n/4):n]);qqline(X[(n/4):n])

###############################################
# to estimate E( exp(Z^16) )
#case: x_0 = 0 ; epslon = 1
x_0<- 0
epslon<- 1
x_2<- x_0
X<- c()

n<-1000
for(i in 1:n){
  x_1<- x_2
  # generate y ~ U( x_1-epslon ,x_1+epslon )
  y<- runif(1, x_1-epslon ,x_1+epslon)
  # acceptance prob.
  A<-function(x,y){ min( exp(-y^2/2 +x^2/2), 1 ) }
  u<- runif(1, 0, 1)
  if( u <= A(x_1,y) ){ x_2<- y }else{ x_2<- x_1 }
  X<- c(X,x_2)
}
windows();par(mfrow=c(3,2))
plot(X ,type = "l" ,xlab = "n=1000" ,ylab = "X_n" ,col = "1",main ="x_0=0 ,epslon=1" )
qqnorm(X);qqline(X)

I.hat<- mean( exp(X^16) )
