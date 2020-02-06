# simulation of mixture normal
# true par
p<-c(p1= 0.2 ,p2= 0.3 ,p3= 0.5)
mu<-c(mu1= -3, mu2=  0, mu3= 2)
sigma.sqr<-c(sigma1.sqr= 0.64, sigma2.sqr= 0.36, sigma3.sqr= 0.25)


## generate data :2000 (400, 1400 ,200)

G<-function(z,n){
  z<-c(z,rep(0,n))
  for(i in 1:n){z[i+1]<-(16807*z[i])%%(2^31-1)}
  u<-z[-1]/(2^31-1) }
#acceptance rejection of N(0,1) with g(x)~exp(1)
windows();plot(x=c(seq(0,3,0.001),seq(0,3,0.001)),y=c(dnorm(seq(0,3,0.001),0,1),2*dexp(seq(0,3,0.001),1)),
               col =c(rep(1,3001),rep(2,3001)) )

envo<-function(x,c){  (sqrt(2/pi)*exp(-x^2/2)) / (c*exp(-x)) }
acc.rej.exp<-function(n,c){
  u1<-G(1,c*n)
  Y<--log(u1)
  u2<-G(u1[c*n],c*n)
  X<-rep(0,n)
  N<-length(which(u2 <= envo(Y,c) ))
  for(i in 1:min(n,N) ){   X[i]<-Y[which( u2 <= envo(Y,c) )][i]    }    #exp(-(Y-1)^2/2)
  uu2<-rep(0,n)
  uu2[n]<-u2[c*n]
  while(N<n){
    uu1<-G(uu2[n],n)
    Y<--log(uu1)
    uu2<-G(uu1[n],n)
    if(length(which(uu2<=envo(Y,c)))>0){
      for(i in 1:min(n,length(which(uu2<=envo(Y,c) ) )) ){ X[N+i]<-Y[which(uu2<=envo(Y,c))][i] } }
    N<-N+length(which(uu2 <= envo(Y,c)))
  }
  u3<-G(uu2[n],n)
  X[which(u3<=0.5)]<--X[which(u3<=0.5)]
  return(X[1:n])
}
A<-acc.rej.exp(2000,2)
windows();qqnorm(A);qqline(A)
ks.test(A ,"pnorm" ,0 ,1)
# X1
X1<- sqrt(sigma.sqr[1])*A[1:(2000*p[1])] +mu[1]
ks.test(X1 ,"pnorm" ,mu[1] ,sqrt(sigma.sqr[1]))
# X2
X2<- sqrt(sigma.sqr[2])*A[(2000*p[1]+1):(2000*(p[1]+p[2]))] +mu[2]
ks.test(X2 ,"pnorm" ,mu[2] ,sqrt(sigma.sqr[2]))
# X3
X3<- sqrt(sigma.sqr[3])*A[(2000*(p[1]+p[2])+1):2000] +mu[3]
ks.test(X3 ,"pnorm" ,mu[3] ,sqrt(sigma.sqr[3]))

# observed data X=(X1,X2,X3)
obs.X<-c(X1,X2,X3)
rm(A,X1,X2,X3)

hist(obs.X ,breaks = 100)
windows();qqnorm(obs.X);qqline(obs.X)

## EM algorithm
# initial value
p1_0<- 1/3 ; mu1_0<- -3 ; sigma1.sqr_0<- 0.64
p2_0<- 1/3 ; mu2_0<- 0 ; sigma2.sqr_0<- 0.36
p3_0<- 1/3 ; mu3_0<- 2 ; sigma3.sqr_0<- 0.25

p1_2<- p1_0
p2_2<- p2_0
p3_2<- p3_0
mu1_2<- mu1_0
mu2_2<- mu2_0
mu3_2<- mu3_0
sigma1.sqr_2<- sigma1.sqr_0
sigma2.sqr_2<- sigma2.sqr_0
sigma3.sqr_2<- sigma3.sqr_0



for(i in 1:1000){
p1_1<- p1_2
p2_1<- p2_2
p3_1<- p3_2
mu1_1<- mu1_2
mu2_1<- mu2_2
mu3_1<- mu3_2
sigma1.sqr_1<- sigma1.sqr_2
sigma2.sqr_1<- sigma2.sqr_2
sigma3.sqr_1<- sigma3.sqr_2

# latent: Z_i1 ,Z_i2 ,Z_i3 update
f<-function(x,mu,sigma.sqr){ 1/sqrt(2*pi*sigma.sqr) *exp( -(x-mu)^2 /sigma.sqr)  }
Z_i1<- ( p1_1 *f(obs.X ,mu1_1 ,sigma1.sqr_1) ) /( p1_1 *f(obs.X ,mu1_1 ,sigma1.sqr_1) + p2_1 *f(obs.X ,mu2_1 ,sigma2.sqr_1) + p3_1 *f(obs.X ,mu3_1 ,sigma3.sqr_1) )
Z_i2<- ( p2_1 *f(obs.X ,mu2_1 ,sigma2.sqr_1) ) /( p1_1 *f(obs.X ,mu1_1 ,sigma1.sqr_1) + p2_1 *f(obs.X ,mu2_1 ,sigma2.sqr_1) + p3_1 *f(obs.X ,mu3_1 ,sigma3.sqr_1) )
Z_i3<- ( p3_1 *f(obs.X ,mu3_1 ,sigma3.sqr_1) ) /( p1_1 *f(obs.X ,mu1_1 ,sigma1.sqr_1) + p2_1 *f(obs.X ,mu2_1 ,sigma2.sqr_1) + p3_1 *f(obs.X ,mu3_1 ,sigma3.sqr_1) )
Z<-cbind(Z_i1,Z_i2,Z_i3)
# p1 ,p2 ,p3 update
p1_2<- sum(Z_i1) /sum(Z_i1+Z_i2+Z_i3)
p2_2<- sum(Z_i2) /sum(Z_i1+Z_i2+Z_i3)
p3_2<- sum(Z_i3) /sum(Z_i1+Z_i2+Z_i3)

# mu1 ,mu2 ,mu3 update
mu1_2<- as.vector(Z_i1 %*% obs.X /sum(Z_i1))
mu2_2<- as.vector(Z_i2 %*% obs.X /sum(Z_i2))
mu3_2<- as.vector(Z_i3 %*% obs.X /sum(Z_i3))

# sigma1^2 ,sigma2^2 ,sigma3^2 update
sigma1.sqr_2<- as.vector((obs.X - mu1_2)^2  %*% Z_i1 /sum(Z_i1))
sigma2.sqr_2<- as.vector((obs.X - mu2_2)^2  %*% Z_i2 /sum(Z_i2))
sigma3.sqr_2<- as.vector((obs.X - mu3_2)^2  %*% Z_i3 /sum(Z_i3))
cat("now p is ",c(p1_2,p2_2,p3_2),",mu is ",c(mu1_2,mu2_2,mu3_2),",sigma is ",c(sigma1.sqr_2,sigma2.sqr_2,sigma3.sqr_2),"\n")
}
