#Exercise 1.
#non-parametric density est of kernal method
#R(f'')
R_f<-function(x){ 1/pi* (x^2-1)^2 *exp(-x^2) }
I_Simpson<-function(a,b,l){
  #a<-0;b<-20;l<-50
  h<-seq(a,b,(b-a)/l)
  if(l%%2==0){
    z1<-0;z2<-0
    for(i in 2:(l/2)){   z1<-z1+2*R_f(h[2*i-1])  }
    for(i in 1:(l/2)){   z2<-z2+4*R_f(h[2*i])  }
    I_S<-(R_f(h[1])+R_f(h[l+1])+z1+z2)*(b-a)/l/3
  } else {print("n need to be even")} }
Rf<-I_Simpson(0,20,50)
#R(k) with Eparechinlcov kernal
Rk<-3/5
#optimal bandwith b_n 
bandwidth<-function(n){ ( Rk / n *(1/25) *Rf )^(1/5) }
#MISE
MISE<-function(n){ bandwidth(n)^4 *(1/25) *Rf /4 + 1/(n*bandwidth(n)) *Rk }
options(digits = 7)
m<-numeric(100)
for(i in 1:100){
  m[i]<-MISE(100+10*i)
}
plot(m ,xlim = c(0,100),ylim = range(m))

#parametric density est 
X<-rnorm(n,0,1)
#sample mean & variance
Xbar<-sum(X)/n
Ssqr<-1/(n-1)*sum( (X-Xbar)^2 )
f_hat<-function(x){  1/sqrt(2*pi) *sqrt(Ssqr) *exp(-(x-Xbar)^2 /2 /Ssqr) }
f<-function(x){ 1/sqrt(2*pi) *exp(-x^2/2) }


#Exercise 2.
#Eparechinlcov kernel
K<-function(x){   if(x>-1 && x<1){ return(3/4*(1-x^2)) }else{ return(0) } }
#Triangle kernel
K<-function(x){ if(x>-1 && x<1){ return(1-abs(x)) }else{ return(0) } }
#Uniform kernel 
K<-function(x){ if(x>-1 && x<1){ return(1/2) }else{ return(0) } }
#Normal kernel
K<-function(x){ return(1/sqrt(2*pi)*exp(-x^2/2)) }
#g1 density
data<-faithful
n<-length(data[,1])
b_n<-c(0.05,0.1,0.5,0.21)
for(k in 1:length(b_n)){
  x<-seq(min(data[,1]),max(data[,1]),0.01)
  f<-rep(0,length(x))
  for(l in 1:length(x)){
    y<-rep(0,n)
    for(i in 1:n){
      y[i]<-K( (x[l]-data[i,1]) /b_n[k] ) /n /b_n[k]
      f[l]<-sum(y)
    }
  }
  if(k==1){ windows();plot(x,f,type = "l",lwd=2) }
  else{ lines(x,f,col= k,lwd=2) }
}
#peseudo likelihood
b_n<-seq(0.01,3,0.01)
peseudo<-rep(0,length(b_n))
for(k in 1:length(b_n)){
  pes<-1
  for(j in 1:length(data[,1])){
    X1<-data[-j,1]
    n<-length(data[,1])-1
    x<-data[j,1]
    y<-rep(0,n)
    for(i in 1:n){                       #f.hat(x_i)
      y[i]<-K( (x-X1[i]) /b_n[k] ) /n /b_n[k]
      f<-sum(y)
    }
    pes<-pes*f
  }
  peseudo[k]<-pes
}
bn.hat<-b_n[which(peseudo==max(peseudo))]
#g2 density
data<-faithful
n<-length(data[,2])
b_n<-c(2.4,2.5,3.5,4)    #c(1,2,3,4)
for(k in 1:length(b_n)){
  x<-seq(min(data[,2]),max(data[,2]),0.1)
  f<-rep(0,length(x))
  for(l in 1:length(x)){
    y<-rep(0,n)
    for(i in 1:n){
      y[i]<-K( (x[l]-data[i,2]) /b_n[k] ) /n /b_n[k]
      f[l]<-sum(y)
    }
  }
  if(k==1){ windows();plot(x,f,type = "l",lwd=2) }
  else{ lines(x,f,col= k,lwd=2) }
}
#peseudo likelihood
data<-faithful
b_n<-seq(0.1,6,0.1)
peseudo<-rep(0,length(b_n))
for(k in 1:length(b_n)){
  pes<-1
  for(j in 1:length(data[,2])){
    X2<-data[-j,2]
    n<-length(data[,2])-1
    x<-data[j,2]
    y<-rep(0,n)
    for(i in 1:n){                       #f.hat(x_i)
      y[i]<-K( (x-X2[i]) /b_n[k] ) /n /b_n[k]
      f<-sum(y)
    }
    pes<-pes+log(f)
    print(log(f))
  }
  peseudo[k]<-pes
}
bn.hat<-b_n[which(peseudo==max(peseudo))]
#g(x1,x2) density
data<-faithful
n<-272
#j=1
b_n1<-0.5
x1<-seq(min(data[,1]),max(data[,1]), (max(data[,1])-min(data[,1]))/50 )
#j=2
x2<-seq(min(data[,2]),max(data[,2]), (max(data[,2])-min(data[,2]))/50 )
b_n2<-2
#g(x1,x2) density 對角線部分
n<-272
x<-cbind(x1,x2)
f<-rep(0,length(x[,1]))
for(l in 1:length(x[,1])){
  y<-rep(0,n)
  for(i in 1:n){
    y[i]<-K( (x[l,1]-data[i,1]) /b_n1 ) /b_n1 /n *K( (x[l,2]-data[i,2]) /b_n2 ) /b_n2 /n
    f[l]<-sum(y)
  }
}
windows();plot(x[,1],f,type = "l",lwd=2)
windows();plot(x[,2],f,col = 2,type = "l",lwd=2)
library(scatterplot3d)
windows();scatterplot3d(x1,x2,f,color = ,type = "h",
                        xlab = "X1",
                        ylab = "X2",
                        zlab = "joint density")
library(plotly)
plot_ly(data.frame(x=x1 ,y=x2 ,z=f)
        , x = ~x ,y = ~y, z = ~z,color = , colors = ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = "X1"),
                      yaxis = list(title = "X2"),
                      zaxis = list(title = "joint density")))
#3維圖全部
n<-272
f<-matrix(0,length(x2),length(x1))
for(k in 1:length(x1)){
  for(l in 1:length(x2)){
    y<-rep(0,n)
    for(i in 1:n){
      y[i]<-K( (x1[k]-data[i,1]) /b_n1 ) /b_n1 /n *K( (x2[l]-data[i,2]) /b_n2 ) /b_n2 /n
      f[k,l]<-sum(y)
    }
  }
}

B<-c()  
for(i in 1:51){
  A<-cbind(x1,x2[i])
  B<-rbind(B,A)
}
windows();scatterplot3d(B[,1],B[,2],f,color = ,type = "p",
                        xlab = "X1",
                        ylab = "X2",
                        zlab = "joint density")
plot_ly(data.frame(x=B[,1] ,y=B[,2] ,z=as.vector(f) )
        , x = ~x ,y = ~y, z = ~z,color = , colors = ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = "X1"),
                      yaxis = list(title = "X2"),
                      zaxis = list(title = "joint density")))

