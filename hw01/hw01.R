rm(list = ls())
f<-function(x){exp(-x^2)}
a<-0;b<-1;n<-50
h<-seq(a,b,(b-a)/n)      #equally partition
#Rectangle rule
right<-0;left<-0
for(i in 2:(n+1)){
  right<-right+f(h[i])*(b-a)/n
  left<-left+f(h[i-1])*(b-a)/n
}
#y1<-(x1+x2)/2;rm(x1,x2)

#Trapezoidal rule
I_T<-0
for(i in 2:(n+1)){
  I_T<-I_T+(f(h[i-1])+f(h[i]))/2
}
I_T<-I_T*(b-a)/n

#Simpsons rule
if(n%%2==0){
  z1<-0;z2<-0
  for(i in 2:(n/2)){   z1<-z1+2*f(h[2*i-1])  }
  for(i in 1:(n/2)){   z2<-z2+4*f(h[2*i])  }
  I_S<-(f(h[1])+f(h[n+1])+z1+z2)*(b-a)/n/3
  rm(z1,z2)
} else {print("n need to be even")
  #z1<-0;z2<-0
  #for(i in 2:((n+1)/2)){   z1<-z1+2*f(h[2*i-1])  }
  #for(i in 1:((n-1)/2)){   z2<-z2+4*f(h[2*i])  }
  #I_S<-(f(h[1])+f(h[n+1])+z1+z2)*(b-a)/n/3
}

I<-integrate(f,a,b)$value
str(integrate(f,a,b))





