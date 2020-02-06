f<-function(x){exp(-x)*(cos(x^2))^2}
a<-0;b<-50
I<-integrate(f,a,Inf)$value

#Trapezoidal rule
I_T<-I+1;n<-1
while(abs(I_T-I)>=10^-3){
  n<-n+1 
  cat("Now is",n-1,"-th iteration","\n",I_T,"")
  h<-seq(a,b,(b-a)/n)
  for(i in 2:(n+1)){ I_T<-I_T+(f(h[i-1])+f(h[i]))/2 }
  I_T<-I_T*(b-a)/n
}
print(n)

#Simpsons rule
I_S<-I+1;n_s<-2
while(abs(I_S-I)>=10^-3){
  n_s<-n_s+2 
  h<-seq(a,b,(b-a)/n_s)
  z1<-0;z2<-0
  for(i in 2:(n_s/2)){   z1<-z1+2*f(h[2*i-1])  }
  for(i in 1:(n_s/2)){   z2<-z2+4*f(h[2*i])  }
  I_S<-(f(h[1])+f(h[n_s+1])+z1+z2)*(b-a)/n_s/3
  cat("Now is",n_s/2-1,"-th iteration","\n",I_S," ")
}
rm(z1,z2)
print(n_s)





