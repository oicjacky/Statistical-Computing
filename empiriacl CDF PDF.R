#empirical CDF & PDF of N(0,1)
CDF<-function(a,b,n){
  x<-a
  cdf<-0 
  while(a<=b){
    I<-0
    for(i in 1:n){ if(rnorm(1,0,1)<=a){ I<-I+1  }  }
    cdf<-c(cdf,I/n)
    a<-a+0.005
  }
  plot(seq(x,b,0.005),cdf[-1]) ;lines(c(0,0),c(0,1)) ;lines(c(x,b),c(0.5,0.5))
  return(cdf)
}
CDF(-5,5,500)

a=-3;b=3;band<-0.001;n=10^4
PDF<-function(a,b,n,band){
  x<-a
  pdf<-rep(0,length(seq(x,b,0.005)))
  j<-1
  while(a<=b){
    I<-0
    for(i in 1:n){ if( abs(a-rnorm(1,0,1))<=band){ I<-I+1  }  }
    pdf[j]<-I/(2*n*band)
    a<-a+0.005
    j<-j+1
  }
  windows();plot(seq(x,b,0.005),pdf) #;lines(c(0,0),c(0,1)) ;lines(c(-x,b),c(0.5,0.5))
  return(pdf)
}
PDF(-3.3,3.3,10^4,0.001)
