# Exercise 1 : multnomial 
theta_0<- 0.5
theta_1<- 0
theta_2<- theta_0

Theta<-c() # storge theta

itera<- 0
while( abs(theta_2-theta_1) >= 10^-7 ){
  cat("now theta is",theta_2,"\n")
  itera<- itera +1
  theta_1<- theta_2
  theta_2<- ( 34 + 125*theta_1 / (2 + theta_1) ) / ( 72 + 125*theta_1 / (2 + theta_1) )
  
  Theta<-c(Theta,theta_2)
}
# observed data (Y1,Y2) ~ multnomial(72 ,p1' ,p2')
log.obser<-function(x){ 38*log( (2-2*x)/(2-x) ) +34*log( x/(2-x) ) }
log.obser(Theta)

# Exercise 2 : 
z<-rpois(4075,0.4)
ks.test(z,x)
Y<-c(3062,587,284,103,33,4,6)
Yplum<- Y/sum(Y)
plot(Yplum ,type = "p")

x <- c(rep(0,3062), rep(1,587), rep(2,284), rep(3,103), rep(4,33), rep(5,4), rep(6,2) )
p_0 <- 0.75
mu_0 <- 0.4

p_2 <- p_0     ;p_1 <-0
mu_2 <- mu_0   ;mu_1 <-0
itera<- 0
while( abs(p_2 - p_1) >= 10^-7 | abs(mu_2 - mu_1) >= 10^-7   ){
  itera<- itera +1
  # Z_i1 Z_i2 update
  p_1 <- p_2
  mu_1 <- mu_2
  Z_i1 <- ( p_1 *(x==0) ) / ( p_1 + (1-p_1) * ( mu_1 ^x /factorial(x) *exp(-mu_1) ) ) 
  Z_i2 <- ( (1-p_1) * ( mu_1 ^x /factorial(x) *exp(-mu_1) ) ) / ( p_1 *(x==0) + (1-p_1) * ( mu_1 ^x /factorial(x) *exp(-mu_1) ) ) 
  # p_2  mu_2 update 
  p_2 <- mean(Z_i1)
  mu_2 <- as.numeric(Z_i2 %*% x /sum(Z_i2))
}
A <- cbind(Z_i1,Z_i2)
