pa <- 12
M <- 1000
times <- 12

fun_1 <- function(M ,pa ,times){
  y <- 0
  for(i in 1:times){
    y <- (y + M) * (1 + pa / 12 /100)
  }
  return(y)
}


