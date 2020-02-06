# test.4  解A的反矩陣
A <- matrix(c(1 ,1 ,-2  ,1 ,3 ,-1,
              2 ,-1 ,1  ,2 ,1 ,-3,
              1 ,3 ,-3  ,-1,2 ,1 ,
              5 ,2 ,-1  ,-1,2 ,1 ,
              -3,-1, 2  ,3 ,1 ,3 , 
              4 ,3 , 1  ,-6,-3,-2 ), 6 ,6 ,byrow = T)
b <- diag(1 , 6 , 6)
Ab <- cbind(A, b)   # det(A) != 0 

# lower-triangle
for(i in 1:dim(A)[1]){
  if(i < dim(A)[1]){
    # pivot
    l <- Ab[i , ]
    change <- which( Ab[ ,i ] ==  max(Ab[-(1:i) ,i ]) )
    Ab[i , ] <- Ab[ change , ]
    Ab[ change , ] <- l
  }
  # identify
  if(Ab[ i , i ] != 1){   Ab[ i , ] <- Ab[ i , ] / ( Ab[ i, i ] )  }
  
  # elimination
  if(i < dim(A)[1]){
    for(k in (i + 1):dim(A)[1]){
      Ab[k , ] <- Ab[k , ] + Ab[i, ] * (- Ab[k , i] / Ab[i , i] ) 
    }
  }
}  
# upper-triangle
for(j in dim(A)[1]:1){
  if(j > 1){
    for(k in 1: (j-1) ){
      Ab[k , ] <- Ab[k , ] + Ab[j, ] * (- Ab[k , j] / Ab[j , j] )
    }
  }
}
cat("for AB = I, using Gaussian-Jordan elimination, we have B is ")
print( Ab[ ,(dim(A)[1]+1) : (2*dim(A)[1]) ] )
solve(A,b) 

# test.1 解A的反矩陣
A <- matrix(c(1 ,2 ,-1,
              2 ,2 ,1,
              3 ,5 ,-2), 3 ,3 ,byrow = T)
b <- diag(1 ,3 ,3)
Ab <- cbind(A, b)

# lower-triangle
for(i in 1:dim(A)[1]){
  
  if(i < dim(A)[1]){
    # pivot
    l <- Ab[i , ]
    change <- which( Ab[ ,i ] ==  max(Ab[-(1:i) ,i ]) )
    Ab[i , ] <- Ab[ change , ]
    Ab[ change , ] <- l
  }
  # identify
  if(Ab[ i , i ] != 1){   Ab[ i , ] <- Ab[ i , ] / ( Ab[ i, i ] )  }
  
  if(i < dim(A)[1]){
    for(k in (i + 1):dim(A)[1]){
      Ab[k , ] <- Ab[k , ] + Ab[i, ] * (- Ab[k , i] / Ab[i , i] ) 
    }
  }
}  
# upper-triangle
for(j in dim(A)[1]:1){
  if(j > 1){
    for(k in 1: (j-1) ){
      Ab[k , ] <- Ab[k , ] + Ab[j, ] * (- Ab[k , j] / Ab[j , j] )
    }
  }
}
cat("for AB = I, using Gaussian-Jordan elimination, we have B is ")
print( Ab[ ,(dim(A)[1]+1) : (2*dim(A)[1]) ] )
solve(A,b)

