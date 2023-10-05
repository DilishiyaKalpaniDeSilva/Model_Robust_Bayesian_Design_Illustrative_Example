expit<-function(a){
  exp(a)/(1+exp(a))
}

logit <- function(a){
  log(a/(1-a))
}

norm_01 <- function(a){
  min <- min(a)
  max <- max(a)
  ans <-lapply(seq_len(length(a)), function(x)
  
      (a[x] - min ) / ( max - min)
  )
  return(unlist(ans))  
}

stand <- function(a){
  return((a-mean(a))/sd(a))
}

norm_01_mm <- function(a,min,max){
  ans <-lapply(seq_len(length(a)), function(x)
    
    (a[x] - min ) / ( max - min)
  )
  return(unlist(ans))  
}

trace <- function(M){
  tr <- sum(diag(M))
  tr
}

