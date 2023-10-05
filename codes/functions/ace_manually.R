library(DiceKriging)

f.pred<-function(x,fit)  #### Change how you predict from GP
{
  nd<-matrix(x,nrow=1)  
  #Predicted values and (marginal of joint) conditional variances based on a km model. 95 % confidence intervals are given, based on strong assumptions: Gaussian process assumption, specific prior distribution on the trend parameters, known covariance parameters. This might be abusive in particular in the case where estimated covariance parameters are plugged in
  ans<-predict(fit,newdata=data.frame(z=nd),type="UK",cov.compute=TRUE,se.compute=TRUE)
  out <- ans$mean
  
  return(out)
}

pval <- function(oldeval, neweval, binary){
  if(binary){
    old_n<-length(oldeval)
    new_n<-length(neweval)
    old_sum<-sum(oldeval)
    new_sum<-sum(neweval)
    new_beta_sam<-rbeta(n = 10000, shape1 = 1 + new_sum, shape2 = 1 + new_n - new_sum)
    out<-mean(pbeta(q = new_beta_sam, shape1 = 1 + old_sum, shape2 = 1 + old_n - old_sum))
    
  } 
  else{
    details<-as.vector(.Call( "pvalcpp", oldeval, neweval, PACKAGE = "acebayes" ))
    out<-1-pt(details[1],df=details[2])
  }
  out
  
}

ACE_krige <- function(utility,start.d,B,Q,lower,upper){

  curr_x <- start.d
  curr_uv <- (utility(curr_x,B[1])) ##current utility, Approximate u(X), B=1000
  curr_u <- mean(curr_uv)
  uz <- as.numeric()
  up <- as.numeric()
  
  z <- round(seq(lower,upper, length.out = Q),3) ## Q
  zp <- seq(lower,upper,by=0.001)
  
  for(i in 1:length(start.d)){
    prop_x <- curr_x
    
    for(j in 1:length(z)){
      prop_x[i] <- z[j]
      uz[j] <- mean(utility(prop_x,B[2]))
    }
    
    fit <- km(formula=~1,design=data.frame(z=z),response=uz,covtype="gauss",control=list(trace=FALSE),nugget.estim=TRUE)
    #### Emulation, #distribution of utility
    
    for(k in 1:length(zp)){
      up[k] <- f.pred(zp[k],fit)
    }
    # plot(z, uz)
    # points(zp, up, type="l", lwd=2, col="red")
    
    ind <- which(up==max(up))#index of maximum utility
    prop_x[i] <- zp[ind[1]] #x value which maximizes the utility, for each x value
    print(prop_x[i])
    prop_uv <- (utility(prop_x,B[1]))
    prop_u <- mean(prop_uv)
    
    if(prop_u >= curr_u){ #compare with cu
      curr_x <- prop_x
      curr_u <- prop_u
    }
    print(curr_x[i])
  }
  
  out <- list(curr_x,curr_u)
  return(out)
  
}

