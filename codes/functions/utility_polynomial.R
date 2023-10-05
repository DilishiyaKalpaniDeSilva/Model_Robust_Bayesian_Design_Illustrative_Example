##  KLD utility function
kld.post <- function(x,y,theta){
  
  C <- x
  com1 <- (1/Sigmae^2)*(t(C)%*%C)
  FIM <- com1+iSp
  Sigma_post <- solve(FIM)
  mu_post <- (1/Sigmae^2)*(Sigma_post%*%t(C)%*%y)
  
  # det.out <- 0.5*(trace(iSp%*%Sigma_post) + t(mp-mu_post)%*%iSp%*%(mp-mu_post) - num_par + 
  #                   log(dSigma_prior/det(Sigma_post)))
  det.out <- 0.5*( trace(iSp%*%Sigma_post) + t(mp-mu_post)%*%iSp%*%(mp-mu_post) - num_par + 
                     ( ldSigma_prior - as.vector(determinant(Sigma_post)$modulus) ) )
  det.out
}

# Approximate expected utility
exp.crit1 <- function(d,B){

  d <- base::rep(d,(n/num_par))
  
  x <- 1
  for(ii in 1:indexB){
    x <- cbind(x,d^(ii))
  }
  x <- as.matrix(x)
  crit  <- foreach(j = 1:B, .packages=c("doParallel","HRW"),
                   .combine = c,.errorhandling = 'stop',
                   .export = ls(globalenv())) %dopar% {
                     theta <- theta_sim[j,]
                     lp <- as.vector(x%*%theta) 
                     y <- qnorm(rand_y[j,],mean=lp,sd=Sigmae)
                     #y <- as.numeric(p > runif(length(x)))
                     kld.post(x,y,theta)             
                   }
  return(crit)
}

# Approximate u(X)
exp.crit2 <- function(d,B){
  
  x <- 1
  for(ii in 1:indexB){
    x <- cbind(x,d^(ii))
  }
  x <- as.matrix(x)
  crit  <- foreach(j = 1:B, .packages=c("doParallel","HRW"),
                   .combine = c,.errorhandling = 'stop',
                   .export = ls(globalenv())) %dopar% {
                     theta <- theta_sim[j,]
                     lp <- as.vector(x%*%theta) 
                     y <- qnorm(rand_y[j,],mean=lp,sd=Sigmae)
                     #y <- as.numeric(p > runif(length(x)))
                     kld.post(x,y,theta)             
                   }
  return(crit)
}
