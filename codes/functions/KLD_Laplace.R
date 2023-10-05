# log posterior 
log.post <- function(theta,x,y,z.spline){
  log.prior <- sum(dnorm(x=theta,mean=prior_mu,sd=prior_sd,log=TRUE))
  u <- qnorm(rand_u,0, Sigmau)
  s <- z.spline %*% t(u)
  lp <- as.vector(theta[1] + theta[2]*x) + s
  log.like <- mean(colSums( dnorm( y, mean =  lp, sd=0.2, log=TRUE )) )
  log.post <- -1*(log.prior + log.like)
  
  if(is.infinite(log.post)){
    return(10000000)
  }else{
    return(log.post)
  }
  
}

#joint log likelihood
joint_like <- function(u,x,y,z.spline,mu_star){
  s <- z.spline %*% u
  mu  <- mu_star[1] + mu_star[2]*x + s
  prior_u <- sum(dnorm(x=u, mean=0,sd=exp(mu_star[4]),log=TRUE))
  log.like <- -(sum( dnorm(y,mean=mu,sd=exp(mu_star[3]),log=TRUE)) + prior_u) 
  
  if(is.infinite(log.like)){
    return(10000000)
  }else{
    return(log.like)
  }
  
}    

# Optimise log posterior density value
gamm_L <- function(x,y,z.spline,theta){
  lp.approx <- optim(par=theta, x=x, y=y,z.spline=z.spline, fn=log.post, hessian=TRUE, method = "L-BFGS-B")
  lp.approx
}


##  Find KLD utility when both prior and posterior are MVN
kld.post <- function(x,y,z.spline,theta){
  
  fitOptim <-gamm_L(x,y,z.spline,theta)
  print(fitOptim$convergence)
  mu_star <- c(fitOptim$par)
  u <- rep(0,num.knots)
  fitOptim_u <- optim(par=u, x=x, y=y,z.spline=z.spline, mu_star=mu_star,
                      fn=joint_like, hessian=TRUE, method = "L-BFGS-B")
  print(fitOptim_u$convergence)
  
  mu_post <- c(fitOptim$par[1:2],fitOptim_u$par)
  Sigma_post <- as.matrix(bdiag(list(solve(fitOptim$hessian[1:2,1:2]), 
                                     solve(fitOptim_u$hessian))))
  
  det.out <- 0.5*(trace(iSp%*%Sigma_post) + t(mp-mu_post)%*%iSp%*%(mp-mu_post) - num_par + log(dSigma_prior/det(Sigma_post)))
  det.out
}


# Approximate expected utility
exp.crit <- function(d,B){
  x <- as.matrix(norm_01_mm(d,-1,1))
  z.spline <- cal_z(num.knots,(x))
  crit  <- foreach(j = 1:B, .packages=c("Matrix","doParallel","HRW"),
                   .combine = c,.errorhandling = 'stop',
                   .export = ls(globalenv())) %dopar% {
                     theta <- theta_sim[j,]
                     s <- z.spline %*% u
                     lp <- as.vector(theta[1] + theta[2]*x) + s
                     y <- qnorm(rand_y,mean=lp,sd=0.2)
                     #y <- as.numeric(p > runif(length(x)))
                     kld.post(x,y,z.spline,theta)             
                   }
  return(crit)
}