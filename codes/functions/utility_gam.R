source("codes/functions/cal_z_knot_loc_given.R")

##  KLD utility function
kld.post <- function(x,y,z.spline,theta){
  
  C <- cbind(1,x,z.spline)
  com1 <- (1/Sigmae^2)*(t(C)%*%C)
  FIM <- com1+iSp
  Sigma_post <- solve(FIM)
  mu_post <- (1/Sigmae^2)*(Sigma_post%*%t(C)%*%y)
  
  det.out <- 0.5*( trace(iSp%*%Sigma_post) + t(mp-mu_post)%*%iSp%*%(mp-mu_post) - num_par + 
                    ( ldSigma_prior - as.vector(determinant(Sigma_post)$modulus) ) )
  #log(dSigma_prior/det(Sigma_post))
  det.out
}

# Approximate expected utility
exp.crit <- function(d,B){
  x <- as.matrix(norm_01_mm(d,-1,1))
  z.spline <- cal_z_knot_loc(num.knots,x, c(0,1), intKnots)
  crit  <- foreach(j = 1:B, .packages=c("Matrix","doParallel","HRW"),
                   .combine = c,.errorhandling = 'stop',
                   .export = ls(globalenv())) %dopar% {
                     theta <- theta_sim[j,]
                     s <- z.spline %*% u[j,]
                     lp <- theta[1] + theta[2]*x + s
                     y <- qnorm(rand_y[j,],mean=lp,sd=Sigmae)
                     kld.post(x,y,z.spline,theta)             
                   }
  return(crit)
}