## ---------------------------------------------------------------------------------------------
## Author: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Obtain Optimal design points for linear additive model 
##           under different values of σu, K, σε, and n to generate  FIGURE 4. 
## ---------------------------------------------------------------------------------------------

################################ HPC-INDEX #####################################################
indexK <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) # index for different number of knots values
indexU <- as.integer(Sys.getenv("indexU"))          # index for different sigma(u) values
indexE <- as.integer(Sys.getenv("indexE"))          # index for different sigma(epsilon) values
indexN <- as.integer(Sys.getenv("indexN"))          # index for different number of design points
indexB <- as.integer(Sys.getenv("indexB"))          # index for different sigma(beta) values

# indexK <- 1
# indexU <- 1
# indexE <- 1
# indexN <- 1
# indexB <- 1

################################ CODE ##########################################################
# library(NCmisc)
# list.functions.in.file("codes/designs/main_gam_designs.R", alphabetic = TRUE)

library(Matrix)
library(foreach)
library(doParallel)
library(here)

source("codes/functions/basic_functions.R")
source("codes/functions/ace_manually.R")
source("codes/functions/utility_gam.R")

## create a cluster for parallel computing in foreach function
name <- paste0("cluster/cl_",indexK,"_",
               indexU,"_",indexE,"_",indexN,"_",indexB)
cl <- makeCluster(12,outfile=paste0(name,".txt"))
registerDoParallel(cl)

set.seed(1332)
n_all <- c(12,24) 
n <- n_all[indexN]
k <- c(3,4,6,12,24)
num.knots <- k[indexK]
intKnots <- seq(0, 1, length.out = num.knots)[-c(1,num.knots)]

## set up for data generation in expected utility evaluation  

Sigmabeta_all <- c(5,10,15)       # different sigma(u) values
Sigmae_all <- c(0.1,0.5,1)        # different sigma(epsilon) values
Sigmau_all <- c(1,5,10,20,30)     # different sigma(u) values

Sigmabeta <- Sigmabeta_all[indexB]
Sigmae <- Sigmae_all[indexE]
Sigmau <- Sigmau_all[indexU]

prior_mu <- c(0,0,rep(0,num.knots))                    # prior mean vector
prior_sd <- c(rep(Sigmabeta,2),rep(Sigmau,num.knots))  # prior sd vector

num_par <- length(prior_mu)

B <- 1000                                                      # Number of Monte Carlo evaluations
theta0 <- rnorm(B,prior_mu[1],prior_sd[1])                     # Draw beta0 from prior
theta1 <- rnorm(B,prior_mu[2],prior_sd[2])                     # Draw beta1 from prior
u <- matrix(rnorm(B*num.knots,0,Sigmau),nrow=B,ncol=num.knots) # Draw u from prior
theta_sim <- cbind(theta0,theta1)

Sigma_prior <- diag(prior_sd^2)
iSp <- solve(Sigma_prior)
# dSigma_prior <- det(Sigma_prior)
ldSigma_prior <- as.vector(determinant(Sigma_prior)$modulus)
mp <- prior_mu

rand_y <- matrix(runif(B*n),nrow=B,ncol=n)   # to make y deterministic 

## Initial designs
set.seed(425)
d <- lapply(c(1:5), function(a){
  (as.matrix(runif(n,-1,1)))
})

d <- d[[1]]
print(d)

rpt <- 1
sub_optim_u <- vector()
sub_optim_u[rpt] <- mean(exp.crit(d,B))
print(sub_optim_u[rpt])

loc <- here("outputs","designs","GAM",paste0("n=",n))

## run ACE algorithm
repeat{
  rpt <- rpt+1
  print(paste("Iteration", rpt-1))
  outCE <- ACE_krige(utility=exp.crit, start.d=d, B = c(B,B/2), Q = 30,
                     lower = -1 , upper = 1
  )
  d <- outCE[[1]]
  print(d)
  sub_optim_u[rpt] <- outCE[[2]]
  util <- sub_optim_u[rpt]
  print(util)
  if((sub_optim_u[rpt]-sub_optim_u[rpt-1])==0){
    v_ut <- paste0(loc,"/OUTPUT_","esd",(Sigmae),"_usd",(Sigmau),"_",num.knots)
    v_fileName <- paste(v_ut,"RData",sep=".")
    save(list=c("outCE","d","util","sub_optim_u"),file=v_fileName)
    break
  }
  
}

stopCluster(cl)
