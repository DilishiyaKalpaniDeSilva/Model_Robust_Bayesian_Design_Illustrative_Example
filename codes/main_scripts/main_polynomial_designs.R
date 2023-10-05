## ---------------------------------------------------------------------------------------------
## Author: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Run polynomial designs for Example 1
## ---------------------------------------------------------------------------------------------

################################ HPC-INDEX #####################################################
indexB <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) # index for degree of the polynomial 
indexN <- as.integer(Sys.getenv("indexN"))          # index for different number of design points
method <- Sys.getenv("log")                         # if number of optimal design points are used
                                                    # as (degree of polynomial +1) method = "fixed"
# indexB <- 1
# indexN <- 1
# method <- 1

################################ CODE ##########################################################
# library(NCmisc)
# list.functions.in.file("codes/main_scripts/main_polynomial_designs.R", alphabetic = TRUE)

library(foreach)
library(doParallel)
library(here)

source("codes/functions/basic_functions.R")
source("codes/functions/ace_manually.R")
source("codes/functions/utility_polynomial.R")

## create a cluster for parallel computing in foreach function
name <- paste0("cluster/parametric_analytical",method,indexB)
cl <- makeCluster(12,outfile=paste0(name,".txt"))
registerDoParallel(cl)

set.seed(1332)
num_par <- indexB+1
n_all <- c(12,24)
n <- n_all[indexN]

## set up for data generation in expected utility evaluation  

prior_mu <- c(rep(0,num_par))   # prior mean vector
prior_sd <- c(rep(10,num_par))  # prior sd vector
Sigmae <- 0.2

B <- 1000                      # Number of Monte Carlo evaluations
theta_sim  <- matrix(NA,nrow=B,ncol=num_par,byrow = TRUE)
for(i in 1:num_par){
  theta_sim[,i] <- rnorm(B,prior_mu[i],prior_sd[i])
}

rand_y <- matrix(runif(B*n),nrow=B,ncol=n) # to make y deterministic

iSp <- solve(diag(prior_sd^2))
# dSigma_prior <- det(diag(prior_sd^2))
ldSigma_prior <- as.vector(determinant(diag(prior_sd^2))$modulus)
mp <- prior_mu

# Initial design
if(method=="fixed"){ ## as there are only (indexB+1) of optimal design points for polynomial of degree (indexB+1)
  set.seed(45)
  d <- lapply(c(1:5), function(a){
    (as.matrix(runif((indexB+1),-1,1)))
  })
  exp.crit <- exp.crit1
}else{
  set.seed(45)
  d <- lapply(c(1:5), function(a){
    (as.matrix(runif(n,-1,1)))
  })
  exp.crit <- exp.crit2
}

d <- d[[1]] ## change to get same result of different initial design 
print(d)

rpt <- 1
sub_optim_u <- vector()
sub_optim_u[rpt] <- mean(exp.crit(d,B))
print(sub_optim_u[rpt])

loc <- here("outputs","designs","polynomial")

repeat{
  rpt <- rpt+1
  print(paste("Iteration", rpt-1))
  outCE <- ACE_krige(utility=exp.crit, start.d=d, B = c(B,B), Q = 30,
                     lower = -1 , upper = 1
  )
  d <- outCE[[1]]
  print(d)
  sub_optim_u[rpt] <- outCE[[2]]
  util <- sub_optim_u[rpt]
  print(util)
  if(sub_optim_u[rpt]-sub_optim_u[rpt-1]<=0){
    v_ut <- paste0(loc,"/OUTPUT_",method,"_N",n,"_deg",indexB)
    v_fileName <- paste(v_ut,"RData",sep=".")
    save(list=c("outCE","d","util","sub_optim_u"),file=v_fileName)
    break
  }
}

stopCluster(cl)
