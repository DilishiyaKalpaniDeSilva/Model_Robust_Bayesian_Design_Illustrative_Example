## ---------------------------------------------------------------------------------------------
## Author: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose: Evaluate expected KLD utility to generate FIG S1.
## ---------------------------------------------------------------------------------------------

################################HPC-INDEX#######################################################
indexK <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))# index for different number of knots
indexU <- as.integer(Sys.getenv("indexU"))         # index for different sigma(u) values
indexE <- as.integer(Sys.getenv("indexE"))         # index for different sigma(epsilon) values
indexN <- as.integer(Sys.getenv("indexN"))         # index for different number of design points
indexB <- as.integer(Sys.getenv("indexB"))         # index for different sigma(beta) values
################################ CODE ##########################################################

# library(NCmisc)
# list.functions.in.file("codes/main_scripts/main_util_corollary.R", alphabetic = TRUE)

library(Matrix)
library(foreach)
library(doParallel)
library(here)

source("codes/functions/cal_z.R")
source("codes/functions/basic_functions.R")
source("codes/functions/utility_gam.R")

name <- paste0("cluster/cl_coll",indexK,"_",
               indexU,"_",indexE,"_",indexN,"_",indexB)
cl <- makeCluster(12,outfile=paste0(name,".txt"))
registerDoParallel(cl)

set.seed(1332)
n_all <- c(12,24) # number of design points
n <- n_all[indexN]
k <- c(3,4,6,12,24)
num.knots <- k[indexK]
intKnots <- seq(0, 1, length.out = num.knots)[-c(1,num.knots)]

# Generate data
Sigmabeta_all <- c(5,10,15)
Sigmae_all <- c(0.1,0.5,1)
Sigmau_all <- c(1,5,10,20,30)

Sigmabeta <- Sigmabeta_all[indexB]
Sigmae <- Sigmae_all[indexE]
Sigmau <- Sigmau_all[indexU]

# Independent normal priors for theta and u
prior_mu <- c(0,0,rep(0,num.knots))
prior_sd <- c(rep(Sigmabeta,2),rep(Sigmau,num.knots))

num_par <- length(prior_mu)

B <- 10000
theta0 <- rnorm(B,prior_mu[1],prior_sd[1])
theta1 <- rnorm(B,prior_mu[2],prior_sd[2])

u <- matrix(rnorm(B*num.knots,0,Sigmau),nrow=B,ncol=num.knots)

theta_sim <- cbind(theta0,theta1)
rand_y <- matrix(runif(B*n),nrow=B,ncol=n)

Sigma_prior <- diag(prior_sd^2)
iSp <- solve(Sigma_prior)
# dSigma_prior <- det(Sigma_prior)
ldSigma_prior <- as.vector(determinant(Sigma_prior)$modulus)
mp <- prior_mu

if(n==12){
  d <- list()
  d[[1]] <- rep(seq(-1,1,length.out = 2),6)
  d[[2]] <- rep(seq(-1,1,length.out = 3),4)
  d[[3]] <- rep(seq(-1,1,length.out = 4),3)
  d[[4]] <- rep(seq(-1,1,length.out = 6),2)
  d[[5]] <- rep(seq(-1,1,length.out = 12),1)
}else{
  d <- list()
  d[[1]] <- rep(seq(-1,1,length.out = 2),12)
  d[[2]] <- rep(seq(-1,1,length.out = 3),8)
  d[[3]] <- rep(seq(-1,1,length.out = 4),6)
  d[[4]] <- rep(seq(-1,1,length.out = 6),4)
  d[[5]] <- rep(seq(-1,1,length.out = 12),2)
  d[[6]] <- rep(seq(-1,1,length.out = 24),1)
  
}

util <- vector()
for(i in 1:length(d)){
  util[i] <- mean(exp.crit(d[[i]],B))
}

loc <- here("outputs","util_coroll",paste0("n=",n))

v_ut <- paste0(loc,"/OUTPUT_","esd",(Sigmae),"_usd",(Sigmau),"_",num.knots)
v_fileName <- paste(v_ut,"RData",sep=".")
save(list=c("d","util"),file=v_fileName)

stopCluster(cl)
