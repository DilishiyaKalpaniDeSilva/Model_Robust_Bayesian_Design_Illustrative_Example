## ---------------------------------------------------------------------------------------------
## Author: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Calculate relative efficiency
## ---------------------------------------------------------------------------------------------

################################ HPC-INDEX #####################################################
indexB <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) # index for degree of the polynomial 
indexN <- as.integer(Sys.getenv("indexN"))          # index for different number of design points
method <- "all"                                     # if number of optimal design points are used
                                                    # as (degree of polynomial +1) method = "fixed"
# indexB <- 2
# indexN <- 2
################################ CODE ##########################################################

library(foreach)
library(doParallel)
library(here)

source("codes/functions/cal_z.R")
source("codes/functions/basic_functions.R")
source("codes/functions/utility_polynomial.R")

name <- paste0("cluster/rel_eff")
cl <- makeCluster(12,outfile=paste0(name,".txt"))
registerDoParallel(cl)

set.seed(1332)
num_par <- indexB+1

# Independent normal priors for theta
prior_mu <- c(rep(0,num_par))
prior_sd <- c(rep(10,num_par))
Sigmae <- 0.2

B <- 20000

theta_sim  <- matrix(NA,nrow=B,ncol=num_par,byrow = TRUE)
for(i in 1:num_par){
  theta_sim[,i] <- rnorm(B,prior_mu[i],prior_sd[i])
}

n_all <- c(12,24)
n <- n_all[indexN] 
rand_y <- matrix(runif(B*n),nrow=B,ncol=n)

iSp <- solve(diag(prior_sd^2))
# dSigma_prior <- det(diag(prior_sd^2))
ldSigma_prior <- as.vector(determinant(diag(prior_sd^2))$modulus)
mp <- prior_mu

# Compute relative utility for parametric designs
utilp <- vector()
loc <- here("outputs","designs","polynomial")
for(i in 1:3){
  v_ut <- paste0(loc,"/OUTPUT_",method,"_N",12,"_deg",i)
  v_fileName <- paste(v_ut,"RData",sep=".")
  load(v_fileName)
  if(n==12){
    print(table(d))
    k <- c(3,4,6,12)
  }else{
    d <- rep(d,2)
    print(table(d))
    k <- c(3,4,6,12,24)
  }
  utilp[i] <- mean(exp.crit2(as.matrix(d),B)) 
}

rel_eff <- utilp/utilp[indexB]

# Compute relative utility for GAMM designs 
Sigmabeta_all <- c(5,10,15)
Sigmae_all <- c(0.1,0.5,1)
Sigmau_all <- c(1,5,10,20,30)
Sigmabeta <- Sigmabeta_all[2]

rel_effe <- list()
loc <- here("outputs","designs","GAM",paste0("n=",n))
for(indexE in 1:3){
  rel_effg <- matrix(NA,ncol = (length(k)+3),nrow = length(Sigmau_all))
  for(indexU in 1:length(Sigmau_all)){
    Sigmau <- Sigmau_all[indexU]
    utilg <- vector()
    for (indexK in 1:length(k)) {
      num.knots <- k[indexK]
      v_ut <- paste0(loc,"/OUTPUT_","esd",(Sigmae_all[indexE]),"_usd",(Sigmau),"_",num.knots)
      v_fileName <- paste(v_ut,"RData",sep=".")
      if(file.exists(v_fileName)){
        load(v_fileName)
        utilg[indexK] <- mean(exp.crit2(as.matrix(round(d,3)),B))
        print(round(as.vector(d),3))
      }else{
        print("not exist")
      }
    }
    rel_effg[indexU,] <-  c(rel_eff,(utilg/utilp[indexB]))
  }
  rel_effe[[indexE]] <- rel_effg
}

loc <- here("outputs","rel_eff")
save(list=c("rel_effe"),file=paste0(loc,"/rel_eff_n",n,"_",indexB,".RData"))

stopCluster(cl)
