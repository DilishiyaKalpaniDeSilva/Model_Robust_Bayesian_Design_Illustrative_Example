## ---------------------------------------------------------------------------------------------
## Author: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Generate Figure 3: Five potential realisations that could be captured by GAM for same 
##           beta values 
## ---------------------------------------------------------------------------------------------

library(ggplot2)

source("codes/functions/cal_z.R")
source("codes/functions/basic_functions.R")

## Function to simulate data from GAM model
wigg_data <- function(n,k,sdu,sde,fixed){
  set.seed(573)
  x1 <- runif(n,-1,1)
  x <- norm_01_mm(x1,-1,1)
  x <- as.matrix((x))
  z.spline <- cal_z(k,x)
  
  prior_mu <- c(0,0,log(sdu),rep(0,k))
  prior_sd <- c(rep(sd_beta,2),0.01,rep(sdu,k))
  
  B <- 5   ## 5 realisations
  theta0 <- rnorm(B,prior_mu[1],prior_sd[1])
  theta1 <- rnorm(B,prior_mu[2],prior_sd[2])
  u <- matrix(rnorm(B*k,0,sdu),nrow=B,ncol=k)
  
  rand_y <- matrix(runif(B*n),nrow=B,ncol=n)
  
  data <- NULL
  for(i in 1:B){
    s <- z.spline %*% u[i,]
    if(fixed){
      theta0 <- -2 ## fixed
      theta1 <- 5  ## fixed
      lp <- as.vector(theta0 + theta1*x) + s
      y <- qnorm(rand_y[i,],mean=lp,sd=sde)
      data <- rbind(data,data.frame(x,Realisation=paste0(i),y=y,n,k=paste0("K = ",k),theta0,theta1,
                                    sdu=paste0("\u03c3\U1D64 = ",sdu),sde))
    }else{
      lp <- as.vector(theta0[i] + theta1[i]*x) + s
      y <- qnorm(rand_y[i,],mean=lp,sd=sde)
      data <- rbind(data,data.frame(x,Realisation=paste0(i),y=y,n,k=paste0("K = ",k),theta0[i],theta1[i],
                                    sdu=paste0("\u03c3\U1D64 = ",sdu),sde))
    }
  }
  
  return(data)
}

## Function to plot the simulated GAMM data
plot_wigg_data_facet <- function(n,knots,sdu_all,sde){
  data_all <- NULL
  for(k in knots){
    for(sdu in sdu_all){
      data_all <- rbind(data_all,wigg_data(n,k,sdu,sde,TRUE))
    }
  }
  
  data_all$k<-factor(data_all$k,levels = unique(data_all$k))
  data_all$sdu<-factor(data_all$sdu,levels = unique(data_all$sdu))
  data_all$Realisation<-factor(data_all$Realisation,levels = unique(data_all$Realisation))
  
  P1 <- ggplot(data_all,aes(x=x,y=y,shape=Realisation,color = Realisation))+
    geom_point(size=1) +
    scale_color_viridis_d() + facet_grid(factor(k)~factor(sdu))+
    theme_bw()+xlab("x")+ylab("y")+
    scale_y_continuous(breaks = seq(-8, 4, by = 4)) +
    theme(panel.grid = element_line(color = "grey",size = 0.3,linetype = 1),
          strip.background = element_rect(fill="grey",color="black"),
          strip.text = element_text(colour = "black",size = 12,face="bold"),
          strip.text.y = element_text(colour = "black",size = 12,face="bold",angle = 270),
          axis.title = element_text(color = "black",size=12,face="bold"),
          axis.text = element_text(color = "black", size=12),
          axis.text.x = element_text(color = "black", size=12,angle = 90),
          legend.title = element_text(size=12,face="bold"),
          legend.text = element_text(color = "black", size=12),
          legend.position = "bottom",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-5,-5,-3,-5))
  
  return(P1)
  
  
}

sdu_all <- c(1,5,10,20,30)
sde <- 0.2
sd_beta <- 10
n <- 96
knots <- c(3,4,6,12,24)

plot_wigg_data_facet(n,knots,sdu_all,sde)
ggsave("plots/Ex1_wigg_data.jpeg",width = 1.8*5,height=1.8*5)
