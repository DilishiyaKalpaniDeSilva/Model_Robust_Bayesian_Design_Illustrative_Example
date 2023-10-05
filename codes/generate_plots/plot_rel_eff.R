## ---------------------------------------------------------------------------------------------
## Author: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Generate FIGURE 5: Relative efficiency (R(d)) of the optimal GAMM designs compared 
##           to optimal designs for polynomial models
## ---------------------------------------------------------------------------------------------

library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(here)
library("viridis")

plot_ref_eff_facet <- function(n,knots){
  data_all <- NULL
  for(i in 1:length(Sigmae_all)){
    Sigmae <- Sigmae_all[i]
    
    data_all_u <- NULL
    for(j in 1:length(Sigmau_all)){
      Sigmau <- Sigmau_all[j]
      Relative_Efficiency <- NULL
      
      for(m in 1:3){ ## true models
        v_ut <- paste0(loc,"/rel_eff_n",n,"_",m)
        v_fileName <- paste(v_ut,"RData",sep=".") 
        load(v_fileName)
        Relative_Efficiency <- c(Relative_Efficiency,round(rel_effe[[i]][j,],3))
      }
      
      data_u <- data.frame(True_Model,Assumed_Model,Relative_Efficiency,
                           sigma_u = paste0("\u03c3\U1D64 = ",Sigmau),
                           sigma_ep = paste0("\u03c3\U03B5 = ",Sigmae))
      data_all_u <- rbind(data_all_u,data_u)
    }
    
    data_all <- rbind(data_all,data_all_u)
  }
  
  data_all$Assumed_Model<-factor(data_all$Assumed_Model,
                                 levels = c("Linear"  , "Quadratic" ,"Cubic" ,paste0("GAM K=",k)),
                                 labels=c("Linear"  , "Quadratic" ,"Cubic" ,paste0("GAM K=",k)))
  
  data_all$True_Model<-factor(data_all$True_Model,
                                 levels = c("Linear"  , "Quadratic" ,"Cubic" ),
                                 labels=c("Linear"  , "Quadratic" ,"Cubic" ))
  data_all$sigma_u <- factor(data_all$sigma_u, levels = unique(data_all$sigma_u) , 
                             labels = unique(data_all$sigma_u) )
  
  dodge <- position_dodge(.25)
  pf <- ggplot(data_all,aes(x=True_Model, y=Relative_Efficiency,
                        colour=Assumed_Model,group=Assumed_Model,
                        shape=Assumed_Model)) +
    ggtitle(paste("n =",n ) ) + theme_bw() + 
    geom_point(size=1,position = dodge)+ 
    geom_line(size=0.35,position = dodge)+
    scale_color_viridis_d(direction = 1)+
    xlab("True Model") + ylab("Relative Efficiency")+
    facet_grid(factor(sigma_ep)~factor(sigma_u)) + 
    theme(strip.text.y = element_text(angle = 270),
          #strip.background = element_rect(fill="white",color="black"),
          strip.text = element_text(colour = "black",size = 12,face="bold"),
          title = element_text(size=12,face="bold"),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(color = "black",size=12,face="bold"),
          axis.text.x = element_text(color = "black", size=8,angle = 0),
          axis.text.y = element_text(color = "black", size=12,angle = 0),
          legend.text = element_text(color = "black", size=10),
          legend.position = "right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-5,-5,-5,-5)) + 
    labs(color="Assumed Model",shape="Assumed Model")+
    rremove("xlab")+rremove("ylab")
  return(pf)
}
  
p <- list()
Sigmabeta_all <- c(5,10,15)
Sigmae_all <- c(0.1,0.5,1)
Sigmau_all <- c(1,5,10,20,30)
Sigmabeta <- Sigmabeta_all[2]
loc <- here("outputs","rel_eff")

n <- 12
k <- c(3,4,6,12)
Assumed_Model <- rep(c("Linear"  , "Quadratic" ,"Cubic" ,paste0("GAM K=",k)),3)
True_Model <- rep(c("Linear"  , "Quadratic" ,"Cubic"),each=7)

p[[1]] <- plot_ref_eff_facet(n,k)

n <- 24
k <- c(3,4,6,12,24)
Assumed_Model <- rep(c("Linear"  , "Quadratic" ,"Cubic" ,paste0("GAM K=",k)),3)
True_Model <- rep(c("Linear"  , "Quadratic" ,"Cubic"),each=8)

p[[2]] <- plot_ref_eff_facet(n,k)

design_plot2 <- ggarrange(p[[1]],p[[2]], #legend.grob = legend_3,
                          legend = "right",
                          nrow = 2, ncol=1,common.legend = FALSE)

annotate_figure(design_plot2, 
                left = textGrob(expression(bold("Relative Efficiency")),
                                rot = 90, vjust = 1),
                bottom = textGrob(expression(bold("True Model")))
)
ggsave("plots/Ex1_rel_eff.jpeg",width = 1.8*5,height=1.8*5,bg="white")
