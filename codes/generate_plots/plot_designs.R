## ---------------------------------------------------------------------------------------------
## Author: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Generate Figure 4: Optimal design points obtained for linear additive model under 
##           different values of σu, K, σε, and n.
## ---------------------------------------------------------------------------------------------

library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(here)
library("viridis")

plot_design_facet <- function(n,knots){
  
  loc <- here("outputs","designs","GAM",paste0("n=",n))
  data_all <- NULL
  for(i in 1:length(Sigmae_all)){
    Sigmae <- Sigmae_all[i]
    for(j in 1:length(Sigmau_all)){
      Sigmau <- Sigmau_all[j]
      for(k in knots){
        v_ut <- paste0(loc,"/OUTPUT_","esd",(Sigmae),"_usd",(Sigmau),"_",k)
        v_fileName <- paste(v_ut,"RData",sep=".")
        if(file.exists(v_fileName)){
          load(v_fileName)
          dataset <- data.frame(design=d,sigma_u = paste0("\u03c3\U1D64 = ",Sigmau),
                                sigma_ep = paste0("\u03c3\U03B5 = ",Sigmae),
                                no_of_knots=k)
          data_all <- rbind(data_all,dataset)
        }else{
          print(v_fileName)
        }
      }
    }
  }
  data_all$no_of_knots <- factor(data_all$no_of_knots, levels = knots)
  data_all$sigma_u <- factor(data_all$sigma_u, levels = unique(data_all$sigma_u) )
  
  pf <- ggplot(data_all, aes(x=design, y=no_of_knots, color=no_of_knots)) +
    xlab("X") + #ylab("Number of Knots (K)") +
    ggtitle(paste("n =",n ) ) +
    geom_point() + theme_bw()+ geom_count() +
    facet_grid(factor(sigma_u)~factor(sigma_ep)) + 
    scale_colour_manual(values = cbbPalette,
                        breaks = knots) +
    theme(strip.text.y = element_text(angle = 270),
        #strip.background = element_rect(fill="white",color="black"),
        strip.text = element_text(colour = "black",size = 12,face="bold"),
        title = element_text(size=12,face="bold"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(color = "black",size=12,face="bold"),
        axis.text.x = element_text(color = "black", size=12,angle = 90),
        axis.text.y = element_text(color = "black", size=12,angle = 0),
        legend.text = element_text(color = "black", size=12),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        #legend.box.margin=margin(-5,-5,-3,-5)
        ) + 
    labs(size="Count",color="K")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE),
           size = guide_legend(order = 1),
           col = guide_legend(order = 0))+
    rremove("ylab") #+ rremove("xlab")+
    
  return(pf)
}

cbbPalette <- viridis(5)

Sigmabeta_all <- c(5,10,15)
Sigmae_all <- c(0.1,0.5,1)
Sigmau_all <- c(1,5,10,20,30)

############## n = 12 ############
p <- list()
n <- 12
knots <- c(3,4,6,12)
Sigmabeta <- Sigmabeta_all[2]
p[[1]] <- plot_design_facet(n,knots) 

############## n = 24 ############
n <- 24
knots <- c(3,4,6,12,24)
Sigmabeta <- Sigmabeta_all[2]
p[[2]] <- plot_design_facet(n,knots)

design_plot2 <- ggarrange(p[[1]],p[[2]], legend = "right",
                          nrow = 1, ncol=2, common.legend = FALSE)

annotate_figure(design_plot2, left = textGrob(expression(bold("Number of Knots (K)")),
                                              rot = 90, vjust = 1)
                )

ggsave("plots/Ex1_design.jpeg",width = 1.8*7,height=1.8*5,bg="white")

