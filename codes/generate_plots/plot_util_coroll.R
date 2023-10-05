## ---------------------------------------------------------------------------------------------
## Author: Dilishiya Kalpani De Silva
## ---------------------------------------------------------------------------------------------
## Purpose:  Generate FIG S1: KLD utility evaluations when n = 12  and KLD
##           utility evaluations when n = 24 
## ---------------------------------------------------------------------------------------------

# library(NCmisc)
# list.functions.in.file("codes/generate_plots/plot_util_coroll.R", alphabetic = TRUE)

library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(here)
library(scales)
library("viridis")

plot_theory_kld <- function(n,Sigmabeta,k){

  loc <- here("outputs","util_coroll",paste0("n=",n))
  data_all <- NULL
  for(i in 1:length(Sigmae_all)){
    Sigmae <- Sigmae_all[i]
    for(j in 1:length(Sigmau_all)){
      Sigmau <- Sigmau_all[j]
      for(knots in 1:length(k)){
        v_ut <- paste0(loc,"/OUTPUT_","esd",(Sigmae),"_usd",(Sigmau),"_",k[knots])
        v_fileName <- paste(v_ut,"RData",sep=".") 
        load(v_fileName)
        dataset <- data.frame(KLD=util,index = c(1:length(d)), 
                              K=paste0("K = ",k[knots]),
                              vu = paste0("\u03c3\U1D64 = ",Sigmau),
                              sigma_ep = paste0("\u03c3\U03B5 = ",Sigmae))
        data_all <- rbind(data_all,dataset)
      }
    }
  }
  
  data_all$K <- factor(data_all$K, levels = unique(data_all$K) , labels = unique(data_all$K) )
  data_all$vu <- factor(data_all$vu, levels = unique(data_all$vu) , labels = unique(data_all$vu) )
  
  p <- ggplot(data_all,aes(x=index,y=KLD,color=sigma_ep,shape=sigma_ep))+
    geom_point(size=1.25) + geom_line(size=0.6)+
    scale_color_viridis_d() + ggtitle(paste("n =",n),) +
    theme_bw()+xlab("Design Index")+
    facet_wrap(factor(vu)~factor(K),scales = "free_y",
               nrow=length(Sigmau_all),ncol=length(k),strip.position = c("top"))+ 
    scale_x_continuous(breaks = seq(1:(length(k)+1)))+
    #scale_y_continuous(breaks= pretty_breaks(),minor_breaks = NULL)
    theme(strip.text.y.right = element_text(angle = 0),
          strip.background = element_rect(color="black"),
          panel.spacing.y = unit(0.1, "lines"),
          title = element_text(size=10,face="bold"),
          strip.text = element_text(colour = "black",size = 8,face="bold"),
          axis.title = element_text(color = "black",size=10,face="bold"),
          axis.text = element_text(color = "black", size=8),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())+
    rremove("ylab") 
  return(p)
}

Sigmabeta_all <- c(5,10,15)
Sigmae_all <- c(0.1,0.5,1)
Sigmau_all <- c(1,5,10,20,30)
Sigmabeta <- Sigmabeta_all[2]
p <- list()

############## n =12 ############
n <- 12
k <- c(3,4,6,12)
p[[1]] <- plot_theory_kld(n,Sigmabeta,k) 

############## n =24 ############
n <- 24
k <- c(3,4,6,12,24)
p[[2]] <- plot_theory_kld(n,Sigmabeta,k) 

design_plot <- ggarrange(plotlist=p, nrow = 1, ncol=2, common.legend = TRUE, widths = c(5,6),
                          legend = "bottom")

annotate_figure(design_plot, 
                      left = textGrob(expression(bold("KLD")),rot = 90, vjust = 1)
                      #bottom = textGrob(expression(bold("Design Index")))
                )
ggsave("plots/Ex1_coroll.jpeg",width = 1.8*5,height=1.8*5,bg="white")
