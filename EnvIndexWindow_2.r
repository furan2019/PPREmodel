#!/usr/bin/env Rscript
#Author: Fran
#Created on: 2023-01-01

##STEP3 Calculate the correlation between index and phenotype

wea_list <- list()
for(i in 1:15)wea_list[[i]] <- read.csv(paste(name_list[i],".csv",sep = ""),
                                        row.names=1,header = T) 
names(wea_list) <- name_list
attach(wea_list)

corr_func_org <- function(col,Day){
  list_cor <- list()
  cor_M_GDD <- data.frame(matrix(data = NA,nrow = Day,ncol = Day))
  cor_M_ptt <- data.frame(matrix(data = NA,nrow = Day,ncol = Day))
  cor_M_ptr <- data.frame(matrix(data = NA,nrow = Day,ncol = Day))
  for(i in 1:Day)
    for(j in i:Day){
      lt_GDD <- c(HN_GDD[i,j],HB_GDD[i,j],BJ_GDD[i,j],LN_GDD[i,j],JL_GDD[i,j])
      lt_ptt <- c(HN_ptt[i,j],HB_ptt[i,j],BJ_ptt[i,j],LN_ptt[i,j],JL_ptt[i,j])
      lt_ptr <- c(HN_ptr[i,j],HB_ptr[i,j],BJ_ptr[i,j],LN_ptr[i,j],JL_ptr[i,j])
      cor_M_GDD[i,j] <- cor(lt_GDD,M_mean[,col])
      cor_M_ptt[i,j] <- cor(lt_ptt,M_mean[,col])
      cor_M_ptr[i,j] <- cor(lt_ptr,M_mean[,col])
    }
  list_cor[[1]] <- cor_M_GDD
  list_cor[[2]] <- cor_M_ptt
  list_cor[[3]] <- cor_M_ptr
  return(list_cor)
}

cor_M_DTT <- corr_func_org(2,120)
cor_M_PH <- corr_func_org(3,120) 
cor_M_EW <- corr_func_org(4,120)
#Jing724F1 population and Zheng58F1 population...#

name_list_cor <- vector()
for(i in c("M","MJ","MZ"))
  for(j in c("FT","PH","EW"))
    for(z in c("GDD","ptt","ptr")){
      name_list_cor <- append(name_list_cor,paste("cor",i,j,z,sep = "_"))
    }

for(i in 1:3){
  write.csv(cor_M_DTT[[i]],paste(name_list_cor[i],".csv",sep=""),row.names=TRUE)
  write.csv(cor_M_PH[[i]],paste(name_list_cor[i+3],".csv",sep=""),row.names=TRUE)
  write.csv(cor_M_EW[[i]],paste(name_list_cor[i+6],".csv",sep=""),row.names=TRUE)
  #Jing724F1 population and Zheng58F1 population...#
}

##STEP4 Draw the heatmap of the correlation between the environment and the phenotypic mean##
##from 0 to 60 days after planting#

library(Rmisc)
library(ggplot2)

corr_func_new <- function(col,Day){
  list_cor <- list()
  cor_M_GDD <- data.frame(matrix(data = NA,nrow = Day,ncol = Day))
  cor_M_ptt <- data.frame(matrix(data = NA,nrow = Day,ncol = Day))
  cor_M_ptr <- data.frame(matrix(data = NA,nrow = Day,ncol = Day))
  for(i in 1:Day)
    for(j in i:Day){
      lt_GDD <- c(HN_GDD[i,j],HB_GDD[i,j],BJ_GDD[i,j],LN_GDD[i,j],JL_GDD[i,j])
      lt_ptt <- c(HN_ptt[i,j],HB_ptt[i,j],BJ_ptt[i,j],LN_ptt[i,j],JL_ptt[i,j])
      lt_ptr <- c(HN_ptr[i,j],HB_ptr[i,j],BJ_ptr[i,j],LN_ptr[i,j],JL_ptr[i,j])
      cor_M_GDD[Day+1-j,i] <- cor(lt_GDD,M_mean[,col])
      cor_M_ptt[Day+1-j,i] <- cor(lt_ptt,M_mean[,col])
      cor_M_ptr[Day+1-j,i] <- cor(lt_ptr,M_mean[,col])
    }
  list_cor[[1]] <- cor_M_GDD
  list_cor[[2]] <- cor_M_ptt
  list_cor[[3]] <- cor_M_ptr
  return(list_cor)
}

cor_M_DTT_60Day <- corr_func_new(2,60)
cor_M_PH_60Day <- corr_func_new(3,60)
cor_M_EW_60Day <- corr_func_new(4,60)

num <- c(60:1)
name_cor_plot <- vector()
for(j in c("Flowering time","Plant height","Ear weight"))
  for(z in c("(GDD)","(PTT)","(PTR)")){
    name_cor_plot <- append(name_cor_plot,paste(j,z,sep = " "))
  }

heatplot_func_judge <- function(X,t,judge=1,cor,mid,q,win_b,win_e){
  cor_60d <- X[[t]]
  cor_60d <- data.frame(cbind(num,cor_60d))
  cor_60d_P <- reshape2::melt(cor_60d,id.vars="num",variable.name="X",value.name="Day")
  if(as.numeric(judge) == -1){
    cor_60d_P$Day[which(cor_60d_P$Day > cor)] <- cor
    cor_60d_P2 <- data.frame(cbind(cor_60d_P,rep(1:60,each=60)))
    colnames(cor_60d_P2) <- c("End","Begin","r","Begin2")
    cor_60d_P2 <- na.omit(cor_60d_P2)
    p <- ggplot(cor_60d_P2,aes(x=Begin2,y=End))+
      geom_tile(aes(fill=r))+
      xlab("Begining of window (days)")+ylab("End of window (days)")+
      scale_fill_gradient2(low = "red", high = "cyan",mid="white",midpoint=mid)+
      theme(panel.background = element_blank(),axis.line = element_line(colour="black"))+
      ggtitle(name_cor_plot[q])+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(family = "sans"),
            legend.position = c(0.9,0.3),
            legend.key.size = unit(6,"pt"))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
      scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
      annotate(geom = "point",x=win_b,y=win_e,colour="black")+
      annotate("segment",x=0,xend=win_b,y=win_e,yend=win_e,colour="black",size=0.5,linetype="dashed")+
      annotate("segment",x=win_b,xend=win_b,y=0,yend=win_e,colour="black",size=0.5,linetype="dashed")
    return(p)
  }
  if(as.numeric(judge) == 1){
    cor_60d_P$Day[which(cor_60d_P$Day < cor)] <- cor
    cor_60d_P2 <- data.frame(cbind(cor_60d_P,rep(1:60,each=60)))
    colnames(cor_60d_P2) <- c("End","Begin","r","Begin2")
    cor_60d_P2 <- na.omit(cor_60d_P2)
    p <- ggplot(cor_60d_P2,aes(x=Begin2,y=End))+
      geom_tile(aes(fill=r))+
      xlab("Begining of window (days)")+ylab("End of window (days)")+
      scale_fill_gradient2(low ="cyan" , high = "red",mid="white",midpoint=mid)+
      theme(panel.background = element_blank(),axis.line = element_line(colour="black"))+
      ggtitle(name_cor_plot[q])+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(family = "sans"),
            legend.position = c(0.9,0.3),
            legend.key.size = unit(6,"pt"))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
      scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
      annotate(geom = "point",x=win_b,y=win_e,colour="black")+
      annotate("segment",x=0,xend=win_b,y=win_e,yend=win_e,colour="black",size=0.5,linetype="dashed")+
      annotate("segment",x=win_b,xend=win_b,y=0,yend=win_e,colour="black",size=0.5,linetype="dashed")
    return(p)
  }
  cor_60d_P2 <- data.frame(cbind(cor_60d_P,rep(1:60,each=60)))
  colnames(cor_60d_P2) <- c("End","Begin","r","Begin2")
  cor_60d_P2 <- na.omit(cor_60d_P2)
  if(as.numeric(judge==0)){
    p <- ggplot(cor_60d_P2,aes(x=Begin2,y=End))+
      geom_tile(aes(fill=r))+
      xlab("Begining of window (days)")+ylab("End of window (days)")+
      scale_fill_gradient2(low = "cyan", high = "red",mid="white",midpoint=mid)+
      theme(panel.background = element_blank(),axis.line = element_line(colour="black"))+
      ggtitle(name_cor_plot[q])+
      theme(plot.title = element_text(hjust = 0.5),
            text = element_text(family = "sans"),
            legend.position = c(0.9,0.3),
            legend.key.size = unit(6,"pt"))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
      scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
      annotate(geom = "point",x=win_b,y=win_e,colour="black")+
      annotate("segment",x=0,xend=win_b,y=win_e,yend=win_e,colour="black",size=0.5,linetype="dashed")+
      annotate("segment",x=win_b,xend=win_b,y=0,yend=win_e,colour="black",size=0.5,linetype="dashed")
    return(p)
  }
  if(as.numeric(judge==2)){
    p <- ggplot(cor_60d_P2,aes(x=Begin2,y=End))+
      geom_tile(aes(fill=r))+
      xlab("Begining of window (days)")+ylab("End of window (days)")+
      scale_fill_gradient2(low = "red", high = "cyan",mid="white",midpoint=mid)+
      theme(panel.background = element_blank(),axis.line = element_line(colour="black"))+
      ggtitle(name_cor_plot[q])+
      theme(plot.title = element_text(hjust = 0.5), 
            text = element_text(family = "sans"),
            legend.position = c(0.9,0.3),
            legend.key.size = unit(6,"pt"))+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
      scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
      annotate(geom = "point",x=win_b,y=win_e,colour="black")+
      annotate("segment",x=0,xend=win_b,y=win_e,yend=win_e,colour="black",size=0.5,linetype="dashed")+
      annotate("segment",x=win_b,xend=win_b,y=0,yend=win_e,colour="black",size=0.5,linetype="dashed")
    return(p)
  }
}

#FT
p1 <- heatplot_func_judge(cor_M_DTT_60Day,1,-1,-0.8,-0.87,1,11,36)
p2 <- heatplot_func_judge(cor_M_DTT_60Day,2,0,NULL,0,2,36,44)
p3 <- heatplot_func_judge(cor_M_DTT_60Day,3,1,0.8,0.85,3,36,44)
#PH 
p4 <- heatplot_func_judge(cor_M_PH_60Day,1,-1,-0.6,-0.75,4,37,46)
p5 <- heatplot_func_judge(cor_M_PH_60Day,2,0,NULL,0,5,13,24)
p6 <- heatplot_func_judge(cor_M_PH_60Day,3,1,0.6,0.75,6,23,34) 
#EW 
p7 <- heatplot_func_judge(cor_M_EW_60Day,1,-1,-0.7,-0.80,7,47,54)
p8 <- heatplot_func_judge(cor_M_EW_60Day,2,2,NULL,0,8,26,35)
p9 <- heatplot_func_judge(cor_M_EW_60Day,3,1,0.4,0.55,9,15,24)

pdf("MaternalWindowHeatmap.pdf",width = 8,height = 8)
multiplot(p1,p4,p7,p2,p5,p8,p3,p6,p9,cols=3)
dev.off()
