#!/usr/bin/env Rscript
#Author: Fran
#Created on: 2023-01-01

##STEP1 Calculate the phenotypic mean of the population##

setwd(".../GITHUB/PPREmodel")

M_all_new <- read.csv("M_pheno.csv",header = TRUE)
for(i in 3:5)M_all_new[,i] <- as.numeric(M_all_new[,i])

loc = c("HN","HeB","BJ","LN","JL")
M_mean <- data.frame(matrix(NA,nrow = 5,ncol = 4))
colnames(M_mean) <- c("loc","M_DTT","M_PH","M_EW")
M_mean$loc <- loc

for(i in 1:5)M_mean$M_DTT[i] <- 
  round(mean(M_all_new$DTT[which(M_all_new$loc==loc[i])],na.rm = TRUE),2)
for(i in 1:5)M_mean$M_PH[i] <- 
  round(mean(M_all_new$PH[which(M_all_new$loc==loc[i])],na.rm=TRUE),2)   
for(i in 1:5)M_mean$M_EW[i] <- 
  round(mean(M_all_new$EW[which(M_all_new$loc==loc[i])],na.rm=TRUE),2) 

##STEP2 Convert weather information to index##

library(readxl)
library("plyr")

weather_func <- function(X){
  X$rDL <- (X$DayLength1 + X$DayLength2) / 2
  X$GDD <- ((X$AverageTemperature1 - 10)+(X$AverageTemperature2 - 10)) / 2
  X$PTT <- X$rDL * X$GDD
  X$PTR <- X$rDL / X$GDD
  X
}
weather_BJ_ori <- read_excel("Weather_5locs.xlsx",sheet=1)
weather_BJ <- weather_func(weather_BJ_ori)
weather_HB_ori <- read_excel("Weather_5locs.xlsx",sheet=2)
weather_HB <- weather_func(weather_HB_ori)
weather_HN_ori <- read_excel("Weather_5locs.xlsx",sheet=3)
weather_HN <- weather_func(weather_HN_ori)
weather_JL_ori <- read_excel("Weather_5locs.xlsx",sheet=4)
weather_JL <- weather_func(weather_JL_ori)
weather_LN_ori <- read_excel("Weather_5locs.xlsx",sheet=5)
weather_LN <- weather_func(weather_LN_ori)

envir_index_func <- function(X){
  list_index <- list()
  index_GDD <- data.frame(matrix(data = NA,nrow = 120,ncol = 120))
  index_ptt <- data.frame(matrix(data = NA,nrow = 120,ncol = 120))
  index_ptr <- data.frame(matrix(data = NA,nrow = 120,ncol = 120))
  for(i in 1:120) 
    for(j in i:120){
      index_GDD[i,j] <- round((sum(X$GDD[i:j]) / (j-i+1)),2)
      index_ptt[i,j] <- round((sum(X$PTT[i:j]) / (j-i+1)),2)
      index_ptr[i,j] <- round((sum(X$PTR[i:j]) / (j-i+1)),2)
    }
  list_index[[1]] <- index_GDD
  list_index[[2]] <- index_ptt
  list_index[[3]] <- index_ptr
  return(list_index)
}

#PlantingTime5.9(9)--HarvestingTime10.19(172)
JL_window <- weather_JL[9:128,]
JL_window$TimeNumber <- c(1:120)
JL_list_index <- envir_index_func(JL_window)

#PlantingTime5.11(11)--HarvestingTime10.2(155)
LN_window <- weather_LN[11:130,]
LN_window$TimeNumber <- c(1:120)
LN_list_index <- envir_index_func(LN_window)

#PlantingTime5.12(12)--HarvestingTime9.29(152)
BJ_window <- weather_BJ[12:131,]
BJ_window$TimeNumber <- c(1:120)
BJ_list_index <- envir_index_func(BJ_window)

#PlantingTime5.29(29)--HarvestingTime10.6(159)
HB_window <- weather_HB[29:148,]
HB_window$TimeNumber <- c(1:120)
HB_list_index <- envir_index_func(HB_window)

#PlantingTime6.10(41)--HarvestingTime10.16(169)
HN_window <- weather_HN[41:160,]
HN_window$TimeNumber <- c(1:120)
HN_list_index <- envir_index_func(HN_window)

name_list <- c("JL_GDD","JL_ptt","JL_ptr",
               "LN_GDD","LN_ptt","LN_ptr",
               "BJ_GDD","BJ_ptt","BJ_ptr",
               "HB_GDD","HB_ptt","HB_ptr",
               "HN_GDD","HN_ptt","HN_ptr")

for(i in 1:3)write.csv(JL_list_index[[i]],paste(name_list[i],".csv",sep=""),
                       row.names = TRUE)
for(i in 4:6)write.csv(LN_list_index[[i-3]],paste(name_list[i],".csv",sep=""),
                       row.names = TRUE)
for(i in 7:9)write.csv(BJ_list_index[[i-6]],paste(name_list[i],".csv",sep=""),
                       row.names = TRUE)
for(i in 10:12)write.csv(HB_list_index[[i-9]],paste(name_list[i],".csv",sep=""),
                         row.names = TRUE)
for(i in 13:15)write.csv(HN_list_index[[i-12]],paste(name_list[i],".csv",sep=""),
                         row.names = TRUE)

