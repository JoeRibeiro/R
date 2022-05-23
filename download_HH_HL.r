rm(list=ls()) #clear environment, start afresh

#useful packages
library("icesDatras")

survey <- "NS-IBTS"
year <- 2019
quarter <- 1

HH <- getHHdata(survey) 


survey_Q_C_S_combinations<-read.csv("R/survey_Q_C_S_combinations.csv")# for QSR


MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR_JR/"
RDIR<- paste0(MAINDIR,"R/")
HHDIR<- paste0(MAINDIR,"HH_HL_download/HH/")
HLDIR<- paste0(MAINDIR,"HH_HL_download/HL/")



