rm(list=ls()) #clear environment, start afresh

#useful packages
library("icesDatras")
years <- 1991:2011
quarters <- 1:4

# Not an ideal approach for downloading the data but the library currently only supports one quarter/year at a time https://github.com/ices-tools-prod/icesDatras/issues/35
i=-1
for(survey in getSurveyList()){
  for(year in years){
    for(quarter in quarters){
    dl <- try(getHHdata(survey,years,quarters); i=i+1)
}



survey_Q_C_S_combinations<-read.csv("R/survey_Q_C_S_combinations.csv")# for QSR


MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR_JR/"
RDIR<- paste0(MAINDIR,"R/")
HHDIR<- paste0(MAINDIR,"HH_HL_download/HH/")
HLDIR<- paste0(MAINDIR,"HH_HL_download/HL/")



