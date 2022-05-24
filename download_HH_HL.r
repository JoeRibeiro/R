rm(list=ls()) #clear environment, start afresh
MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR_JR/"


#useful packages
library("icesDatras")
years <- 1991:2011
quarters <- 1:4

# Whilst not an ideal approach for downloading the data, it works. The library currently only supports one quarter/year at a time https://github.com/ices-tools-prod/icesDatras/issues/35
for(survey in getSurveyList()){
  i=-1
  for(year in years){
    for(quarter in quarters){
        skip_to_next <- FALSE
        tryCatch(dl <- getFlexFile(survey,year,quarter), error = function(e) { skip_to_next <<- TRUE})
        if(skip_to_next | class(dl)!="data.frame") { next } else { 
          i=i+1 
          if(i==0){ #first file
            HH = dl} else { # append if not first successful download
              HH = rbind(HH,dl)
            }
        }
    }
  }
  write.csv(HH,paste0(MAINDIR,"HH_HL_download/HH/HH-",survey,".csv"))
}

