rm(list=ls()) #clear environment, start afresh
MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR_JR/"


#useful packages
library("icesDatras")
years <- 1991:2011
quarters <- 1:4

# no function to get sab file yet

getSweptFile <- function(survey, year, quarter) {

  # check survey name
  if (!checkSurveyOK(survey)) return(FALSE)

  # check year
  if (!checkSurveyYearOK(survey, year, checksurvey = FALSE)) return(FALSE)

  # check quarter
  if (!checkSurveyYearQuarterOK(survey, year, quarter, checksurvey = FALSE, checkyear = FALSE)) return(FALSE)

  # read url and parse to data frame
  url <-
    sprintf(
      "https://datras.ices.dk/WebServices/DATRASWebService.asmx/getSAdata?survey=%s&year=%i&quarter=%i",
      survey, year, quarter)
  print(url)
  out <- readDatras(url)
  out <- parseDatras(out)

  out
}

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

