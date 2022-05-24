rm(list=ls()) #clear environment, start afresh
MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR_JR/"


#useful packages
library("icesDatras")
years <- 1999:2011
quarters <- 1:4

# Original getflexfile function not working as it points to HH not ff
getFlexFile_fixed <- function(survey, year, quarter) {
library("icesDatras")

  # check survey name
  if (!checkSurveyOK(survey)) return(FALSE)

  # check year
  if (!checkSurveyYearOK(survey, year, checksurvey = FALSE)) return(FALSE)

  # check quarter
  if (!checkSurveyYearQuarterOK(survey, year, quarter, checksurvey = FALSE, checkyear = FALSE)) return(FALSE)

  # read url and parse to data frame
  url <-
    sprintf(
      "https://datras.ices.dk/WebServices/DATRASWebService.asmx/getFlexFiledata?survey=%s&year=%i&quarter=%i",
      survey, year, quarter)
  out <- readDatras(url)
  out <- parseDatras(out)

  out
}


# Whilst not an ideal approach for downloading the data, it works. The library currently only supports one quarter/year at a time https://github.com/ices-tools-prod/icesDatras/issues/35
for(survey in getSurveyList()){
  i=0
  dl = NULL
  for(year in years){
    for(quarter in quarters){
        skip_to_next <- FALSE
        tryCatch(dl <- getFlexFile_fixed(survey,year,quarter), error = function(e) { skip_to_next <<- TRUE})
        if(skip_to_next | class(dl)!="data.frame") { next } else { 
          i=i+1 
          if(i==1){ #first file
            HH = dl} else { # append if not first successful download
              HH = rbind(HH,dl)
            }
        }
    }
  }
  if(i>0){write.csv(HH,paste0(MAINDIR,"HH_HL_download/HH/HH-",survey,".csv"))}
}

