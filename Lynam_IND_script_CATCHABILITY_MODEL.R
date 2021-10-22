# Author: C Lynam, Cefas
# Contact: chris.lynam@cefas.co.uk
# Version: 1 
# Date: May 2020 

# add catchability by length group
if(survey == "CSFraOT4"){
  #Q <- read.table("O:/C5716_FizzyFish/Working_Area/WP3 Food web indicators/lenmodels/CSeaModSpecies.csv",header=T,sep=',',as.is=T);         SPECIES_SUBSET <- Q[,1]
  Q <- read.table("//lowfilecds/CDP/C5716_FizzyFish/Working_Area/WP3 Food web indicators/lenmodels/CSea RATIO SURVEY_MODEL.csv",header=T,sep=',',as.is=T)        
  
  SPECIES_SUBSET <- names(Q[1,-c(1,2)])
  SPECIES_SUBSET <-sub("[.]", " ", SPECIES_SUBSET) #remove the R created .
  Q[1,-c(1,2)] <- SPECIES_SUBSET
}
if(survey == "GNSIntOT1"){
  Q <- read.table("//lowfilecds/CDP/C5716_FizzyFish/Working_Area/WP3 Food web indicators/lenmodels/NSea RATIO SURVEY_MODEL.csv",header=T,sep=',',as.is=T)
  SPECIES_SUBSET <- names(Q[1,-c(1,2)])
  SPECIES_SUBSET <-sub("[.]", " ", SPECIES_SUBSET) #remove the R created .
  Q[1,-c(1,2)] <- SPECIES_SUBSET
}
