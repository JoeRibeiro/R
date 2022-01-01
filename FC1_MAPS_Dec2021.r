# ChrisLynam@Cefas.co.uk
# Dec 2021
#shapefiles for assessment subdivisions
#useful packages
library(tidyverse)
library(mgcv) # to investigate indicator trends with gam
library(lattice)
library(maps) # for eyeball of data
library(Hmisc)# for errbars
library(maptools)
require(rgdal) # to read shapefiles for subdivisions
#to add shapefile info to data in Lynam_OSPARsubdiv.r
library(sp)
library(spatstat)
library(rgeos)    
library(stringr) 
library(sf)
library(data.table)

#create one file with all haul locations
READSURV <- F
if(READSURV){
  setwd("C:/Users/cl06/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/out/latlon_haul/hauls/")
  HAULS <- list.files()
  
  ALLHAULS <- bind_rows(lapply(HAULS, 
                               fread, select=c("Year","HaulID","ShootLat_degdec","ShootLong_degdec") ) 
  )
  with(ALLHAULS, plot(ShootLong_degdec,ShootLat_degdec, pch=10,cex=0.01,col="grey",xlab="",ylab="") ); map(add=T)
  # write.csv(ALLHAULS,"ALLHAULS.csv")
  
  #create one file with all locations of FC1 records
  setwd("C:/Users/cl06/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/out/latlon_haul/")
  RECORDS <- list.files()
  RECORDS <- RECORDS[RECORDS!="hauls"]
  ALLRECORDS <- bind_rows(lapply(RECORDS, 
                               fread, select=c("Year","HaulID","SpeciesSciName","lon","lat","N") ) 
  )
  # write.csv(ALLRECORDS,"ALLRECORDS.csv")
} else {
  setwd("C:/Users/cl06/OneDrive - CEFAS/Fish_dataproduct_QSR/FC1_26Nov2021/final")
  ALLHAULS <- read.csv("ALLHAULS.csv")
  ALLRECORDS <- read.csv("ALLRECORDS.csv")
}
#plot by species
setwd("C:/Users/cl06/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/out/")
SPP<-sort(unique(ALLRECORDS$SpeciesSciName))
for( i in SPP){
  png(paste(i,".png",sep=""),width = 480, height = 800)
  with(ALLHAULS, plot(ShootLong_degdec,ShootLat_degdec, pch=19,cex=0.01,col="light grey",xlab="",ylab="", main=i) )
  with(ALLRECORDS[ALLRECORDS$SpeciesSciName==i & ALLRECORDS$Year<2015,], points(lon,lat, pch=19,cex=1,col="black") )
  with(ALLRECORDS[ALLRECORDS$SpeciesSciName==i & ALLRECORDS$Year>=2015,], points(lon,lat, pch=19,cex=0.75,col="red") )
  map(add=T)
  dev.off()
}


setwd("C:/Users/cl06/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/")
SHAPEPATH<-"Strata/"

SUBDIV1 <- rgdal::readOGR(paste(SHAPEPATH,"GNS_rectstrat/GNSIntOT/GNSstrat_Atlantis.shp",sep='') )
SUBDIV2 <- rgdal::readOGR(paste(SHAPEPATH,"GNS_rectstrat/GNSGerBT3/GNSstrat_Atlantis.shp",sep='') ) 
SUBDIV4 <- rgdal::readOGR(paste(SHAPEPATH,"EChanEngBeamSimple/EChanBT3Simple.shp",sep=''))
SUBDIV5<-rgdal::readOGR(paste(SHAPEPATH,"WChanEngBeam/WChanBT4.shp",sep='')) #no rect "CSEngBT4" carhelmar discontinued
SUBDIV6 <- rgdal::readOGR(paste(SHAPEPATH,"BChanEngBeam/BChanBT3.shp",sep=''))
SUBDIV7<-rgdal::readOGR(paste(SHAPEPATH,"irish_seaBT//NI_IBTS.WGS84.shp",sep=''))
SUBDIV8<-rgdal::readOGR(paste(SHAPEPATH,"SWCQ1.WGS84//SWC_Q1.shp",sep=''))
SUBDIV9<-rgdal::readOGR(paste(SHAPEPATH,"SWCQ4.WGS84//SWC_Q4.shp",sep=''))
SUBDIV10<-rgdal::readOGR(paste(SHAPEPATH,"SWC-RockQ3.WGS84//SWC_Q3.shp",sep=''))
SUBDIV11<-rgdal::readOGR(paste(SHAPEPATH,"Sp-NGFS.WGS84//Sp_North.WGS84.shp",sep=''))
SUBDIV12 <- rgdal::readOGR(paste(SHAPEPATH,"Sp-Cadiz.WGS84/Sp_Cadiz.WGS84.shp",sep=''))
SUBDIV13 <- rgdal::readOGR(paste(SHAPEPATH,"Sp-Cadiz.WGS84/Sp_Cadiz.WGS84.shp",sep=''))
SUBDIV14 <- rgdal::readOGR(paste(SHAPEPATH,"Sp-PorcGFS.WGS84/Porcupine.WGS84.shp",sep=''))
SUBDIV15<-rgdal::readOGR(paste(SHAPEPATH,"IGFS.WGS84//IGFS.WGS84.shp",sep=''))
SUBDIV16<-rgdal::readOGR(paste(SHAPEPATH,"Fr-EVHOE.WGS84//EVHOE.WGS84.shp",sep=''))
SUBDIV17<-rgdal::readOGR(paste(SHAPEPATH,"Fr-EVHOE.WGS84_original//EVHOE.WGS84.shp",sep=''))
SUBDIV18<-rgdal::readOGR(paste(SHAPEPATH,"BBICPorOT4/Contour3strata_sector.shp",sep='')); 
SUBDIV19 <- rgdal::readOGR(paste(SHAPEPATH,"GNSFraOT4/GNSFraOT4_EngBT3Simple.shp",sep='') ) # EChannel



with(ALLHAULS, plot(ShootLong_degdec,ShootLat_degdec, pch=19,cex=0.01,col="light grey",xlab="",ylab="", main=i) )
plot(SUBDIV1,col=2,add=T)
plot(SUBDIV2,col=2,add=T)
plot(SUBDIV3,col=2,add=T)
plot(SUBDIV4,col=2,add=T)
plot(SUBDIV5,col=2,add=T)
plot(SUBDIV6,col=2,add=T)
plot(SUBDIV7,col=2,add=T)
plot(SUBDIV8,col=2,add=T)
plot(SUBDIV9,col=2,add=T)
plot(SUBDIV10,col=2,add=T)
plot(SUBDIV11,col=2,add=T)
plot(SUBDIV12,col=2,add=T)
plot(SUBDIV13,col=2,add=T)
plot(SUBDIV14,col=2,add=T)
plot(SUBDIV15,col=2,add=T)
plot(SUBDIV16,col=2,add=T)
plot(SUBDIV17,col=2,add=T)
plot(SUBDIV18,col=2,add=T)
plot(SUBDIV19,col=2,add=T)
with(ALLHAULS, points(ShootLong_degdec,ShootLat_degdec, pch=1,col="grey",cex=0.01) )
map(add=T)



#data
dhspp <- dhspp[dhspp$ShootLat_degdec<= 62,]
if(substr(survey,1,2) == "GN"){
  #OSPAR GREATER NORTH SEA
  dhspp <- dhspp[dhspp$ShootLong_degdec>= neg(5),]
  dhspp <- dhspp[dhspp$ShootLat_degdec> 48,]
} else { #CS and WA
  ##OSPAR CELTIC SEAS - do this by shapefile, split CSScot and EVHOE survey by strata not line of latitude, assign strata to sea that most encompasses
  #if(substr(survey,1,2) == "CS" & !survey %in% c("CSScoOT1","CSScoOT4") ) dhspp <- dhspp[!(dhspp$ShootLong_degdec>neg(5) & dhspp$ShootLat_degdec>58.25),]
  if(survey=="CSFraOT4") dhspp <- dhspp[dhspp$ShootLat_degdec> 48,] # CSFraOT4 like CSBBFraOT4 but only north
  if(substr(survey,1,2) == "BB" | survey=="CSBBFraOT4") dhspp <- dhspp[dhspp$ShootLat_degdec<= 48,]# "CSBBFraOT4" only south
}


#SUBDIV 
  dhspp <- dhspp[,-which(names(dhspp)=="SurvStratum")]
  dhspp0<- dhspp #altered 17jul2017
  coordinates(dhspp0) <- ~ ShootLong_degdec +ShootLat_degdec #altered 17jul2017
  suppressWarnings(proj4string(dhspp0) <- CRS("+init=epsg:4326") ) #Warning message: In proj4string(obj) : CRS object has comment, which is lost in output
  suppressWarnings(proj4string(SUBDIV) <- CRS("+init=epsg:4326") ) #Warning message: In proj4string(obj) : CRS object has comment, which is lost in output
  

  #"GNSBelBT3","GNSEngBT3","GNSIntOT1_channel"
  if((SAMP_STRAT & (survey %in% c("GNSGerBT3", "GNSNetBT3", "GNSIntOT1", "GNSIntOT3"))) ){
    dhspp <- dhspp[,-which(names(dhspp) == "sampstrat")]
    dhspp$sampstrat <- dhspp$ICESStSq
  }

  #exclude poorly sampled
  if(BYSUBDIV){ dhspp$SurvStratum <- ac(dhspp$SurvStratum)
                #exclude non-strata
                dhspp <- dhspp[!is.na(dhspp$SurvStratum),]
  }
  #survey specific excludes
  if(survey=="CSScoOT4" | surveyread=="CSScoOT4"){
    dhspp <- dhspp[ dhspp$Year>1996,]
  }
  if(survey=="BBICPorOT4" | surveyread=="BBICPorOT4"){    #exclude poorly sampled strata
    EXCLUDES<- c("17","21","37","42")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
    dhspp <- dhspp[ dhspp$Year>2004,]
  }
  
  if(survey=="CSBBFraOT4" | surveyread=="CSBBFraOT4"){
    #exclude northern strata
    EXCLUDES <- c("Cs3","Cs4","Cs5","Cs6","Cs7", "Cc5","Cc6","Cc7", "Cn3","Cn2", "Cc3e","Cc4e","Cc4w","Cc3w")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
    }
  if(survey=="CSFraOT4" | surveyread=="CSFraOT4"){
    #exclude southern strata and "Cs3" not sampled
    EXCLUDES<- c("Gn1","Gn2","Gn3","Gn4","Gn5","Gn6", "Gn7","Gs1","Gs2","Gs3","Gs4","Gs5","Gs6","Gs7", "Cs3","Cs6","Cs7","Cc6","Cc7","Cn2e")#deep and missing yrs ,"Cc3w"
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
  }
  
  if(survey=="CSScoOT1" | surveyread=="CSScoOT1"){
    EXCLUDES<-c("red1_lam")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
  }
  if(survey=="WAScoOT3" | surveyread=="WAScoOT3"){
    EXCLUDES <- "mylightblue_lam" 
    dhspp <- dhspp[dhspp$SurvStratum!=EXCLUDES,]
  }
  if(survey == "CSIreOT4"){
    #exclude slope as poorly sampled and not technically CSeas ecoregion
    EXCLUDES <- c("VIIj_Slope","VIIb_Slope","VIa_Slope")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
  }
  
  if(survey == "CSEngBT1" | surveyread=="CSEngBT1"){
    EXCLUDES <- c("StratumA","StratumI","StratumJ")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
  }
  
  if(survey == "GNSIntOT1_channel"){
    EXCLUDES <- c("UKcoast <25m")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
  }
  
  #outside survey strata?
  PLOTIT<-T #sometimes have hauls outside of strata! not if exclude NA above
  if(PLOTIT){
    #par(bg='light grey')
    x11()
    plot(SUBDIV,col=af(SUBDIV$Name))#with(dhspp,plot(ShootLong_degdec,ShootLat_degdec))
    with(dhspp,points(ShootLong_degdec,ShootLat_degdec,col="white",pch=19))
    for(i in 1:length(unique(dhspp$SurvStratum))){
      with(dhspp[dhspp$SurvStratum== unique(dhspp$SurvStratum)[i],],points(ShootLong_degdec,ShootLat_degdec,col=i+1,cex=0.7))
    }
    map(add=T)
    savePlot(filename= paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),"samp_subdiv",".bmp",sep=''),type="bmp")
    dev.off()
  }

  
  #outside sampstrata?
  if(SAMP_STRAT){  PLOTIT<-F #sometimes have hauls outside of strata! not if exclude NA above
  if(PLOTIT){
    plot(SUBDIV,col=af(SUBDIV$Name))#with(dhspp,plot(ShootLong_degdec,ShootLat_degdec))
    with(dhspp,points(ShootLong_degdec,ShootLat_degdec,col="white",pch=19))
    for(i in 1:length(unique(dhspp$sampstrat))){
      with(dhspp[dhspp$sampstrat== unique(dhspp$sampstrat)[i],],points(ShootLong_degdec,ShootLat_degdec,col=i,cex=0.7))
    }
    map(add=T)
    savePlot(filename= paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),"samp_strat",".bmp",sep=''),type="bmp")
    dev.off()
  }
  }
  

 