# ChrisLynam@Cefas.co.uk
# update Feb 2019 to calc num hauls by quadrants and area by quad (see Lynam_INDfn_Feb2019_guild.r for code to smooth over 60km radius)
# update 19 Dec 2016 - new shapes for North Sea and differ for BT and OT surveys
# which hauls are in the subregions
#ospar<-readShapeSpatial(paste(SHAPEPATH,"ospar//OSPAR_inner_Boundary.shp",sep=''))#not polys
#coordinates(dhspp) <- ~ ShootLong_degdec +ShootLat_degdec
#ox <- over(dhspp, ospar) #load("C:/Users/cl06/Desktop/biodiv19 ref/RUNuptoOSPARsubdiv.rdata.RData")
#march 2017 add EHDs as used for PP work excluding PMS and UNC_EAST >57.7; UNC_WEST < -2E
#jan 2018 for GNS surveys do not lose inherited samp$ICSStSq with sampstat value from Attributes
#may 2020 added IBTS so can use exchange data
if( !exists("EHDS_PP") ){ EHDS_PP<-F } else { if(EHDS_PP) SAMP_STRAT<-F } # not for standard indicators #mar2017 test for cf to PP 
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

#shapefiles for assessment subdivisions
require(maptools)

if(survey %in% c("GNSIntOT1","GNSIntOT3")){
  SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSIntOT/GNSstrat_Atlantis.shp",sep='') )
  if(EHDS_PP) SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNS_EHDPP/ehu_polygons.shp",sep='') ) 
  if(BYSUBDIV) NAMsubdiv <- "NAME"     #new areas as used for FC/FW3
  #if(BYSUBDIV) NAMsubdiv <- "LFIregion" #old spatial areas - 25 year plan
  NAMsampstrat<-"ICESNAME"
} 


if(survey %in% c("GNSGerBT3")){
  SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSGerBT3/GNSstrat_Atlantis.shp",sep='') ) 
  if(BYSUBDIV) NAMsubdiv <- "NAME"
  NAMsampstrat<-"ICESNAME"
} 

if(survey %in% c("GNSNetBT3", "GNSBelBT3")){
  SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSNetBT3/GNSstrat_Atlantis.shp",sep='') ) 
  if(BYSUBDIV) NAMsubdiv <- "NAME"
  NAMsampstrat<-"ICESNAME"
}

if(survey %in% c("GNSEngBT3","GNSIntOT1_channel") ){ 
  SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"EChanEngBeamSimple/EChanBT3Simple.shp",sep=''))
  NAMsubdiv <- "name"
}
if(survey%in% c("CSEngBT4","CSEngBT1")) { 
  SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"WChanEngBeam/WChanBT4.shp",sep='')) #no rect "CSEngBT4" carhelmar discontinued
  #SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"Q1SWBEAM/2016onwards/CSEngBT1_rect.shp",sep=''))#does not work with over  #strata split by rect
  #NAMsampstrat <- "ICESNAME"
  #Q1SWECOS_Strata3WGS84.shp   #svens#names(SUBDIV)[which(names(SUBDIV)=="Stratum_Ar")]<-"KM2_LAM"#1e6
  #CelticSea_StrataSurveyGrid_20151015_final_Albers.shp #minigrid
  NAMsubdiv <- "Stratum"
}
if(survey=="CSEngBT3") { 
  SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"irish_seaBT//NI_IBTS.WGS84.shp",sep=''))
  NAMsubdiv <- "Features"
}
if(survey=="CSEngBT3_Bchannel"){ 
  SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"BChanEngBeam/BChanBT3.shp",sep=''))
  NAMsubdiv <- "Stratum"
}
if(survey %in% c("CSNIrOT4","CSNIrOT1")){ #SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"irish_seaGFS//NI_IBTS.WGS84.shp",sep=''))#not using as not sampled deepwater
  SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"irish_seaBT//NI_IBTS.WGS84.shp",sep=''))
  NAMsubdiv <- "Features"
}
if(survey=="CSScoOT1"){ SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"SWCQ1.WGS84//SWC_Q1.shp",sep=''))
  NAMsubdiv <- "Name"
}
if(survey=="CSScoOT4"){ SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"SWCQ4.WGS84//SWC_Q4.shp",sep=''))
NAMsubdiv <- "Name"
}
if(survey=="WAScoOT3"){ SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"SWC-RockQ3.WGS84//SWC_Q3.shp",sep=''))
NAMsubdiv <- "Name"
}
if(survey %in% c("BBICnSpaOT1","BBICnSpaOT4") ){ SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"Sp-NGFS.WGS84//Sp_North.WGS84.shp",sep=''))
NAMsubdiv <- "Primary"
}
if(survey=="BBICsSpaOT1"){ SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"Sp-Cadiz.WGS84/Sp_Cadiz.WGS84.shp",sep=''))
NAMsubdiv <- "ESTRATO"
}
if(survey=="BBICsSpaOT4"){ SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"Sp-Cadiz.WGS84/Sp_Cadiz.WGS84.shp",sep=''))
NAMsubdiv <- "ESTRATO"
}
if(survey=="WASpaOT3"){ SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"Sp-PorcGFS.WGS84/Porcupine.WGS84.shp",sep=''))
NAMsubdiv <- "Primary"
}
if(survey=="CSIreOT4"){ SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"IGFS.WGS84//IGFS.WGS84.shp",sep=''))
NAMsubdiv <- "Primary"
}
if(survey == "CSBBFraOT4"){ SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"Fr-EVHOE.WGS84//EVHOE.WGS84.shp",sep=''))
NAMsubdiv <- "STRATE"
}

if(survey %in% c("CSFraOT4","BBICFraBT4")){ 
  SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"Fr-EVHOE.WGS84_original//EVHOE.WGS84.shp",sep=''))
NAMsubdiv <- "STRATE"
}
if(survey=="BBICPorOT4"){ SUBDIV<-readShapeSpatial(paste(SHAPEPATH,"BBICPorOT4/Contour3strata_sector.shp",sep='')); 
NAMsubdiv <- "ID_1" }

if(survey=="GNSFraOT4"){
  SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNSFraOT4/GNSFraOT4_EngBT3Simple.shp",sep='') ) # EChannel
  if(BYSUBDIV) NAMsubdiv <- "name" # combined with depth layers 
  NAMsampstrat<-"FID_GNSFra" # mini grid 0.25deg # sampstrat
}   

#read in area for strata
#if(OVERWITE_SUBDIV){ BYSUBDIV<-F; if(!SAMP_STRAT){ NAMsampstrat<-NAMsubdiv; SAMP_STRAT<-T }}###01Feb2017
surveyread<-survey
if(survey=="GNSBelBT3"){ surveyread<-"GNSIntOT1" } 
if(survey=="GNSIntOT1_channel"){ surveyread<-"GNSEngBT3" } 
if(survey %in% c("BBICFraBT4")){ surveyread<-"CSFraOT4" } 
if(survey %in% c("CSEngBT1")){ surveyread<-"CSEngBT1_attrib_subdivkm2" }
ATTRIB <- read.csv(paste(SHAPEPATH,"attributes/",surveyread,".csv",sep=''))
SAMP_FACT <- "KM2_LAM"
if(SAMP_STRAT){ names(ATTRIB)[which(names(ATTRIB) %in% NAMsampstrat)] <- "sampstrat"; SAMP_FACT <- c(SAMP_FACT, "sampstrat") }
if(BYSUBDIV){ names(ATTRIB)[which(names(ATTRIB) %in% NAMsubdiv)] <- "SurvStratum"; SAMP_FACT <- c(SAMP_FACT, "SurvStratum")} 
if(EHDS_PP){  ATTRIB <- read.csv(paste(SHAPEPATH,"attributes/GNS_EHDPP.csv",sep='') ) 
              names(ATTRIB)[which(names(ATTRIB) %in% NAMsubdiv)] <- "SurvStratum"; SAMP_FACT <- c("KM2_LAM", "SurvStratum") } 
#if(OVERWITE_SUBDIV) dhspp$SurvStratum<-dhspp$sampstrat###01Feb2017
  
ATTRIB <- ATTRIB[,which(names(ATTRIB) %in% SAMP_FACT )]
#area relates to lowest sampling strata (i.e. rects, minigrid or survey strata poly)
#subdiv area - if using by rectangle sampstrat need to sum area for SUBDIV
if(survey %in% c("GNSIntOT1","GNSIntOT3","GNSNetBT3","GNSGerBT3","GNSBelBT3","GNSEngBT3","GNSIntOT1_channel")){
  ATTRIB_SUBDIV <- aggregate(x=ATTRIB$KM2_LAM,by=list(SurvStratum=ATTRIB$SurvStratum), FUN=sum)
  names(ATTRIB_SUBDIV)[2] <- "KM2_LAM"
}

#SUBDIV 
  dhspp <- dhspp[,-which(names(dhspp)=="SurvStratum")]
  dhspp0<- dhspp #altered 17jul2017
  coordinates(dhspp0) <- ~ ShootLong_degdec +ShootLat_degdec #altered 17jul2017
  ox <- over(dhspp0, SUBDIV) #bring in all attributes of location i..e both sampstrat and subdiv if applicable #head(ox)
  ## which are the subdivisions and sampling stratification units
  if(BYSUBDIV) names(ox)[which(names(ox)==NAMsubdiv)] <- "SurvStratum" 
  if(SAMP_STRAT) names(ox)[which(names(ox)==NAMsampstrat)] <- "sampstrat"
  if(EHDS_PP){  names(ox)[which(names(ox)=="area_1")] <- "KM2_LAM" #rename as not in shp correct
                dhspp <- dhspp[!is.na(ox$SurvStratum) & ox$SurvStratum!="Other", ] # some areas were cut for PP
                ox <- ox[!is.na(ox$SurvStratum) & ox$SurvStratum!="Other", ]
  }
  dhspp <- cbind(dhspp,ox[,SAMP_FACT])  ### issue here for GNSGerBT3 41F4 and 41F5 do not agree with samp file for 2001 since long shoot == 5 exactly

#if(QUAD) #replace sampstrat with quads
  ###### feb2019 quadrants
  if(QUAD){
    if(survey %in% c("GNSIntOT1","GNSIntOT3")) SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSIntOT/GNSstrat_Atlantis.shp",sep='') ) 
    if(survey %in% "GNSGerBT3") SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSGerBT3/GNSstrat_Atlantis.shp",sep='') ) 
    if(survey %in% "GNSNetBT3") SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSNetBT3/GNSstrat_Atlantis.shp",sep='') ) 
      
    QUAD_DAT <- SUBDIV@data
    #split rect into 2 x 2 quadrants
    QUAD_DAT$QUADNAME <- QUAD_DAT$ICESNAME
    QUAD_DAT$NS <- QUAD_DAT$NORTH - QUAD_DAT$SOUTH
    QUAD_DAT$EW <- QUAD_DAT$EAST - QUAD_DAT$WEST
    QUAD_DAT$KM2_LAM<- QUAD_DAT$KM2_LAM/4
    
    quads<-expand.grid(1:2,1:2)    
    QUAD_DATA<-list()
    QNAMES<-c("SE","SW","NE","NW")
    for(v in 1:4){ 
      x<-quads[v,1];  y<-quads[v,2]
      QUAD_DATA[[v]] <- QUAD_DAT[,c("QUADNAME","KM2_LAM","SOUTH","NORTH","EAST","WEST","NS","EW")]
      
      QUAD_DATA[[v]]$SOUTH<- QUAD_DATA[[v]]$SOUTH + (QUAD_DATA[[v]]$NS)*(y-1)/2
      QUAD_DATA[[v]]$NORTH<- QUAD_DATA[[v]]$SOUTH + (QUAD_DATA[[v]]$NS)*(1/2)
      
      QUAD_DATA[[v]]$WEST<- QUAD_DATA[[v]]$WEST + (QUAD_DATA[[v]]$EW)*(x-1)/2 
      QUAD_DATA[[v]]$EAST<- QUAD_DATA[[v]]$WEST + (QUAD_DATA[[v]]$EW)*(1/2) 
      
      
      QUAD_DATA[[v]]$QUADNAME<- paste(QUAD_DATA[[v]]$QUADNAME,QNAMES[v],sep="_")
      QUAD_DATA[[v]] <- QUAD_DATA[[v]][,c("QUADNAME","KM2_LAM","SOUTH","NORTH","WEST","EAST")]
    }
    QUADS <- NULL
    for(v in 1:4) QUADS <- rbind(QUADS,QUAD_DATA[[v]])
    rm(QUAD_DATA,QUAD_DAT) #QUADS[QUADS[,1] %in% c(paste("40E7",QNAMES[1:16],sep="")),]
    #determine centre points
    QUADS$cent_lat <- QUADS$SOUTH + (QUADS$NORTH - QUADS$SOUTH)/2 
    QUADS$cent_lon <- QUADS$WEST + (QUADS$EAST - QUADS$WEST)/2 
    
    #replace sampstrat with quads
    #QUADS[ QUADS$SOUTH== unique(QUADS$SOUTH[QUADS$SOUTH==max(QUADS$SOUTH[QUADS$SOUTH<=dhspp[1,]$ShootLat_degdec])]) & 
     #   QUADS$WEST== unique( QUADS$WEST[QUADS$WEST ==max( QUADS$WEST[QUADS$WEST<=dhspp[1,]$ShootLong_degdec])]), ]
    #
    #test<-QUADS[ (QUADS$SOUTH<=dhspp[1,]$ShootLat_degdec & QUADS$NORTH>dhspp[1,]$ShootLat_degdec) & (QUADS$WEST<=dhspp[1,]$ShootLong_degdec & QUADS$EAST>dhspp[1,]$ShootLong_degdec),]
    gridq<-QUADS
    coordinates(gridq) = c("cent_lon", "cent_lat") 
     #coordinates(gridq) = c("EAST", "SOUTH") 
    gridded(gridq) <- TRUE # promote to SpatialPixelsDataFrame #
    gridQ = as(gridq, "SpatialGridDataFrame")
    
    ox <- over(dhspp0, gridQ) #bring in all attributes of location i.e. both sampstrat and subdiv if applicable 
    names(ox)[which(names(ox)=="QUADNAME")] <- "sampstrat"
    #dhspp2 <- data.frame(dhspp0,sampstrat=ox[,"sampstrat"])
    
    dhspp <- cbind(dhspp[,-which(names(dhspp)=="sampstrat")],sampstrat=ox[,"sampstrat"]) 
    #determine distance of haul to centre points
    #should use simon's areas as the SUBDIV
  }
  

  if(!QUAD & (SAMP_STRAT & (survey %in% c("GNSGerBT3", "GNSNetBT3", "GNSIntOT1", "GNSIntOT3"))) ){
    dhspp <- dhspp[,-which(names(dhspp) == "sampstrat")]
    dhspp$sampstrat <- dhspp$ICESStSq
  }
  rm(ox,dhspp0)   #altered 17jul2017#dhspp <- dhspp[,-which(names(dhspp)=="optional")]
  
  #exclude poorly sampled
  if(BYSUBDIV){ dhspp$SurvStratum <- ac(dhspp$SurvStratum)
                #exclude non-strata
                dhspp <- dhspp[!is.na(dhspp$SurvStratum),]
  }
  #survey specific excludes
  if(survey=="CSScoOT4"){
    dhspp <- dhspp[ dhspp$Year>1996,]
    print("Excluded pre 1997 poor sampling")
  }
  if(survey=="BBICPorOT4"){
    #exclude poorly sampled strata
    EXCLUDES<- c("17","21","37","42")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB$SurvStratum %in% EXCLUDES,]
    print("Excluded strata 17, 21, 37 and 42")
    dhspp <- dhspp[ dhspp$Year>2004,]
    print("Excluded pre 2005 poor sampling")
  }
  
  if(survey=="CSBBFraOT4"){
    #exclude northern strata
    EXCLUDES <- c("Cs3","Cs4","Cs5","Cs6","Cs7", "Cc5","Cc6","Cc7", "Cn3","Cn2", "Cc3e","Cc4e","Cc4w","Cc3w")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB$SurvStratum %in% EXCLUDES,]
    print("Excluded northern strata")
    }
  if(survey=="CSFraOT4"){
    #exclude southern strata and "Cs3" not samppled
    EXCLUDES<- c("Gn1","Gn2","Gn3","Gn4","Gn5","Gn6", "Gn7","Gs1","Gs2","Gs3","Gs4","Gs5","Gs6","Gs7", "Cs3"
    ,"Cs6","Cs7","Cc6","Cc7","Cn2e","Cc3w")#deep and missing yrs
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB$SurvStratum %in% EXCLUDES,]
    print("Excluded southern strata, deep and poorly sampled")
    #dhspp <- dhspp[ dhspp$Year>2000,]; print("Excluded pre 2001 poor sampling")
    #dhspp <- dhspp[ dhspp$Year>2001,]
    #print("Excluded pre 2002 poor sampling")
  }
  
  if(survey=="CSScoOT1"){
    #exclude red1_lam	 windsock_lam  edge as poorly sampled
    EXCLUDES<-c("red1_lam")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB$SurvStratum %in% EXCLUDES,]
    print("Excluded red1_lam")
  }
  if(survey=="WAScoOT3"){
    #exclude mylightblue_lam edge as poorly sampled
    EXCLUDES <- "mylightblue_lam" 
    dhspp <- dhspp[dhspp$SurvStratum!=EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB$SurvStratum %in% EXCLUDES,]
    print("Excluded mylightblue_lam")
  }
  if(survey %in% c("CSNIrOT4","CSNIrOT1")){
    #exclude st georges channel as poorly sampled
    EXCLUDES <- "St George's Channel <100m"
    dhspp <- dhspp[dhspp$SurvStratum!=EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB$SurvStratum %in% EXCLUDES,]
    print("Excluded St George's Channel <100m")
  }
  if(survey == "CSIreOT4"){
    #exclude slope as poorly sampled and not technically CSeas ecoregion
    EXCLUDES <- c("VIIj_Slope","VIIb_Slope","VIa_Slope")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB$SurvStratum %in% EXCLUDES,]
    print("Excluded VIIj_Slope, VIIb_Slope and VIa_Slope")
  }
  
  if(survey == "CSEngBT1"){
    #exclude A J and l as poorly sampled 
    EXCLUDES <- c("StratumA","StratumI","StratumJ")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB$SurvStratum %in% EXCLUDES,]
    print("Excluded StratumA, StratumI and StratumJ")
  }
  
  if(survey == "GNSIntOT1_channel"){
    #exclude "UKcoast <25m" as poorly sampled 
    EXCLUDES <- c("UKcoast <25m")
    dhspp <- dhspp[!dhspp$SurvStratum %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB$SurvStratum %in% EXCLUDES,]
    print("Excluded UKcoast <25m")
  }
  
  #outside survey strata?
  if(BYSUBDIV){ PLOTIT<-F #sometimes have hauls outside of strata! not if exclude NA above
  if(PLOTIT){
    #par(bg='light grey')
    with(dhspp,plot(ShootLong_degdec,ShootLat_degdec))
    for(i in 1:length(unique(dhspp$SurvStratum))){
      with(dhspp[dhspp$SurvStratum== unique(dhspp$SurvStratum)[i],],points(ShootLong_degdec,ShootLat_degdec,col=i+1))
    }
  }
  }
  
  #outside sampstrata?
  if(SAMP_STRAT){  PLOTIT<-F #sometimes have hauls outside of strata! not if exclude NA above
  if(PLOTIT){
    with(dhspp,plot(ShootLong_degdec,ShootLat_degdec))
    for(i in 1:length(unique(dhspp$sampstrat))){
      with(dhspp[dhspp$sampstrat== unique(dhspp$sampstrat)[i],],points(ShootLong_degdec,ShootLat_degdec,col=i))
    }
  }
  }
  

 