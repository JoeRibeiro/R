# Joe rewrite to get from one shapefile

surveyread<-survey #for ATTRIBUTES
if( substr(survey,nchar(survey)-4,nchar(survey))=="_hist" ){ surveyread <- substr(survey,1, nchar(survey)-5) } else { surveyread <- survey}

subdiv<-rgdal::readOGR(paste(SHAPEPATH,"combined_strata/combined_simple.shp",sep='')); 
subdiv = subdiv[subdiv$Survey_Acr==survey,]

SAMP_FACT <- c("KM2_LAM", "S_REG","L_REG")

# get areas in km2 when in lambert azimuthal equal area projection
subdiv$KM2_LAM=raster::area(spTransform(subdiv, CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))) / 1000000


 
if( !exists("EHDS_PP") ){ EHDS_PP<-F } else { if(EHDS_PP) BY_SREG<-F } # not for standard indicators #mar2017 test for cf to PP 
dhspp <- dhspp[dhspp$ShootLat_degdec<= 62,]
if(substr(survey,1,2) == "GN"){
  #OSPAR GREATER NORTH SEA
  dhspp <- dhspp[dhspp$ShootLong_degdec>= neg(5),]
  dhspp <- dhspp[dhspp$ShootLat_degdec> 48,]
} else { #CS and WA
  ##OSPAR CELTIC SEAS - do this by shapefile, split CSScot and EVHOE survey by strata not line of latitude, assign strata to sea that most encompasses
  #if(substr(survey,1,2) == "CS" & !survey %in% c("CSScoOT1","CSScoOT4") ) dhspp <- dhspp[!(dhspp$ShootLong_degdec>neg(5) & dhspp$ShootLat_degdec>58.25),]
  if(survey=="CSFraOT4") dhspp <- dhspp[dhspp$ShootLat_degdec> 48,] # CSFraOT4 like BBICFraOT4 but only north
  if(substr(survey,1,2) == "BB") dhspp <- dhspp[dhspp$ShootLat_degdec<= 48,]# "BBICFraOT4" only south
}

#shapefiles for assessment subdivisions
require(maptools)


ATTRIB <- subdiv@data[,SAMP_FACT]
#area relates to lowest sampling strata (i.e. rects, minigrid or survey strata poly)
#subdiv area - if using by rectangle S_REG need to sum area for subdiv
if(survey %in% c("GNSIntOT1","GNSIntOT3","GNSNetBT3","GNSGerBT3","GNSBelBT3","GNSEngBT3","GNSIntOT1_channel", "GNSNetBi3", "GNSIntBi3")){
  ATTRIB_L_REG <- aggregate(x=ATTRIB$KM2_LAM,by=list(ATTRIB[,"L_REG"]), FUN=sum)
  names(ATTRIB_L_REG) <- c("L_REG","KM2_LAM")
}

#subdiv 
  dhspp <- dhspp[,names(dhspp)!="S_REG"]
  dhspp0<- dhspp #altered 17jul2017
  coordinates(dhspp0) <- ~ ShootLong_degdec +ShootLat_degdec #altered 17jul2017
  suppressWarnings(proj4string(dhspp0) <- CRS("+init=epsg:4326") ) #Warning message: In proj4string(obj) : CRS object has comment, which is lost in output
  suppressWarnings(proj4string(subdiv) <- CRS("+init=epsg:4326") ) #Warning message: In proj4string(obj) : CRS object has comment, which is lost in output
  
  ox <- over(dhspp0, subdiv) #bring in all attributes of location i..e both S_REG and subdiv if applicable #head(ox)
  ## which are the subdivisions and sampling stratification units
  # if(BY_LREG) names(ox)[which(names(ox)=="L_REG")] <- "S_REG" 
  # if(BY_SREG) names(ox)[which(names(ox)=="S_REG")] <- "S_REG"
  if(EHDS_PP){  names(ox)[which(names(ox)=="area_1")] <- "KM2_LAM" #rename as not in shp correct
                dhspp <- dhspp[!is.na(ox["S_REG"]) & ox["S_REG"]!="Other", ] # some areas were cut for PP
                ox <- ox[!is.na(ox["S_REG"]) & ox["S_REG"]!="Other", ]
  }
  names(dhspp)[which(names(dhspp)=="L_REG")] <- "duplicatedsurvstratum"

  dhspp <- cbind(dhspp,ox[,SAMP_FACT])  ### issue here for GNSGerBT3 41F4 and 41F5 do not agree with samp file for 2001 since long shoot == 5 exactly

#if(QUAD) #replace S_REG with quads
  ###### feb2019 quadrants
  if(QUAD){
    if(survey %in% c("GNSIntOT1","GNSIntOT3")) subdiv <- rgdal::readOGR(paste(SHAPEPATH,"GNS_rectstrat/GNSIntOT/GNSstrat_Atlantis.shp",sep='') ) 
    if(survey %in% "GNSGerBT3") subdiv <- rgdal::readOGR(paste(SHAPEPATH,"GNS_rectstrat/GNSGerBT3/GNSstrat_Atlantis.shp",sep='') ) 
    if(survey %in% "GNSNetBT3") subdiv <- rgdal::readOGR(paste(SHAPEPATH,"GNS_rectstrat/GNSNetBT3/GNSstrat_Atlantis.shp",sep='') ) 
      
    QUAD_DAT <- subdiv@data
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
    
    #replace S_REG with quads
    #QUADS[ QUADS$SOUTH== unique(QUADS$SOUTH[QUADS$SOUTH==max(QUADS$SOUTH[QUADS$SOUTH<=dhspp[1,]$ShootLat_degdec])]) & 
     #   QUADS$WEST== unique( QUADS$WEST[QUADS$WEST ==max( QUADS$WEST[QUADS$WEST<=dhspp[1,]$ShootLong_degdec])]), ]
    #
    #test<-QUADS[ (QUADS$SOUTH<=dhspp[1,]$ShootLat_degdec & QUADS$NORTH>dhspp[1,]$ShootLat_degdec) & (QUADS$WEST<=dhspp[1,]$ShootLong_degdec & QUADS$EAST>dhspp[1,]$ShootLong_degdec),]
    gridq<-QUADS
    coordinates(gridq) = c("cent_lon", "cent_lat") 
     #coordinates(gridq) = c("EAST", "SOUTH") 
    gridded(gridq) <- TRUE # promote to SpatialPixelsDataFrame #
    gridQ = as(gridq, "SpatialGridDataFrame")
    
    ox <- over(dhspp0, gridQ) #bring in all attributes of location i.e. both S_REG and subdiv if applicable 
    names(ox)[which(names(ox)=="QUADNAME")] <- "S_REG"
    #dhspp2 <- data.frame(dhspp0,S_REG=ox[,"S_REG"])
    
    dhspp <- cbind(dhspp[,-which(names(dhspp)=="S_REG")],S_REG=ox[,"S_REG"]) 
    #determine distance of haul to centre points
    #should use simon's areas as the subdiv
  }

  # #c("GNSIntOT1","GNSIntOT1_channel","GNSIntOT3","GNSNetBT3","GNSGerBT3","GNSBelBT3", "GNSNetBi3", "GNSIntBi3")
  # if(!QUAD & (BY_SREG & (survey %in% c("GNSGerBT3", "GNSNetBT3", "GNSIntOT1", "GNSIntOT3", "GNSNetBi3", "GNSIntBi3"))) ){
  #   dhspp$S_REG <- dhspp$ICESStSq
  # }
  # rm(ox,dhspp0)   #altered 17jul2017#dhspp <- dhspp[,-which(names(dhspp)=="optional")]
  # 
  # #exclude poorly sampled
  # dhspp["S_REG"] <- ac(dhspp["S_REG"])
  #exclude non-strata
  dhspp <- dhspp[!is.na(dhspp["S_REG"]),]
  #survey specific excludes
  if(survey=="CSScoOT4" | surveyread=="CSScoOT4"){
    dhspp <- dhspp[ dhspp$Year>1996,]
    print("Excluded pre 1997 poor sampling")
  }
  if(survey=="BBICPorOT4" | surveyread=="BBICPorOT4"){
    #exclude poorly sampled strata
    EXCLUDES<- c("17","21","37","42")
    
    excluded <- dhspp[dhspp["S_REG"] %in% EXCLUDES,]
    if(nrow(excluded)>0){ print(paste0("losing ",nrow(excluded)," HL rows from ",nrow(dhspp)," when excluding strata")) } else { print("no HL rows lost when exclude strata"); }
    
    dhspp <- dhspp[!dhspp["S_REG"] %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB["S_REG"] %in% EXCLUDES,]
    print("Excluded strata 17, 21, 37 and 42")
    dhspp <- dhspp[ dhspp$Year>2004,]
    print("Excluded pre 2005 poor sampling")
  }
  
  if(survey=="BBICFraOT4" | survey=="BBICFraBT4"){
    #exclude northern strata
    EXCLUDES <- c("Cs3","Cs4","Cs5","Cs6","Cs7", "Cc5","Cc6","Cc7", "Cn3","Cn2", "Cc3e","Cc4e","Cc4w","Cc3w")
    
    excluded <- dhspp[dhspp["S_REG"] %in% EXCLUDES,]
    if(nrow(excluded)>0){ print(paste0("losing ",nrow(excluded)," HL rows from ",nrow(dhspp)," when  when excluding strata")) } else { print("no HL rows lost when exclude strata"); }
    
    dhspp <- dhspp[!dhspp["S_REG"] %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB["S_REG"] %in% EXCLUDES,]
    print("Excluded northern strata")
    }
  if(survey=="CSFraOT4"){
    #exclude southern strata and "Cs3" not sampled
    EXCLUDES<- c("Gn1","Gn2","Gn3","Gn4","Gn5","Gn6", "Gn7","Gs1","Gs2","Gs3","Gs4","Gs5","Gs6","Gs7", "Cs3"
    ,"Cs6","Cs7","Cc6","Cc7","Cn2e")#deep and missing yrs ,"Cc3w"
    
    excluded <- dhspp[dhspp["S_REG"] %in% EXCLUDES,]
    if(nrow(excluded)>0){ print(paste0("losing ",nrow(excluded)," HL rows from ",nrow(dhspp)," when  when excluding strata")) } else { print("no HL rows lost when exclude strata"); }
    
    dhspp <- dhspp[!dhspp["S_REG"] %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB["S_REG"] %in% EXCLUDES,]
    print("Excluded southern strata, deep and poorly sampled")
    #dhspp <- dhspp[ dhspp$Year>2000,]; print("Excluded pre 2001 poor sampling")
    #dhspp <- dhspp[ dhspp$Year>2001,]
    #print("Excluded pre 2002 poor sampling")
  }
  
  if(survey=="CSScoOT1" | surveyread=="CSScoOT1"){
    #exclude red1_lam	 windsock_lam  edge as poorly sampled
    EXCLUDES<-c("red1_lam")
    
    excluded <- dhspp[dhspp["S_REG"] %in% EXCLUDES,]
    if(nrow(excluded)>0){ print(paste0("losing ",nrow(excluded)," HL rows from ",nrow(dhspp)," when  when excluding strata")) } else { print("no HL rows lost when exclude strata"); }
    
    dhspp <- dhspp[!dhspp["S_REG"] %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB["S_REG"] %in% EXCLUDES,]
    print("Excluded red1_lam")
  }
  if(survey=="WAScoOT3" | surveyread=="WAScoOT3"){
    #exclude mylightblue_lam edge as poorly sampled
    EXCLUDES <- "mylightblue_lam" 
    
    excluded <- dhspp[dhspp["S_REG"] %in% EXCLUDES,]
    if(nrow(excluded)>0){ print(paste0("losing ",nrow(excluded)," HL rows from ",nrow(dhspp)," when  when excluding strata")) } else { print("no HL rows lost when exclude strata"); }
    
    dhspp <- dhspp[dhspp["S_REG"]!=EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB["S_REG"] %in% EXCLUDES,]
    print("Excluded mylightblue_lam")
  }
  #if(survey %in% c("CSNIrOT4","CSNIrOT1")){
  #  #exclude st georges channel as poorly sampled
  #  EXCLUDES <- "St George's Channel <100m"
  #  dhspp <- dhspp[dhspp["S_REG"]!=EXCLUDES,]
  #  ATTRIB <- ATTRIB[!ATTRIB["S_REG"] %in% EXCLUDES,]
  #  print("Excluded St George's Channel <100m")
  #}
  if(survey == "CSIreOT4" | surveyread=="CSIreOT4"){
    #exclude slope as poorly sampled and not technically CSeas ecoregion
    EXCLUDES <- c("VIIj_Slope","VIIb_Slope","VIa_Slope")
    
    excluded <- dhspp[dhspp["S_REG"] %in% EXCLUDES,]
    if(nrow(excluded)>0){ print(paste0("losing ",nrow(excluded)," HL rows from ",nrow(dhspp)," when  when excluding strata")) } else { print("no HL rows lost when exclude strata"); }
    
    dhspp <- dhspp[!dhspp["S_REG"] %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB["S_REG"] %in% EXCLUDES,]
    print("Excluded VIIj_Slope, VIIb_Slope and VIa_Slope")
  }
  
  if(survey == "CSEngBT1" | surveyread=="CSEngBT1"){
    #exclude A J and l as poorly sampled 
    EXCLUDES <- c("StratumA","StratumI","StratumJ")
    
    excluded <- dhspp[dhspp["S_REG"] %in% EXCLUDES,]
    if(nrow(excluded)>0){ print(paste0("losing ",nrow(excluded)," HL rows from ",nrow(dhspp)," when  when excluding strata")) } else { print("no HL rows lost when exclude strata"); }
    
    dhspp <- dhspp[!dhspp["S_REG"] %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB["S_REG"] %in% EXCLUDES,]
    print("Excluded StratumA, StratumI and StratumJ")
  }
  
  if(survey == "GNSIntOT1_channel" | surveyread=="GNSIntOT1_channel"){
    #exclude "UKcoast <25m" as poorly sampled 
    EXCLUDES <- c("UKcoast <25m")
    
    excluded <- dhspp[dhspp["S_REG"] %in% EXCLUDES,]
    if(nrow(excluded)>0){ print(paste0("losing ",nrow(excluded)," HL rows from ",nrow(dhspp)," when  when excluding strata")) } else { print("no HL rows lost when exclude strata"); }
    
    dhspp <- dhspp[!dhspp["S_REG"] %in% EXCLUDES,]
    ATTRIB <- ATTRIB[!ATTRIB["S_REG"] %in% EXCLUDES,]
    print("Excluded UKcoast <25m")
  }
  
  #outside survey strata?
  if(BY_LREG){ PLOTIT<-T #sometimes have hauls outside of strata! not if exclude NA above
  if(PLOTIT){
    #par(bg='light grey')
    x11()
    plot(subdiv,col=af(subdiv$Name))#with(dhspp,plot(ShootLong_degdec,ShootLat_degdec))
    with(dhspp,points(ShootLong_degdec,ShootLat_degdec,col="white",pch=19))
    for(i in 1:length(unique(dhspp["S_REG"]))){
      with(dhspp[dhspp["S_REG"]== unique(dhspp["S_REG",i]),],points(ShootLong_degdec,ShootLat_degdec,col=i+1,cex=0.7))
    }
    map(add=T)
    savePlot(filename= paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),"samp_subdiv",".bmp",sep=''),type="bmp")
    dev.off()
  }
  }
  
  #outside S_REGa?
  if(BY_SREG){  PLOTIT<-F #sometimes have hauls outside of strata! not if exclude NA above
  if(PLOTIT){
    plot(subdiv,col=af(subdiv$Name))#with(dhspp,plot(ShootLong_degdec,ShootLat_degdec))
    with(dhspp,points(ShootLong_degdec,ShootLat_degdec,col="white",pch=19))
    for(i in 1:length(unique(dhspp$S_REG))){
      with(dhspp[dhspp$S_REG== unique(dhspp$S_REG)[i],],points(ShootLong_degdec,ShootLat_degdec,col=i,cex=0.7))
    }
    map(add=T)
    savePlot(filename= paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),"samp_strat",".bmp",sep=''),type="bmp")
    dev.off()
  }
  }
  