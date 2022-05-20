# Author: C Lynam, Cefas
# Contact: chris.lynam@cefas.co.uk
# Version: 6 (exchange inc) 
# Date: May 2020 new guilds and spatial smooth inc
# reads shapefiles for subdivisions and includes spatial area of strata in calcs
# can edit to run old ns areas by changing subdiv name in subdiv.r
# link to loess ploting code for subdivs complete
# with option of bootstapping by haul and subdiv and trendline fitting
# for pelagics can use volume correction
# calc average densities across quadrants - to do: hauls in a 60km radius of quadrant mid-points 
# May 2020 use exchange data once again
# OCT 2021 - new HH output from SWEPT AREA ASSESSMENT OUTPUT [Klaas] to combine with HL data from exchange
# NOV 2021 - define SSA (code by Joseph Ribeiro)
# DEC 2021 - modify application of 'durraise' for 'S' subsampled data, 'R' raw and C' CPUE data in Lynam_IND_script_process_exchange_HL2021.R to produce numbers per hour and numbers per swept area

rm(list=ls()) #clear environment, start afresh

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
library(mapplots)

#some shorthand
an <- as.numeric
af <- as.factor
ac <- as.character
h <- head
u <- unique
count<- function(DATA){length(unique(DATA))}
neg<-function(x) -x 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#location of the data and subscripts

#RDIR<- dirname(parent.frame(2)$ofile) #Where you have saved the folder called R. Note this will only work if the file is SOURCED, not if it is run in the console. Alternatively, please define your WD
#MAINDIR<- paste0(strsplit(RDIR,"/R")[[1]],"/")
#RDIR = paste0(RDIR,"//")

MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR_JR/"
RDIR<- paste0(MAINDIR,"R/")
#definedSSA = sf::st_read(paste0(RDIR,"rectanglesICESdl29oct2021/shp_dir_for_original/SSAspatial.shp")) # read.csv(paste0(RDIR,"/defined_SSA.csv"))
definedSSA = sf::st_read(paste0(RDIR,"rectanglesICESdl29oct2021/shp/SSAspatial.shp")) # read.csv(paste0(RDIR,"/defined_SSA.csv"))
definedSSA <- sf:::as_Spatial(st_zm(definedSSA))
# plot(definedSSA)
# plot(definedSSA[definedSSA$survey=="CSScoOT4_hist",],col=4);map(add = T)

#Trophic Level data also requires update for MTL analyses

#location of the subscripts in function area
 PROC_SCRIPT<- "//lowfilecds/function/Eco Indicators/DATRASoutput/" #for HH and HL processing scripts incl strata by survey
 SHAPEPATH<-paste("Strata/",sep="") #used in Lynam_OSPARsubdiv.r
 #SUBSCRIPTS_TRAITS<-paste(PROC_SCRIPT,"MarScot/INDscripts2021/",sep="")#max length
 
 setwd(MAINDIR)
 #read subscripts and traits
 # PRODDAT<-read.csv(paste(SUBSCRIPTS_TRAITS,"SpeciesAtLength_Productivity.csv",sep=""))
 SPPFILE<-"R/Species_List_Final_19Jan2022.csv" # originally prepared by Meadhbh Moriarty and Simon Greenstreet for IA02017 shared in ICES WGBIODIV. FC1 identified added by CL based on WKABSENS 2021
 
 LW<-read.csv(SPPFILE)
 names(LW)[2]<-"ScientificName_WoRMS"
 LW$SpeciesSciName <- LW$ScientificName_WoRMS
 LW$a <- LW$LWRa
 LW$b <- LW$LWRb
 LW$"Max.L..cm." <- LW$MaxL
 
 FC1Sp<-LW[LW$SensFC1=="FC1",2] #which species are used for sensitive species indicator
 SPPLIST <- LW$ScientificName_WoRMS
 trait_MAXL<-LW[,which(names(LW) %in% c("ScientificName_WoRMS","SpeciesSciName","Loo","Max.L..cm.","MaxL","Lm","Order","Group"))]
 
#where save output?
OUTPATHstem<-paste(MAINDIR,"out/",sep="")


## choices for analyses upfront

#do you want to write outputs along the way and save workspace
WRITE <- T #save csvs as we go?
  WRITE_LDs <- T #write Length distributions by species/year/subdiv? 
BOOTSTRAP <- F # invoke slow slow code? if F next 3 lines redundant
  NBOOT<- 9; 
  B<-0 # restart counter as only output LD once before boostrap starts i.e. when B=0 
  WRITE_BOOT <- F # every bootstrap dataset and indicator output
  SAVE <- F # save workspace (after bootstrap)
SSA_WRITE_NEW <- T #create standard sampling area from rects
SSA_WRITE_ONLY <- F #skip analysis and output SSA shp only
FILTER_COUNTRY <- F

# Catchability for general species groups
CATCHABILITY_COR_MOD<-SPECIES_IN_MOD_ONLY <-F # for comparison to ewe or lemans'
CATCHABILITY_COR_WALKER<- F # read estimates from nsea paper for q##problem somewhere looking for sweptbefore when this false

#which indicators?
  TyL_GeoM <- F # OSPAR FW3
  TyL_SPECIES <- F #ALSO PROVIDE BY SPECIES MEAN LENGTH
  MaxL <- F # OSPAR FC3
  Loo <- F # alt for OSPAR FC3
  Lm <- F # 
  MEANTL <- F #similar to OSPAR FW4 # not will not be calc'd for WAsurveys as no data file for TL
  MeanL <- F #not OSPAR but simple
  LFI_NULL <- F # for no group/guild calc 25% biom thresh
  LFI <- T
  FINALPLOTS<- T#create indicator plots with smooths
  IEO_FW4 <- F #keep inverts
  
  #Deciding on the species to be included in the analysis is the first step in calculating each survey LFI time series.
  #The FishBase website (www.fishbase.org) provided an 'ecotype' classification for all species encountered. 
  #here in 2021...
  #The LFI has been designed as an indicator of size composition within demersal fish communities,
  #so species assigned to the Pelagic and Bathypelagic ecotypes were automatically excluded, 
  #Species belonging to the Benthopelagic ecotype were generally included, but with the exceptions of 
   #Clupeidae+Clupea harengus, Dysomma brivirostre, 
   #Hyperoplys, Anguilla, Ammodytidae (family-level ID code), Gymnammodytes
   #Salmo (genus-level ID code), Salmo salar and Sarpa salpa, which were all excluded
  #on the basis that these species are relatively poorly sampled by the survey gear.
  
  #SPPFILE is read just after running Lynam_IND_script_process_exchange_HL2021.r 
  #and removes rare species beyound the 550
  
  #subscripts for indicators 
  BYGROUP<-F; if(BYGROUP){ BYGUILD<-F } # elasmos etc 
  BYICESGROUP<-F; if(BYICESGROUP){ BYGUILD<-F } # apex predators etc
  BYGUILD <- F; if(BYGUILD){ BYICESGROUP<-F; BYGROUP<-F }#not run BOOTSTRAP
  FILL_GUILD<-F #if false have a no-guild group
  
# Spatial analyses
AREASCALE <- T # #raise LD data and CPUE and Catch_Biomass by area e.g rectangle KM2_LAM of lowest resolution of sampling strategy/subdiv in Lynam_INDfn_Nov2021.r then x'scale' for missing subdiv area
EHDS_PP <-F #ecohydrodynamic zones
STRATA <- F  #for output
QUAD_NS <- F #only north sea
QUAD_SMOOTH_NS<- F # repeats per guild at present must change - SLOW!
QUADS<-NULL #created by Lynam_OSPARsubdiv

##NOTE additional choices below for plotting
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd(RDIR)

#additional functions
source("required.funcs.r")              # using tapply.ID
source("Lynam_INDfn_Dec2021.r") #Feb2019 now use QUAD to average biomass prior to indicators # Dec2017 update to include TAXA_GROUPINGS, Jan to calc Loo and Lm through Lynam_INDfn_Jan2018_Mtrait.r; Mar to include BX020_guilds
source("Lynam_INDPLOTFN_Nov2018.r")             # plot options 
if(BOOTSTRAP) source("Lynam_IND_BOOTfn_Aug_2017 - OSPAR.r")  # bootstrap the hauls by STSQ and subdiv

if(LFI) source("Lynam_INDfn_Dec2021_LFI.r")   # Large Fish Indicator # Nov2016 no longer fall over if no fish above LFI_THRESHOLD
 if(TyL_GeoM) source("Lynam_INDfn_Dec2021_TyL_GeoM.r")#Geometric mean length (Typical Length cm)
if(Lm & !MaxL) source("Lynam_INDfn_Dec2021_GeoMtrait.r")  # Geometric MaxL, Geometric L.infinity, Geometric L.maturity
if(MaxL | Loo) source("Lynam_INDfn_Dec2021_Mtrait.r")  # mean MaxL, mean L.infinity, mean L.maturity
if(MEANTL) source("Lynam_INDfn_Dec2017_MeanTL.r")# TL output by rectangle and year 
 if(MeanL) source("Lynam_INDfn_Dec2017_MeanL.r") # Mean Length cm
 #if(MaxL | Loo | Lm) source("Lynam_INDfn_Jan2018_Mtrait.r")  # MaxL, Mean L.infinity, Mean L.maturity
#also sourced below:
## Lynam_OSPARsubdiv_Jan2022.r
## Lynam_IND_script_process_exchange_HL2021.R
## and optionally:
## Lynam_IND_script_CATCHABILITY_MODEL.R
## Lynam_IND_script_MAKEGUILDIND.R
## Lynam_IND_script_FINALPLOTS_Dec2021.R 
## Lynam_IND_script_BOOTSTRAP.r

#max length observed in ospar dataproduct by species
#Hyperoplus immaculatus 'Greater sand-eel' were missing -> recorded in GNS IBTS Q1 but without length  treat as Ammodytidae 'sandeel'
#trait_file <- paste("traits_by_species_Mar2019.csv",sep='')#incl elasmo taxa
#trait_MAXL <- read.csv(trait_file)
#can get DEMPEL from SPP_FILE
#trait_MAXL<- trait_MAXL[,-which(names(trait_MAXL)=="DEMPEL")] 
#trait_MAXL<- merge(trait_MAXL, LW, by.x="SpeciesSciName", by.y="ScientificName_WoRMS")
#rm(trait_file)

if(BYGUILD){#2020 paper
  guild_file <- paste("feeding.guilds_05.05.19_processed.csv",sep='')#incl trophic guild ID
  guild_dat <- read.csv(guild_file)
  guild_dat$fguild <- (guild_dat$F.guild)  #fguild_7_soras.numeric
  #guild_dat<- guild_dat[guild_dat$data=="fguild6",]#1:6
}
#trophic level info
TLceltic <- TLnorth <- TLwest <- NULL
if(MEANTL==T){
  #cnrs mtl fw4
  TLceltic <- read.csv("TL_complete14102016-1.csv")
  names(TLceltic)[1] <- "SpeciesSciName"
  TLceltic <- TLceltic[,c("SpeciesSciName","TL")]
  
  # west of scotland
  TLwest <- NULL
  
  # north sea
  TLnorth <- read.csv("refspp_LWmerged_edit06Nov2015.csv")
  names(TLnorth)[2] <- "SpeciesSciName"
  ECOPATH_TL <- F; if(ECOPATH_TL)  TLnorth$TL <- TLnorth$TL_Ecopath #for comparison with ewe
  TLnorth <- TLnorth[,c("SpeciesSciName","TL")]
}

# general species groups # Jan 2018 added 2 stingray Dasyatis spp to lookup and Cetorhinus maximus	Basking shark and Alopias vulpinus to SciName2GroupEff to stop error
 SciName2GroupEff <- read.csv("SciName2GroupEff_Jan2018.csv")
 SciName2GroupEff$Group <- paste("GRP",SciName2GroupEff$Group, sep="")
 SciName2GroupEff$sciName<-as.character(SciName2GroupEff$sciName)
 substr(SciName2GroupEff$sciName,1,1) <- toupper(substr(SciName2GroupEff$sciName,1,1))#correct case
##select only those species from LeMans model and correct relative biomasses to match: Q
if(CATCHABILITY_COR_MOD | SPECIES_IN_MOD_ONLY) source("Lynam_IND_script_CATCHABILITY_MODEL.R")

 # for making SSA
SSAdf=data.frame('rectangle'='','survey'='')

#### indicators survey loop ####
setwd(MAINDIR)
survey_Q_C_S_combinations<-read.csv("R/survey_Q_C_S_combinations.csv")# for QSR
survey_Q_C_S_combinations[,1:8]
for(combrow in 1:nrow(survey_Q_C_S_combinations)){#skipping the inshore surveys
#for(combrow in 27){  # combrow <- 1
  combs=survey_Q_C_S_combinations[combrow,]
  #combs
  QUARTER=combs$Quarter
  COUNTRY=combs$Country
  SEA=combs$Sea
  survey=combs$Surveynam1
  survey_alt_name=combs$Surveynam2
  LFI_THRESHOLD = combs$LFI_threshold
  GEAR= combs$Gear
  BY_SREG= combs$BY_SREG # SMALLER REGIONS?
  BY_LREG = combs$BY_LREG    # AND THEN AGAIN BY LARGER REGIONS?
  
  SPECIES= combs$Species
  if(IEO_FW4) SPECIES<- "ALL"
  
  FIRSTYEAR = combs$First_year
  LASTYEAR = combs$Last_year #note line 401: if(survey=="BBICnSpaOT4" & FIRSTYEAR< 2017 & LASTYEAR > 2017) bio <- bio[bio$Year>=2017,] #no valid data 2010-2016
  STDGEAR = combs$Std_gear
  GEARSUBSCRIPTS = combs$Std_gear_subscripts
  if(is.na(GEARSUBSCRIPTS)) GEARSUBSCRIPTS<-""
  GEARSUBSCRIPTS = unlist(str_split(GEARSUBSCRIPTS, ", "))  
  MINTDUR=combs$MinTrawlDuration
  MAXTDUR=combs$MaxTrawlDuration
  if(BYGUILD) SPECIES <- c("ALL")#,"DEM"); 

  SAMP_FILE=paste0("HH//",combs$hhfilename)
  BIOL_FILE=paste0("HL//",combs$hlfilename)
  
  print(survey) #new OSPAR name as in Mar Scot dataproduct
  if(survey=="CSEngBT4"){print("commercial FV Carhelmar survey stopped - skipping"); next} 
  #add a directory for files out
  if(!file.exists(paste(OUTPATHstem,survey,sep=''))) dir.create(file.path(OUTPATHstem, survey))
  OUTPATH <- paste(OUTPATHstem,survey,"/",sep='')
 
   
  #split ICES Rect into Quadrants? with/out smoothing?
  if(SEA == "GNS"){#,
    QUAD<- QUAD_NS #only north sea
    QUAD_SMOOTH<-QUAD_SMOOTH_NS
  } else { #CS, BB and WA
    QUAD<- F #only north sea
    QUAD_SMOOTH<-F
  }
  #Trophic Level FW4 check
  if(MEANTL==T & !(substr(survey,1,2) %in% c("CS","BB", "GN","IB")) ) { MEANTL<-F; print("no data for TL for WA surveys, deleting MTL selection") }
  
  
  ##load data
  #### sampling data #### 
  samp <- read.table(SAMP_FILE ,as.is = c(1,2,4,5,6,10,11,23),header = TRUE,sep=",") 
  samp$StatRec <- ices.rect2(samp$ShootLon,samp$ShootLat) #sometime StatRec has been excel-ed into a number
  if(substr(survey,7,8)!="Bi") samp <- samp[samp$Quarter==QUARTER,] #young fish surveys are q3-4
  # where two different surveys use the same input data files the correct preprocessing needs to be done to define SSA
  # need to split SEA DATA FOR BTS Q3 ENG'
  #otter surveys
  if(survey %in% "BBICFraOT4") samp<-samp[samp$ShootLat<= 48,]
  if(survey %in% "CSFraOT4") samp<-samp[samp$ShootLat> 48,]
  if(survey %in% "GNSIntOT1_channel") samp<-samp[samp$ShootLat<= 51,]
  #beam surveys
  if(survey %in% "CSEngBT3_Bchannel") samp<-samp[samp$ShootLat<= 52 & samp$ShootLong< -3,]
  if(survey %in% "GNSEngBT3") samp<-samp[samp$ShootLong>= -2 & samp$ShootLat>= 49.5,]
  if(survey %in% "CSEngBT3") samp<-samp[samp$ShootLong< -3 & samp$ShootLat>= 52 & samp$ShootLat< 56,]
  # BTS is downloaded as all countries in one file
  if(GEAR == "BEAM"){
    samp<-samp[samp$Country==COUNTRY,]
    samp$DoorSpread <- samp$WingSpread <- samp$BeamWidth
    samp$SweptAreaDSKM2 <- samp$SweptAreaWSKM2 <- samp$SweptAreaBWKM2 
    
    if(survey %in% c("CSEngBT3","CSEngBT1") ){
      samp$DoorSpread[is.na(samp$DoorSpread) | samp$DoorSpread<2] <- 4
      samp$WingSpread[is.na(samp$WingSpread) | samp$WingSpread<2] <- 4
    }
    
    if(survey_alt_name == "BTS"){ #not french BT in BoB/IC
      samp$SweptAreaWSKM2[is.na(samp$SweptAreaWSKM2)] <- samp$WingSpread[is.na(samp$SweptAreaWSKM2)]*samp$Distance[is.na(samp$SweptAreaWSKM2)]
      samp$SweptAreaDSKM2[is.na(samp$SweptAreaDSKM2)] <- samp$DoorSpread[is.na(samp$SweptAreaDSKM2)]*samp$Distance[is.na(samp$SweptAreaDSKM2)]
      samp$WingSwpArea_sqkm[is.na(samp$WingSwpArea_sqkm)] <- samp$SweptAreaWSKM2[is.na(samp$WingSwpArea_sqkm)]*samp$Distance[is.na(samp$WingSwpArea_sqkm)]
    }
  }
  if(substr(survey,7,8)=="Bi"){
    samp$WingSpread <- an(substr(STDGEAR,nchar(STDGEAR),nchar(STDGEAR)))
    samp$DoorSpread <-  an(substr(STDGEAR,nchar(STDGEAR),nchar(STDGEAR)))
    samp$SweptAreaWSKM2 <-  samp$SweptAreaDSKM2 <-  samp$SweptAreaBWKM2 
    
    samp = samp[samp$Gear %in% paste0(STDGEAR,GEARSUBSCRIPTS),]  
  }  
 
  
  #  SMFS 0816 Derivation report Step 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  samp<- samp[samp$Year >= FIRSTYEAR,] 
  samp<- samp[samp$Year <= LASTYEAR,] 
  samp<- samp[samp$HaulDur >= MINTDUR,] 
  samp<- samp[samp$HaulDur <= MAXTDUR,] 
  #  SMFS 0816 Derivation report Step 1 END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  if(nrow(samp)==0){ print(paste("no data for",survey, COUNTRY)); next; }
  #add a directory for files out
  if(!file.exists(paste(OUTPATHstem,survey,sep=''))) dir.create(file.path(OUTPATHstem, survey))
  
  samp<-samp[samp$HaulVal=="V",] ###
  samp<-samp[samp$StdSpecRecCode==1,]
  samp$StNo[is.na(samp$StNo)] <- -9 #to match bio
  names(samp)[which(names(samp) == "HaulDur")] <- "HaulDur_min";
  names(samp)[which(names(samp) == "Depth_gam")] <- "Depth_m";
  names(samp)[which(names(samp) == "WingSpread")] <- "WingSpread_m";
  names(samp)[which(names(samp) == "DoorSpread")] <- "DoorSpread_m";
  names(samp)[which(names(samp) == "Netopening")] <- "NetOpen_m";
  names(samp)[which(names(samp) == "SweptAreaWSKM2")] <- "WingSwpArea_sqkm"
  names(samp)[which(names(samp) == "Distance")] <- "Distance_km"
  samp$NetOpen_m[samp$NetOpen_m== -9] <- mean(samp$NetOpen_m[samp$NetOpen_m!= -9])
  samp$WingSwpVol_CorF<- 1/samp$NetOpen_m
  samp$DoorSwptArea_CorF <- samp$WingSwpArea_sqkm / samp$SweptAreaDSKM2
  samp$DoorSwptVol_CorF <- samp$DoorSwptArea_CorF*samp$WingSwpVol_CorF
  names(samp)[which(names(samp) == "ShootLong")] <- "ShootLong_degdec";
  names(samp)[which(names(samp) == "ShootLat")] <- "ShootLat_degdec";
  names(samp)[which(names(samp) == "Year")] <- "YearShot"
  names(samp)[which(names(samp) == "Month")] <- "MonthShot";
  names(samp)[which(names(samp) == "StatRec")] <- "ICESStSq"
  names(samp)[which(names(samp) == "DepthStratum")] <- "L_REG"

  if( nrow(samp[is.na(samp$ShootLat_degdec),])>0) samp<-samp[!is.na(samp$ShootLat_degdec),]
  if( nrow(samp[is.na(samp$ShootLong_degdec),])>0) samp<-samp[!is.na(samp$ShootLong_degdec),]
  #correction needed!
  samp$Ship[samp$Ship=="7.40E+10"]<- "74E9" #correction! excel error pre-upload!
  
  #samp<-  merge(samp,flex,by="HaulID", all.x=F,all.y=F)
  #with(samp[samp$ShootLong_degdec>7,], xyplot(ShootLat_degdec~ShootLong_degdec | ac(YearShot) ))
  #with(samp, xyplot(ShootLat_degdec~ShootLong_degdec | ac(YearShot) ))
  
  
  #  SMFS 0816 Derivation report Step 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  samp$Month <- sprintf("%02s",samp$Month)
  samp$Day <- sprintf("%02s",samp$Day)
  samp$date =  paste0(samp$Year,"-",samp$Month,"-",samp$Day)
  samp$date = lubridate::as_date(samp$date)
  rects <- sf::st_read(paste0(RDIR,"rectanglesICESdl29oct2021/ICESRECT.shp"))
  rects <- sf:::as_Spatial(st_zm(rects))
  coordinates(samp) <- ~ ShootLong_degdec + ShootLat_degdec
  suppressWarnings(proj4string(samp) <- CRS("+init=epsg:4326")) #Warning message: In proj4string(obj) : CRS object has comment, which is lost in output
  suppressWarnings(proj4string(rects) <- CRS("+init=epsg:4326"))#Warning message: In proj4string(obj) : CRS object has comment, which is lost in output
  #plot(rects)
  #SSA_WRITE_NEW<-T
  if(SSA_WRITE_NEW){ # & !( substr(survey,nchar(survey)-4,nchar(survey))=="_hist" ) 
  SSAlist=list()
    for(re in 1:nrow(rects)){
      rect=rects[re,]
      samples_in=samp[rect,]
      if(length(samples_in)>0){

          # loop trough rectangle
          # subset spatial data by rect
          # get years with data in samp

          # 50 percent rule
          yearssampled=unique(samples_in$YearShot)
          nyearssampled=length(yearssampled)
          nyears = length(FIRSTYEAR:LASTYEAR)          
          over_50pct <- (nyearssampled/nyears)>=0.5

          # start & end 20% rule
          firstdate=lubridate::date_decimal(FIRSTYEAR)
          lastdate=lubridate::date_decimal(LASTYEAR)
          datessampled=unique(samples_in$date)
          mintimesampled = min(samples_in$date)
          maxtimesampled = max(samples_in$date)
          speriod = lubridate::date_decimal(FIRSTYEAR+(0.2*nyears))
          eperiod = lubridate::date_decimal(LASTYEAR-(0.2*nyears))
          inendperiod =  any( lastdate > datessampled & datessampled > eperiod)
          instartperiod =  any( firstdate < datessampled & datessampled < speriod)
          sampled_20pct = inendperiod & instartperiod

          # Accepted y/n?
          accepted_as_SSA = over_50pct & sampled_20pct
          if(accepted_as_SSA){SSAlist = c(SSAlist,rect$ICESNAME)}
        }
      }
    SSAdfs=data.frame('rectangle'=unlist(SSAlist))
    if(length(SSAdfs)>0){
      SSAdfs$survey=survey
    }
    #  SMFS 0816 Derivation report Step 2 END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  } else {
    # Filter down to SSA
    if( substr(survey,nchar(survey)-4,nchar(survey))=="_hist" ){ surveyread <- substr(survey,1, nchar(survey)-5) } else { surveyread <- survey}
    surveyread <- survey
    surveySSA = definedSSA[definedSSA$survey==surveyread,]
    # plot(surveySSA,col=3);map(add = T)
    #samp = samp[surveySSA,] #do later to keep track of lost hauls# now in "process_exchange_HL.r"
  }
  if(SSA_WRITE_ONLY & SSA_WRITE_NEW){ SSAdf=rbind(SSAdf,SSAdfs); next } 
  samp@data$ShootLong_degdec = samp$ShootLong_degdec
  samp@data$ShootLat_degdec = samp$ShootLat_degdec
  samp=samp@data
  
  
  
  #### biological data ####
  bio <- read.csv(BIOL_FILE,as.is = c(1,2,4,5,6,10,11) )  #avoid conversions of rect to e+07 etc #,as.is=1 ) #
  bio <- bio[bio$Quarter==QUARTER,]
  bio <- bio[!is.na(bio$ValidAphiaID),] #remove invalids
  if(GEAR == "BEAM" | FILTER_COUNTRY==T)  bio<-bio[bio$Country==COUNTRY,] #beam trawl surveys in one file
  if(nrow(bio[!is.na(bio$LenMeasType) & bio$LenMeasType!=1 & bio$LenMeasType!= -9,])>0){
    #write.csv(file=paste("out/",survey,"/",survey,"_bio_LenMeasType.csv",sep=""), x=bio[!is.na(bio$LenMeasType) & bio$LenMeasType!=1 & bio$LenMeasType!= -9,])
    write.csv(file=paste("out/",survey,"/bio_LenMeasType.csv",sep=""), x=bio[!is.na(bio$LenMeasType) & bio$LenMeasType!=1 & bio$LenMeasType!= -9,])
  }
  #in SPPLIST? 1 is Total Length	
  if(nrow(bio[!is.na(bio$LenMeasType) & bio$LenMeasType!=1 & bio$LenMeasType!= -9 & (bio$ScientificName_WoRMS %in% SPPLIST),])>0){
    write.csv(file=paste("out/",survey,"_biospplist_LenMeasType.csv",sep=""), x=bio[!is.na(bio$LenMeasType) & bio$LenMeasType!=1 & bio$LenMeasType!= -9 & (bio$ScientificName_WoRMS %in% SPPLIST),])
    write.csv(file=paste("out/",survey,"/biospplist_LenMeasType.csv",sep=""), x=bio[!is.na(bio$LenMeasType) & bio$LenMeasType!=1 & bio$LenMeasType!= -9 & (bio$ScientificName_WoRMS %in% SPPLIST),])
  }
  #next
  #correction needed to remove all leading zeroes
  if(length(bio$HaulNo[substr(bio$HaulNo,1,1)%in%"0"])>0) bio$HaulNo <- str_replace(bio$HaulNo, "^0+" ,"") 
  if(length(bio$StNo[substr(bio$StNo,1,1)%in%"0"])>0) bio$StNo <- str_replace(bio$StNo, "^0+" ,"") 
  if(length(samp$HaulNo[substr(samp$HaulNo,1,1)%in%"0"])>0) samp$HaulNo <- str_replace(samp$HaulNo, "^0+" ,"") 
  if(length(samp$StNo[substr(samp$StNo,1,1)%in%"0"])>0) samp$StNo <- str_replace(samp$StNo, "^0+" ,"") 

  ###### Corrections due to known DATRAS input errors #################################################################################################################################################################################
  #write.csv(bio[bio$ScientificName_WoRMS=="Gadus morhua",],"bio_cod.csv")
  #write.csv(bio[bio$SpeciesSciName=="Gadus morhua",],"C:/Users/cl06/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/bio_cod_SSASMP_kNN.csv")
  #write.csv((FISHDATA[FISHDATA$SpeciesSciName=="Gadus morhua",]),"biotest_cod.csv")
  
  ##### Dec 2021 corrections made here AWAITING CORRECTIONS ON DATRAS #####
  #check LenMeasType is total length
  
  #LngtCode 1 is cm; 0 is 0.5 cm reported in mm; '.' is mm
  # SPAIN - SP-NORTH
  # corrections to SpecVal 10 Nov 2021
  if(survey=="BBICnSpaOT4"){
     
    if(FIRSTYEAR< 2017 & nrow(bio[bio$Year%in% 2011:2016 & bio$SpecVal==0,])>0 ) bio[bio$Year%in% 2011:2016 & bio$SpecVal==0,]$SpecVal <- 1 #IEO confirm mistake
   
    # Capros aper boarfish cm -> mm
    if(LASTYEAR>2015 & nrow(bio[bio$Year%in% 2016:2020 & bio$ValidAphiaID==127419 & bio$LngtCode==1,])){ 
      bio[bio$Year%in% 2016:2020 & bio$ValidAphiaID==127419 & bio$LngtCode==1,]$LngtClass <- 
        bio[bio$Year%in% 2016:2020 & bio$ValidAphiaID==127419 & bio$LngtCode==1,]$LngtClass/10 #IEO confirm mistake
    }
    
    #Maurolicus muelleri 11  
    if(nrow(bio[bio$ValidAphiaID==127312 & bio$LngtCode==1 & bio$LngtClass>14,])){ 
      bio[bio$ValidAphiaID==127312 & bio$LngtCode==1 & bio$LngtClass>14,]$LngtClass <- 
        bio[bio$ValidAphiaID==127312 & bio$LngtCode==1 & bio$LngtClass>14,]$LngtClass/10
    }
    
    #Trachyrincus scabrus 57
    #if LenMeasType = 4 (preanal) so for Total Length = LngtClass*3.1 [Mindel et al 2016; SMFS 0818]
    if(nrow(bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>60,])){ 
      bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>60,]$LngtClass <- 
        bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>60,]$LngtClass/10 }
    
    #Engraulis encrasicolus and >400 and up to 1800 mm 2017 AB high - solution lengthclass/10 ?
    if(nrow(bio[bio$ValidAphiaID==126426 & bio$LngtCode==0 & bio$LngtClass>400,])){ 
      bio[bio$ValidAphiaID==126426 & bio$LngtCode==0 & bio$LngtClass>400,]$LngtClass <- 
        bio[bio$ValidAphiaID==126426 & bio$LngtCode==0 & bio$LngtClass>400,]$LngtClass/10 }
    
    #LngtCode 1# should it be . or 0 for  
    #Coelorinchus caelorhincus
    #if LenMeasType should be 4 (preanal)  Total Length = LngtClass*2.82 [Mindel et al 2016; SMFS 0818]
    if(nrow(bio[bio$ValidAphiaID==398381 & bio$LngtCode==1 & bio$LngtClass>40,])){ 
      bio[bio$ValidAphiaID==398381 & bio$LngtCode==1 & bio$LngtClass>40,]$LngtClass <- 
        bio[bio$ValidAphiaID==398381 & bio$LngtCode==1 & bio$LngtClass>40,]$LngtClass/10 }
    
    #Coelorinchus labiatus all in A hauls #IEO confirm mistake in LngtCode should be 0 (measured to nearest 0.5 cm; units in mm) 
    #AND LenMeasType should be 4 (preanal) so for Total Length = LngtClass*2.5 [Mindel et al 2016; SMFS 0818]
    if(nrow(bio[bio$ValidAphiaID==280299 & bio$LngtCode==1 & bio$LngtClass>40,])){ 
      bio[bio$ValidAphiaID==280299 & bio$LngtCode==1 & bio$LngtClass>40,]$LngtClass <- 
        bio[bio$ValidAphiaID==280299 & bio$LngtCode==1 & bio$LngtClass>40,]$LngtClass*2.5/10 }
    
    # Coryphaenoides rupestris 121
    #if LenMeasType = 4 (preanal) so for Total Length = LngtClass*4.7399 [Atkinson 1981; SMFS 0818]
    if(nrow(bio[bio$ValidAphiaID==158960 & bio$LngtCode==1 & bio$LngtClass>105,])){ 
      bio[bio$ValidAphiaID==158960 & bio$LngtCode==1 & bio$LngtClass>105,]$LngtClass <- 
        bio[bio$ValidAphiaID==158960 & bio$LngtCode==1 & bio$LngtClass>105,]$LngtClass/10 }
    
    #	Hymenocephalus italicus 22  #
    #if LenMeasType = 4 (preanal) so for Total Length = LngtClass*1 [Mindel et al 2016; SMFS 0818]
    if(nrow(bio[bio$ValidAphiaID==158961 & bio$LngtCode==1 & bio$LngtClass>22,])){ 
      bio[bio$ValidAphiaID==158961 & bio$LngtCode==1 & bio$LngtClass>22,]$LngtClass <- 
        bio[bio$ValidAphiaID==158961 & bio$LngtCode==1 & bio$LngtClass>22,]$LngtClass/10 }
     
  }

  if(survey=="WASpaOT3"){
    #Q) SP-PORC. issues with SpecVal=0 in 98% of HL records 2002:2005
    #A) Dec 2021 - solution - set first year set to 2006
    
    #Q) LngtCode 1# should it be . ? for   #Coelorinchus caelorhincus 51  #Hymenocephalus italicus 22  #Maurolicus muelleri 11  #Trachyrincus scabrus 57 #Trachyrincus scabrus 57
    #A) IEO: Yes, I've checked our data and yes we have some problems with data in Porcupine survey for both Trachyrincus and Coelorinchus 
    #during the first years were a few inconsistencies on how the fish were measured, sometimes to the close cm [1] and sometimes to the .5 cm [0] and it is not correctly coded, 
    #we will check those and correct them will not be able to have all the data in the same units, since it is not possible recover the .5 cm from the data taken as 1 cm.
    
    #Coelorinchus caelorhincus #bimodal distribution so setting threshold to 40 not 51 to be conservative
    if(nrow(bio[bio$ValidAphiaID==398381 & bio$LngtCode==1 & bio$LngtClass>=40,])){ 
      bio[bio$ValidAphiaID==398381 & bio$LngtCode==1 & bio$LngtClass>=40,]$LngtClass <- 
        bio[bio$ValidAphiaID==398381 & bio$LngtCode==1 & bio$LngtClass>=40,]$LngtClass/10 }
    
    #Trachyrincus scabrus 57
    if(nrow(bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>57,])){ 
      bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>57,]$LngtClass <- 
        bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>57,]$LngtClass/10 }
  

    # Capros aper boarfish cm->mm
    if(nrow(bio[bio$ValidAphiaID==127419 & bio$LngtCode==1 & bio$LngtClass>20,])){ 
        bio[bio$ValidAphiaID==127419 & bio$LngtCode==1 & bio$LngtClass>20,]$LngtClass <- 
          bio[bio$ValidAphiaID==127419 & bio$LngtCode==1 & bio$LngtClass>20,]$LngtClass/10 
    }

    # Coryphaenoides rupestris 121 #bimodal distribution so setting threshold to 105 not 121 to be conservative
    if(nrow(bio[bio$ValidAphiaID==158960 & bio$LngtCode==1 & bio$LngtClass>105,])){ 
      bio[bio$ValidAphiaID==158960 & bio$LngtCode==1 & bio$LngtClass>105,]$LngtClass <- 
        bio[bio$ValidAphiaID==158960 & bio$LngtCode==1 & bio$LngtClass>105,]$LngtClass/10 }
    
    #	Hymenocephalus  italicus 22  #
    if(nrow(bio[bio$ValidAphiaID==158961 & bio$LngtCode==1 & bio$LngtClass>22,])){ 
      bio[bio$ValidAphiaID==158961 & bio$LngtCode==1 & bio$LngtClass>22,]$LngtClass <- 
        bio[bio$ValidAphiaID==158961 & bio$LngtCode==1 & bio$LngtClass>22,]$LngtClass/10 }
    
    #Maurolicus muelleri 11  
    if(nrow(bio[bio$ValidAphiaID==127312 & bio$LngtCode==1 & bio$LngtClass>11,])){ 
      bio[bio$ValidAphiaID==127312 & bio$LngtCode==1 & bio$LngtClass>11,]$LngtClass <- 
        bio[bio$ValidAphiaID==127312 & bio$LngtCode==1 & bio$LngtClass>11,]$LngtClass/10 }

    #Malacocephalus laevis 79  
    if(nrow(bio[bio$ValidAphiaID==272392 & bio$LngtCode==1,])){ 
      bio[bio$ValidAphiaID==272392 & bio$LngtCode==1,]$LngtClass <- 
        bio[bio$ValidAphiaID==272392 & bio$LngtCode==1,]$LngtClass/10 }
    
    #Nezumia aequalis 49
    if(nrow(bio[bio$ValidAphiaID==126473 & bio$LngtCode==1,])){ 
      bio[bio$ValidAphiaID==126473 & bio$LngtCode==1,]$LngtClass <- 
        bio[bio$ValidAphiaID==126473 & bio$LngtCode==1,]$LngtClass/10 }
    
  } 
  
  ##SCOTLAND
  #length code for smooth sandeel in the SWC-IBTS for years 2006-2010 - in these years they are recorded as 'cm' rather than 'mm'. # MSS confirm mistake
  if( (COUNTRY=="GB-SCT" | COUNTRY=="Int") & (nrow(bio[bio$ValidAphiaID==126754 & bio$LngtCode==1 & bio$LngtClass>50,])>0) ){ #approx corrections
      bio[bio$ValidAphiaID==126754 & bio$LngtCode==1 & bio$LngtClass>50,]$LngtClass <- 
      bio[bio$ValidAphiaID==126754 & bio$LngtCode==1 & bio$LngtClass>50,]$LngtClass/10 #however in MSS script 170-172 cm should be 170/175/180 mm
  }
  
  #Nir-gfs Q1
  #subfactor for 1 record of Scyliorhinus stellaris is non zero - corrected in HL. AFBI confirmed 20/01/2022
  
  #FRANCE
  # error in a single record for Galeorhinus galeus in CGFS in 2018 (haul no 6)
  # the length class should be 1193 rather than 11930 [where LngtCode .] - corrected in HL and confirmaed by ifremer 20/01/2022
  
  #length code for TWO ALOSA SPECIES in the EVHOE DATA - in these years they are recorded as 'cm' rather than 'mm'.
  #126415 |   126413
  if( (survey=="BBICFraOT4" | survey=="CSFraOT4") & (nrow(bio[bio$ValidAphiaID==126413 | bio$ValidAphiaID==126415 & bio$LngtCode==1 & bio$LngtClass>100,])>0) ){ 
                             bio[bio$ValidAphiaID==126413 | bio$ValidAphiaID==126415 & bio$LngtCode==1 & bio$LngtClass>100,]$LngtClass <- 
                             bio[bio$ValidAphiaID==126413 | bio$ValidAphiaID==126415 & bio$LngtCode==1 & bio$LngtClass>100,]$LngtClass/10 #confirmed by ifremer 17/12/2021
  }
  #### select SpecVal ############################################################################################################################################################################################
  sort(unique(bio$SpecVal))# IBTS Q1 1983on has: 0 1 4 5 6 7  #data from 2016 download has only 0 1 4 5
  #6 No length measurements, only category catch weight	# could be ok for FC1 but not used by WKABSENS
  #5 Observed only, not measured, not counted, but only presence/absence is registered	# could be ok for FC1 but not used by WKABSENS
  # bioPA <- bio[(bio$SpecVal %in% c(5,6) ),]
  # sort(unique(PAworms<-bioPA$ScientificName_WoRMS)) # inverts mostly benthos/squid/shrimps
  # FC1Sp[FC1Sp %in% PAworms] # not losing any FC1 species in IBTS Q1 by losing 5 and 6
  # however if lose 7 would lose some Hippocampus hippocampus and Raja clavata
  # and if lose 4 would lose 18 species from IBTS Q1 data
  bio <- bio[(bio$SpecVal %in% c(1,4,7,10) ),]#remove invalids 0=Invalid information	 2=Partly valid information	   https://vocab.ices.dk/?ref=5
  #7 No length measurements, only total number and category catch weight	-> ok for FC1
  #4 No length measurements only total number		-> ok for FC1
  #10 No category catch weight, only total numbers and length composition	
  #1 Valid information for use in DATRAS data products	
  ## Note Meadhbh script no.7 includes 1 and V only - are there Vs in DATRAS? should not be as not a valid character

  #### make HaulID ############################################################################################################################################################################################
  #make ID and look for mismatches
  bio$HaulID	<- paste(paste(survey,QUARTER,sep=""), bio$Country, #bio$Gear, 
                      bio$Ship,  bio$StNo, bio$HaulNo, bio$Year,sep=":")
  
  samp$HaulID <- paste(paste(survey,QUARTER,sep=""),samp$Country,#samp$Gear,
                       samp$Ship, samp$StNo, samp$HaulNo, samp$YearShot,sep=":" )
  
    haulidhl <- sort(unique(samp$HaulID))
    haulidhh <- sort(unique(bio$HaulID))
    
    haulidhhm <- subset(samp, !samp$HaulID %in% bio$HaulID)
    haulidhhm <- sort( unique(haulidhhm$HaulID) ); length(haulidhhm) 
    100*length(haulidhhm)/nrow(samp) # 2samp hauls not match in bio
    
    haulidhlm <- subset(bio, !bio$HaulID %in% samp$HaulID)# bio hauls not match in samp
    haulidhlm <- sort(unique(haulidhlm$HaulID) ); length(haulidhlm) 
    100*length(haulidhlm)/nrow(bio) #samp hauls not match in bio
  
  print("process_exchange_HL") #QC of data and eyeball of sampling and bio data
   source('R/Lynam_IND_script_process_exchange_HL2021.R') #inc CATCHABILITY_COR_WALKER
  #now have merged samp and bio to create 'dhspp'
  
  #### link subdivisional shapefiles to data ####
  # and read attributes of shapefiles create table ATTRIB with names L_REG & KM2_LAM 
  ##### strata ##### 
  print("now add strata")
  source(paste(MAINDIR,"R/Lynam_OSPARsubdiv_Jan2022.r",sep=""))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Calc ALL indicators ####  
  if(EHDS_PP) BY_SREG<-F
  # if(BY_SREG){ # dependent on survey see Lynam_OSPARsubdiv.r above# if(survey=="GNSFraOT4"){ dhspp$S_REG <- dhspp$S_REG # minigrid of sqs
  #   dhspp <- dhspp[!is.na(dhspp$L_REG),] #rm outside area## lots of KS
  # }# else { dhspp$L_REG <- NA; ATTRIB$L_REG <- NA } 
 # if(BY_LREG){   dhspp$L_REG <- dhspp$S_REG; ATTRIB$L_REG <- ATTRIB$S_REG;
 # #dhspp$L_REG <-substr(dhspp$L_REG,1,2); # dependent on survey
 # #ATTRIB$L_REG <-substr(ATTRIB$L_REG,1,2); # dependent on survey
 #   dhspp$S_L_REG <- paste(dhspp$L_REG, dhspp$L_REG,sep="_"); ATTRIB$S_L_REG <- paste(ATTRIB$L_REG, ATTRIB$L_REG,sep="_")
 # } else { dhspp$L_REG <- dhspp$S_L_REG <- NA;  ATTRIB$L_REG <- ATTRIB$S_L_REG <- NA; }
 # #could be hauls(L_REG) within subdiv or only 'NA_subdiv'

  dhspp$S_L_REG <- paste(dhspp$S_REG, dhspp$L_REG,sep="_"); # STRAT_DIV renamed to S_L_REG
  ATTRIB$S_L_REG <- paste(ATTRIB$S_REG, ATTRIB$L_REG,sep="_") 

  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #save(dhspp,file=paste("DHSPP_DAT_",Sys.Date(),".RData",sep=""))
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #dhspp<-trawldf
  YRS <- sort(unique(dhspp$Year))
  
  #all indicators calc here and biomass by rectangle
  print("Calc Biomass and Indicators")
  time_of_run = format(Sys.time(), "%d%b%Y")
  FILENAM<-paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),sep="")
  if(!IEO_FW4){
    
    if(BYGROUP){
      IND_OUT_BYGROUP<-list()
      for(SPGROUP in levels(dhspp$Group)){ #SPGROUP<-"Elasmobranchii"
        print(paste(survey,SPGROUP,sep=""))
        if(SPGROUP=="Elasmobranchii"){ LFI_THRESHOLD_APPLY <- 60 } 
        if(SPGROUP=="Pleuronectiformes"){ LFI_THRESHOLD_APPLY <- 20 } 
        if(SPGROUP=="Scorpaeniformes"){ LFI_THRESHOLD_APPLY <- 20 } 
        if(SPGROUP=="Other"){ LFI_THRESHOLD_APPLY <- 30 } 
        if(SPGROUP=="Gadiformes"){ LFI_THRESHOLD_APPLY <- LFI_THRESHOLD} 
        if(nrow(dhspp[dhspp$Group==SPGROUP,])>3){
          try(
          IND_OUT_BYGROUP[[SPGROUP]] <- INDfn( DATA=dhspp, WRITE=T, SPECIES=SPECIES, GROUP=SPGROUP,
                        LFI_THRESHOLD=LFI_THRESHOLD_APPLY, 
                        BY_SREG=BY_SREG, BY_LREG=BY_LREG, FILENAM=FILENAM,
                        LFI=LFI, MEANTL=F, MaxL=T, Loo=T, Lm=T, MeanL=F, TyL_GeoM=T,
                        QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS)
          ,silent=TRUE)
        }
      }
    }#BYGROUP
    
    if(BYICESGROUP){
      IND_OUT_BYICESGROUP<-list()
      trait_MAXL <- merge(x=trait_MAXL, y=SciName2GroupEff[,c(1,9)], by.x="SpeciesSciName", by.y="sciName",all.x=T)#INDfn need to re-merge this
      trait_MAXL$Group <- ac(trait_MAXL$habitat.guild )
      for(SPGROUP in u(dhspp$habitat.guild)){
        if(nrow(dhspp[dhspp$Group==SPGROUP,])>3){
        print(paste(survey,SPGROUP,sep=""))  
          try(
            IND_OUT_BYICESGROUP[[SPGROUP]] <- INDfn( DATA=dhspp, WRITE=WRITE, SPECIES=SPECIES, GROUP=SPGROUP,
                        LFI=LFI, LFI_THRESHOLD=NULL, FILENAM=FILENAM,
                        BY_SREG=BY_SREG, BY_LREG=BY_LREG, 
                        MEANTL=F, MaxL=T, Loo=T, Lm=T, MeanL=F, TyL_GeoM=T,
                        QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS)
            ,silent=TRUE)
        }
      }
    }#BYICESGROUP
    
    #trophic GUILD - this subsets dhspp so copy as _raw and the replace at end
    if(BYGUILD){
      
      print("make GUILD GROUPING")
      IND_OUT_BYGUILD<-list()
      source(paste(PROC_SCRIPT,"Lynam_IND_script_MAKEGUILDIND.R",sep="")) 
      print("indicator for Trophic Guilds")
      for(SPGROUP in sort(u(dhspp$fguild)) ){
        if(SPGROUP=="5") next
        print(paste(survey," Guild ",SPGROUP,sep=""))  
        try(
          IND_OUT_BYGUILD[[SPGROUP]] <- INDfn( DATA=dhspp, WRITE=WRITE, SPECIES=SPECIES, GROUP=SPGROUP,
                                               LFI=F, LFI_THRESHOLD=NULL, FILENAM=FILENAM,
                                               BY_SREG=BY_SREG, BY_LREG=BY_LREG, 
                                               MEANTL=F, MaxL=T, Loo=F, Lm=F, MeanL=F, TyL_GeoM=T,BYGUILD=BYGUILD,
                                               QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS)
          ,silent=TRUE)
      }
      
    }#BYGUILD
  }
  ##### NO TAXA GROUPING ##### 
  #do this last to make sure have all species in final biomass plots
  print("NO TAXA GROUPING")
  FISHDATA <- dhspp[dhspp$SpeciesSciName%in%SPPLIST,] #if IEO_FW4 FALSE then this is already done in processing script
  #if IEO_FW4 true make sure not including inverts in LFI etc
  #INDfn creates haul_by_spp and hauls.csv i.e. the sampling and biological data used for indicator assessments
  ghj
  if(LFI_NULL) LFI_THRESHOLD<-NULL
   try(
    IND_OUT <- INDfn( DATA=FISHDATA, WRITE=WRITE, BOOTSTRAP=BOOTSTRAP, LFI=LFI, LFI_THRESHOLD=LFI_THRESHOLD,
                      FILENAM=FILENAM,BY_SREG=BY_SREG, BY_LREG=BY_LREG, 
                    MEANTL=MEANTL, MaxL=MaxL,Loo=Loo, Lm=Lm, MeanL=MeanL, TyL_GeoM=TyL_GeoM, SPECIES=SPECIES, 
                    GROUP=NULL, TyL_SPECIES=TyL_SPECIES, BYGUILD=F, QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS,ATTRIB=ATTRIB)
   ,silent=F)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### BOOTSTRAP ####
  #lists to collate multiple bootstrapped indicators
  BOOTDATA2PLOT <- NULL
  LFI_regional<- TL_regional <- TyLrect_regional <- TyL_regional <- TyL_reg_var <- Len_regional <- MaxL_regional <- Loo_regional <- Lm_regional <- NULL
  if(BOOTSTRAP){ print("Bootstrap dataset"); source(paste(PROC_SCRIPT,"Lynam_IND_script_BOOTSTRAP.r",sep="")) }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  #investigate an indicator with GAM    library(mgcv)
  #and plot some (with bootstrap?)
  #IND_OUT <- IND_OUT_BYGROUP[["Elasmobranchii"]]; BOOT_OUT <- BOOT_OUT_BYGROUP[["Elasmobranchii"]]
  if(FINALPLOTS){ 
    print("Final plots"); 
    SPECIES= combs$Species #EVEN IF HAVE IEO_SPECIES == T, ie FULL LIST, DO ONLY THOSE SELECTED FOR INDS HERE
    try( source(paste(MAINDIR,"R/Lynam_IND_script_FINALPLOTS_Dec2021.R",sep="")) ,silent=F); 
    dev.off()
  }
  if(SSA_WRITE_NEW) { SSAdf=rbind(SSAdf,SSAdfs) }
  #dhspp<-dhspp[,-which(names(dhspp) %in% c("sciName","SubFactor","DataType","DurRaise","LogLngtClass","LogLngtBio","LFI_Fung_list","LFI_OSPAR_list","S_REG","L_REG","S_L_REG")),]
  #SciName		Number
  dhspp$Survey_Acronym <- survey
  names(dhspp)[which(names(dhspp)=="sciName")] <- "SciName"
  dhspp$Gear <- GEAR 
  dhspp$GearType <- STDGEAR 
  names(dhspp)[which(names(dhspp)=="MonthShot")] <- "Month"
  names(dhspp)[which(names(dhspp)=="WingSwpArea_sqkm")] <- "SweptArea_KM2"
  dhspp <- dhspp[, which(names(dhspp) %in% c("HaulID","Survey_Acronym","ICESStSq",
                                             "Year","HaulDur_min","SweptArea_KM2","L_REG","S_REG",
                                             "SciName","SensFC1","DEMPEL","FishLength_cm",
                                             "DensAbund_N_perhr","DensBiom_kg_perhr",
                                             "DensBiom_kg_Sqkm","DensAbund_N_Sqkm",
                                             "WingSwpArea_sqkm",
                                             "NetOpen_m","Gear","Ship",
                                             "Month","Day","TimeShot",
                                             "ShootLong_degdec","ShootLat_degdec"
  ) )]
  if(WRITE_LDs | IEO_FW4) write.csv(dhspp,paste0(FILENAM,"_haul_by_spp.txt"),row.names = F)
  
  print(paste("Finished",survey, "survey",sep=" "))
  rm(dhspp)
} #next survey


if(SSA_WRITE_NEW){
  SSAdf=SSAdf[-1,] # first row is null from setup
  write.csv(SSAdf,paste0(RDIR,"/SSA.csv"),row.names = F)
  rects$rectangle = rects$ICESNAME
  spatialSSA=sp::merge(rects,SSAdf,on='rectangle',all.x=F,all.y=T, duplicateGeoms = TRUE)
  setwd(paste0(RDIR,"/rectanglesICESdl29oct2021/"))
  writeOGR(spatialSSA, "shp", "SSAspatial" , driver = "ESRI Shapefile", overwrite_layer = T) 
  # SSA saved to csv and shapefile now.
  
  
  # Split rects into 4 using a grid defined in QGIS.
  # This grid was made using 'create_grid' tool with a spacing of 0.5 degrees longitude and 0.25 degrees latitude in wgs84. Any other settings were default.
  quarter_rects_grid = sf::st_read(paste0(RDIR,"rectanglesICESdl29oct2021/qgis_create_grid_spacing_05h_025v_wgs84.shp")) # read.csv(paste0(RDIR,"/defined_SSA.csv"))
  quarter_rects_grid <- sf:::as_Spatial(st_zm(quarter_rects_grid))

  # For some reason I can't find an r function which intersects polygons by a grid to get smaller polygons, hence this long approach which involves merging again afterwards:
  lpi = gIntersection(spatialSSA,quarter_rects_grid)
  blpi <- gBuffer(lpi, width = 0.00000001, byid = T)  # create a very thin polygon 
  dpi <- gDifference(spatialSSA, blpi, byid = T)
  
  # bring desired attributes over  
  dpi$rectangle <- spatialSSA$rectangle
  dpi$Ecoregion<- spatialSSA$Ecoregion
  dpi$survey <- spatialSSA$survey
  dpi$dummy <- NULL
  
  # cast from multipolygons to single with disaggregate
  dpi = sp::disaggregate(dpi)
  dpi$quadid <- 1:length(dpi)
  # write
  dpi = as(dpi, "SpatialPolygonsDataFrame" )
  writeOGR(dpi, "shp", "SSAspatial_quarters" , driver = "ESRI Shapefile", overwrite_layer = T) 
}

WRITE_to_DB = F
source(paste(MAINDIR,"R/combine_all_hauls.r",sep=""))  
source(paste(MAINDIR,"R/database_upload.r",sep=""))  

print("script complete")