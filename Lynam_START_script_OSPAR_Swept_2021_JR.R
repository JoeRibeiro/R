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
MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/"
RDIR<- paste0(MAINDIR,"R/")
definedSSA = sf::st_read(paste0(RDIR,"rectanglesICESdl29oct2021/shp_dir_for_original/SSAspatial.shp")) # read.csv(paste0(RDIR,"/defined_SSA.csv"))
definedSSA <- sf:::as_Spatial(st_zm(definedSSA))


#Trophic Level data also requires update for MTL analyses

#location of the subscripts in function area
 PROC_SCRIPT<- "//lowfilecds/function/Eco Indicators/DATRASoutput/" #for HH and HL processing scripts incl strata by survey
 SHAPEPATH<-paste("Strata/",sep="") #used in Lynam_OSPARsubdiv.r
 #SUBSCRIPTS_TRAITS<-paste(PROC_SCRIPT,"MarScot/INDscripts2021/",sep="")#max length
 
 setwd(MAINDIR)
 #read subscripts and traits
 # PRODDAT<-read.csv(paste(SUBSCRIPTS_TRAITS,"SpeciesAtLength_Productivity.csv",sep=""))
 SPPFILE<-"R/Species_List_Final_06Dec2021.csv" # prepared by Meadhbh Moriarty and Simon Greenstreet for IA02017 shared in ICES WGBIODIV. FC1 identified added by CL based on WKABSENS 2021
 
 LW<-read.csv(SPPFILE)
 names(LW)[2]<-"ScientificName_WoRMS"
 LW$SpeciesSciName <- LW$ScientificName_WoRMS
 LW$a <- LW$LWRa
 LW$b <- LW$LWRb
 LW$"Max.L..cm." <- LW$MaxL
 
 FC1Sp<-LW[LW$SensFC1=="FC1",2]
 SPPLIST <- LW$ScientificName_WoRMS
 trait_MAXL<-LW[,which(names(LW) %in% c("ScientificName_WoRMS","SpeciesSciName","Loo","Max.L..cm.","MaxL","Lm","Order","Group"))]
 
#where save output?
OUTPATHstem<-paste(MAINDIR,"out_for_chibuzor/",sep="")


## choices for analyses upfront

#do you want to write outputs along the way and save workspace
WRITE <- T #save csvs as we go?
  WRITE_LDs <- T #write Length distributions by species/year/subdiv? 
BOOTSTRAP <- F # invoke slow slow code? if F next 3 lines redundant
  NBOOT<- 9; 
  B<-0 # restart counter as only output LD once before boostrap starts i.e. when B=0 
  WRITE_BOOT <- F # every bootstrap dataset and indicator output
  SAVE <- F # save workspace (after bootstrap)
SSA_WRITE_NEW <- F
FILTER_COUNTRY <- F
WRITE_to_DB <- F

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
  LFI <- F
  FINALPLOTS<- F#create indicator plots with smooths
  IEO_FW4<-T #keep inverts
  
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
  USE_GUILD_COVARIATE_SITES<- F #'to match Murray et al'
  FILL_GUILD<-F #if false have a no-guild group
  
# Spatial analyses
AREASCALE <- T # #raise LD data and CPUE and Catch_Biomass by area e.g rectangle KM2_LAM of lowest resolution of sampling strategy/subdiv in Lynam_INDfn_Nov2021.r then x'scale' for missing subdiv area
EHDS_PP <-F #ecohydrodynamic zones
STRATA <- F  #for output
QUAD_NS <- F #only north sea
QUAD_SMOOTH_NS<- F # repeats per guild at present must change - SLOW!
QUADS<-NULL #created by Lynam_OSPARsubdiv_Feb2019.r

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
## Lynam_OSPARsubdiv_Nov2021.r
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
for(combrow in 1:nrow(survey_Q_C_S_combinations) ){
  #for(combrow in c(13,26,1,14,18,19,20,21) ){
  #for(combrow in c(26,21) ){
  #### combrow<-9 #### combrow<-nrow(survey_Q_C_S_combinations) # combrow<-13
  combs=survey_Q_C_S_combinations[combrow,]
  #combs
  QUARTER=combs$Quarter
  COUNTRY=combs$Country
  SEA=combs$Sea
  survey=combs$Surveynam1
  survey_alt_name=combs$Surveynam2
  LFI_THRESHOLD = combs$LFI_threshold
  GEAR= combs$Gear
  SAMP_STRAT= combs$SAMP_STRAT # average hauls by ICESStSq rect in north sea #if set to FALSE need to update Attibutes table as area only given by rect
  BYSUBDIV = combs$BYSUBDIV    # average indicator by LFI-subdivision
  
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
  #sampling data 
  samp <- read.table(SAMP_FILE ,as.is = c(1,2,4,5,6,10,11,23),header = TRUE,sep=",") 
  samp$StatRec <- ices.rect2(ices.rect(samp$StatRec)$lon,ices.rect(samp$StatRec)$lat) #sometime StatRec has been excel-ed into a number
  samp <- samp[samp$Quarter==QUARTER,]
  # where two different surveys use the same input data files the correct preprocessing needs to be done to define SSA
  # need to split SEA DATA FOR BTS Q3 ENG'
  #otter surveys
  if(survey %in% "CSBBFraOT4") samp<-samp[samp$ShootLat<= 48,]
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

 
  
  #  SMFS 0816 Derivation report Step 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  samp = samp[samp$Gear %in% paste0(STDGEAR,GEARSUBSCRIPTS),]  
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
  names(samp)[which(names(samp) == "DepthStratum")] <- "SurvStratum"

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
  
  if(SSA_WRITE_NEW & !( substr(survey,nchar(survey)-4,nchar(survey))=="_hist" ) ){
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
    surveySSA = definedSSA[definedSSA$survey==surveyread,]
    #samp = samp[surveySSA,] #do later to keep track of lost hauls# now in "process_exchange_HL.r"
  }
  samp@data$ShootLong_degdec = samp$ShootLong_degdec
  samp@data$ShootLat_degdec = samp$ShootLat_degdec
  samp=samp@data
  
  
  
  # biological data
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

  ################################################################################################################################################################################################
  ##### Dec 2021 corrections made here AWAITING CORRECTIONS ON DATRAS #####
  #check LenMeasType is total length
  
  #LngtCode 1 is cm; 0 is 0.5 cm reported in mm; '.' is mm
  # SPAIN - SP-NORTH
  # corrections to SpecVal 10 Nov 2021
  if(survey=="BBICnSpaOT4"){
    if(FIRSTYEAR< 2017 & nrow(bio[bio$Year%in% 2011:2016 & bio$SpecVal==0,])>0 ) bio[bio$Year%in% 2011:2016 & bio$SpecVal==0,]$SpecVal <- 1 #IEO confirm mistake
    
    # Capros aper boarfish cm->mm
    if(LASTYEAR>2015 & nrow(bio[bio$Year%in% 2016:2020 & bio$ValidAphiaID==127419 & bio$LngtCode==1,])){ 
      bio[bio$Year%in% 2016:2020 & bio$ValidAphiaID==127419 & bio$LngtCode==1,]$LngtClass <- 
        bio[bio$Year%in% 2016:2020 & bio$ValidAphiaID==127419 & bio$LngtCode==1,]$LngtClass/10 #IEO confirm mistake
    }
    
    #Engraulis encrasicolus and >400 and up to 1800 mm 2017 AB high - solution lengthclass/10 ?
    if(nrow(bio[bio$ValidAphiaID==126426 & bio$LngtCode==0 & bio$LngtClass>400,])){ 
      bio[bio$ValidAphiaID==126426 & bio$LngtCode==0 & bio$LngtClass>400,]$LngtClass <- 
        bio[bio$ValidAphiaID==126426 & bio$LngtCode==0 & bio$LngtClass>400,]$LngtClass/10 }
    
    #Maurolicus muelleri 11  
    if(nrow(bio[bio$ValidAphiaID==127312 & bio$LngtCode==1 & bio$LngtClass>14,])){ 
      bio[bio$ValidAphiaID==127312 & bio$LngtCode==1 & bio$LngtClass>14,]$LngtClass <- 
        bio[bio$ValidAphiaID==127312 & bio$LngtCode==1 & bio$LngtClass>14,]$LngtClass/10 }
    
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
      
    #Trachyrincus scabrus 57
    #if LenMeasType = 4 (preanal) so for Total Length = LngtClass*3.1 [Mindel et al 2016; SMFS 0818]
    if(nrow(bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>60,])){ 
      bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>60,]$LngtClass <- 
        bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>60,]$LngtClass/10 }
    
    ##note outstanding issue
    #incorrect matching of HL to A hauls since Day and Time not in HL data!!!!!!!!!!!!!!!!!!!!!!!!
    
  }
  if(survey=="WASpaOT3"){#rerun
    #Q) SP-PORC. issues with SpecVal=0 in 98% of HL records 2002:2005
    #A) Dec 2021 - solution - set first year set to 2006
    
    #Q) LngtCode 1# should it be . ? for   #Coelorinchus caelorhincus 51  #Hymenocephalus italicus 22  #Maurolicus muelleri 11  #Trachyrincus scabrus 57 #Trachyrincus scabrus 57
    #A) IEO: Yes, I've checked our data and yes we have some problems with data in Porcupine survey for both Trachyrincus and Coelorinchus 
    #during the first years were a few inconsistencies on how the fish were measured, sometimes to the close cm [1] and sometimes to the .5 cm [0] and it is not correctly coded, 
    #we will check those and correct them will not be able to have all the data in the same units, since it is not possible recover the .5 cm from the data taken as 1 cm.
    
    
    # Capros aper boarfish cm->mm
    if(nrow(bio[bio$ValidAphiaID==127419 & bio$LngtCode==1 & bio$LngtClass>20,])){ 
        bio[bio$ValidAphiaID==127419 & bio$LngtCode==1 & bio$LngtClass>20,]$LngtClass <- 
          bio[bio$ValidAphiaID==127419 & bio$LngtCode==1 & bio$LngtClass>20,]$LngtClass/10 
    }
    
    #Coelorinchus caelorhincus #bimodal distribution so setting threshold to 40 not 51 to be conservative
    if(nrow(bio[bio$ValidAphiaID==398381 & bio$LngtCode==1 & bio$LngtClass>=40,])){ 
      bio[bio$ValidAphiaID==398381 & bio$LngtCode==1 & bio$LngtClass>=40,]$LngtClass <- 
        bio[bio$ValidAphiaID==398381 & bio$LngtCode==1 & bio$LngtClass>=40,]$LngtClass/10 }

    
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
    
    #Trachyrincus scabrus 57
    if(nrow(bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>57,])){ 
      bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>57,]$LngtClass <- 
        bio[bio$ValidAphiaID==126482 & bio$LngtCode==1 & bio$LngtClass>57,]$LngtClass/10 }
    
  } 
  
  ##SCOTLAND
  #length code for smooth sandeel in the SWC-IBTS for years 2006-2010 - in these years they are recorded as 'cm' rather than 'mm'. #Q-emailed to MSS are checking.
  if( (survey=="CSScoOT4" | survey=="CSScoOT1" | survey=="CSScoOT4_hist" | survey=="CSScoOT1_hist" ) & (nrow(bio[bio$ValidAphiaID==126754 & bio$LngtCode==1 & bio$LngtClass>50,])>0) ){ 
      bio[bio$ValidAphiaID==126754 & bio$LngtCode==1 & bio$LngtClass>50,]$LngtClass <- 
      bio[bio$ValidAphiaID==126754 & bio$LngtCode==1 & bio$LngtClass>50,]$LngtClass/10 #MSS?
  }
  
  #FRANCE
  #length code for TWO ALOSA SPECIES in the EVHOE DATA - in these years they are recorded as 'cm' rather than 'mm'. MSS are checking.
  #126415 |   126413
  if( (survey=="CSBBFraOT4" | survey=="CSFraOT4") & (nrow(bio[bio$ValidAphiaID==126413 | bio$ValidAphiaID==126415 & bio$LngtCode==1 & bio$LngtClass>100,])>0) ){ 
                             bio[bio$ValidAphiaID==126413 | bio$ValidAphiaID==126415 & bio$LngtCode==1 & bio$LngtClass>100,]$LngtClass <- 
                             bio[bio$ValidAphiaID==126413 | bio$ValidAphiaID==126415 & bio$LngtCode==1 & bio$LngtClass>100,]$LngtClass/10 #confirmed by ifremer 17/12/2021
  }
  ################################################################################################################################################################################################
  
  #make ID and look for mismatches
  bio <- bio[(bio$SpecVal %in% c(1,4,7,10) ),]#remove invalids 0=Invalid information	 2=Partly valid information	   https://vocab.ices.dk/?ref=5
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
  # and read attributes of shapefiles create table ATTRIB with names SurvStratum & KM2_LAM 
  ##### strata ##### 
  print("now add strata")
  #source("//lowfilecds/Function/Eco Indicators/DATRASoutput/MarScot/INDscriptsForV3/Lynam_OSPARsubdiv_Feb2019.r")
  #source(paste(MAINDIR,"R/Lynam_OSPARsubdiv_Oct2021.r",sep=""))
  source(paste(MAINDIR,"R/Lynam_OSPARsubdiv_Nov2021.r",sep=""))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### Calc ALL indicators ####  
  if(EHDS_PP) SAMP_STRAT<-F
  if(SAMP_STRAT){ # dependent on survey see Lynam_OSPARsubdiv.r above# if(survey=="GNSFraOT4"){ dhspp$sampstrat <- dhspp$SurvStratum # minigrid of sqs
    dhspp <- dhspp[!is.na(dhspp$sampstrat),] #rm outside area## lots of KS
  } else { dhspp$sampstrat <- NA } 
  if(BYSUBDIV){   dhspp$subdiv <- dhspp$SurvStratum; #dhspp$subdiv <-substr(dhspp$subdiv,1,2); # dependent on survey
    dhspp$STRAT_DIV <- paste(dhspp$sampstrat, dhspp$subdiv,sep="_") 
  } else { dhspp$subdiv <- dhspp$STRAT_DIV <- NA; } 
  #could be hauls(sampstrat) within subdiv or only 'NA_subdiv'
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #save(dhspp,file=paste("DHSPP_DAT_",Sys.Date(),".RData",sep=""))
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #dhspp<-trawldf
  YRS <- sort(unique(dhspp$Year))
  
  #all indicators calc here and biomass by rectangle
  print("Calc Biomass and Indicators")
  time_of_run = format(Sys.time(), "%d%b%Y")
  FILENAM<-paste(OUTPATH,survey,"_",time_of_run,sep="")
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
                      SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, FILENAM=FILENAM,
                      LFI=LFI, MEANTL=F, MaxL=T, Loo=T, Lm=T, MeanL=F, TyL_GeoM=T,
                      QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS)
        ,silent=TRUE)
      }
    }
  }
  
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
                      SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, 
                      MEANTL=F, MaxL=T, Loo=T, Lm=T, MeanL=F, TyL_GeoM=T,
                      QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS)
          ,silent=TRUE)
      }
    }
  }
  
  #trophic GUILD - this subsets dhspp so copy as _raw and the replace at end
  if(BYGUILD){#dhspp<-dhspp_raw
    dhspp_raw <- dhspp
    if(survey == "GNSIntOT1" & USE_GUILD_COVARIATE_SITES | survey == "GNSIntOT1_channel" & USE_GUILD_COVARIATE_SITES){
      print("reading GUILD_COVARIATE_SITES")
      load(file=paste(PROC_SCRIPT,'Processed_data_for_models_13.11.18.RData',sep="")) #fulldat
      sites = unique(fulldat$sampstrat)
      dhspp<-dhspp[dhspp$sampstrat %in% sites,]
    }
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
                                             SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, 
                                             MEANTL=F, MaxL=T, Loo=F, Lm=F, MeanL=F, TyL_GeoM=T,BYGUILD=BYGUILD,
                                             QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS)
        ,silent=TRUE)
    }
    dhspp <- dhspp_raw
    }
  }
  ##### NO TAXA GROUPING ##### 
  #do this last to make sure have all species in final biomass plots
  print("NO TAXA GROUPING")
  FISHDATA <- dhspp[dhspp$SpeciesSciName%in%SPPLIST,] #if IEO_FW4 FALSE then this is already done in processing script
  #if IEO_FW4 true make sure not including inverts in LFI etc
  if(LFI_NULL) LFI_THRESHOLD<-NULL
   try(
    IND_OUT <- INDfn( DATA=FISHDATA, WRITE=WRITE, BOOTSTRAP=BOOTSTRAP, LFI=LFI, LFI_THRESHOLD=LFI_THRESHOLD,
                      FILENAM=FILENAM,SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, 
                    MEANTL=MEANTL, MaxL=MaxL,Loo=Loo, Lm=Lm, MeanL=MeanL, TyL_GeoM=TyL_GeoM, SPECIES=SPECIES, 
                    GROUP=NULL, TyL_SPECIES=TyL_SPECIES, BYGUILD=F, QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS)
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
  if(FINALPLOTS){ print("Final plots"); 
    try( source(paste(MAINDIR,"R/Lynam_IND_script_FINALPLOTS_Dec2021.R",sep="")) ,silent=F); 
    dev.off()
  }
  if(SSA_WRITE_NEW) { SSAdf=rbind(SSAdf,SSAdfs) }
  dhspp<-dhspp[,-which(names(dhspp) %in% c("sciName","SubFactor","DataType","DurRaise",
                                         "LogLngtClass","LogLngtBio","LFI_Fung_list","LFI_OSPAR_list",
                                         "sampstrat","subdiv","STRAT_DIV")),]
  if(WRITE_LDs | IEO_FW4) write.csv(dhspp,paste0(FILENAM,"_haul_by_spp_",survey,".csv"),row.names = F)
  
  print(paste("Finished",survey, "survey",sep=" "))
  rm(dhspp)
} #next survey
#
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



if(WRITE_LDs | IEO_FW4){
  # append all the data that are required for (a) Chibuzor's modelling work an (b) an output to go into the BX5 app
  fileslisted=list.files(OUTPATHstem,'haul_by_spp',full.names=T,recursive = T);
  # In case of reruns done on previous date
  fileslisted = fileslisted[grepl(time_of_run, fileslisted)]
  file1=read.csv(fileslisted[1]) 
  # Also get survey name from filename, crude approach
  file1$survey_name = strsplit(fileslisted[1],'/')[[1]][length(strsplit(fileslisted[1],'/')[[1]])-1]
  for(file in 2:length(fileslisted)){file2=read.csv(fileslisted[file]); 
    file2$survey_name = strsplit(fileslisted[file],'/')[[1]][length(strsplit(fileslisted[file],'/')[[1]])-1]
    # SurvStratum is required but not present for some surveys, these must be GNS surveys
    if(!"SurvStratum" %in% colnames(file2)){file2$SurvStratum = file2$ICESStSq}
    #file2_names = colnames(file2); file1_names = colnames(file1); common_names = intersect(file2_names, file1_names); file1 = rbind(file2[common_names], file1[common_names])
    file1 = rbind(file2, file1)
    }
  
  # Add species ID from lookup table originating from database table bx005.public.species
  specieslookup=read.csv("C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/public_species_table.csv")
  specieslookup$SpeciesSciName = specieslookup$latin_name
  specieslookup$latin_name <- NULL
  specieslookup = specieslookup[,c("SpeciesSciName","species_id")]
  file1=merge(file1,specieslookup,all.x=T)
  
  # File is too big - aggregate and remove some columns. And rename some to be more consistent with the expected file
  desiredcols=colnames(cod_final)
  desiredcols[!desiredcols %in% colnames(file1)]
  file1$Survey_Acronym = file1$survey_name
  file1$survey_name <- NULL
  file1$GearType = file1$Gear
  file1$Gear <- NULL
  file1$YearShot = file1$Year
  file1$Year <- NULL
  file1$DayShot = file1$Day
  file1$Day <- NULL
  file1=file1[,!colnames(file1) %in% c("KM2_LAM","checkLen","Order","Group","Loo","Lm","MaxL","Max.L..cm.")]
  file1$numhauls=1
  file1=aggregate(cbind(DensBiom_kg_Sqkm,DensAbund_N_Sqkm,DensAbund_N_perhr,DensBiom_kg_perhr,numhauls) ~ SpeciesSciName + HaulID + SensFC1 + DEMPEL + HaulDur_min + ICESStSq + WingSwpArea_sqkm + WingSwpVol_CorF + NetOpen_m + Ship + MonthShot + TimeShot + ShootLong_degdec + ShootLat_degdec + habitat.guild + ScientificName_WoRMS + SurvStratum + species_id + Survey_Acronym + GearType + YearShot + DayShot , data = file1, FUN = sum, na.rm = TRUE)
  file1$sumDensBiom_kg_Sqkm = file1$DensBiom_kg_Sqkm
  file1$sumDensBiom_kg_perhr = file1$DensBiom_kg_perhr
  file1$sumDensAbund_N_Sqkm = file1$DensAbund_N_Sqkm
  file1$sumDensAbund_N_perhr = file1$DensAbund_N_perhr
  file1$DensBiom_kg_Sqkm <- NULL
  file1$DensBiom_kg_perhr <- NULL
  file1$DensAbund_N_Sqkm <- NULL
  file1$DensAbund_N_perhr <- NULL
    
  # Write file for Chibuzor's modelling work with lats and longs
  write.csv(file1,paste0(OUTPATHstem,"hauls_by_spp_all.csv"))
}


if(WRITE_to_DB){
  # upload to bx5 postgres database
  library(DBI)

  # append all the data that are required for an output to go into the BX5 app
  fileslisted=list.files(OUTPATHstem,'LD_tonnes_Year_W.by.subdiv',full.names=T,recursive = T);
  # In case of reruns done on previous date
  fileslisted = fileslisted[grepl(time_of_run, fileslisted)]
  file1=read.csv(fileslisted[1]) 
  # Also get survey name from filename, crude approach
  file1$survey_name = strsplit(fileslisted[1],'/')[[1]][length(strsplit(fileslisted[1],'/')[[1]])-1]
  for(file in 2:length(fileslisted)){file2=read.csv(fileslisted[file]); 
    file2$survey_name = strsplit(fileslisted[file],'/')[[1]][length(strsplit(fileslisted[file],'/')[[1]])-1]
    file2_names = colnames(file2); file1_names = colnames(file1); common_names = intersect(file2_names, file1_names); file1 = rbind(file2[common_names], file1[common_names])
  }
  
  # Add species ID from lookup table originating from database table bx005.public.species
  specieslookup=read.csv("C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/public_species_table.csv")
  specieslookup$SpeciesSciName = specieslookup$latin_name
  specieslookup$latin_name <- NULL
  specieslookup = specieslookup[,c("SpeciesSciName","species_id")]
  file1=merge(file1,specieslookup,all.x=T)

  # chosen column names to lower, this next file will be going into a postgresql database for the app:
  colnames(file1) <-tolower(colnames(file1))
  fileout = file1[,c("survstratum","year","catcatchwgtswept","fishlength_cm","dempel","order","group","survey_name","species_id")]
#   fileout=fileout[c(!is.na(fileout$species_id) & !is.na(fileout$catcatchwgtswept) & !is.na(fileout$year) & !is.na(fileout$survey_name) & !is.na(fileout$wingswparea_sqkm)) ,]
  fileout$fishlength_cm=as.integer(fileout$fishlength_cm)

  connect_t0_db <- function(connection_info_text="W:/WPx_Application/appBX005IK/ZAPP beta/bxappserver.txt"){
    connection_info <- read.csv(connection_info_text, header=T, stringsAsFactors=FALSE, sep="\t")
    #establish connection to the database we want to make the table in
    drv <- RPostgres::Postgres()
    con <- RPostgres::dbConnect( #before there was DBI
      drv,
      dbname = connection_info$res[connection_info$var == "db"],
      host = connection_info$res[connection_info$var == "host"],
      port = connection_info$res[connection_info$var == "port"],
      user = connection_info$res[connection_info$var == "user"],
      password = connection_info$res[connection_info$var == "pword"]
    )
    return(c(con, drv,connection_info))
  }
  conF<-connect_t0_db()
  con1 <-conF[[1]]
  drv <- conF[[2]]
  con_info <- conF[[3]]

  # There seem to be surveys that don't exist in the ICES data product. I don't want the app to lose functionality by losing these surveys, so download them and merge before reuploading
  cpua_original_2021 <- DBI::dbReadTable(con1, Id(schema = "wp2", table = "cpua_original_2021"))
  cpua_original_2021 = cpua_original_2021[cpua_original_2021$survey_name %in% c("heras","peltic","Q1SWOTTER","Q1SWBEAM"),]
  fileout$cpua_id = rownames(fileout)
  fileout=rbind(fileout,cpua_original_2021)

  # Also make a cpua_id
  fileout$cpua_id = rownames(fileout)

  # Write to table so field types are consistent with other table called cpua
  DBI::dbWriteTable(con1, Id(schema = "wp2", table = "cpua"), value = fileout, field.types = c(cpua_id="INTEGER PRIMARY KEY",survstratum="varchar",year="int",catcatchwgtswept="real",fishlength_cm="int",dempel="varchar",order="varchar",group="varchar",survey_name="varchar",species_id="int"), row.names=FALSE)
}

print("script complete")


# ## Code for some comparisons vs the old cpua table. This is a dataproduct from quite a different derivation so will be different due to different correction factors and SSA. But there should at least be a correlation
# cpua_original_2021 <- DBI::dbReadTable(con1, Id(schema = "wp2", table = "cpua_original_2021"))
# fileout = fileout[complete.cases(fileout), ]; cpua_original_2021 = cpua_original_2021[complete.cases(cpua_original_2021), ]
# cpua_original_2021$mergeon = paste0(cpua_original_2021$survstratum,cpua_original_2021$year,cpua_original_2021$fishlength_cm,cpua_original_2021$dempel,cpua_original_2021$order,cpua_original_2021$group,cpua_original_2021$survey_name,cpua_original_2021$species_id)
# fileout$mergeon = paste0(fileout$survstratum,fileout$year,fileout$fishlength_cm,fileout$dempel,fileout$order,fileout$group,fileout$survey_name,fileout$species_id)
# fileout$newcpua <- fileout$catcatchwgtswept; fileout$catcatchwgtswept <- NULL
# plotme=merge(fileout,cpua_original_2021,by="mergeon")
# plotme=plotme[!duplicated(plotme$mergeon),] # unique joins only, otherwise something has gone wrong in the setup of the join
# 
# # Reasonable correlation with previous
# plot(plotme$catcatchwgtswept[plotme$survey_name.x!='GNSIntOT1'],plotme$newcpua[plotme$survey_name.x!='GNSIntOT1'],ylim=c(0,50000), xlim=c(0,50000))
# 
# # GNSINTOT1 survey has a poorer agreement but still seems consistent with a corrected dataproduct
# # Drill into gnsintot1 - problem seems to be specific to 2018 and 2019 in this survey. Original cpua estimates for 2018/19 were too low, they must have been erroneous. The replacement is more consistent with other years.
# plot(plotme$catcatchwgtswept[plotme$survey_name.x=='GNSIntOT1'],plotme$newcpua[plotme$survey_name.x=='GNSIntOT1'],ylim=c(0,50000), xlim=c(0,50000))
# for(yr in unique(plotme$year.x)){plot(plotme$catcatchwgtswept[plotme$survey_name.x=='GNSIntOT1' & plotme$year.x==yr],plotme$newcpua[plotme$survey_name.x=='GNSIntOT1' & plotme$year.x==yr],ylim=c(0,10000), xlim=c(0,10000))}