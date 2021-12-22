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
# OCT 2021 - new HH output from SWEPT AREA ASSESSMENT OUTPUT [KLAAS] to combine with HL data from exchange
# NOV 2021 - define SSA (code by Joseph Ribeiro)

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

# Known Warning messages to sort out
#1: readShapeSpatial is deprecated; use rgdal::readOGR or sf::st_read 
#2: readShapePoly is deprecated; use rgdal::readOGR or sf::st_read 

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

RDIR<- dirname(parent.frame(2)$ofile) #Where you have saved the folder called R. Note this will only work if the file is SOURCED, not if it is run in the console. Alternatively, please define your WD
MAINDIR<- paste0(strsplit(RDIR,"/R")[[1]],"//")
RDIR = paste0(RDIR,"//")

definedSSA = sf::st_read(paste0(RDIR,"rectanglesICESdl29oct2021/shp_dir_for_original/SSAspatial.shp")) # read.csv(paste0(RDIR,"/defined_SSA.csv"))
definedSSA <- sf:::as_Spatial(st_zm(definedSSA))


#Trophic Level data also requires update for MTL analyses
 GOV.SPECIES.OVERWRITE<-F #"DEM" #apr2020 update to speed up and avoid PEL species and ALL in loop
 #SPECIES <- c("ALL","DEM","PEL");  if(GOV.SPECIES.OVERWRITE!=F) SPECIES <- GOV.SPECIES.OVERWRITE

#location of the subscripts in function area
 PROC_SCRIPT<- "//lowfilecds/function/Eco Indicators/DATRASoutput/" #for HH and HL processing scripts incl strata by survey
 SHAPEPATH<-paste("Strata/",sep="") #used in Lynam_OSPARsubdiv.r
 SUBSCRIPTS_TRAITS<-paste(PROC_SCRIPT,"MarScot/INDscriptsForV3/",sep="")#max length
  #read subscripts and traits
 # PRODDAT<-read.csv(paste(SUBSCRIPTS_TRAITS,"SpeciesAtLength_Productivity.csv",sep=""))
 FC1Sp<-read.csv(paste(SUBSCRIPTS_TRAITS,"FC1Specieslist.csv",sep=""))# for IA2017
 #
 SPPFILE<-"R/Species_List_Final.csv" # prepared by Meadhbh Moriarty and Simon Greenstreet for IA02017 shared in ICES WGBIODIV
 #LWFILE<-"TakFungLW - plusHakan2020.csv" #prepared by Tak Fung for LFI paper with some update by Haken Wennhage 
 #still need to include LW for all species
setwd(MAINDIR)
LW<-read.csv(SPPFILE)
 names(LW)[2]<-"ScientificName_WoRMS"
 LW$a <- LW$LWRa
 LW$b <- LW$LWRb
 LW$"Max.L..cm." <- LW$LmaxFB
 #LW_COLS<-c("ScientificName_WoRMS","a","b","Max.L..cm.") #this is used in Lynam_IND_script_process_exchange_HL2021.r
#where save output?
OUTPATHstem<-paste(MAINDIR,"out/",sep="")

## choices for analyses upfront

#do you want to write outputs along the way and save workspace
WRITE <- T #save csvs as we go?
  WRITE_LDs <- T #write Length distributions by species/year/subdiv? 
BOOTSTRAP <- F # invoke slow slow code? if F next 3 lines redundant
  NBOOT<- 9; 
  B<-0 # restart counter as only output LD once before boostrap starts i.e. when B=0 
  WRITE_BOOT <- T # every bootstrap dataset and indicator output
SAVE <- T # save workspace (after bootstrap)
FINALPLOTS<-T #create indicator plots with smooths
FILTER_COUNTRY <- T
SSA_WRITE_NEW <- F

# Catchability for general species groups
CATCHABILITY_COR_MOD<-SPECIES_IN_MOD_ONLY <-F # for comparison to ewe or lemans'
CATCHABILITY_COR_WALKER<- F # read estimates from nsea paper for q##problem somewhere looking for sweptbefore when this false

#which indicators?
  TyL_GeoM <- T # OSPAR FW3
  TyL_SPECIES<-F
  MaxL <- T # OSPAR FC3
  Loo <-F # similar to OSPAR FC3
  Lm <- F # similar to OSPAR FC3
  MEANTLs <- F #similar to OSPAR FW4 # not will not be calc'd for WAsurveys as no data file for TL
  MeanL <- F #not OSPAR but simple
  LFI_NULL<-F#T # for no group/guild calc 25% biom thresh
  LFI <- T; 
  #remove species??
  LFI_SP <-F # similar to OSPAR FC2. LFI_SP alters demersal sp list for all indicators selected
  FC1_SP<-F #if T remove others not considered in sens/resil/inter spreadsheet for OSPAR FC1
  #SPPFILE is read just after running Lynam_IND_script_process_exchange_HL2021.r 
  #and removes rare species
  
  #subscripts for indicators 
  BYGROUP<-F; if(BYGROUP){ BYGUILD<-F } # elasmos etc 
  BYICESGROUP<-F; if(BYICESGROUP){ BYGUILD<-F } # apex predators etc
  BYGUILD <- F; if(BYGUILD){ BYICESGROUP<-F; BYGROUP<-F }#not run BOOTSTRAP
  USE_GUILD_COVARIATE_SITES<- F #'to match Murray et al'
  FILL_GUILD<-F #if false have a no-guild group
  
# Spatial analyses
AREASCALE <- T # #raise by area of lowest resolution of sampling strategy/subdiv
EHDS_PP <-F #ecohydrodynamic zones
STRATA <- F  #for output
QUAD_NS <- F #only north sea
QUAD_SMOOTH_NS<- F # repeats per guild at present must change - SLOW!
QUADS<-NULL #created by Lynam_OSPARsubdiv_Feb2019.r

##NOTE additional choices below for plotting
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd(RDIR)

#LW<-   SPPLIST<-read.csv(SPPFILE)
 #read.csv(LWFILE)#LWfile needs to be complete - TakFungLW - plusHakan2020.csv is not


#additional functions
source("required.funcs.r")              # using tapply.ID
source("Lynam_INDfn_Nov2021.r") #Feb2019 now use QUAD to average biomass prior to indicators # Dec2017 update to include TAXA_GROUPINGS, Jan to calc Loo and Lm through Lynam_INDfn_Jan2018_Mtrait.r; Mar to include BX020_guilds
source("Lynam_INDPLOTFN_Nov2018.r")             # plot options 
if(BOOTSTRAP) source("Lynam_IND_BOOTfn_Aug_2017 - OSPAR.r")  # bootstrap the hauls by STSQ and subdiv

if(LFI) source("Lynam_INDfn_Dec2017_LFI.r")   # Large Fish Indicator # Nov2016 no longer fall over if no fish above LFI_THRESHOLD
 if(TyL_GeoM) source("Lynam_INDfn_Dec2017_TyL_GeoM.r")#Geometric mean length (Typical Length cm)
if(Lm & !MaxL) source("Lynam_INDfn_Jan2018_GeoMtrait.r")  # Geometric MaxL, Geometric L.infinity, Geometric L.maturity
if(MaxL | Loo) source("Lynam_INDfn_Jan2018_Mtrait.r")  # mean MaxL, mean L.infinity, mean L.maturity
if(MEANTLs) source("Lynam_INDfn_Dec2017_MeanTL.r")# TL output by rectangle and year 
 if(MeanL) source("Lynam_INDfn_Dec2017_MeanL.r") # Mean Length cm
 #if(MaxL | Loo | Lm) source("Lynam_INDfn_Jan2018_Mtrait.r")  # MaxL, Mean L.infinity, Mean L.maturity
# info on trophic level, demersal/pelagic,  maxLength, etc
#trait_file <- paste("fishspecieslifehistorytraits_v3.csv",sep='')
# should have some similar file in final dataproduct!
# dataset with traits by species for indicator script
#trait <- read.csv(trait_file)
#trait$SpeciesSciName <- paste(trait$Genus,trait$Species,sep='')

#max length observed in ospar dataproduct by species
#Hyperoplus immaculatus 'Greater sand-eel' were missing -> recorded in GNS IBTS Q1 but without length  treat as Ammodytidae 'sandeel'
trait_file <- paste("traits_by_species_Mar2019.csv",sep='')#inclelasmo taxa
trait_MAXL <- read.csv(trait_file)

if(BYGUILD){#2020 paper
  #guild_file <-  paste("feeding.guilds_12.11.18.csv",sep='')#incl trophic guild ID
guild_file <- paste("feeding.guilds_05.05.19_processed.csv",sep='')#incl trophic guild ID
guild_dat <- read.csv(guild_file)
  guild_dat$fguild <- (guild_dat$F.guild)  #fguild_7_soras.numeric
  #guild_dat<- guild_dat[guild_dat$data=="fguild6",]#1:6
}
#trophic level info
TLceltic <- TLnorth <- TLwest <- NULL
if(MEANTLs==T){
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
 #"WASpaOT3" "BBICnSpaOT4" "BBICPorOT4"  "BBICsSpaOT1" "BBICsSpaOT4" 
#"CSBBFraOT4" "CSEngBT3"    "CSIreOT4"    "CSNIrOT1"  "CSNIrOT4"    "CSScoOT1"    "CSScoOT4"    
#"GNSEngBT3"   "GNSFraOT4"  "GNSGerBT3"   "GNSIntOT1"   "GNSIntOT3"  "GNSNetBT3"   
#"WAScoOT3"     "CSFraOT4"
#SURVEY_LOOP <- c("IBTS."GNSIntOT1","CSScoOT1","CSFraOT4","CSNIrOT1")# "CSScoOT1","GNSEngBT3","CSFraOT4","CSNIrOT1")#
 
setwd(MAINDIR) #
survey_Q_C_S_combinations<-read.csv("R/survey_Q_C_S_combinations.csv")# for IA2017
for(combrow in 1:nrow(survey_Q_C_S_combinations)){#16
  #### combrow<-9 #### combrow<-nrow(survey_Q_C_S_combinations) # combrow<-2
  combs=survey_Q_C_S_combinations[combrow,]
  
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
  FIRSTYEAR = combs$First_year
  LASTYEAR = combs$Last_year #note line 401: if(survey=="BBICnSpaOT4" & FIRSTYEAR< 2017 & LASTYEAR > 2017) bio <- bio[bio$Year>=2017,] #no valid data 2010-2016
  STDGEAR = combs$Std_gear
  GEARSUBSCRIPTS = combs$Std_gear_subscripts
  GEARSUBSCRIPTS = unlist(str_split(GEARSUBSCRIPTS, ", "))  

  if(GOV.SPECIES.OVERWRITE!=F) SPECIES <- GOV.SPECIES.OVERWRITE
  if(BYGUILD) SPECIES <- c("ALL","DEM"); 
  if(FC1_SP) SPECIES <- c("DEM"); 
  
  SAMP_FILE=paste0("HH//",combs$hhfilename)
  BIOL_FILE=paste0("HL//",combs$hlfilename)
  # survey <-"IBTS"; QUARTER<-1; COUNTRY<-"Int"; SEA<-"GNS"
  # survey <-"BTS"; QUARTER<-3; COUNTRY<-"GB"; SEA<-"GNS"
  
  print(survey) #new OSPAR name as in Mar Scot dataproduct
  if(survey=="CSEngBT4") next
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
  if(MEANTLs==T & !(substr(survey,1,2) %in% c("CS","BB", "GN","IB")) ) { MEANTLs<-F; print("no data for TL for WA surveys, deleting MTL selection") }
  
  
  ##load data
  #sampling data 
  samp <- read.table(SAMP_FILE ,as.is = c(1,2,4,5,6,10,11,23),header = TRUE,sep=",") 
  samp$StatRec <- ices.rect2(ices.rect(samp$StatRec)$lon,ices.rect(samp$StatRec)$lat)

  print(max(samp$Ship))
  
  samp <- samp[samp$Quarter==QUARTER,]

  # BTS is downloaded as all countries in one file
  if(GEAR == "BEAM"){ samp<-samp[samp$Country==COUNTRY,]
    samp$DoorSpread <- samp$WingSpread <- samp$BeamWidth
    samp$SweptAreaDSKM2 <- samp$SweptAreaWSKM2 <- samp$SweptAreaBWKM2 
  }

  # where two different surveys use the same input data files the correct preprocessing needs to be done to define SSA
  # need to split SEA DATA FOR BTS Q3 ENG'
    if(survey %in% "CSBBFraOT4") samp<-samp[samp$ShootLat<= 48,]
    if(survey %in% "CSFraOT4") samp<-samp[samp$ShootLat> 48,]
    if(survey %in% "GNSIntOT1_channel") samp<-samp[samp$ShootLat<= 51,]
    if(survey %in% "CSEngBT3_Bchannel") samp<-samp[samp$ShootLat<= 52 & samp$ShootLong< -3,]
    if(survey %in% "GNSEngBT3") samp<-samp[samp$ShootLong>= -2 & samp$ShootLat>= 49.5,]
    if(survey %in% "CSEngBT3") samp<-samp[samp$ShootLong< -3 & samp$ShootLat>= 52 & samp$ShootLat< 56,]
    if(survey %in% c("CSEngBT3","CSEngBT1") ){
    samp$DoorSpread[is.na(samp$DoorSpread) | samp$DoorSpread<2] <- 4
    samp$WingSpread[is.na(samp$WingSpread) | samp$WingSpread<2] <- 4
    }
  
    #  SMFS 0816 Derivation report Step 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  samp = samp[samp$Gear %in% paste0(STDGEAR,GEARSUBSCRIPTS),]  
  samp<- samp[samp$Year >= FIRSTYEAR,] 
  samp<- samp[samp$Year <= LASTYEAR,] 
  samp<- samp[samp$HaulDur >= 13,] 
  samp<- samp[samp$HaulDur <= 66,] 
  #  SMFS 0816 Derivation report Step 1 END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  if(survey_alt_name == "BTS"){ 
    samp$SweptAreaWSKM2[is.na(samp$SweptAreaWSKM2)] <- samp$WingSpread[is.na(samp$SweptAreaWSKM2)]*samp$Distance[is.na(samp$SweptAreaWSKM2)]
    samp$SweptAreaDSKM2[is.na(samp$SweptAreaDSKM2)] <- samp$DoorSpread[is.na(samp$SweptAreaDSKM2)]*samp$Distance[is.na(samp$SweptAreaDSKM2)]
    samp$WingSwpArea_sqkm[is.na(samp$WingSwpArea_sqkm)] <- samp$SweptAreaWSKM2[is.na(samp$WingSwpArea_sqkm)]*samp$Distance[is.na(samp$WingSwpArea_sqkm)]
  }
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
  proj4string(samp) <- CRS("+init=epsg:4326")
  suppressWarnings(proj4string(rects) <- CRS("+init=epsg:4326"))
  if(SSA_WRITE_NEW){
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
    surveySSA = definedSSA[definedSSA$survey==survey,]
  }
  samp@data$ShootLong_degdec = samp$ShootLong_degdec
  samp@data$ShootLat_degdec = samp$ShootLat_degdec
  samp=samp@data
  
  
  
  # biological data
  bio <- read.csv(BIOL_FILE,as.is = c(1,2,4,5,6,10,11) )  #avoid conversions of rect to e+07 etc #,as.is=1 ) #
  
  #correction needed to remove all leading zeroes
  if(length(bio$HaulNo[substr(bio$HaulNo,1,1)%in%"0"])>0) bio$HaulNo <- str_replace(bio$HaulNo, "^0+" ,"") 
  if(length(bio$StNo[substr(bio$StNo,1,1)%in%"0"])>0) bio$StNo <- str_replace(bio$StNo, "^0+" ,"") 
  if(length(samp$HaulNo[substr(samp$HaulNo,1,1)%in%"0"])>0) samp$HaulNo <- str_replace(samp$HaulNo, "^0+" ,"") 
  if(length(samp$StNo[substr(samp$StNo,1,1)%in%"0"])>0) samp$StNo <- str_replace(samp$StNo, "^0+" ,"") 
  
  bio <- bio[bio$Quarter==QUARTER,]
  bio <- bio[!is.na(bio$ValidAphiaID),] #remove invalids
  bio <- bio[(bio$SpecVal %in% c(1,4,7,10) ),]#remove invalids 0=Invalid information	 2=Partly valid information	   https://vocab.ices.dk/?ref=5
  if(survey=="BBICnSpaOT4" & FIRSTYEAR< 2017 & LASTYEAR > 2017) bio <- bio[bio$Year>=2017,] #no valid data 2010-2016
  #summary( bio$HaulNo)
  #unique( bio$StNo)
  bio$HaulID	<- paste(paste(survey,QUARTER,sep=""), bio$Country, #bio$Gear, 
                      bio$Ship,  bio$StNo, bio$HaulNo, bio$Year,sep=":")
  
  samp$HaulID <- paste(paste(survey,QUARTER,sep=""),samp$Country,#samp$Gear,
                       samp$Ship, samp$StNo, samp$HaulNo, samp$YearShot,sep=":" )#HaulNo and StNo ,samp$Month,samp$Day
  # KLASS/RUTH use:  Survey,Country,Quarter,Gear,Ship,StNo,HaulNo,Year,Month,Day
  
  haulidhl <- sort(unique(samp$HaulID))
  haulidhh <- sort(unique(bio$HaulID))
  
  haulidhhm <- subset(samp, !samp$HaulID %in% bio$HaulID)
  haulidhhm <- sort( unique(haulidhhm$HaulID) ); length(haulidhhm) 
  100*length(haulidhhm)/nrow(samp) # 2samp hauls not match in bio
  
  haulidhlm <- subset(bio, !bio$HaulID %in% samp$HaulID)# bio hauls not match in samp
  haulidhlm <- sort(unique(haulidhlm$HaulID) ); length(haulidhlm) 
  100*length(haulidhlm)/nrow(bio) #samp hauls not match in bio
  
  
   SPPLIST <- LW$ScientificName_WoRMS
   #check which species would be removed by SPPLIST
   bio$SpeciesSciName <- ac(bio$ScientificName_WoRMS)#edit
   if(nrow(bio[bio$SpeciesSciName == "Aspitrigla cuculus",]) >0) bio[bio$SpeciesSciName == "Aspitrigla cuculus",]$SpeciesSciName  <- "Chelidonichthys cuculus"
   if(nrow(bio[bio$SpeciesSciName == "Chelidonichthys lucernus",]) >0) bio[bio$SpeciesSciName == "Chelidonichthys lucernus",]$SpeciesSciName <- "Chelidonichthys lucerna"
   if(nrow(bio[bio$SpeciesSciName == "Phrynorhombus norvegicus",]) >0) bio[bio$SpeciesSciName == "Phrynorhombus norvegicus",]$SpeciesSciName <- "Zeugopterus norvegicus"
   if(nrow(bio[bio$SpeciesSciName == "Psetta maxima",]) >0) bio[bio$SpeciesSciName == "Psetta maxima",]$SpeciesSciName <- "Scophthalmus maximus"
   if(nrow(bio[bio$SpeciesSciName == "Mullus barbatus",]) >0) bio[bio$SpeciesSciName == "Mullus barbatus",]$SpeciesSciName <- "Mullus barbatus barbatus"
   if(nrow(bio[bio$SpeciesSciName == "Trisopterus esmarki",]) >0) bio[bio$SpeciesSciName == "Trisopterus esmarki",]$SpeciesSciName <- "Trisopterus esmarkii"
   if(nrow(bio[bio$SpeciesSciName == "Solea",]) >0) bio[bio$SpeciesSciName == "Solea",]$SpeciesSciName <- "Solea solea"
   if(nrow(bio[bio$SpeciesSciName == "Solea vulgaris",]) >0) bio[bio$SpeciesSciName == "Solea vulgaris",]$SpeciesSciName <- "Solea solea"
   if(nrow(bio[bio$SpeciesSciName == "Dipturus flossada",]) >0) bio[bio$SpeciesSciName == "Dipturus flossada",]$SpeciesSciName <- "Dipturus"
   if(nrow(bio[bio$SpeciesSciName == "Dipturus intermedia",]) >0) bio[bio$SpeciesSciName == "Dipturus intermedia",]$SpeciesSciName <- "Dipturus"
   if(nrow(bio[bio$SpeciesSciName == "Microchirus (Microchirus) variegatus",]) >0) bio[bio$SpeciesSciName == "Microchirus (Microchirus) variegatus",]$SpeciesSciName <- "Microchirus variegatus"
   if(nrow(bio[bio$SpeciesSciName == "Malacocephalus (Malacocephalus) laevis",]) >0) bio[bio$SpeciesSciName == "Malacocephalus (Malacocephalus) laevis",]$SpeciesSciName <- "Malacocephalus laevis"
   if(nrow(bio[bio$SpeciesSciName == "Engraulis",]) >0) bio[bio$SpeciesSciName == "Engraulis",]$SpeciesSciName <- "Engraulis albidus"
   if(nrow(bio[bio$SpeciesSciName == "Trisopterus esmarki",]) >0) bio[bio$SpeciesSciName == "Trisopterus esmarki",]$SpeciesSciName <- "Trisopterus esmarkii"
   bio$ScientificName_WoRMS <- (bio$SpeciesSciName)#edit
   
   #happy to lose snail "Liparis liparis"   #wo-spotted clingfish "Diplecogaster bimaculata"   #three-spined stickleback "Gasterosteus aculeatus"   etc
   
   BIONAM<-ac(unique(bio$SpeciesSciName))
   sort(BIONAM[!BIONAM %in% SPPLIST])#species excluded 
   write.table(BIONAM[!BIONAM %in% SPPLIST],paste(OUTPATH,"lostspp_",survey,"_Q",QUARTER,".txt",sep=""))
   
   bio <- bio[bio$SpeciesSciName %in% SPPLIST,]#remove rare spp keep 550 species
   bio$SpeciesSciName<-af(ac(bio$SpeciesSciName))  #relevel 131 from >700
   
   print("process_exchange_HL") #moved down
   source('R/Lynam_IND_script_process_exchange_HL2021.R')
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#quick eyeball of sampling data 
  # quick check   unique(bio[bio$LngtCode==1,"SubFactor"])
  #par(mfrow=c(1,1))
  x11(); 
  if(substr(survey,7,8)!="BT" & substr(survey,6,7)!="BT" & substr(survey,8,9)!="BT"){ 
    par(mfrow=c(4,4)) } else { par(mfrow=c(2,4)) }
  plot(samp$ShootLong_degdec, samp$ShootLat_degdec); map(add=T)#,xlim=c(4,14)
  abline(v=c(4,8));  abline(h=55.5)
  with(samp, hist(MonthShot))
  with(samp, hist(HaulDur_min))
  with(samp, hist(Depth_m))
  if(substr(survey,7,8)!="BT" & substr(survey,6,7)!="BT" & substr(survey,8,9)!="BT"){ 
   with(samp, hist(WingSpread_m))
   with(samp, hist(DoorSpread_m))
   with(samp, hist(NetOpen_m))
   with(samp, hist(WingSwpArea_sqkm) )
   with(samp, hist(WingSwpVol_CorF ))
   with(samp, hist(DoorSwptArea_CorF) )
   with(samp, hist(DoorSwptVol_CorF))
  }
  with(samp, hist(Distance_km))
  #samp[(samp$Distance_km)==max(samp$Distance_km),]
  #samp[(samp$Distance_km)>5,] # big value
  ##x11(); par(mfrow=c(2,3))
  with(bio, hist(FishLength_cm))
  #with(bio, hist(IndivFishWgt_g))
  #bio[!is.na(bio$IndivFishWgt_g) & (bio$IndivFishWgt_g)==max(bio$IndivFishWgt_g,na.rm=T),]
  with(bio, hist(DensAbund_N_Sqkm))
  #bio[!is.na(bio$DensAbund_N_Sqkm) & (bio$DensAbund_N_Sqkm)==max(bio$DensAbund_N_Sqkm,na.rm=T),]
  # samp[samp$HaulID == bio[(bio$DensAbund_N_Sqkm)==max(bio$DensAbund_N_Sqkm),'HaulID'], ] # sometimes get multiple matches e.g sprat on north sea IBTS Q1
  with(bio, hist(DensBiom_kg_Sqkm))
  #bio[!is.na(bio$DensBiom_kg_Sqkm) & (bio$DensBiom_kg_Sqkm)==max(bio$DensBiom_kg_Sqkm,na.rm=T),]
  #with(bio, hist(Number))
  #eyeball end
  savePlot(filename= paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),"eyeball",".bmp",sep=''),type="bmp")
  dev.off()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #now merge datasets on haul and species for indicator script
  #bio <- read.csv(BIOL_FILE) #reduce memory usage
  bio <- bio[,c("HaulID","SpeciesSciName","FishLength_cm","DensBiom_kg_Sqkm","DensAbund_N_Sqkm","SubFactor")]
  bio$sciName<-as.character(bio$SpeciesSciName)
  
  
  samp <-samp[,c("HaulDur_min","DataType","HaulID","YearShot","ShootLat_degdec","ShootLong_degdec",
                 "ICESStSq","SurvStratum","WingSwpArea_sqkm","WingSwpVol_CorF","NetOpen_m","Gear","Ship","YearShot","MonthShot","Day","TimeShot")] 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  #add efficiency of E=GOV gear
  #use gear efficiency to correct catches
  #the probability that fish in the path of a trawl will be caught and retained
  if(CATCHABILITY_COR_WALKER) print("now CATCHABILITY correct data with WALKER")
  if(CATCHABILITY_COR_WALKER) source(paste(PROC_SCRIPT,"Lynam_IND_script_CATCHABILITY_COR_WALKER.R",sep="/"))
  
  ave_NetOpen_m<-mean(samp$NetOpen_m)
  #samp$WingSwpVol_CorF <- ave_NetOpen_m / samp$NetOpen_m # scale down if larger than usual net opening
  #"SubFactor",  
  dhspp <- merge(bio,samp,by="HaulID",all = FALSE)
  # plot(dhspp$ShootLong_degdec, dhspp$ShootLat_degdec); map(add=T)#,xlim=c(4,14)
  #with(dhspp, xyplot(ShootLat_degdec~ShootLong_degdec | ac(YearShot) ))
  lostID_from_lack_of_matching_haulID<-unique(bio[!(bio$HaulID %in% dhspp$HaulID),"HaulID"])#all bio matched with StnNo
   # JR edit - dhspp contains NA lats and longs
  dhspp = dhspp[is.finite(dhspp$ShootLat_degdec) & is.finite(dhspp$ShootLong_degdec),]
   
  if(length(lostID_from_lack_of_matching_haulID) >0){print(paste("losing", length(lostID_from_lack_of_matching_haulID) ,"hauls from",length(unique(bio$HaulID)),"when merge bio and samp")) } else { print("successful merge HL and HH to create dhspp")}
  write.table(lostID_from_lack_of_matching_haulID,paste(OUTPATH,"lostID_from_lack_of_matching_haulID_",survey,"_Q",QUARTER,".txt",sep=""))
  rm(bio,samp)
  
  # Make spatial, filter by SSA, then drop spatial aspect
  predhspp <- dhspp
  coordinates(dhspp) <- ~ ShootLong_degdec + ShootLat_degdec
  proj4string(dhspp) <- CRS("+init=epsg:4326")
  dhspp = dhspp[surveySSA,]
  dhspp@data$ShootLong_degdec = dhspp$ShootLong_degdec
  dhspp@data$ShootLat_degdec = dhspp$ShootLat_degdec
  dhspp=dhspp@data

  lostID_from_SSA_filtering<-unique(predhspp[!(predhspp$HaulID %in% dhspp$HaulID),"HaulID"])#all predhspp matched with StnNo
  if(length(lostID_from_SSA_filtering) >0){print(paste("losing", length(lostID_from_SSA_filtering) ,"hauls from",length(unique(predhspp$HaulID)),"when merge predhspp and samp")) } else { print("successful merge HL and HH to create dhspp")}
  write.table(lostID_from_SSA_filtering,paste(OUTPATH,"lostID_from_SSA_filtering_",survey,"_Q",QUARTER,".txt",sep=""))
  rm(predhspp)
  
  #bioraw$HaulID<- paste(paste(survey,QUARTER,sep=""), bioraw$Ship, bioraw$YearShot, bioraw$HaulNo,sep="/")#
  #bioraw[(bioraw$HaulID%in% lostID),]

    #LOSE RECORDS MISSING LENGTH (already done above) and MISSING DATATYPE
    #dhspp <- dhspp[!(dhspp$FishLength_cm== -9 & !is.na(dhspp$FishLength_cm)),] 
    #dhspp <- dhspp[!is.na(dhspp$DataType),]

    # subfactor to raise measured num at length to total num fish in haul 
    # datatype C is per hour #dhspp$DataType == C DATA ARE CPUE, SO
    # if R or S: HLNoAtLngt*SubFactor   
    dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensAbund_N_Sqkm),]$DensAbund_N_Sqkm <-
    dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensAbund_N_Sqkm),]$DensAbund_N_Sqkm *
    dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensAbund_N_Sqkm),]$SubFactor
    
    dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensBiom_kg_Sqkm),]$DensBiom_kg_Sqkm <-
      dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensBiom_kg_Sqkm),]$DensBiom_kg_Sqkm *
      dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensBiom_kg_Sqkm),]$SubFactor
    
    #standardise per 60 min tow with DurRaise
    dhspp$DurRaise <- dhspp$HaulDur_min/60;  # per hour # hist(dhspp$DurRaise)  
    if(nrow(dhspp[dhspp$DataType=="C" ,])>0) dhspp[dhspp$DataType=="C" ,]$DurRaise <- 1 # already per hour
    dhspp$DensAbund_N_Sqkm <- dhspp$DensAbund_N_Sqkm/dhspp$DurRaise  #num per 60 min as in 'C'
    dhspp$DensBiom_kg_Sqkm <- dhspp$DensBiom_kg_Sqkm/dhspp$DurRaise  #CPUE per 60 min
    
    ##add swept area from ICES?
    dhspp$DensAbund_N_Sqkm<-( dhspp$DensAbund_N_Sqkm*dhspp$DurRaise)/ dhspp$WingSwpArea_sqkm
    dhspp$DensBiom_kg_Sqkm<-( dhspp$DensBiom_kg_Sqkm*dhspp$DurRaise)/ dhspp$WingSwpArea_sqkm
    
  #rm(samp,bio)
  dhspp <- merge(x=dhspp, y=SciName2GroupEff[,c(1,9)], by.x="SpeciesSciName", by.y="sciName",all.x=T)
  if(SPECIES_IN_MOD_ONLY) dhspp <- dhspp[dhspp$SpeciesSciName %in% SPECIES_SUBSET,] #use only those from thorpe model
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #simplify names
  names(dhspp)[which(names(dhspp) =="YearShot")] <-"Year"
  #for tyl
  dhspp$LogLngtClass <- log(dhspp$FishLength_cm)
  dhspp$LogLngtBio <- dhspp$LogLngtClass * dhspp$DensBiom_kg_Sqkm
  
  #### link subdivisional shapefiles to data ####
  # and read attributes of shapefiles create table ATTRIB with names SurvStratum & KM2_LAM 
  #source("//lowfilecds/Function/Eco Indicators/DATRASoutput/MarScot/INDscriptsForV3/Lynam_OSPARsubdiv_Jan2018.r")  
  
  #now use quadrants## load("C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR/biodiv19 ref/RUNuptoOSPARsubdiv.rdata.RData")
  #if(survey=="IBTS"){ survey<-"GNSIntOT1" }
  ##### strata ##### 
  print("now add strata")
  #source("//lowfilecds/Function/Eco Indicators/DATRASoutput/MarScot/INDscriptsForV3/Lynam_OSPARsubdiv_Feb2019.r")
  source(paste(MAINDIR,"R/Lynam_OSPARsubdiv_Oct2021.r",sep=""))
  
  #replace sampstrat with  quadrants if QUAD==T
  #if(SURVEY=="IBTS"){ survey<-SURVEY }
  

  #MaxLength and demersal or pelagic
  #MaxLength observed for MML
  names(dhspp)[names(dhspp)=="Group"] <- "QGrouping"
  dhspp <- merge(dhspp,trait_MAXL[,c("SpeciesSciName","maximum.length", "DEMPEL","Order","Group", "LFI_Fung_list", "Loo","K","Winfinity","tmax","tm","M","Lm","LFI_OSPAR_list")],by="SpeciesSciName",all.x=T)
  names(dhspp)[which( names(dhspp)=="maximum.length")] <- "MaxL"
  #DEM are ecotypes Demersal + Bathydemersal + Bathypelagic  (except Micromesistius poutassou Blue whiting) + Benthopelagic (except Clupea harengus herring)
  #PEL are Pelagic  #sprat, mac, hmx, her, whb etc
  dhspp$DEMPEL <-as.character(dhspp$DEMPEL)
  dhspp$DEMPEL[dhspp$DEMPEL=="Demersal"] <- "DEM"
  dhspp$DEMPEL[dhspp$DEMPEL=="Pelagic"] <- "PEL"
  dhspp$DEMPEL <-as.character(dhspp$DEMPEL)
  if(LFI_SP){ dhspp$DEMPEL[dhspp$LFI_OSPAR_list=="Demersal"] <- "DEM"
              dhspp$DEMPEL[dhspp$LFI_OSPAR_list=="Other"] <- "Other"
              dhspp$DEMPEL[dhspp$LFI_OSPAR_list=="Pelagic"] <- "PEL"
  }#unique( dhspp[dhspp$LFI_OSPAR_list=="Demersal",c("DEMPEL","SpeciesSciName")])
  #also check line 256-269 of INDfn
  if(LFI_SP) dhspp$DEMPEL[dhspp$DEMPEL=="DEM" & dhspp$LFI_Fung_list!="Demersal"] <- "Other"; 
  if(FC1_SP) dhspp <- dhspp[dhspp$SpeciesSciName %in% FC1Sp$SpeciesSciName,]
  #for pelagics use volume correction
  # density*samp$WingSwptVol_CorF
  #if(survey!="IBTS") dhspp[dhspp$DEMPEL %in% "PEL",]$DensBiom_kg_Sqkm <- dhspp[dhspp$DEMPEL %in% "PEL",]$DensBiom_kg_Sqkm*dhspp[dhspp$DEMPEL %in% "PEL",]$WingSwpVol_CorF
  if(CATCHABILITY_COR_WALKER) dhspp[dhspp$DEMPEL %in% "PEL",]$DensBiom_kg_Sqkm_beforeQmult <- dhspp[dhspp$DEMPEL %in% "PEL",]$DensBiom_kg_Sqkm_beforeQmult*dhspp[dhspp$DEMPEL %in% "PEL",]$WingSwpVol_CorF
  if(CATCHABILITY_COR_WALKER) dhspp[dhspp$DEMPEL %in% "PEL",]$DensAbund_N_Sqkm_beforeQmult <- dhspp[dhspp$DEMPEL %in% "PEL",]$DensAbund_N_Sqkm_beforeQmult*dhspp[dhspp$DEMPEL %in% "PEL",]$WingSwpVol_CorF

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
  FILENAM<-paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),sep="")
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
                      LFI_THRESHOLD=LFI_THRESHOLD_APPLY, LFI_SP = F,
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
          IND_OUT_BYICESGROUP[[SPGROUP]] <- INDfn( DATA=dhspp, WRITE=T, SPECIES=SPECIES, GROUP=SPGROUP,
                      LFI=LFI, LFI_THRESHOLD=NULL, LFI_SP = F, FILENAM=FILENAM,
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
      #load("Z:/Foodweb Models/Feeding guilds_BX020/Script/results/Ices_rectangles_for_analysis.RData")
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
        IND_OUT_BYGUILD[[SPGROUP]] <- INDfn( DATA=dhspp, WRITE=T, SPECIES=SPECIES, GROUP=SPGROUP,
                                             LFI=F, LFI_THRESHOLD=NULL, LFI_SP = F, FILENAM=FILENAM,
                                             SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, 
                                             MEANTL=F, MaxL=T, Loo=F, Lm=F, MeanL=F, TyL_GeoM=T,BYGUILD=BYGUILD,
                                             QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS)
        ,silent=TRUE)
    }
    dhspp <- dhspp_raw
    }
  
  ##### NO TAXA GROUPING ##### 
  #do this last to make sure have all species in final biomass plots
  print("NO TAXA GROUPING")
  if(LFI_NULL) LFI_THRESHOLD<-NULL
   try(
    IND_OUT <- INDfn( DATA=dhspp, WRITE=T, BOOT=F, LFI=LFI, LFI_THRESHOLD=LFI_THRESHOLD, LFI_SP = LFI_SP,
                      FILENAM=FILENAM,SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, 
                    MEANTL=MEANTLs, MaxL=MaxL,Loo=Loo, Lm=Lm, MeanL=MeanL, TyL_GeoM=TyL_GeoM, SPECIES=SPECIES, 
                    GROUP=NULL, TyL_SPECIES=TyL_SPECIES, BYGUILD=F, QUAD=QUAD,QUAD_SMOOTH=QUAD_SMOOTH,QUADS=QUADS)
   ,silent=F)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##### BOOTSTRAP ####
  #lists to collate multiple bootstrapped indicators
  BOOTDATA2PLOT <- NULL
  LFI_regional<- TL_regional <- TyLrect_regional <- TyL_regional <- TyL_reg_var <- Len_regional <- MaxL_regional <- Loo_regional <- Lm_regional <- NULL
  if(BOOTSTRAP){ print("Bootstrap dataset"); source(paste(PROC_SCRIPT,"Lynam_IND_script_BOOTSTRAP.r",sep="")) }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  #invesigate an indicator with GAM    library(mgcv)
  #and plot some (with bootstrap?)
  #IND_OUT <- IND_OUT_BYGROUP[["Elasmobranchii"]]; BOOT_OUT <- BOOT_OUT_BYGROUP[["Elasmobranchii"]]
  if(FINALPLOTS){ print("Final plots"); try( source(paste(MAINDIR,"R/Lynam_IND_script_FINALPLOTS.R",sep="")) ,silent=F) }
  print(paste("Finished",survey, "survey",sep=" "))
  dev.off()
  if(SSA_WRITE_NEW) { SSAdf=rbind(SSAdf,SSAdfs)}
  write.csv(dhspp,paste0(MAINDIR,"out/haul_by_spp_",survey,".csv"),row.names = F)
  rm(dhspp)
} #next survey

if(SSA_WRITE_NEW) { SSAdf=SSAdf[-1,] # first row is null from setup
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




# append all the outputs that relate to Chibuzor modelling work. Column for sampstrat and SurvStratum is inconsistent across files, so dropping these. 
fileslisted=list.files(OUTPATHstem,'.csv',full.names=T); file1=read.csv(fileslisted[1]); file1 = file1[ , -which(names(file1) %in% c("sampstrat","SurvStratum"))]
for(file in 2:length(fileslisted)){file2=read.csv(fileslisted[file]); file2 = file2[ , -which(names(file2) %in% c("sampstrat","SurvStratum"))]; file1=rbind(file1,file2)}
write.csv(file1,paste0(OUTPATHstem,"haul_by_spp_all.csv"))
print("script complete")


