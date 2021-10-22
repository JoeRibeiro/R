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

rm(list=ls()) #start afresh

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
MAINDIR<- "C:/Users/cl06/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/"
#FLEX<-T  #?use ICES swept areas

# sampling data and biological data 
SURVEY_LOOP <- "BTS"#"SP-ARSA","SCOROC","IBTS"  )
QUARTER_LOOP <-c(1,2,3,4) #c(1,2,3,4) # 3#c(1,3) 
COUNTRY_LOOP <- c("Int","GB","NL","GE","BE"); #important for BTS where data for seperate surveys combined in downloads
SEA_LOOP <- c("GNS")#,"CS")#  #important for EVHOE as split across OSPAR assessment regions and for BTS where data for seperate surveys combined in downloads

#Trophic Level data also requires update for MTL analyses
 GOV.SPECIES.OVERWRITE<-"DEM" #apr2020 update to speed up and avoid PEL species and ALL in loop
 #SPECIES <- c("ALL","DEM","PEL");  if(GOV.SPECIES.OVERWRITE!=F) SPECIES <- GOV.SPECIES.OVERWRITE

#location of the subscripts in function area
 PROC_SCRIPT<- "//lowfilecds/function/Eco Indicators/DATRASoutput/" #for HH and HL processing scripts incl strata by survey
 SHAPEPATH<-paste(PROC_SCRIPT,"Strata/",sep="") #used in Lynam_OSPARsubdiv.r
 SUBSCRIPTS_TRAITS<-paste(PROC_SCRIPT,"MarScot/INDscriptsForV3/",sep="")#max length #location of the shapefiles and where is the 'attributes' folder
 #read subscripts and traits
 PRODDAT<-read.csv(paste(SUBSCRIPTS_TRAITS,"SpeciesAtLength_Productivity.csv",sep=""))
 FC1Sp<-read.csv(paste(SUBSCRIPTS_TRAITS,"FC1Specieslist.csv",sep=""))# for IA2017
 #
 LWFILE<-"TakFungLW - plusHakan2020.csv" #prepared by Tak Fung for LFI paper with some update by Haken Wennhage 
 #still need to include LW for all species
 
#where save output?
OUTPATHstem<-paste(MAINDIR,"out/",sep="")


#do you want to write outputs along the way and save workspace
WRITE <- T #save csvs as we go?
  WRITE_LDs <- F #write Length distributions by species/year/subdiv? 
BOOTSTRAP <- F # invoke slow slow code? if F next 3 lines redundant
  NBOOT<- 9; 
  B<-0 # restart counter as only output LD once before boostrap starts i.e. when B=0 
  WRITE_BOOT <- T # every bootstrap dataset and indicator output
SAVE <- T # save workspace (after bootstrap)
FINALPLOTS<-T #create indicator plots with smooths

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
  LFI_NULL<-T # for no group/guild calc 25% biom thresh
  LFI <- T; 
  #remove species??
  LFI_SP <-F # similar to OSPAR FC2. LFI_SP alters demersal sp list for all indicators selected
  FC1_SP<-F #if T remove others not considered in sens/resil/inter spreadsheet for OSPAR FC1
  SPPFILE<-"R/SpeciesInfoSG.csv" # prepared by Simon Greenstreet for IA02017 shared in ICES WGBIODIV
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
setwd(MAINDIR)
setwd("R")

LW<-read.csv(LWFILE)#LWfile needs to be complete - TakFungLW - plusHakan2020.csv is not

#additional functions
source("required.funcs.r")              # using tapply.ID
source("Lynam_INDfn_Oct2021_guild.r") #Feb2019 now use QUAD to average biomass prior to indicators # Dec2017 update to include TAXA_GROUPINGS, Jan to calc Loo and Lm through Lynam_INDfn_Jan2018_Mtrait.r; Mar to include BX020_guilds
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
trait_file <- paste(SUBSCRIPTS_TRAITS,"traits_by_species_Mar2019.csv",sep='')#inclelasmo taxa
trait_MAXL <- read.csv(trait_file)

if(BYGUILD){#2020 paper
  #guild_file <-  paste("feeding.guilds_12.11.18.csv",sep='')#incl trophic guild ID
  guild_file <-  paste(SUBSCRIPTS_TRAITS,"feeding.guilds_05.05.19_processed.csv",sep='')#incl trophic guild ID
  guild_dat <- read.csv(guild_file)
  guild_dat$fguild <- (guild_dat$F.guild)  #fguild_7_soras.numeric
  #guild_dat<- guild_dat[guild_dat$data=="fguild6",]#1:6
}
#trophic level info
TLceltic <- TLnorth <- TLwest <- NULL
if(MEANTLs==T){
  #cnrs mtl fw4
  TLceltic <- read.csv(paste(SUBSCRIPTS_TRAITS,"TL_complete14102016-1.csv",sep=""))
  names(TLceltic)[1] <- "SpeciesSciName"
  TLceltic <- TLceltic[,c("SpeciesSciName","TL")]
  
  # west of scotland
  TLwest <- NULL
  
  # north sea
  TLnorth <- read.csv(paste(SUBSCRIPTS_TRAITS,"refspp_LWmerged_edit06Nov2015.csv",sep=""))
  names(TLnorth)[2] <- "SpeciesSciName"
  ECOPATH_TL <- F; if(ECOPATH_TL)  TLnorth$TL <- TLnorth$TL_Ecopath #for comparison with ewe
  TLnorth <- TLnorth[,c("SpeciesSciName","TL")]
}

# general species groups # Jan 2018 added 2 stingray Dasyatis spp to lookup and Cetorhinus maximus	Basking shark and Alopias vulpinus to SciName2GroupEff to stop error
 SciName2GroupEff <- read.csv(paste(SUBSCRIPTS_TRAITS,"SciName2GroupEff_Jan2018.csv",sep=""))
 SciName2GroupEff$Group <- paste("GRP",SciName2GroupEff$Group, sep="")
 SciName2GroupEff$sciName<-as.character(SciName2GroupEff$sciName)
 substr(SciName2GroupEff$sciName,1,1) <- toupper(substr(SciName2GroupEff$sciName,1,1))#correct case
##select only those species from LeMans model and correct relative biomasses to match: Q
if(CATCHABILITY_COR_MOD | SPECIES_IN_MOD_ONLY) source("Lynam_IND_script_CATCHABILITY_MODEL.R")
  

#### indicators survey loop ####"WASpaOT3" "BBICnSpaOT4" "BBICPorOT4"  "BBICsSpaOT1" "BBICsSpaOT4" 
#"CSBBFraOT4" "CSEngBT3"    "CSIreOT4"    "CSNIrOT1"  "CSNIrOT4"    "CSScoOT1"    "CSScoOT4"    
#"GNSEngBT3"   "GNSFraOT4"  "GNSGerBT3"   "GNSIntOT1"   "GNSIntOT3"  "GNSNetBT3"   
#"WAScoOT3"     "CSFraOT4"
#SURVEY_LOOP <- c("IBTS."GNSIntOT1","CSScoOT1","CSFraOT4","CSNIrOT1")# "CSScoOT1","GNSEngBT3","CSFraOT4","CSNIrOT1")#
setwd(MAINDIR) #
for(QUARTER in QUARTER_LOOP){ 
for(COUNTRY in COUNTRY_LOOP){
for(SEA in SEA_LOOP){
  for(survey in SURVEY_LOOP){ #need survey last so that once overwritten below resets for next Quarter/Country/Sea
    # survey <-"IBTS"; QUARTER<-1; COUNTRY<-"Int"; SEA<-"GNS"
    # survey <-"BTS"; QUARTER<-3; COUNTRY<-"GB"; SEA<-"GNS"
    if(survey!="BTS" & COUNTRY!=COUNTRY_LOOP[1]) next
    if(survey=="BTS" & COUNTRY=="Int") next
    
    print(paste(survey," Q",QUARTER,sep=""))
    SAMP_FILE   <- paste(MAINDIR,"HH/HH-",survey,".csv",sep="")  
    BIOL_FILE   <- paste(MAINDIR,"HL/HL-",survey,".csv",sep="")
    #rename to OSPAR survey names:
    if(survey=="IBTS"){  if(COUNTRY!="Int") next
      survey<-"GNSIntOT"
      SSA <- read.csv(paste(MAINDIR,"R/IBTS_SSA.csv",sep="")) # to reduce output for standard survey area as in IA2017 
      #SSA is not currently implemented and datafiles are required for other surveys beyond GNSIntOT1 and GNSIntOT3
      }   
    
    if(survey=="SP-PORC"){  survey<-"WASpaOT" }
    if(survey=="SP-NORTH"){  survey<-"BBICnSpaOT"} 
    if(survey=="SP-ARSA"){  survey<-"BBICsSpaOT"} 
    if(survey=="SP-ARSA"){  survey<-"BBICsSpaOT"}  
    
    if(survey=="PT-IBTS"){  survey<-"BBICPorOT"}  
   
    if(survey=="IE-IGFS"){  survey<-"CSIreOT"}    
    
    if(survey=="NIGFS"){  survey<-"CSNIrOT"}   
    
    if(survey=="SCOWCGFS"){  survey<-"CSScoOT"} #from 2010
    if(survey=="SWC-IBTS"){  survey<-"CSScoOT"}
    if(survey=="SCOROC"){  survey<-"WAScoOT"}  #from 2011 ROCKALL previous years
    
    
    if(survey=="FR-CGFS" ){  survey<-"GNSFraOT"}
    if(survey=="BTS-VIII"){  survey<-"CSFraBT"}   
    if(survey=="EVHOE" & SEA=="CS"){  survey<-"CSFraOT"}
    if(survey=="EVHOE" & SEA!="CS"){  survey<-"CSBBFraOT"} # this is split from above for OSPAR assessment - need to edit
    
    if(survey=="BTS" & SEA=="GNS"& COUNTRY=="GB"){  survey<-"GNSEngBT"}#north sea no data for Q1 as this is CS
    if(survey=="BTS" & SEA=="CS" & COUNTRY=="GB"){  survey<-"CSEngBT"}#new Eng Q4 data not in IA2017 for Celtic sea beam trawl for SW  
    if(survey=="BTS" & COUNTRY=="GE"){  survey<-"GNSGerBT" }  
    if(survey=="BTS" & COUNTRY=="NL"){  survey<-"GNSNetBT" }  
    if(survey=="BTS" & COUNTRY=="BE"){  survey<-"GNSBelBT" }#new data not in IA2017 as not previously in DATRAS
    
    survey <- paste(survey,QUARTER,sep="")
  print(survey) #new OSPAR name as in Mar Scot dataproduct
  if(survey=="GNSEngBT1") next ## BTS ENG Q1 is CS
  
  OUTPATH <- paste(OUTPATHstem,survey,"/",sep='')
  LFI_THRESHOLD <- 40 
  if(survey %in% c("GNSGerBT3","GNSNetBT3","GNSBelBT3","GNSEngBT3","CSEngBT1","CSEngBT3")){
    SPECIES <- "DEM"; GEAR<-"BEAM"
    if(survey %in% c("GNSNetBT3") ) LFI_THRESHOLD <- 45
    if(survey %in% c("CSEngBT3") ) LFI_THRESHOLD <- 35 
    if(survey %in% c("GNSGerBT3") ) LFI_THRESHOLD <- 50 
  } else { 
    GEAR<-"GOV" #or at least gov like e.g bak 
    SPECIES <- c("ALL","DEM","PEL");  if(GOV.SPECIES.OVERWRITE!=F) SPECIES <- GOV.SPECIES.OVERWRITE
    if(BYGUILD) SPECIES <- c("ALL","DEM"); 
    if(FC1_SP) SPECIES <- c("DEM"); 
    if(survey %in% c("BBICsSpaOT1","GNSFraOT4","GNSEngBT3",
                     "GNSIntOT1","GNSIntOT3","IBTS"))  LFI_THRESHOLD <- 50  #GNSIntOT1 50 UKBI #note use of LFIregion SUBDIV in ...Lynam_OSPARsubdiv.r
    if(survey %in% c("BBICsSpaOT4","CSNIrOT1"))  LFI_THRESHOLD <- 45 #CSNIrOT1 45 UKBI
    if(survey %in% c("BBICPorOT4","CSBBFraOT4","CSNIrOT4",
                     "CSScoOT4","WASpaOT3"))  LFI_THRESHOLD <- 40 
    if(survey %in% c("CSScoOT1","WAScoOT3"))  LFI_THRESHOLD <- 35  #CSScoOT1 35 UKBI
    if(survey %in% c("CSIreOT4") ) LFI_THRESHOLD <- 30 
    if(survey %in% c("BBICnSpaOT4"))  LFI_THRESHOLD <- 25
    if(survey %in% c("CSFraOT4"))  LFI_THRESHOLD <- 40 #CSFraOT4 40 UKBI
  } 
  
  if(survey %in% c("GNSGerBT3","GNSBelBT3","GNSNetBT3","GNSIntOT1","GNSIntOT3","GNSEngBT3","IBTS")){#,
    SAMP_STRAT <- T # average hauls by ICESStSq rect in north sea #if set to FALSE need to update Attibutes table as area only given by rect
    BYSUBDIV <-   T # average indicator by LFI-subdivision
    QUAD<- QUAD_NS #only north sea
    QUAD_SMOOTH<-QUAD_SMOOTH_NS
    
  } else { #CS, BB and WA
    SAMP_STRAT <- F # do not average hauls by rect
    BYSUBDIV <-   T # average indicator by survey strata from shapefiles
    QUAD<- F #only north sea
    QUAD_SMOOTH<-F
  }
  #exceptions
  if(survey=="GNSFraOT4"){
    SAMP_STRAT <- T # average hauls over surveyed minigrid0.25deg  #if set to FALSE need to update Attibutes table as area only given by minigrid
    BYSUBDIV <-   F # no higher survey strata - could use depth strata from GNSEngBT3 - not implemented
  }
  
  #Trophic Level FW4 check
  if(MEANTLs==T & !(substr(survey,1,2) %in% c("CS","BB", "GN","IB")) ) { MEANTLs<-F; print("no data for TL for WA surveys, deleting MTL selection") }
  
  
  ##load data
  #sampling data 
  #flex <- read.csv(FLEXFILE) 
  #samp0 <- read.table("C:/Users/cl06/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/HL/HH-NS-IBTS.csv" ,as.is = c(1,2,4,5,6,10,11,23),header = TRUE,sep=",") 
  #samp0<-samp0[samp0$Quarter==QUARTER,]
  #samp0<-samp0[samp0$HaulVal=="V",] ###
  #samp0<-samp0[samp0$StdSpecRecCode==1,]
  #samp0$HaulID <- paste(paste(survey,QUARTER,sep=""),samp0$Country,#samp$Gear,samp0$Ship, samp0$StNo, samp0$HaulNo, samp0$Year,sep=":")
  #samp0 matches samp exactly! :)
  
  samp <- read.table(SAMP_FILE ,as.is = c(1,2,4,5,6,10,11,23),header = TRUE,sep=",") 
  
  samp<-samp[samp$Quarter==QUARTER,]
  if(survey %in% c("GNSEngBT1","GNSEngBT2","GNSEngBT3","GNSEngBT4",
                   "CSEngBT1","CSEngBT2","CSEngBT3","CSEngBT4",
                   "GNSBelBT1","GNSBelBT2","GNSBelBT3","GNSBelBT4",
                   "GNSGerBT1","GNSGerBT2","GNSGerBT3","GNSGerBT4",
                   "GNSNetBT1","GNSNetBT2","GNSNetBT3","GNSNetBT4")) samp<-samp[samp$Country==COUNTRY,]
  #need split SEA DATA FOR BTS Q3 ENG' 
  if(survey %in% "GNSEngBT3") samp<-samp[samp$ShootLong_degdec>= -2 & samp$ShootLat_degdec>= 49.5,]
  if(survey %in% "CSEngBT3") samp<-samp[samp$ShootLong_degdec< -3  & samp$ShootLat_degdec>= 50.5 & samp$ShootLat_degdec< 56,]
  # plot(samp$ShootLong_degdec, samp$ShootLat_degdec); map(add=T)#,xlim=c(4,14)
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

  #samp<-  merge(samp,flex,by="HaulID", all.x=F,all.y=F)
  #with(samp[samp$ShootLong_degdec>7,], xyplot(ShootLat_degdec~ShootLong_degdec | ac(YearShot) ))
  with(samp, xyplot(ShootLat_degdec~ShootLong_degdec | ac(YearShot) ))
  # biological data
  bio <- read.csv(BIOL_FILE,as.is = c(1,2,4,5,6,10,11) )  #avoid conversions of rect to e+07 etc #,as.is=1 ) #
  
  #correction needed!
  samp$Ship[samp$Ship=="7.40E+10"]<- "74E9" #correction! excel error pre-upload!
  #if(nrow(bio[nchar(bio$StNo)==6 & substr(bio$StNo,1,3)=="000",]) >0) bio[nchar(bio$StNo)==6 & substr(bio$StNo,1,3)=="000",]$StNo <- substr(bio[nchar(bio$StNo)==6 & substr(bio$StNo,1,3)=="000",]$StNo,4,6) #correction!
  #if(nrow(samp[nchar(samp$StNo)==6 & substr(samp$StNo,1,3)=="000",]) >0) samp[nchar(samp$StNo)==6 & substr(samp$StNo,1,3)=="000",]$StNo <- substr(samp[nchar(samp$StNo)==6 & substr(samp$StNo,1,3)=="000",]$StNo,4,6) #correction!
  #correction!
  #sometimes missing zero in front. add one in everytime just in case
  #if(nrow(bio[nchar(bio$StNo)==1,]) >0)    bio[nchar(bio$StNo)==1, ]$StNo <-  paste("0",bio[nchar(bio$StNo)==1,]$StNo,sep="") #correction!
  #if(nrow(samp[nchar(samp$StNo)==1,]) >0) samp[nchar(samp$StNo)==1,]$StNo <- paste("0",samp[nchar(samp$StNo)==1,]$StNo,sep="")#correction!
  #sometimes missing zero in front. add one in everytime just in case
  #if(nrow(bio[nchar(bio$HaulNo)==1,]) >0)    bio[nchar(bio$HaulNo)==1, ]$HaulNo <-  paste("0",bio[nchar(bio$HaulNo)==1,]$HaulNo,sep="") #correction!
  #if(nrow(samp[nchar(samp$HaulNo)==1,]) >0) samp[nchar(samp$HaulNo)==1,]$HaulNo <- paste("0",samp[nchar(samp$HaulNo)==1,]$HaulNo,sep="") #correction!
  #corrections complete!
  #alternatively remove all leading zeroes
  if(length(bio$HaulNo[substr(bio$HaulNo,1,1)%in%"0"])>0) bio$HaulNo <- str_replace(bio$HaulNo, "^0+" ,"") 
  if(length(bio$StNo[substr(bio$StNo,1,1)%in%"0"])>0) bio$StNo <- str_replace(bio$StNo, "^0+" ,"") 
  if(length(samp$HaulNo[substr(samp$HaulNo,1,1)%in%"0"])>0) samp$HaulNo <- str_replace(samp$HaulNo, "^0+" ,"") 
  if(length(samp$StNo[substr(samp$StNo,1,1)%in%"0"])>0) samp$StNo <- str_replace(samp$StNo, "^0+" ,"") 
  
  bio <- bio[bio$Quarter==QUARTER,]
  bio <- bio[!is.na(bio$ValidAphiaID),] #remove invalids
  bio <- bio[(bio$SpecVal %in% c(1,4,7,10) ),]#remove invalids
  #summary( bio$HaulNo)
  unique( bio$StNo)
  bio$HaulID	<- paste(paste(survey,QUARTER,sep=""), bio$Country, #bio$Gear, 
                      bio$Ship,  bio$StNo, bio$HaulNo, bio$Year,sep=":")
  
  samp$HaulID <- paste(paste(survey,QUARTER,sep=""),samp$Country,#samp$Gear,
                       samp$Ship, samp$StNo, samp$HaulNo, samp$YearShot,sep=":" )#HaulNo and StNo ,samp$Month,samp$Day
  # KLASS/RUTH use:  Survey,Country,Quarter,Gear,Ship,StNo,HaulNo,Year,Month,Day
  
  
    #sort(unique(samp$Gear)) # "GOV" "H18"
    #sort(unique(bio$Gear))  # "ABD" "DHT" "GOV" "GRT" "H18" "HOB"
    #bio$Gear[bio$Gear=="ABD"] <- "GOV"
  
    #library(openxlsx)
    #gears <- read.xlsx("Fishing_gears_input_RKynoch_FBurns.xlsx")
    #gear_cats <- gears[,c(1,2,6,7)]
    #colnames(gear_cats) <- c("Survey","Gear","GearCat","Comment")
  
    #bio <- bio[!bio$Gear%in% c("HOB","H18"),]
    #samp <-samp[!samp$Gear%in% c("HOB","H18"),]
  
  haulidhl <- sort(unique(samp$HaulID))
  haulidhh <- sort(unique(bio$HaulID))
  
  haulidhhm <- subset(samp, !samp$HaulID %in% bio$HaulID)
  haulidhhm <- sort( unique(haulidhhm$HaulID) ); length(haulidhhm) 
  100*length(haulidhhm)/nrow(samp) # 2samp hauls not match in bio
  
  haulidhlm <- subset(bio, !bio$HaulID %in% samp$HaulID)# bio hauls not match in samp
  haulidhlm <- sort(unique(haulidhlm$HaulID) ); length(haulidhlm) 
  100*length(haulidhlm)/nrow(bio) #samp hauls not match in bio
  
  #haulidm <- sort(c(haulidhlm,haulidhhm))
  #haulidm <- haulidm[!duplicated(haulidm)]
  #
  #unique(samp$HaulNo[!samp$HaulNo %in% bio$HaulNo])#all match
  #sort(unique(samp$StNo[!samp$StNo %in% bio$StNo]))
  #sort(unique(bio$StNo[!bio$StNo %in% samp$StNo]))
  
  #samp[!samp$StNo %in% bio$StNo,c("Country","Ship","YearShot","StNo","HaulNo")]
  #view(bio[!(bio$StNo %in% samp$StNo) ,c("Country","Ship","Year","StNo","HaulNo")])
  
   print("process_exchange_HL")
   source('R/Lynam_IND_script_process_exchange_HL2021.R')
   SPPLIST<-read.csv(SPPFILE)
   #SPPLIST<- SPPLIST[SPPLIST$Habit=="DEMERSAL",]
   bio <- bio[bio$SpeciesSciName %in% SPPLIST[,1],]#remove rare spp
   bio$SpeciesSciName<-af(ac(bio$SpeciesSciName))#relevel 131 from >700
   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#quick eyeball of sampling data 
  # quick check   unique(bio[bio$LngtCode==1,"SubFactor"])
  #par(mfrow=c(1,1))
  x11(); par(mfrow=c(4,4))
  plot(samp$ShootLong_degdec, samp$ShootLat_degdec); map(add=T)#,xlim=c(4,14)
  abline(v=c(4,8));  abline(h=55.5)
  with(samp, hist(MonthShot))
  with(samp, hist(HaulDur_min))
  with(samp, hist(Depth_m))
  if(substr(survey,7,8)!="BT" & substr(survey,6,7)!="BT"){ 
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
  bio <- bio[,c("HaulID","SpeciesSciName","FishLength_cm","DensBiom_kg_Sqkm","SubFactor")]
  bio$sciName<-as.character(bio$SpeciesSciName)
  
  
  samp <-samp[,c("HaulDur_min","DataType","HaulID","YearShot","ShootLat_degdec","ShootLong_degdec",
                 "ICESStSq","SurvStratum","WingSwpArea_sqkm","WingSwpVol_CorF","NetOpen_m")] 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  #add efficiency of E=GOV gear
  #use gear efficiency to correct catches
  #the probability that fish in the path of a trawl will be caught and retained
  if(CATCHABILITY_COR_WALKER) print("now CATCHABILITY correct data with WALKER")
  if(CATCHABILITY_COR_WALKER) source(PROC_SCRIPT,"Lynam_IND_script_CATCHABILITY_COR_WALKER.R",sep="/")

  ave_NetOpen_m<-mean(samp$NetOpen_m)
  #samp$WingSwpVol_CorF <- ave_NetOpen_m / samp$NetOpen_m # scale down if larger than usual net opening
  #"SubFactor",  
  dhspp <- merge(bio,samp,by="HaulID")
  with(dhspp, xyplot(ShootLat_degdec~ShootLong_degdec | ac(YearShot) ))
  lostID<-unique(bio[!(bio$HaulID %in% dhspp$HaulID),"HaulID"])#all matched with StnNo
  if(length(lostID) >0){print(paste("losing", length(lostID) ,"hauls from",length(unique(bio$HaulID)),"when merge bio and samp")) } else { print("successful merge HL and HH to create dhspp")}
  write.table(lostID,paste(OUTPATH,"lostID_",survey,"_Q",QUARTER,".txt",sep=""))
  #bioraw$HaulID<- paste(paste(survey,QUARTER,sep=""), bioraw$Ship, bioraw$YearShot, bioraw$HaulNo,sep="/")#
  #bioraw[(bioraw$HaulID%in% lostID),]

    #LOSE RECORDS MISSING LENGTH (already done above) and MISSING DATATYPE
    #dhspp <- dhspp[!(dhspp$FishLength_cm== -9 & !is.na(dhspp$FishLength_cm)),] 
    #dhspp <- dhspp[!is.na(dhspp$DataType),]

    # subfactor to raise measured num at length to total num fish in haul 
    # datatype C is per hour #dhspp$DataType == C DATA ARE CPUE, SO
    # if R or S: HLNoAtLngt*SubFactor   
    #dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensAbund_N_Sqkm),]$DensAbund_N_Sqkm <-
     # dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensAbund_N_Sqkm),]$DensAbund_N_Sqkm *
      #dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensAbund_N_Sqkm),]$SubFactor
    
    dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensBiom_kg_Sqkm),]$DensBiom_kg_Sqkm <-
      dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensBiom_kg_Sqkm),]$DensBiom_kg_Sqkm *
      dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensBiom_kg_Sqkm),]$SubFactor
    
    #standardise per 60 min tow with DurRaise
    dhspp$DurRaise <- dhspp$HaulDur_min/60;  # per hour # hist(dhspp$DurRaise)  
    if(nrow(dhspp[dhspp$DataType=="C" ,])>0) dhspp[dhspp$DataType=="C" ,]$DurRaise <- 1 # already per hour
    #dhspp$DensAbund_N_Sqkm <- dhspp$DensAbund_N_Sqkm/dhspp$DurRaise  #num per 60 min as in 'C'
    dhspp$DensBiom_kg_Sqkm <- dhspp$DensBiom_kg_Sqkm/dhspp$DurRaise  #CPUE per 60 min
    
    ##add swept area from ICES?
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
  
  #now use quadrants## load("C:/Users/cl06/OneDrive - CEFAS/Fish_dataproduct_QSR/biodiv19 ref/RUNuptoOSPARsubdiv.rdata.RData")
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
    if(survey == "GNSIntOT1" & USE_GUILD_COVARIATE_SITES){
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
  if(FINALPLOTS){ print("Final plots"); source(paste(MAINDIR,"R/Lynam_IND_script_FINALPLOTS.R",sep="")) }
  print(paste("Finished",survey, "survey",sep=" "))
  dev.off()
  }
  }
  }#QUARTER # COUNTRY # SEA loops
} # next survey
print("script complete")
  