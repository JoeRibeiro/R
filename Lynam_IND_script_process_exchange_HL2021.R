# Author: C Lynam, Cefas
# Contact: chris.lynam@cefas.co.uk
# Version: 6
# Date: Dec 2021 

# bioraw<-bio # bio<-bioraw


SPPLIST <- LW$ScientificName_WoRMS
#check which species would be removed by SPPLIST
bio$SpeciesSciName <- ac(bio$ScientificName_WoRMS)#edit

if(nrow(bio[bio$SpeciesSciName == "Dipturus lintea",]) >0) bio[bio$SpeciesSciName == "Dipturus lintea",]$SpeciesSciName  <- "Rajella lintea"
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
#trout Salmo trutta trutta
if(nrow(bio[bio$SpeciesSciName == "Salmo trutta",]) >0) bio[bio$SpeciesSciName == "Salmo trutta",]$SpeciesSciName <- "Salmo trutta trutta"
#Synaphobranchus kaupii
if(nrow(bio[bio$SpeciesSciName == "Synaphobranchus kaupi",]) >0) bio[bio$SpeciesSciName == "Synaphobranchus kaupi",]$SpeciesSciName <- "Synaphobranchus kaupii"
#Raja brachyura		include: Bathyraja brachyurops
if(nrow(bio[bio$SpeciesSciName == "Bathyraja brachyurops",]) >0) bio[bio$SpeciesSciName == "Bathyraja brachyurops",]$SpeciesSciName <- "Raja brachyura"

bio$ScientificName_WoRMS <- (bio$SpeciesSciName)


BIONAM<-ac(unique(bio$SpeciesSciName))
sort(BIONAM[!BIONAM %in% SPPLIST])#species excluded 
write.table(BIONAM[!BIONAM %in% SPPLIST],paste(OUTPATH,"lostspp_",survey,"_Q",QUARTER,".txt",sep=""))

#remove -9 entries  #unique(bio[bio$LngtCode== -9,'Number'])#all -9
bio<- bio[bio$HLNoAtLngt!= -9 & bio$HLNoAtLngt!= 0,]
#rename to match MSS names
names(bio)[which(names(bio) == "HLNoAtLngt")] <- "Number";
names(bio)[which(names(bio) == "Year")] <- "YearShot";
names(bio)[which(names(bio) == "LngtClass")] <- "FishLength_cm";
#bio$HaulID	<- paste(paste(survey,QUARTER,sep=""), bio$Ship, bio$YearShot, bio$StNo,sep="/")
bio$SpeciesSciName<-af(ac(bio$SpeciesSciName))  #relevel 131 from >700



## pelagic species in mm and others in cm
# herring, sprat, sandeel, mackerel and horse mackerel in mm  
# to nearest 0.1cm for shellfish, 0.5 cm for herring and sprat, and 1 cm for all other species
## ALL SPECIES data from ENG in mm 
#LngtCode=Length class code for the given species.
#0.1cm for shellfish (reporting units mm)
#0.5 cm for herring and sprat (reporting units mm)
#1 cm for all other species (reporting units cm)

#check LenMeasType is total length

#NS-IBTS    North Sea International Bottom Trawl Survey    .    1 mm   [mm]
#NS-IBTS    North Sea International Bottom Trawl Survey    0    0.5 cm [mm]
#NS-IBTS    North Sea International Bottom Trawl Survey    1    1 cm   [cm]
##NS-IBTS    North Sea International Bottom Trawl Survey - before 2004    -9    Missing Value   #before 2004 only
##NS-IBTS    North Sea International Bottom Trawl Survey - before 2004    5    5 cm             #before 2004 only
##NS-IBTS    North Sea International Bottom Trawl Survey - before 2004    9    + group          #before 2004 only
# check the LngtCode/FishLength_cm

#scottish error SWC-IBTS and NS-IBTS for 126754 "Gymnammodytes semisquamatus" entered as cm when in fact mm
if(nrow(bio[bio$ValidAphiaID==126754 & bio$LngtCode==1 & bio$LngtClass>30,])>0) bio[bio$ValidAphiaID==126754 & bio$LngtCode==1 & bio$LngtClass>30,]$LngtCode <- 0

if( !is.null(bio[bio$LngtCode=='.' | bio$LngtCode=='0','FishLength_cm']) ){ 
  bio[bio$LngtCode=='.' | bio$LngtCode=='0','FishLength_cm'] <- 
    bio[bio$LngtCode=='.' | bio$LngtCode=='0','FishLength_cm']/10; }
## correct to cm unique(bio$LngtCode)#Levels: . 0 1
#& Lophius budegassa, the blackbellied angler,
#CORRECT SPECIES CODES
if(nrow(bio[bio$SpecCode==170992 & bio$SpecCodeType=="t",])>0) bio[bio$SpecCode==170992 & bio$LngtClass<20 & bio$SpecCodeType=="t",'SpecCode'] <- 170991      #'#'#'# lesser weevers 170991 misidentified as greater weeevers 170992
if(nrow(bio[bio$SpecCode==127082 & bio$SpecCodeType=="w",])>0) bio[bio$SpecCode==127082 & bio$LngtClass<20 & bio$SpecCodeType=="t",'SpecCode'] <- 150630      #'#'#'# lesser weevers 170991 misidentified as greater weeevers 170992

if(nrow(bio[bio$SpecCode==170991 & bio$SpecCodeType=="t",])>0) bio[bio$SpecCode==170991 & bio$LngtClass>=20 & bio$SpecCodeType=="t",'SpecCode']<- 170992      #'#'#'# lesser weevers 170991 misidentified as greater weeevers 170992
if(nrow(bio[bio$SpecCode==150630 & bio$SpecCodeType=="w",])>0) bio[bio$SpecCode==150630 & bio$LngtClass<20 & bio$SpecCodeType=="t",'SpecCode'] <- 127082      #'#'#'# lesser weevers 170991 misidentified as greater weeevers 170992
 
PLOTcheck<-T
if(PLOTcheck){
  pdf(file=paste(OUTPATH,survey,"_len_hist.pdf",sep=""))
  par(mfrow=c(4,4))
  for(i in sort(unique(bio$ScientificName_WoRMS)) ){ hist(bio[bio$ScientificName_WoRMS==i,]$FishLength_cm, main=i,xlab="cm")}
  dev.off() 
  }


### length cm to weight g 
#### LW #########################################################################################################################
#separate and keep all records for IEO-MTL-FW4
if(IEO_FW4){
  bio_oth <- bio[!(bio$SpeciesSciName %in% SPPLIST),]
  bio_oth$a <- 1; bio_oth$b <- 1; bio_oth$SensFC1 <- F; bio_oth$DEMPEL <- "Other"; bio_oth$Max.L..cm. <- F; bio_oth$checkLen<-F
  
  #trait_MAXL_other<- data.frame( matrix(NA, ncol=8, nrow=length(unique(bio_oth$SpeciesSciName)) ) )
  #colnames(trait_MAXL_other) <- c("ScientificName_WoRMS","Loo","Lm","Order","Group","MaxL","SpeciesSciName", "Max.L..cm.") 
  #trait_MAXL_other$ScientificName_WoRMS <- unique(bio_oth$SpeciesSciName)
  #trait_MAXL_other$SpeciesSciName <- trait_MAXL_other$ScientificName_WoRMS
  #trait_MAXL_other$Order <- "Other"
  #trait_MAXL_other$Group <- "Other"
  
  bio_oth$DensBiom_kg_Sqkm <- NA #bio_oth$Number*(bio_oth$a*(bio_oth$FishLength_cm^bio_oth$b))/1000
  names(bio_oth)[which(names(bio_oth) == "ScientificName_WoRMS")] <- "SpeciesSciName" 
  #still need to correct DensAbund_N_Sqkm for haul duration or swept area - in main script
  bio_oth$DensAbund_N_Sqkm <- bio_oth$Number
}
#### SPPLIST #########################################################################################################################
#now select the fish/elasmo species only 
bio <- bio[bio$SpeciesSciName %in% SPPLIST,]#remove rare spp keep 550 species
  #sort( unique(bio_oth$SpeciesSciName) )
  #sort( unique(bio$SpeciesSciName) )

LW_COLS<-c("ScientificName_WoRMS","a","b","Max.L..cm.","SensFC1","DEMPEL")
bio <- merge(x=bio,y=LW[,names(LW) %in% LW_COLS],by="ScientificName_WoRMS",all.x=F,all.y=F)
#any still missing?
unique(bio$ScientificName_WoRMS[is.na(bio$a)|is.na(bio$b)] )
bio$a[is.na(bio$a)|is.na(bio$b)] <- 0.00635
bio$b[is.na(bio$a)|is.na(bio$b)] <- 3.116


#check lens rel to max
bio$checkLen<- bio$FishLength_cm/bio$Max.L..cm.
summary(bio$checkLen)
#hist((bio$checkLen))
NAM<-unique(bio$ScientificName_WoRMS[bio$checkLen>1]) # The lesser weever (Echiichthys vipera) corrected
NAM<-ac(NAM[!is.na(NAM)])
if(PLOTcheck){
  pdf(file=paste(OUTPATH,survey,"_maxlen.pdf",sep=""))
  par(mfrow=c(3,3));
  for(i in NAM){
    #if(which(NAM==i) == c(1,10,19,28) ){ x11(); par(mfrow=c(3,3));}
    hist(bio[bio$ScientificName_WoRMS==i,"FishLength_cm"],main = paste(i, 
                                                                       unique(bio[bio$ScientificName_WoRMS==i,"Max.L..cm."])
    ))
    abline(v=unique(bio[bio$ScientificName_WoRMS==i,"Max.L..cm."]),col=2,lwd=2) 
  }
  dev.off()
}#ok
#calc bio at len [in g so /1000 -> kg] from LW
bio$DensBiom_kg_Sqkm <- bio$Number*(bio$a*(bio$FishLength_cm^bio$b))/1000
names(bio)[which(names(bio) == "ScientificName_WoRMS")] <- "SpeciesSciName" 
#still need to correct DensAbund_N_Sqkm for haul duration or swept area - in main script
bio$DensAbund_N_Sqkm <- bio$Number



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#quick eyeball of sampling and bio data 
# quick check   unique(bio[bio$LngtCode==1,"SubFactor"])
#par(mfrow=c(1,1))
x11(); 
if(substr(survey,7,8)!="BT" & substr(survey,6,7)!="BT" & substr(survey,8,9)!="BT" & substr(survey,7,8)!="Bi"){ 
  par(mfrow=c(4,4)) } else { par(mfrow=c(2,4)) }
plot(samp$ShootLong_degdec, samp$ShootLat_degdec); map(add=T)#,xlim=c(4,14)
abline(v=c(4,8));  abline(h=55.5)
with(samp, hist(MonthShot))
with(samp, hist(HaulDur_min))
with(samp, hist(Depth_m))
if(substr(survey,7,8)!="BT" & substr(survey,6,7)!="BT" & substr(survey,8,9)!="BT" & substr(survey,7,8)!="Bi"){ 
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
if(IEO_FW4) bio <- rbind(bio,bio_oth)

bio <- bio[,c("HaulID","SpeciesSciName","FishLength_cm","DensBiom_kg_Sqkm","DensAbund_N_Sqkm","SubFactor","SensFC1","DEMPEL")]
bio$sciName<-as.character(bio$SpeciesSciName)


samp <-samp[,c("HaulDur_min","DataType","HaulID","YearShot","ShootLat_degdec","ShootLong_degdec",
               "ICESStSq","SurvStratum","WingSwpArea_sqkm","WingSwpVol_CorF","NetOpen_m"
               ,"Gear","Ship","MonthShot","Day","TimeShot")] 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  #add efficiency of E=GOV gear
#use gear efficiency to correct catches
#the probability that fish in the path of a trawl will be caught and retained
if(CATCHABILITY_COR_WALKER){
  print("now CATCHABILITY correct data with WALKER")
  source(paste(PROC_SCRIPT,"Lynam_IND_script_CATCHABILITY_COR_WALKER.R",sep="/"))
} 
ave_NetOpen_m<-mean(samp$NetOpen_m)
#samp$WingSwpVol_CorF <- ave_NetOpen_m / samp$NetOpen_m # scale down if larger than usual net openinge
#"SubFactor",  
dhspp <- merge(bio,samp,by="HaulID",all = FALSE)
# plot(dhspp$ShootLong_degdec, dhspp$ShootLat_degdec); map(add=T)#,xlim=c(4,14)
#with(dhspp, xyplot(ShootLat_degdec~ShootLong_degdec | ac(YearShot) ))
# plot(samp[!(samp$HaulID %in% dhspp$HaulID),]$ShootLong_degdec,samp[!(samp$HaulID %in% dhspp$HaulID),]$ShootLat_degdec)
lostID_from_lack_of_matching_haulID<-unique(bio[!(bio$HaulID %in% dhspp$HaulID),"HaulID"])#all bio matched with StnNo
# JR edit - dhspp contains NA lats and longs
dhspp = dhspp[is.finite(dhspp$ShootLat_degdec) & is.finite(dhspp$ShootLong_degdec),]
# head(samp[!(samp$HaulID %in% dhspp$HaulID),])
# view(bio[!(bio$HaulID %in% dhspp$HaulID),])
if(length(lostID_from_lack_of_matching_haulID) >0){
  print(paste("losing", length(lostID_from_lack_of_matching_haulID) ,"hauls from",length(unique(bio$HaulID)),"when merge bio and samp")) 
} else { 
  print("successful merge HL and HH to create dhspp")
}
write.table(lostID_from_lack_of_matching_haulID,paste(OUTPATH,"lostID_from_lack_of_matching_haulID_",survey,"_Q",QUARTER,".txt",sep=""))
rm(bio,samp)
if(IEO_FW4) rm(bio_oth)

# Make spatial, filter by SSA, then drop spatial aspect
predhspp <- dhspp #pre SSA
#drop hauls not in SSA
if(SSA_WRITE_NEW) dhspp = dhspp[dhspp$ICESStSq %in% SSAdfs$rectangle,] #cadiz should be 02E2 03E3 etc  as in SSAdfs$rectangle
coordinates(dhspp) <- ~ ShootLong_degdec + ShootLat_degdec
suppressWarnings(proj4string(dhspp) <- CRS("+init=epsg:4326"))
if(!SSA_WRITE_NEW) dhspp = dhspp[surveySSA,] #drop hauls not in SSA

dhspp@data$ShootLong_degdec = dhspp$ShootLong_degdec
dhspp@data$ShootLat_degdec = dhspp$ShootLat_degdec
dhspp=dhspp@data

lostID_from_SSA_filtering<-unique(predhspp[!(predhspp$HaulID %in% dhspp$HaulID),"HaulID"])#all predhspp matched with StnNo
if(length(lostID_from_SSA_filtering) >0){print(paste("losing", length(lostID_from_SSA_filtering) ,"hauls from",length(unique(predhspp$HaulID)),"when merge predhspp and samp i.e. from SSA")) } else { print("no hauls lost due to SSA")}
write.table(lostID_from_SSA_filtering,paste(OUTPATH,"lostID_from_SSA_filtering_",survey,"_Q",QUARTER,".txt",sep=""))
rm(predhspp)

# subfactor to raise measured num at length to total num fish in haul 
# datatype C is per hour #dhspp$DataType == C DATA ARE CPUE, SO
# if R or S: HLNoAtLngt*SubFactor   
# also done in Meadhbh script no.7
# note the Total number (not used here) is already Raised no need to multiply by the Subfactor 
dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensAbund_N_Sqkm),]$DensAbund_N_Sqkm <-
  dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensAbund_N_Sqkm),]$DensAbund_N_Sqkm *
  dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensAbund_N_Sqkm),]$SubFactor

dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensBiom_kg_Sqkm),]$DensBiom_kg_Sqkm <-
  dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensBiom_kg_Sqkm),]$DensBiom_kg_Sqkm *
  dhspp[dhspp$DataType!="C" & dhspp$SubFactor!=1 & !is.na(dhspp$DensBiom_kg_Sqkm),]$SubFactor

#### standardise per 60 min tow with DurRaise #### 
dhspp$DurRaise <- dhspp$HaulDur_min/60;  # per hour # hist(dhspp$DurRaise)  
if(nrow(dhspp[dhspp$DataType=="C" ,])>0) dhspp[dhspp$DataType=="C" ,]$DurRaise <- 1 # already per hour
dhspp$DensAbund_N_perhr <- dhspp$DensAbund_N_Sqkm/dhspp$DurRaise  #num per 60 min as in 'C'
dhspp$DensBiom_kg_perhr <- dhspp$DensBiom_kg_Sqkm/dhspp$DurRaise  #CPUE per 60 min

#### use swept area from ICES to give per sq km #### 
dhspp$DurRaise <- 1
if(nrow(dhspp[dhspp$DataType=="C" ,])>0) dhspp[dhspp$DataType=="C" ,]$DurRaise <- dhspp[dhspp$DataType=="C" ,]$HaulDur_min/60 # do not want per hour #also done in Meadhbh scripts no.7
dhspp$DensAbund_N_Sqkm<-( dhspp$DensAbund_N_Sqkm*dhspp$DurRaise) / dhspp$WingSwpArea_sqkm
dhspp$DensBiom_kg_Sqkm<-( dhspp$DensBiom_kg_Sqkm*dhspp$DurRaise) / dhspp$WingSwpArea_sqkm
#dhspp$DurRaise for C data

#rm(samp,bio)
dhspp <- merge(x=dhspp, y=SciName2GroupEff[,c(1,9)], by.x="SpeciesSciName", by.y="sciName",all.x=T)
if(SPECIES_IN_MOD_ONLY) dhspp <- dhspp[dhspp$SpeciesSciName %in% SPECIES_SUBSET,] #use only those from thorpe model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#simplify names
names(dhspp)[which(names(dhspp) =="YearShot")] <-"Year"
#for tyl
dhspp$LogLngtClass <- log(dhspp$FishLength_cm)
dhspp$LogLngtBio <- dhspp$LogLngtClass * dhspp$DensBiom_kg_Sqkm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#MaxLength and demersal or pelagic
#MaxLength observed for MML
names(dhspp)[names(dhspp)=="Group"] <- "QGrouping"
#dhspp <- merge(dhspp,trait_MAXL[,c("SpeciesSciName","maximum.length","Order","Group", "LFI_Fung_list", "Loo","K","Winfinity","tmax","tm","M","Lm","LFI_OSPAR_list")],by="SpeciesSciName",all.x=T)#"DEMPEL",
#names(dhspp)[which( names(dhspp)=="maximum.length")] <- "MaxL"
dhspp <- merge(dhspp,trait_MAXL[,c("ScientificName_WoRMS","SpeciesSciName","Loo","Lm","MaxL","Max.L..cm.","Order","Group")],by="SpeciesSciName",all.x=T)#"DEMPEL",
#if(IEO_FW4){   dhspp <- merge(dhspp,trait_MAXL_other,by="SpeciesSciName",all.x=T) }

#DEM are ecotypes Demersal + Bathydemersal + Bathypelagic  (except Micromesistius poutassou Blue whiting) + Benthopelagic (except Clupea harengus herring)
#PEL are Pelagic  #sprat, mac, hmx, her, whb etc
dhspp$DEMPEL <-as.character(dhspp$DEMPEL)
dhspp$DEMPEL[dhspp$DEMPEL=="Demersal"] <- "DEM"
dhspp$DEMPEL[dhspp$DEMPEL=="Pelagic"] <- "PEL"
dhspp$DEMPEL <-as.character(dhspp$DEMPEL)

#notes in ICES WKABSENS:
#Alosa spp: remove data if Longitude < 7 
# as Alosa spp in the Baltic/Kattegat/Skagerrak were spatially separated from
# more western Alosa spp and there were few western Alosa spp caught 
# as judged from the spatial distribution of survey catches.
#if(nrow(dhspp[dhspp$SpeciesSciName == "Alosa fallax" & dhspp$ShootLong_degdec< 7,]) >0) dhspp<- dhspp[!(dhspp$SpeciesSciName == "Alosa fallax" & dhspp$ShootLong_degdec< 7),]
#not implemented above as better to keep records and remove from FC1 analysis if too few

#Amblyraja radiata: correct data if 
# Longitude < -6, 
# Longitude < -4 and Latitude < 60, 
# Longitude < -3 and Latitude < 55 
# as WGEF considered that these specimens were misidentified as Raja clavata.
if(nrow(dhspp[dhspp$SpeciesSciName == "Amblyraja radiata" & dhspp$ShootLong_degdec< -6,]) >0) dhspp[dhspp$SpeciesSciName == "Amblyraja radiata" & dhspp$ShootLong< -6,]$SpeciesSciName <- "Raja clavata"
if(nrow(dhspp[dhspp$SpeciesSciName == "Amblyraja radiata" & dhspp$ShootLong_degdec< -4 & dhspp$ShootLat_degdec< 60,]) >0) dhspp[dhspp$SpeciesSciName == "Amblyraja radiata" & dhspp$ShootLong< -4 & dhspp$ShootLat_degdec< 60,]$SpeciesSciName <- "Raja clavata"
if(nrow(dhspp[dhspp$SpeciesSciName == "Amblyraja radiata" & dhspp$ShootLong_degdec< -3 & dhspp$ShootLat_degdec< 55,]) >0) dhspp[dhspp$SpeciesSciName == "Amblyraja radiata" & dhspp$ShootLong< -3 & dhspp$ShootLat_degdec< 55,]$SpeciesSciName <- "Raja clavata"
#
#ices wkabsens
#species groups needed		
#Alosa		Alosa alosa and Alosa fallax
if(nrow(dhspp[dhspp$SpeciesSciName == "Alosa alosa",]) >0) dhspp[dhspp$SpeciesSciName == "Alosa alosa",]$SpeciesSciName <- "Alosa"
if(nrow(dhspp[dhspp$SpeciesSciName == "Alosa fallax",]) >0) dhspp[dhspp$SpeciesSciName == "Alosa fallax",]$SpeciesSciName <- "Alosa"
#Galeus		Galeus melastomus
  if(nrow(dhspp[dhspp$SpeciesSciName == "Galeus atlanticus",]) >0) dhspp[dhspp$SpeciesSciName == "Galeus atlanticus",]$SpeciesSciName <- "Galeus"
if(nrow(dhspp[dhspp$SpeciesSciName == "Galeus melastomus",]) >0) dhspp[dhspp$SpeciesSciName == "Galeus melastomus",]$SpeciesSciName <- "Galeus"
#Hippocampus spp.		Hippocampus hippocampus with Hippocampus guttulatus
if(nrow(dhspp[dhspp$SpeciesSciName == "Hippocampus hippocampus",]) >0) dhspp[dhspp$SpeciesSciName == "Hippocampus hippocampus",]$SpeciesSciName <- "Hippocampus"
if(nrow(dhspp[dhspp$SpeciesSciName == "Hippocampus guttulatus",]) >0) dhspp[dhspp$SpeciesSciName == "Hippocampus guttulatus",]$SpeciesSciName <- "Hippocampus"
#Mustelus		Mustelus spp. and M. mustelus and M. asterias.
if(nrow(dhspp[dhspp$SpeciesSciName == "Mustelus mustelus",]) >0) dhspp[dhspp$SpeciesSciName == "Mustelus mustelus",]$SpeciesSciName <- "Mustelus"
if(nrow(dhspp[dhspp$SpeciesSciName == "Mustelus asterias",]) >0) dhspp[dhspp$SpeciesSciName == "Mustelus asterias",]$SpeciesSciName <- "Mustelus"
#Sebastes		S. marinus, S. mentella, S. norvegicus
if(nrow(dhspp[dhspp$SpeciesSciName == "Sebastes marinus",]) >0) dhspp[dhspp$SpeciesSciName == "Sebastes marinus",]$SpeciesSciName <- "Sebastes" #Sebastes marinus misspelling of Sebastes norvegicus
if(nrow(dhspp[dhspp$SpeciesSciName == "Sebastes mentella",]) >0) dhspp[dhspp$SpeciesSciName == "Sebastes mentella",]$SpeciesSciName <- "Sebastes"
if(nrow(dhspp[dhspp$SpeciesSciName == "Sebastes norvegicus",]) >0) dhspp[dhspp$SpeciesSciName == "Sebastes norvegicus",]$SpeciesSciName <- "Sebastes"
#Coregonus spp.	Coregonus oxyrinchus	Coregonus maraena# only 8 records in total
if(nrow(dhspp[dhspp$SpeciesSciName == "Coregonus maraena",]) >0) dhspp[dhspp$SpeciesSciName == "Coregonus maraena",]$SpeciesSciName <- "Coregonus"
if(nrow(dhspp[dhspp$SpeciesSciName == "Coregonus oxyrinchus",]) >0) dhspp[dhspp$SpeciesSciName == "Coregonus oxyrinchus",]$SpeciesSciName <- "Coregonus"
#Deania calcea		
#includes from D. profundorum (105905) prior to 2002 (Portuguese surveys) #QSR not going back this far
#/2012 (Spanish surveys) #problem here though
if(nrow(dhspp[dhspp$SpeciesSciName == "Deania profundorum",]) >0) dhspp[dhspp$SpeciesSciName == "Deania profundorum",]$SpeciesSciName <- "Deania calcea"
if(nrow(dhspp[dhspp$SpeciesSciName == "Deania profundorum",]) >0) dhspp[dhspp$SpeciesSciName == "Deania profundorum",]$SpeciesSciName <- "Deania calcea"
#Dipturus batis complex		Uncertain ID, combined Dipturus, D. batis, D. flossada and D. intermedia in abundance indices.
if(nrow(dhspp[dhspp$SpeciesSciName == "Dipturus batis",]) >0) dhspp[dhspp$SpeciesSciName == "Dipturus batis",]$SpeciesSciName <- "Dipturus"
if(nrow(dhspp[dhspp$SpeciesSciName == "Dipturus flossada",]) >0) dhspp[dhspp$SpeciesSciName == "Dipturus flossada",]$SpeciesSciName <- "Dipturus"
if(nrow(dhspp[dhspp$SpeciesSciName == "Dipturus intermedia",]) >0) dhspp[dhspp$SpeciesSciName == "Dipturus intermedia",]$SpeciesSciName <- "Dipturus"
if(nrow(dhspp[dhspp$SpeciesSciName == "Dipturus spp",]) >0) dhspp[dhspp$SpeciesSciName == "Dipturus spp",]$SpeciesSciName <- "Dipturus"
#mediterranean species found in the channel in GNSFraOT4 - mis ID  between these two species so group together:
if(survey=="GNSFraOT4" & nrow(dhspp[dhspp$SpeciesSciName == "Dasyatis tortonesei",]) >0) dhspp[dhspp$SpeciesSciName == "Dasyatis tortonesei",]$SpeciesSciName <- "Dasyatis pastinaca"



#check lens rel to max for those with Valid sampling data (excluding A and I)
dhspp$checkLen<- dhspp$FishLength_cm/dhspp$Max.L..cm.
summary(dhspp$checkLen)
#hist((dhspp$checkLen))
NAM<-unique(dhspp$sciName[dhspp$checkLen>1]) # The lesser weever (Echiichthys vipera) corrected
NAM<-ac(NAM[!is.na(NAM)])
if(PLOTcheck){
  pdf(file=paste(OUTPATH,survey,"_maxlen_Validsamp.pdf",sep=""))
  par(mfrow=c(3,3));
  for(i in NAM){
    #if(which(NAM==i) == c(1,10,19,28) ){ x11(); par(mfrow=c(3,3));}
    hist(dhspp[dhspp$sciName==i,"FishLength_cm"],main = paste(i, 
                                                                       unique(dhspp[dhspp$sciName==i,"Max.L..cm."])
    ))
    abline(v=unique(dhspp[dhspp$sciName==i,"Max.L..cm."]),col=2,lwd=2) 
  }
  dev.off()
}#ok
