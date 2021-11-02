# Author: C Lynam, Cefas
# Contact: chris.lynam@cefas.co.uk
# Version: 3
# Date: Oct 2021 

bioraw<-bio # bio<-bioraw
bio<-bio[bio$Quarter==QUARTER,]
#remove -9 entries  #unique(bio[bio$LngtCode== -9,'Number'])#all -9
bio<- bio[bio$HLNoAtLngt!= -9 & bio$HLNoAtLngt!= 0,]
#rename to match MSS names
names(bio)[which(names(bio) == "HLNoAtLngt")] <- "Number";
names(bio)[which(names(bio) == "Year")] <- "YearShot";
names(bio)[which(names(bio) == "LngtClass")] <- "FishLength_cm";
#bio$HaulID	<- paste(paste(survey,QUARTER,sep=""), bio$Ship, bio$YearShot, bio$StNo,sep="/")


## pelagic species in mm and others in cm
# herring, sprat, sandeel, mackerel and horse mackerel in mm  
# to nearest 0.1cm for shellfish, 0.5 cm for herring and sprat, and 1 cm for all other species
## ALL SPECIES data from ENG in mm 
#LngtCode=Length class code for the given species.
#0.1cm for shellfish (reporting units mm)
#0.5 cm for herring and sprat (reporting units mm)
#1 cm for all other species (reporting units cm)

#NS-IBTS    North Sea International Bottom Trawl Survey    .    1 mm   [mm]
#NS-IBTS    North Sea International Bottom Trawl Survey    0    0.5 cm [mm]
#NS-IBTS    North Sea International Bottom Trawl Survey    1    1 cm   [cm]
##NS-IBTS    North Sea International Bottom Trawl Survey - before 2004    -9    Missing Value   #before 2004 only
##NS-IBTS    North Sea International Bottom Trawl Survey - before 2004    5    5 cm             #before 2004 only
##NS-IBTS    North Sea International Bottom Trawl Survey - before 2004    9    + group          #before 2004 only
# check the LngtCode/FishLength_cm
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

PLOTcheck<-F
if(PLOTcheck){
  pdf(file=paste(OUTPATH,survey,QUARTER,"_len_hist.pdf",sep=""))
  par(mfrow=c(5,5))
  for(i in sort(unique(bio$ScientificName_WoRMS)) ){ hist(bio[bio$ScientificName_WoRMS==i,]$FishLength_cm, main=i,xlab="cm")}
  dev.off() 
  }

### length cm to weight g 
#LW<-read.csv(LWFILE)
#names(LW)[1]<-"ScientificName_WoRMS"
LW_COLS<-c("ScientificName_WoRMS","a","b","Max.L..cm.")
#LW[LW$ScientificName_WoRMS=="Psetta maxima","a"]<-LW[LW$ScientificName_WoRMS=="Scophthalmus maximus","a"]
#LW[LW$ScientificName_WoRMS=="Psetta maxima","b"]<-LW[LW$ScientificName_WoRMS=="Scophthalmus maximus","b"]

length( ISECT<-intersect(LW$ScientificName_WoRMS,bio$ScientificName_WoRMS) )
sort(unique(bio[bio$ScientificName_WoRMS %in% ISECT,]$ScientificName_WoRMS))
sort(unique(bio[!bio$ScientificName_WoRMS %in% ISECT,]$ScientificName_WoRMS))


bio <- merge(x=bio,y=LW[,names(LW) %in% LW_COLS],by="ScientificName_WoRMS",all.x=F,all.y=F)
#any still missing?
unique(bio$ScientificName_WoRMS[is.na(bio$a)|is.na(bio$b)] )
bio$a[is.na(bio$a)|is.na(bio$b)] <- 0.00635
bio$b[is.na(bio$a)|is.na(bio$b)] <- 3.116

#check lens rel to max
bio$checkLen<- bio$FishLength_cm/bio$Max.L..cm.
summary(bio$checkLen); hist((bio$checkLen))
NAM<-unique(bio$ScientificName_WoRMS[bio$checkLen>1]) # The lesser weever (Echiichthys vipera) corrected
NAM<-ac(NAM[!is.na(NAM)])
if(PLOTcheck){
  pdf(file=paste(OUTPATH,survey,QUARTER,"_maxlen_ibts.pdf",sep=""))
  par(mfrow=c(3,3));
  for(i in NAM){
    #if(which(NAM==i) == c(1,10,19,28) ){ x11(); par(mfrow=c(3,3));}
    hist(bio[bio$ScientificName_WoRMS==i,"FishLength_cm"],main = paste(i, 
                                                                       unique(bio[bio$ScientificName_WoRMS==i,"Max.L..cm."])
    ))
    abline(v=unique(bio[bio$ScientificName_WoRMS==i,"Max.L..cm."])) 
  }
  dev.off()
}#ok
#calc bio at len [kg]
bio$DensBiom_kg_Sqkm <- bio$Number*(bio$a*(bio$FishLength_cm^bio$b))/1000
names(bio)[which(names(bio) == "ScientificName_WoRMS")] <- "SpeciesSciName" 
#still need to correct DensAbund_N_Sqkm for haul duration or swept area
bio$DensAbund_N_Sqkm <- bio$Number



