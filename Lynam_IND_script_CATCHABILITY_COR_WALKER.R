#use gear efficiency to correct catches
#the probability that fish in the path of a trawl will be caught and retained

# add catchability by length group
Q <- read.csv("//lowfilecds/Function/Eco Indicators/DATRASoutput/Walker_GearEfficiency_ICESJMarSci17_Supp/EfficiencyTab.csv")
Q$Code <- ac(Q$Code)
# add species names
Qsp <- read.csv("//lowfilecds/Function/Eco Indicators/DATRASoutput/Walker_GearEfficiency_ICESJMarSci17_Supp/specieslist.csv");
names(Qsp)[1] <- "SciName"
Qsp$SpCode <- ac(Qsp$SpCode)

Qsp$SciName <- ac(Qsp$SciName)#spaces at the end!
Qsp$SciName <- substr(Qsp$SciName, start=1,stop=(nchar(Qsp$SciName))-1)

Qsp$Common <- ac(Qsp$Common)#spaces at the end!
Qsp$Common <- substr(Qsp$Common, start=1,stop=(nchar(Qsp$Common))-1)

# merge species name in Q file
Qnam<- merge(x=Q, y=Qsp, by.x=c("Species","Code"), by.y=c("Common","SpCode"), all.x=T)
##Qnam[Qnam$Species=="Sandeel","SciName"]<-"Ammodytidae"#only to 26cm
#rm(Qsp)

#biopreQ<-bio;
#if(CATCHABILITY_COR_WALKER){ #bio<-bioraw; 
QnamG <- Qnam[Qnam$Gear==GEAR,]
QnamG[!is.na(QnamG$Efficiency) & QnamG$Efficiency<0.01,"Efficiency"] <- 0.01
QnamG[,'Group'] <- ac(QnamG[,'Group']);  QnamG[,'Species'] <- ac(QnamG[,'Species'])
QnamG[,'Gear'] <- ac(QnamG[,'Gear']);    QnamG[,'Ref'] <- ac(QnamG[,'Ref'])   
bio <- merge(bio, QnamG,  by.x=c("SpeciesSciName","FishLength_cm"), by.y=c("SciName", "Length"), all.x=T, all.y=F)
bio$mult <- 1/bio$Efficiency
# missing for###  
u(bio[is.infinite(bio$mult),'SpeciesSciName'])
INFS<- (bio[is.infinite(bio$mult),c('SpeciesSciName',"FishLength_cm")])
u(bio[is.na(bio$mult),'SpeciesSciName'])
NAS<- (bio[is.na(bio$mult),c('SpeciesSciName',"FishLength_cm")])
#lose factors #bio[,"Survey_Acronym"] <- ac(bio[,"Survey_Acronym"])
bio[,"SpeciesSciName"] <- ac(bio[,"SpeciesSciName"])
bio[,"HaulID"] <- ac(bio[,"HaulID"])


if( nrow(NAS) > 0 | nrow(INFS) > 0){
  #for lengths not in the walker study use closest values
  for(NA_SP in u(bio[is.na(bio$mult),'SpeciesSciName'])){
    #are there other lengths for this species? # NA_SP<-u(bio[is.na(bio$mult),'SpeciesSciName'])[2]
    OTH_LEN <- QnamG[!is.na(QnamG$SciName) & !is.na(QnamG$Efficiency) & QnamG$SciName == NA_SP,]
    names(OTH_LEN)[which(names(OTH_LEN)=="SciName")] <- "SpeciesSciName"
    if(nrow(OTH_LEN)>0){ 
      print(NA_SP)
      NA_LENS <- u(bio[is.na(bio$mult) & bio$SpeciesSciName == NA_SP, ]$FishLength_cm)
      for(NA_LEN in NA_LENS){
        #min abs difference # NA_LEN<-NA_LENS[1]
        USE_LEN <- OTH_LEN$Length[which(abs((OTH_LEN$Length - NA_LEN)) == min(abs((OTH_LEN$Length - NA_LEN))))]
        #  
        #replace missing eff for some lens
        ROWN <-an(row.names(bio[is.na(bio$mult) & bio$FishLength_cm == NA_LEN & bio$SpeciesSciName == NA_SP,]) )#which row?
        bio[ROWN, which(names(bio) %in% names(OTH_LEN))] <- 
          OTH_LEN[ OTH_LEN$Length==USE_LEN,names(bio)[which(names(bio) %in% names(OTH_LEN))] ]
        
        bio[ROWN,"mult"] <- 1/ bio[ROWN,"Efficiency"]
      }
    }
  } 
  #else for species not in Walker study use my Groupings SciName2GroupEff
  bio.na.Eff <- bio[is.na(bio$Efficiency) | is.infinite(bio$Efficiency),c(1:4)]
  bio.na.Eff <- merge(x=bio.na.Eff, y=SciName2GroupEff[,c(1,7)], by.x="SpeciesSciName", by.y="sciName",all.x=T)
  #add effic
  
  bio.na.Eff <- merge(bio.na.Eff, QnamG[,-which(names(QnamG)=="Group")],  
                      by.x=c("FishLength_cm","Group"),by.y=c("Length","Species"), all.x=T, all.y=F)
  names(bio.na.Eff)[which(names(bio.na.Eff)=="SciName")] <- "Species"
  bio.na.Eff$mult <- 1/bio.na.Eff$Efficiency
  for(NA_SP in u(bio.na.Eff[,'Group']) ){
    #are there other lengths for this GRP? # NA_SP <- u(bio.na.Eff[,'Group'])[1]
    OTH_LEN <- QnamG[!is.na(QnamG$Efficiency) & QnamG$Species == NA_SP,]
    OTH_LEN$Species <- ac(OTH_LEN$Species );  OTH_LEN$Gear <- ac(OTH_LEN$Gear )
    OTH_LEN$Group <- ac(OTH_LEN$Group );      OTH_LEN$Ref <- ac(OTH_LEN$Ref )
    if(nrow(OTH_LEN)>0){ 
      print(NA_SP)
      NA_LENS <- u(bio.na.Eff[bio.na.Eff$Group == NA_SP, ]$FishLength_cm)
      for(NA_LEN in NA_LENS){
        #min abs difference NA_LEN<-NA_LENS[1]
        USE_LEN <- OTH_LEN$Length[which(abs((OTH_LEN$Length - NA_LEN)) == min(abs((OTH_LEN$Length - NA_LEN))))]
        #replace missing eff for some lens            #bio.na.Eff[,c("Code","Gear","Ref")] <- 
        ROWN <-an(row.names(bio.na.Eff[is.na(bio.na.Eff$mult) & bio.na.Eff$FishLength_cm == 79 & bio.na.Eff$Group == NA_SP,]) )#which row?
        bio.na.Eff[ROWN, which(names(bio.na.Eff) %in% names(OTH_LEN))] <- 
          OTH_LEN[ OTH_LEN$Length==USE_LEN, names(bio.na.Eff)[which(names(bio.na.Eff) %in% names(OTH_LEN))] ]
        bio.na.Eff[ROWN,"mult"] <- 1/ bio.na.Eff[ROWN,"Efficiency"]
      }
    }
    
  } 
  bio.na.Eff$sciName<-ac(bio.na.Eff$SpeciesSciName)
  bio <- rbind(bio[!(is.na(bio$Efficiency) | is.infinite(bio$Efficiency)),], bio.na.Eff)
}    #u(bio[is.infinite(bio$mult),'SpeciesSciName']); u(bio[is.na(bio$mult),'SpeciesSciName']); #bio.na.Eff<-bio[is.na(bio$mult),]

bio[is.na(bio$mult),'mult'] <- 1/bio[is.na(bio$mult),'Efficiency']
#? missing for###  check Q for large TOPE Galeorhinus galeus, COMMON SKATE Dipturus batis and Hexanchus griseus

INFS<- (bio.na.Eff[is.infinite(bio.na.Eff$mult),c('SpeciesSciName',"FishLength_cm")])
NAS<- (bio.na.Eff[is.na(bio.na.Eff$mult),c('SpeciesSciName',"FishLength_cm")])
if( nrow(NAS) > 0 | nrow(INFS)>0 ){ 
  print(paste("missing Qmultiplier for", u( bio[ (is.na(bio$mult) | is.infinite(bio$mult)),'SpeciesSciName']) ))
  write(file=paste(OUTPATHstem,"missing Q.csv",sep=""), x=paste("missing Qmultiplier for", u(bio[(is.na(bio$mult) | is.infinite(bio$mult)),'SpeciesSciName'])), sep=",")
  write.csv(file=paste(OUTPATHstem,survey,"_missing Q_bio.csv",sep=""), x=bio[(is.na(bio$mult) | is.infinite(bio$mult)),],append=T)
  bio[is.na(bio$mult),"mult"]<-1
}
bio$DensBiom_kg_Sqkm_beforeQmult <- bio$DensBiom_kg_Sqkm
bio$DensBiom_kg_Sqkm <- bio$DensBiom_kg_Sqkm*bio$mult
