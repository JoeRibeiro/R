#ChrisLynam@Cefas.co.uk
# update 21 Oct 2021 AREASCALE is now an argument with a TRUE default. also added Bel Beam Trawl to lookups for NSea 
# update 07 Feb 2019 to use quadrants and average withing 60km radius on centers
# update 15 Nov 2018 to stop issues with missing data in guilds 
# update 19 Dec 2017 to give IND by taxa Grouping
# update 20 Dec 2016 to give LD by species and total biomass by strata
# update 09 Jan 2017 to pass SP to indicators and enable plotting LOESS fns
# update 09 Jan 2017 now biomass out and plots with and without Q correction
# species_bio_by_area includes numhauls by subdivision if no sampling strata applied (i.e. no rectangles or minigrid)
# average over hauls then raise up by area of S_REG or subdivision
##second level raising if needed, here weighting is correct for coverage of L_REG (sum of S_REG areas should be L_REG area - but prone to issues with missing S_REG)
# check no zeroes in LD output
# Dec 2021. changed looping when species ALL - no longer create combined group
INDfn <- function(DATA, WRITE=F, BOOTSTRAP=F, LFI=T, LFI_THRESHOLD=NULL, FILENAM="", BY_SREG=T, BY_LREG=F, AREASCALE=T,
                  MEANTL=T, MaxL=T, Loo=T, Lm=T, MeanL=T, TyL_GeoM=F, SPECIES = c("DEM"), 
                  GROUP=NULL, TyL_SPECIES=F, BYGUILD=F,QUAD=F,QUAD_SMOOTH=F,QUADS=QUADS, ATTRIB=ATTRIB){
    #BY_SREG<-T; BY_LREG<-T; SPECIES<-c("DEM"); #
    #BOOTSTRAP<-F;  WRITE=T; LFI=T; MEANTL=F; MaxL=T; Loo=F; Lm=F; MeanL=F; TyL_GeoM=T; LFI_THRESHOLD<-50; QUAD<-F; QUAD_SMOOTH<-F
  # DATA<-dhspp; GROUP<-NULL;
  #AREASCALE <- T; #similar to CefMAT if F
  #DATA<-FISHDATA
  LFIout    <- LFI_by_sub    <- FishLength_cmsea    <- MaxLsea    <- Loosea    <- Lmsea    <- TLsea    <- TyL.cm.sea    <- NULL
  LFIout$LFIregional <- NULL
  LFIoutpel <- LFI_by_subpel <- FishLength_cmseapel <- MaxLseapel <- Looseapel <- Lmseapel <- TLseapel <- TyL.cm.seapel <- LFIout$LFIregionalpel <- NULL
  LFIoutdem <- LFI_by_subdem <- FishLength_cmseadem <- MaxLseadem <- Looseadem <- Lmseadem <- TLseadem <- TyL.cm.seadem <- LFIout$LFIregionaldem <- NULL
  
  TyL.cm.sea_pel<- FishLength_cmsea_pel <- MaxLsea_pel <- Loosea_pel <- Lmsea_pel <- TLsea_pel <-  LFIbind_pel <- NULL  
  TyL.cm.sea_dem<- FishLength_cmsea_dem <- MaxLsea_dem <- Loosea_dem <- Lmsea_dem <- TLsea_dem <-  LFIbind_dem <- NULL  
  TyL.cm.sea_all<- FishLength_cmsea_all <- MaxLsea_all <- Loosea_all <- Lmsea_all <- TLsea_all <-  LFIbind_all <- NULL  
  
  
  # how many hauls?    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  numhauls <- DATA;
  
  numhauls$ones <- 0 # 1 val as a marker to pivot around
  #HaulID and StNo are lost so now using HaulID#

  FACTHAUL <-  c("HaulID","Year","Ship","MonthShot","Day","TimeShot", "HaulDur_min","ShootLat_degdec","ShootLong_degdec","ICESStSq",
                 "NetOpen_m", "WingSwpArea_sqkm")
  if(BY_SREG) FACTHAUL <-  c(FACTHAUL,"S_REG")
  if(BY_LREG)   FACTHAUL <-  c(FACTHAUL,"L_REG","S_L_REG")
  
  numhauls <- tapply.ID(df=numhauls, datacols=c("ones"),factorcols=FACTHAUL,sum,c("ones"))
  numhauls$ones <- 1  # now 1 val per haul    
  numhauls$Survey_Acronym <- survey 
  numhauls$Gear <- GEAR 
  numhauls$GearType <- STDGEAR 
  names(numhauls)[which(names(numhauls)=="MonthShot")] <- "Month"
  names(numhauls)[which(names(numhauls)=="WingSwpArea_sqkm")] <- "SweptArea_KM2"
  #save as txt since we have ICES rect codes here and excel will corrupt otherwise
  if(WRITE & (!BOOTSTRAP | (BOOTSTRAP & B==0)  ) ) write.table(numhauls[,-which(names(numhauls)=="ones")],paste(FILENAM,"hauls.txt",sep="_"),row.names =F,sep=',')
  
  FACTHAUL <-  c("ones","HaulID","Year","ShootLat_degdec","ShootLong_degdec")
  if(BY_SREG) FACTHAUL <-  c(FACTHAUL,"S_REG")
  if(BY_LREG)   FACTHAUL <-  c(FACTHAUL,"L_REG","S_L_REG")
  numhauls<- numhauls[,FACTHAUL]
  numhaulsyr<-aggregate(x=numhauls$HaulID,by=list(numhauls$Year), FUN=length)
  names(numhaulsyr)[1]<-"Year"
  names(numhaulsyr)[ncol(numhaulsyr)]<-"H"
  if(WRITE & (!BOOTSTRAP | (BOOTSTRAP & B==0)  ) ) write.table(numhaulsyr,paste(FILENAM,"numhauls_yr.csv",sep="_"),row.names =F,sep=',')
  
  
  #add centlat centlon i.e. centre points of ICES stsqs###
  numhauls$centlon <- floor(numhauls$ShootLong_degdec)+.5
  numhauls$centlat <- round(numhauls$ShootLat_degdec)
  cor <- ifelse(numhauls$centlat < numhauls$ShootLat_degdec, 0.25,  
                ifelse(numhauls$centlat > numhauls$ShootLat_degdec, -0.25,  
                       + 0.25) #if x.00
  )
  numhauls$centlat <- (numhauls$centlat + cor) 
  if(QUAD){
    #Feb2019 correct centers here as have already identified which hauls are in which quads
    numhauls[substr(numhauls$S_REG,6,7)=="SE",]$centlon<- numhauls[substr(numhauls$S_REG,6,7)=="SE",]$centlon-0.25; 
    numhauls[substr(numhauls$S_REG,6,7)=="SE",]$centlat<- numhauls[substr(numhauls$S_REG,6,7)=="SE",]$centlat-0.125
    numhauls[substr(numhauls$S_REG,6,7)=="SW",]$centlon<- numhauls[substr(numhauls$S_REG,6,7)=="SW",]$centlon+0.25
    numhauls[substr(numhauls$S_REG,6,7)=="SW",]$centlat<- numhauls[substr(numhauls$S_REG,6,7)=="SW",]$centlat+0.125 
    numhauls[substr(numhauls$S_REG,6,7)=="NE",]$centlon<- numhauls[substr(numhauls$S_REG,6,7)=="NE",]$centlon-0.25
    numhauls[substr(numhauls$S_REG,6,7)=="NE",]$centlat<- numhauls[substr(numhauls$S_REG,6,7)=="NE",]$centlat-0.125 
    numhauls[substr(numhauls$S_REG,6,7)=="NW",]$centlon<- numhauls[substr(numhauls$S_REG,6,7)=="NW",]$centlon+0.25
    numhauls[substr(numhauls$S_REG,6,7)=="NW",]$centlat<- numhauls[substr(numhauls$S_REG,6,7)=="NW",]$centlat+0.125 
    
    fd<-NULL
    if(QUAD_SMOOTH){
      ##if(BYGUILD & GROUP!="1") next
      source("//lowfilecds/Function/Eco Indicators/DATRASoutput/MarScot/INDscriptsForV3/smoothquad.r") #uses fd lat/lon 
      #need to work out dist to nearest quads and find quadrants within 60km to smooth run by year
      fdhspp <- DATA[,c("Year","S_REG","ICESStSq","HaulID","ShootLat_degdec","ShootLong_degdec")] 
      YRS<-sort(unique(DATA$Year))
      dhspp_match_yrs<-NULL
      for(YR in YRS){   # YR<-1998
        print(paste("smooth data for",survey,YR,sep=" "))
        #fd <- fdhspp[fdhspp$Year==YR,]
        fd <- aggregate(x=fdhspp[fdhspp$Year==YR,c("ShootLat_degdec","ShootLong_degdec")] ,
                        by=list(fdhspp[fdhspp$Year==YR,"HaulID"]), margin=1,FUN=mean)
        
        XYmatch <- smoothquad(FD=fd)
        #XYmatch has the HaulID to average over# link QUADS$cent_lon and QUADS$cent_lat to Xmatch and Ymatch 
        DATCOL<-c("DensBiom_kg_Sqkm","DensBiom_kg_Sqkm_beforeQmult","WingSwpVol_CorF","NetOpen_m")
        FACCOL<-c("SpeciesSciName","FishLength_cm","sciName","Species","Code","Gear","Group","Year","habitat.guild")
        
        dhspp_match<-NULL
        for(m in 1:length(XYmatch)){  # m <- 55
          ##          if(!XYmatch[[m]][3] %in%  DATA[DATA$Year==YR,"S_REG"]) next #skip quads where no original sample within i.e. all from smooth
          dmatch <- DATA[DATA$Year==YR & DATA$HaulID %in% XYmatch[[m]][-c(1:4)], c(DATCOL,FACCOL)] #<- # catch at len per species from matches
          
          if( ncol(XYmatch[[m]])>5 ){  #greater than 5 otherwise only 1 haul in the quad and no need to average (here find sum later /numhauls)
            dmatch <- tapply.ID(df=dmatch, datacols=DATCOL,factorcols=FACCOL, func=sum, newnames=DATCOL,na.stuff=T)
          } ##lose: HaulID,"mult","Ref","Absolute","Abs.l.95","Abs.u.95","Efficiency","Eff.l.95","Eff.u.95","QGroup","LogLngtClass","LogLngtBio","KM2_LAM","L_REG","S_REG"
          dhspp_match<-rbind(dhspp_match,data.frame(ICESStSq=substr(XYmatch[[m]][3],1,4), S_REG=XYmatch[[m]][3],numhauls=XYmatch[[m]][4], dmatch[,c(DATCOL,FACCOL)]))
        }
        dhspp_match$"LogLngtClass" <- log(dhspp_match$FishLength_cm)
        dhspp_match$"LogLngtBio" <- dhspp_match$"LogLngtClass"*dhspp_match$DensBiom_kg_Sqkm
        #dhspp_match$"LogLngtBio_beforeQmult" <- dhspp_match$"LogLngtClass"*dhspp_match$DensBiom_kg_Sqkm_beforeQmult
        #nrow(dhspp_match) > nrow(dhspp[dhspp$Year==YR,])# have more samples as replicated data across quads through smoothing
        #length(unique(dhspp_match$FishLength_cm))==length(unique(dhspp[dhspp$Year==YR,]$FishLength_cm))
        #length(unique(dhspp_match$sciName))==length(unique(dhspp[dhspp$Year==YR,]$sciName))
        #length(unique(dhspp_match$S_REG))==length(unique(dhspp[dhspp$Year==YR,]$S_REG))
        dhspp_match_yrs<- rbind(dhspp_match_yrs, dhspp_match)
      } 
      #add back L_REG S_L_REG centlon centlat fguild
      DATA <- dhspp_match_yrs
      rm(dhspp_match_yrs,fd,fdhspp)
      ##if(BYGUILD & GROUP=="1") dhspp <<- DATA #<<- to make sure this is updated in the global env #problem as lose some cols #note will overwrite with original dhspp_raw after guild loop
    }#end smooth
    
    DATA <- merge(DATA,QUADS[,c("QUADNAME","KM2_LAM","cent_lat","cent_lon")],by.x=("S_REG"),by.y=("QUADNAME"))
    SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNS_rectstrat/GNSIntOT/GNSstrat_Atlantis.shp",sep='') ) 
    if(EHDS_PP) SUBDIV <- readShapeSpatial(paste(SHAPEPATH,"GNS_EHDPP/ehu_polygons.shp",sep='') ) 
    #if(BY_LREG) "L_REG" <- "LFIregion" #old spatial areas - 25 year plan
    #NAMS_REG<-"ICESNAME"
    #ATTRIB <- read.csv(paste(SHAPEPATH,"attributes/",survey,".csv",sep=''))
    SAMP_FACT <- "KM2_LAM"
    if(BY_SREG){ SAMP_FACT <- c(SAMP_FACT, "S_REG") }
    if(BY_LREG){ SAMP_FACT <- c(SAMP_FACT, "L_REG")} 
    if(EHDS_PP){  
      ATTRIB <- read.csv(paste(SHAPEPATH,"attributes/GNS_EHDPP.csv",sep='') ) 
      SAMP_FACT <- c("KM2_LAM", "L_REG") 
    } 
    #if(OVERWITE_SUBDIV) DATA$L_REG<-DATA$S_REG###01Feb2017
    
    ATTRIB <- ATTRIB[,which(names(ATTRIB) %in% SAMP_FACT )]
    #area relates to lowest sampling strata (i.e. rects, minigrid or survey strata poly)
    #L_REG area - if using by rectangle S_REG need to sum area for SUBDIV
    # if(survey %in% c("GNSIntOT1","GNSIntOT1_channel","GNSIntOT3","GNSNetBT3","GNSGerBT3","GNSBelBT3", "GNSNetBi3", "GNSIntBi3")){
    #   ATTRIB_SUBDIV <- aggregate(x=ATTRIB$KM2_LAM,by=list(L_REG=ATTRIB$L_REG), FUN=sum)
    #   names(ATTRIB_SUBDIV) <- c("L_REG","KM2_LAM")
    # }
    #add ShootLong_degdec ShootLat_degdec as average of hauls? no need use cent_lat and cent_lon
    
    #SUBDIV 
    dhspp0<- DATA 
    coordinates(dhspp0) <- ~ cent_lon +cent_lat
    ox <- over(dhspp0, SUBDIV) #bring in all attributes of location i..e both S_REG and L_REG if applicable 
    ## which are the subdivisions and sampling stratification units
    if(BY_LREG) names(ox)[which(names(ox)=="L_REG")] <- "L_REG" 
    if(EHDS_PP){  
      names(ox)[which(names(ox)=="area_1")] <- "KM2_LAM" #rename as not in shp correct
      DATA <- DATA[!is.na(ox$L_REG) & ox$L_REG!="Other", ] # some areas were cut for PP
      ox <- ox[!is.na(ox$L_REG) & ox$L_REG!="Other", ]
    }
    DATA <- cbind(DATA,ox[,c("L_REG")])  ### issue here for GNSGerBT3 41F4 and 41F5 do not agree with samp file for 2001 since long shoot == 5 exactly
    names(DATA)[ncol(DATA)]<-"L_REG"
    if(QUAD) DATA<- DATA[!is.na(DATA$L_REG),]#plot(sdr); map(add=T); with(DATA[is.na(DATA$L_REG),],points(cent_lon,cent_lat,pch=19,col=4)) # a couple of odd points - poss on land - can be smoothed in
  } #end quad 
  ##DATA is updated so have smooth total catch from hauls <60km from q_center by quadrant
  #DATA[DATA$L_REG!=DATA$L_REG,]
  ##load("C:/Users/cl06/Desktop/biodiv19 ref/allyrs_DATA_QUADmatch.RData")
  if(BY_SREG){ #might be STSQ, QUADrants or minigrid see S_REG
    FACTHAUL <-  c("Year","centlon","centlat","S_REG")
    if(BY_LREG) FACTHAUL <-  c(FACTHAUL,"L_REG","S_L_REG")
    
    if(!QUAD | !QUAD_SMOOTH) numhaulsBYS_REG <- tapply.ID(df=numhauls, datacols=c("ones"), factorcols=FACTHAUL, sum,c("numhauls"));
    if(QUAD & QUAD_SMOOTH){ 
      DATA$numhauls<- as.numeric(as.character(DATA$numhauls))#numbers read in as factors
      numhaulsBYS_REG <- aggregate( x=DATA$numhauls,by=list(Year=DATA$Year,S_REG=DATA$S_REG,L_REG=DATA$L_REG),FUN=mean)
      names(numhaulsBYS_REG)[which(names(numhaulsBYS_REG)=="x")]<-"numhauls"
    }
    #if(WRITE & (!BOOTSTRAP | (BOOTSTRAP & B==0)  ) ) write.table(numhaulsBYS_REG,paste(FILENAM,"numhaulsBYS_REG.csv",sep="_"),row.names =F,sep=',')
    #and reshape since have one value per year and L_REG combination
    numhaulsBYS_REGout <- (tapply(numhaulsBYS_REG$numhauls,list(numhaulsBYS_REG$Year, numhaulsBYS_REG$S_REG), FUN=sum, na.rm=T))
    if(WRITE  & (!BOOTSTRAP | (BOOTSTRAP & B==0)  ) ) write.table(numhaulsBYS_REGout,paste(FILENAM,"numhaulsBYS_REG.csv",sep="_"),row.names =T,sep=',')
    rm(numhaulsBYS_REGout)
    
  } else { numhaulsBYS_REG <- NULL }
  
  if(BY_LREG){#user_defined or survey poly
    if(!QUAD) numhaulsBYsubdiv <- tapply.ID(df=numhauls, datacols=c("ones"),
                                            factorcols=c("Year","L_REG"), sum,c("numhauls"));
    if(QUAD)  numhaulsBYsubdiv <- tapply.ID(df=numhaulsBYS_REG, datacols=c("numhauls"),
                                            factorcols=c("Year","L_REG"),sum,c("numhauls"));
    #and reshape since have one value per year and L_REG combination
    numhaulsBYsubdivout <- (tapply(numhaulsBYsubdiv$numhauls,list(numhaulsBYsubdiv$Year, numhaulsBYsubdiv$L_REG), FUN=sum, na.rm=T))
    if(WRITE  & (!BOOTSTRAP | (BOOTSTRAP & B==0)  ) ) write.table(numhaulsBYsubdivout,paste(FILENAM,"numhaulsBYsubdiv.csv",sep="_"),row.names =T,sep=',')
    rm(numhaulsBYsubdivout)
  } else { numhaulsBYsubdiv <- NULL }
  if(!QUAD){  
    numhaulsyr <- tapply.ID(df=numhauls, datacols=c("ones"),factorcols=c("Year"),sum,c("numhauls"));  # now 1 val per STSQ
    numhauls<-numhauls[,-1] 
  } #now just a list of hauls
  if(QUAD){ 
    numhauls<-numhaulsBYS_REG
    numhaulsBYS_REG$ones <- 1
    numhaulsyr <- tapply.ID(df=numhaulsBYS_REG, datacols=c("ones"),factorcols=c("Year"),sum,c("numhauls"));  # now 1 val per STSQ
  }
  #browser()
  #plot(numhaulsBYsubdiv[numhaulsBYsubdiv$BOX_ID==11,2:1])
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #### species_occurrence_area ####
  # species records by rect (strata1) and subdivision (strata2) from haul data
  if(WRITE  & (!BOOTSTRAP | (BOOTSTRAP & B==0) ) ){
    
    dhfc1 <- DATA[DATA$SpeciesSciName %in% FC1Sp,]
    dhfc1 <- aggregate(x=dhfc1$DensAbund_N_Sqkm, 
                       by=list(Year=dhfc1$Year, HaulID=dhfc1$HaulID, SpeciesSciName=dhfc1$SpeciesSciName,
                               dhfc1$ShootLong_degdec, dhfc1$ShootLat_degdec), FUN=sum) 
    names(dhfc1)[(length(names(dhfc1))-2):length(names(dhfc1))]<-c("lon","lat","N")
    write.csv(dhfc1,paste0(FILENAM,"_haul_by_FC1_",survey,".csv"),row.names = F)
    
    if(!BY_LREG){ 
      if(!BY_SREG){
        DATANOLEN<-aggregate(x=DATA$DensBiom_kg_Sqkm, by=list(Year=DATA$Year, SpeciesSciName=DATA$SpeciesSciName, Group=DATA$Group, HaulID=DATA$HaulID), FUN=sum)  
        species_rec_by_area <- aggregate(x=DATANOLEN$x, by=list(Year=DATANOLEN$Year, SpeciesSciName=DATANOLEN$SpeciesSciName), FUN=length)  
        #SPECIESSTRATRECS <- xtabs(N ~ Year + SpeciesSciName,data=species_rec_by_area)
        #write.table( SPECIESSTRATRECS,paste(FILENAM,SPECIES,"species_records_yr.csv",sep="_"),row.names =T,sep=',')
      } else { 
        DATANOLEN<-aggregate(x=DATA$DensBiom_kg_Sqkm, by=list(Year=DATA$Year, SpeciesSciName=DATA$SpeciesSciName, Group=DATA$Group, HaulID=DATA$HaulID, S_REG=DATA$S_REG), FUN=sum)  
        species_rec_by_area <- aggregate(x=DATANOLEN$x, by=list(Year=DATANOLEN$Year, SpeciesSciName=DATANOLEN$SpeciesSciName, S_REG=DATANOLEN$S_REG), FUN=length)  
        names(species_rec_by_area)[length(names(species_rec_by_area))]<-"N"
        SPECIESSTRATRECS <- xtabs(N ~ Year + SpeciesSciName + S_REG,data=species_rec_by_area)
        #if(survey != "GNSFraOT4") for(i in 1:dim(SPECIESSTRATRECS)[3]) write.table( SPECIESSTRATRECS[,,i],paste(FILENAM,SPECIES,"species_records_S_REG",i,"_yr.csv",sep="_"),row.names =T,sep=',')
      }
    }
    
    if(BY_LREG){        
      if(!BY_SREG){ 
        DATANOLEN<-aggregate(x=DATA$DensBiom_kg_Sqkm, by=list(Year=DATA$Year, SpeciesSciName=DATA$SpeciesSciName, Group=DATA$Group, HaulID=DATA$HaulID, L_REG=DATA$L_REG), FUN=sum)  
        species_rec_by_area <- aggregate(x=DATANOLEN$x, by=list(Year=DATANOLEN$Year, SpeciesSciName=DATANOLEN$SpeciesSciName, L_REG=DATANOLEN$L_REG), FUN=length)  
        names(species_rec_by_area)[length(names(species_rec_by_area))]<-"N"
        SPECIESSTRATRECS <- xtabs(N ~ Year + SpeciesSciName + L_REG,data=species_rec_by_area)
        for(i in 1:dim(SPECIESSTRATRECS)[3]){
          write.table( matrix(c( dimnames(SPECIESSTRATRECS)$L_REG[i], colnames(SPECIESSTRATRECS[,,i]) ),nrow=1),paste(FILENAM,SPECIES,"species_records_subdiv",i,"_yr.csv",sep="_"),row.names =F,col.names =F,sep=',',append = FALSE)
          write.table( SPECIESSTRATRECS[,,i],paste(FILENAM,SPECIES,"species_records_subdiv",i,"_yr.csv",sep="_"),row.names =T,col.names =F,sep=',',append = TRUE)
        } 
      }
      if(BY_SREG){ 
        DATANOLEN<-aggregate(x=DATA$DensBiom_kg_Sqkm, by=list(Year=DATA$Year, SpeciesSciName=DATA$SpeciesSciName, Group=DATA$Group, HaulID=DATA$HaulID, S_REG=DATA$S_REG,L_REG=DATA$L_REG), FUN=sum)  
        species_rec_by_area <- aggregate(x=DATANOLEN$x, by=list(Year=DATANOLEN$Year, SpeciesSciName=DATANOLEN$SpeciesSciName, S_REG=DATANOLEN$S_REG, L_REG=DATANOLEN$L_REG), FUN=length)  
        names(species_rec_by_area)[length(names(species_rec_by_area))]<-"N"
        SPECIESSTRATRECS <- xtabs(N ~ Year + SpeciesSciName + L_REG + S_REG,data=species_rec_by_area)
        #for(i in 1:dim(SPECIESSTRATRECS)[3]) write.table( SPECIESSTRATRECS[,,i,],paste(FILENAM,SPECIES,"species_records_strat_div",i,"_yr.csv",sep="_"),row.names =T,sep=',')
      }
      
      species_records_BY_reg<- aggregate(species_rec_by_area$N, by=list(Year=species_rec_by_area$Year, SpeciesSciName=species_rec_by_area$SpeciesSciName, 
                                                                        L_REG=species_rec_by_area$L_REG), FUN=sum)
      names(species_records_BY_reg)[length(names(species_records_BY_reg))]<-"N"
      species_records_BY_reg <- xtabs(N ~ Year + SpeciesSciName + L_REG,data=species_records_BY_reg)
      write.table(species_records_BY_reg,paste(FILENAM,SPECIES,"species_records_BY_reg_yr.csv",sep="_"),row.names =F,sep=',')
    }
    
    #species_rec_by_area <- cbind(numhaulsyr,species_rec_by_area)
    
    FC1spp_rec_by_area <- xtabs(N ~ Year + SpeciesSciName,data=species_rec_by_area[species_rec_by_area$SpeciesSciName %in% FC1Sp,])
    write.table(matrix(c(survey,colnames(FC1spp_rec_by_area)),nrow=1),paste(FILENAM,SPECIES,"_FC1_records_BY_yr.csv",sep="_"),row.names =F,col.names =F,sep=',',append = FALSE)
    write.table(FC1spp_rec_by_area,paste(FILENAM,SPECIES,"_FC1_records_BY_yr.csv",sep="_"),row.names =T,col.names =F,sep=',',append = TRUE)
    
    species_rec_by_area <- xtabs(N ~ Year + SpeciesSciName,data=species_rec_by_area)
    write.table(matrix(c(survey,colnames(species_rec_by_area)),nrow=1),paste(FILENAM,SPECIES,"species_records_BY_yr.csv",sep="_"),row.names =F,col.names =F,sep=',',append = FALSE)
    write.table(species_rec_by_area,paste(FILENAM,SPECIES,"species_records_BY_yr.csv",sep="_"),row.names =T,col.names =F,sep=',',append = TRUE)
    rm(DATANOLEN,species_rec_by_area)
  }#end WRITE
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   species_bio_by_area    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # species bioLD by rect (strata1) and subdivision (strata2) from haul data
  
  # run tapply.id over these factorcols FACT
  FACT <- c("Year","FishLength_cm","SpeciesSciName")
  if(BYGUILD) FACT <- c(FACT,"Group")
  if(BY_SREG) FACT <- c(FACT,"S_REG")
  if(BY_LREG)  FACT <- c(FACT,"L_REG")
  if(QUAD)  FACT <- c(FACT,"ICESStSq")
  #length(dhspp$DensBiom_kg_Sqkm[is.na(dhspp$DensBiom_kg_Sqkm)]) #check
  #length(DATA$DensBiom_kg_Sqkm[is.na(DATA$DensBiom_kg_Sqkm)]) #check
  #length(species_bio_by_area$DensBiom_kg_Sqkm[is.na(species_bio_by_area$DensBiom_kg_Sqkm)]) #check
  
  #sum by species by length cat by S_REG (e.g ICES STSQ)
  DATACOLS<-c("DensBiom_kg_Sqkm","DensAbund_N_Sqkm")
  if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) DATACOLS<- c(DATACOLS, "DensBiom_kg_Sqkm_beforeQmult")
  NEWDATANAM<-c("CatCatchWgtSwept","Abund_N_Swept")
  if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) NEWDATANAM<-c(NEWDATANAM,"CatCatchWgtSwept_beforeQmult")
  # and subdivisional strata (e.g. NE North Sea or 'survstrata') #DATA<-dhspp
  suppressWarnings( #NAs introduced by coercion since some MaxL and TL are NA
    species_bio_by_area <- tapply.ID(df=DATA, datacols=DATACOLS, factorcols=FACT, sum,NEWDATANAM)  
  )  ##Feb2019 poss to  create NA subdivs!!!!! 
  
  if(LFI & is.null(LFI_THRESHOLD)){
    #work out LF threshold for 20%biomass
    total_biomass <- sum(species_bio_by_area$CatCatchWgtSwept,na.rm=T)  
    bio_cum <- cumsum(x= species_bio_by_area$CatCatchWgtSwept[ order(species_bio_by_area$FishLength_cm) ])
    plot( species_bio_by_area$FishLength_cm[ order(species_bio_by_area$FishLength_cm) ], bio_cum)
    abline(h=total_biomass,col=2)
    abline(h=total_biomass*.8,col=3)
    ABSDIFF <- abs(bio_cum - (total_biomass*.8))
    LFI_THRESHOLD <- 
      species_bio_by_area$FishLength_cm[ order(species_bio_by_area$FishLength_cm) ][which(ABSDIFF== min(abs(ABSDIFF),na.rm=T) )]
    abline(v=LFI_THRESHOLD,col=3)
  }
  
  # species_bio_by_area[species_bio_by_area$L_REG!=species_bio_by_area$L_REG,]
  # to average LD must add num hauls to bio data
  if(BY_SREG){
    species_bio_by_area <- merge(x = species_bio_by_area,
                                 y = numhaulsBYS_REG[,which(names(numhaulsBYS_REG) != "S_L_REG" & names(numhaulsBYS_REG) != "L_REG")], #avoid replicating names and creating .x .y
                                 by = c("Year","S_REG"),all.x=T)
    
    #average species cpue over hauls by rectangle-strata for MaxL, TL , Len, TyL                                                                      
    species_bio_by_area$CatCatchWgtSwept <- species_bio_by_area$CatCatchWgtSwept / species_bio_by_area$numhauls                                        
    species_bio_by_area$Abund_N_Swept <- species_bio_by_area$Abund_N_Swept / species_bio_by_area$numhauls                                        
    if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) species_bio_by_area$CatCatchWgtSwept_beforeQmult <- species_bio_by_area$CatCatchWgtSwept_beforeQmult / species_bio_by_area$numhauls
    
  } else {
    if(BY_LREG){ # only do here if not using rects/minigrid/etc
      species_bio_by_area <- merge(x = species_bio_by_area,
                                   y = numhaulsBYsubdiv[,which(names(numhaulsBYsubdiv) != "S_L_REG")], #avoid replicating names and creating .x .y
                                   by = c("Year","L_REG"),all.x=T)
      
      #average species cpue over hauls by rectangle-strata for MaxL, TL , Len, TyL                                                                      
      species_bio_by_area$CatCatchWgtSwept <- species_bio_by_area$CatCatchWgtSwept / species_bio_by_area$numhauls                                      
      species_bio_by_area$Abund_N_Swept <- species_bio_by_area$Abund_N_Swept / species_bio_by_area$numhauls
      if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) species_bio_by_area$CatCatchWgtSwept_beforeQmult <- species_bio_by_area$CatCatchWgtSwept_beforeQmult / species_bio_by_area$numhauls
    }
  }
  #if both BY_SREG and BY_LREG are true will have to sum up catch by BY_SREG within SUBDIV later to avoid change between years due to change in relative sampling of strata..
  #length(species_bio_by_area$DensBiom_kg_Sqkm[is.na(species_bio_by_area$DensBiom_kg_Sqkm)]) #check
  #length(species_bio_by_area$L_REG[is.na(species_bio_by_area$L_REG)]) #check
  
  #re-introduce MAXL, Loo, Lm, TL, DEMPEL - warning here is an opportunity for NAs to appear!
  species_bio_by_area <- merge(species_bio_by_area,LW[,c("ScientificName_WoRMS","SensFC1","DEMPEL")], by.x="SpeciesSciName",by.y="ScientificName_WoRMS",all.x=T)
  #if(!BYGUILD) species_bio_by_area <- merge(species_bio_by_area,trait_MAXL[,c("SpeciesSciName","maximum.length","Loo","Lm","Order","Group","LFI_Fung_list","LFI_OSPAR_list")],by="SpeciesSciName",all.x=T)
  #if(BYGUILD) species_bio_by_area <- merge(species_bio_by_area,trait_MAXL[,c("SpeciesSciName","maximum.length","Loo","Lm","Order","LFI_Fung_list","LFI_OSPAR_list")],by="SpeciesSciName",all.x=T)
  if(!BYGUILD) species_bio_by_area <- merge(species_bio_by_area,trait_MAXL[,c("SpeciesSciName","MaxL","Loo","Lm","Order","Group")],by="SpeciesSciName",all.x=T)
  if(BYGUILD)  species_bio_by_area <- merge(species_bio_by_area,trait_MAXL[,c("SpeciesSciName","MaxL","Loo","Lm","Order")],by="SpeciesSciName",all.x=T)
  
  species_bio_by_area$Group <- ac(species_bio_by_area$Group)
  # species_bio_by_area[is.na(species_bio_by_area$DEMPEL),] 
  #names(species_bio_by_area)[which( names(species_bio_by_area)=="maximum.length")] <- "MaxL"
  species_bio_by_area$DEMPEL <-as.character(species_bio_by_area$DEMPEL)
  species_bio_by_area$DEMPEL[species_bio_by_area$DEMPEL=="Demersal"] <- "DEM"
  species_bio_by_area$DEMPEL[species_bio_by_area$DEMPEL=="Pelagic"] <- "PEL"
  species_bio_by_area$DEMPEL <-as.character(species_bio_by_area$DEMPEL)
  
  #Trophic Level FW4
  if(MEANTL){
    if(substr(survey,1,2) == "GN") species_bio_by_area <- merge(x=species_bio_by_area,y=TLnorth,by="SpeciesSciName",all.x=T,all.y=F)
    if(substr(survey,1,2) %in% c("CS","BB")) species_bio_by_area <- merge(x=species_bio_by_area,y=TLceltic,by="SpeciesSciName",all.x=T,all.y=F)
    if(substr(survey,1,2) == "WA") MEANTLs <- F
  }
  #save a copy 'species_bio_by_area_DEMPEL' so can loop through DEM or PEL etc
  species_bio_by_area_DEMPEL <- species_bio_by_area
  #species_bio_by_area[species_bio_by_area$L_REG!=species_bio_by_area$L_REG,]
  #length(species_bio_by_area$DensBiom_kg_Sqkm[is.na(species_bio_by_area$DensBiom_kg_Sqkm)]) #check
  #length(species_bio_by_area$L_REG[is.na(species_bio_by_area$L_REG)]) #check
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # record sampling effort for indicators
  # e.g. num rects sampled by L_REG  sumS_REG_by_sub   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(BY_SREG){
    FACT<-c("Year","S_REG")
    if(BY_LREG) FACT<-c(FACT,"L_REG")
    numS_REG_by_sea <- tapply.ID(df=species_bio_by_area, datacols=c("CatCatchWgtSwept"), 
                                     factorcols=FACT, sum,c("CatCatchWgtSwept"))
    numS_REG_by_sea$numS_REG <- 1
    if(BY_LREG){ 
      numS_REG_by_sub <- tapply.ID(df=numS_REG_by_sea, datacols=c("numS_REG"), factorcols=c("Year","L_REG"), sum,c("numS_REG"))
      sumS_REG_by_sub <- xtabs(numS_REG ~ Year + L_REG, numS_REG_by_sub)
    }
    numS_REG_by_sea <- tapply.ID(df=numS_REG_by_sea, datacols=c("numS_REG"), factorcols=c("Year"), sum,c("numS_REG"))
    if(BY_LREG){ 
      numS_REG_by_sea <- cbind(sumS_REG_by_sub,sea=numS_REG_by_sea[,1])
      rm(sumS_REG_by_sub)
    }
    if(WRITE  & (!BOOTSTRAP | (BOOTSTRAP & B==0)  ) ) write.table(numS_REG_by_sea,paste(FILENAM,"num_rects_sampled_BY_reg_yr.csv",sep="_"),row.names =T,sep=',')
  }
  if(BY_LREG){ #if no BY_SREG and S_REG=NA, then above gives same as this
    num_by_sub <- tapply.ID(df=species_bio_by_area, datacols=c("CatCatchWgtSwept"), 
                            factorcols=c("Year","L_REG"), sum,c("CatCatchWgtSwept")) 
    num_by_sub$numsamp <- 1
    sum_by_sub <- xtabs(numsamp ~ Year + L_REG, num_by_sub)
    num_by_sea <- tapply.ID(df=num_by_sub, datacols=c("numsamp"), factorcols=c("Year"), sum,c("numsamp"))
    if(BY_LREG) num_by_sea <- cbind(sum_by_sub,sea=num_by_sea[,1])
    rm(sum_by_sub)
    if(WRITE  & (!BOOTSTRAP | (BOOTSTRAP & B==0)  ) ) write.table(num_by_sea,paste(FILENAM,"num_subdiv_sampled_BY_yr.csv",sep="_"),row.names =T,sep=',')
    
    #correction for regional sea sampling area required if missing part of SUBDIV
    #area sampled
    num_by_sub <- merge(x=num_by_sub,y=ATTRIB,by.x="L_REG",by.y="L_REG",all=T)
    areasurveyed_by_sub <- tapply.ID(df=num_by_sub, datacols=c("KM2_LAM"), 
                                     factorcols=c("Year"), sum,c("KM2_LAM")) 
    #proportion of regional sea area sampled #ATTRIB_SUBDIV is same as totalarea for GNS 'BY_SREG+BY_LREG'
    if(survey %in% c("GNSIntOT1","GNSIntOT1_channel","GNSIntOT3","GNSNetBT3","GNSGerBT3","GNSBelBT3", "GNSNetBi3", "GNSIntBi3")) {
      totalarea <- sum(ATTRIB_SUBDIV$KM2_LAM)
    } else { totalarea <- sum(ATTRIB$KM2_LAM) }
    areasurveyed_by_sub$scale <- totalarea/areasurveyed_by_sub$KM2_LAM
    if(length(areasurveyed_by_sub[areasurveyed_by_sub$scale>1,'Year']) >0) print( paste("survey area not fully covered in", areasurveyed_by_sub[areasurveyed_by_sub$scale>1,'Year'] ) )
    if(WRITE  & (!BOOTSTRAP | (BOOTSTRAP & B==0)  ) ) write.table(areasurveyed_by_sub,paste(FILENAM,"areasurveyed_by_sub.csv",sep="_"),row.names =F,sep=',')
    #not used further
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #LD+catch raised by kmsq for species    meanLD_bio_by_area[meanLD_bio_by_area$L_REG!=meanLD_bio_by_area$L_REG,]
  #summarise by S_REG+L_REG raised by spatial area in kmsq and then kg->tonnes
  if(BY_LREG | BY_SREG){
    #find sampling areas from ATTRIB and merge with species data
    if(BY_SREG){ 
      if(QUAD){  
        ATTRIB$KM2_LAM <- ATTRIB$KM2_LAM/4
        ATTRIB1<-ATTRIB; ATTRIB1$S_REG<- paste(ATTRIB1$S_REG,"SW",sep="_")
        ATTRIB2<-ATTRIB; ATTRIB2$S_REG<- paste(ATTRIB2$S_REG,"SE",sep="_")
        ATTRIB3<-ATTRIB; ATTRIB3$S_REG<- paste(ATTRIB3$S_REG,"NW",sep="_")
        ATTRIB4<-ATTRIB; ATTRIB4$S_REG<- paste(ATTRIB4$S_REG,"NE",sep="_")
        ATTRIB <-rbind(ATTRIB1,ATTRIB2,ATTRIB3,ATTRIB4)
        rm(ATTRIB1,ATTRIB2,ATTRIB3,ATTRIB4)
      }
      meanLD_bio_by_area <- merge(x=ATTRIB, 
                                  y=species_bio_by_area[,!colnames(species_bio_by_area) %in% "L_REG"], 
                                  all.y=TRUE,
                                  by=c("S_REG") ) 
    } else {
      meanLD_bio_by_area <- merge(x=ATTRIB, y=species_bio_by_area, all.y=TRUE, by.x=c("L_REG") , by.y=c("L_REG") ) 
      if( any(names(meanLD_bio_by_area) %in% "L_REG") & !any(names(meanLD_bio_by_area) %in% "L_REG") ) meanLD_bio_by_area$L_REG <-meanLD_bio_by_area$L_REG
      #meanLD_bio_by_area <- meanLD_bio_by_area[,-which(names(meanLD_bio_by_area) == "L_REG")] #rm duplicate col
    }
    #raise by area of lowest resolution of sampling strategy/L_REG
    if(AREASCALE){
      meanLD_bio_by_area$CatCatchWgtSwept <- meanLD_bio_by_area$CatCatchWgtSwept*meanLD_bio_by_area$KM2_LAM
      meanLD_bio_by_area$Abund_N_Swept <- meanLD_bio_by_area$Abund_N_Swept*meanLD_bio_by_area$KM2_LAM
      if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) meanLD_bio_by_area$CatCatchWgtSwept_beforeQmult <- meanLD_bio_by_area$CatCatchWgtSwept_beforeQmult*meanLD_bio_by_area$KM2_LAM
    } 
    #return in tonnes
    meanLD_bio_by_area$CatCatchWgtSwept <- meanLD_bio_by_area$CatCatchWgtSwept/1000
    if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) meanLD_bio_by_area$CatCatchWgtSwept_beforeQmult <- meanLD_bio_by_area$CatCatchWgtSwept_beforeQmult/1000
    
    #output
    if(BY_SREG & WRITE_LDs & (!BOOTSTRAP | (BOOTSTRAP & B==0) )) write.csv(meanLD_bio_by_area,file = paste(paste(FILENAM,sep='_'),"LD_tonnes_Year_W.by.S_REG.csv",sep=".") )
    if(BY_LREG & !BY_SREG & WRITE_LDs & (!BOOTSTRAP | (BOOTSTRAP & B==0)) ) write.csv(meanLD_bio_by_area,file = paste(paste(FILENAM,sep='_'),"LD_tonnes_Year_W.by.L_REG.csv",sep=".") )
    
    #second level raising if both levels applied i.e. c("GNSGerBT3","GNSBelBT3","GNSNetBT3","GNSIntOT1","GNSIntOT3")
    if(BY_LREG & BY_SREG){
      #use catches scaled by size of grid (rects not constant over sea area)
      # and scale to SUBDIV (beware GNSGerBT3 only sampled a small part of NE so should not do this)
      if(!survey %in% c("GNSIntOT1","GNSIntOT1_channel","GNSIntOT3","GNSNetBT3","GNSGerBT3","GNSBelBT3", "GNSNetBi3", "GNSIntBi3")) print(paste(survey,"survey does not have two level stratification"))
      #work out value to scale up L_REG by
      #lose species and length (otherwise inflate sum of areas)                           factorcols=c("S_REG","Year","L_REG")
      area_by_subdiv <- tapply.ID(df=species_bio_by_area, datacols=c("CatCatchWgtSwept"), factorcols=c("S_REG","Year"), sum,c("CatCatchWgtSweptsum"))  
      #merge in area for sample coverage
      area_by_subdiv <- merge(x=ATTRIB, y=area_by_subdiv, all.y=TRUE, by=c("S_REG") )#area_by_subdiv can be from sum of ICESStSq or quadrant here 
      #sum area by L_REG sampled
      #area_by_subdiv1 <- tapply.ID(df=area_by_subdiv, datacols=c("KM2_LAM"), factorcols=c("Year","L_REG"), sum,c("KM2_LAMsum"))  
      area_by_subdiv <- tapply.ID(df=area_by_subdiv, datacols=c("KM2_LAM"), factorcols=c("Year","L_REG"), sum,c("KM2_LAMsum"))  
      
      #compare to area of subdivision for survey (all years)
      area_by_subdiv <- merge(x=ATTRIB, y=area_by_subdiv, all.y=TRUE, by.x=c("L_REG") , by.y=c("L_REG") ) 
      #ratio to scale up to L_REG estimate
      area_by_subdiv$scale <- area_by_subdiv$KM2_LAM / area_by_subdiv$KM2_LAMsum
      if(WRITE  & (!BOOTSTRAP | (BOOTSTRAP & B==0)  ) ) write.csv(area_by_subdiv, file = paste(FILENAM,"_area_by_subdiv.csv",sep=''))
      #big raising factors?#            #area_by_subdiv[area_by_subdiv$scale>2,]
      #problems?   # area_by_subdiv[area_by_subdiv$scale<1,]
      #now sum catch by S_REG to L_REG area and scale for missing area (each year)
      FACT <- c("Year","FishLength_cm","SpeciesSciName","L_REG","DEMPEL","Order","Group")
      DATACOLS<-c("CatCatchWgtSwept","Abund_N_Swept")
      if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) DATACOLS<- c(DATACOLS, "CatCatchWgtSwept_beforeQmult")
      
      meanLD_bio_by_subdiv <- tapply.ID(df=meanLD_bio_by_area, datacols=DATACOLS, factorcols=FACT,sum,DATACOLS)  
      meanLD_bio_by_subdiv <- merge(x=area_by_subdiv[,c("L_REG","Year","scale")], y=meanLD_bio_by_subdiv, all.y=TRUE, by.x=c("L_REG","Year") , by.y=c("L_REG","Year") )  
      if(AREASCALE){
        meanLD_bio_by_subdiv$CatCatchWgtSwept <- meanLD_bio_by_subdiv$scale*meanLD_bio_by_subdiv$CatCatchWgtSwept
        meanLD_bio_by_subdiv$Abund_N_Swept <- meanLD_bio_by_subdiv$scale*meanLD_bio_by_subdiv$Abund_N_Swept
        if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) meanLD_bio_by_subdiv$CatCatchWgtSwept_beforeQmult <- meanLD_bio_by_subdiv$scale*meanLD_bio_by_subdiv$CatCatchWgtSwept_beforeQmult
      }
      if(WRITE_LDs & (!BOOTSTRAP | (BOOTSTRAP & B==0))) write.csv(meanLD_bio_by_subdiv,file = paste(paste(FILENAM,sep='_'), "LD_tonnes_Year_W.by.L_REG.csv",sep=".") )
      #lost DEMPEL, MAXL and TL again
    }
    
    #if missing completely a survey stratum from sampling will have underestimate here
    DATACOLS<-c("CatCatchWgtSwept","Abund_N_Swept")
    if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) DATACOLS<- c(DATACOLS, "CatCatchWgtSwept_beforeQmult")
    #mean length distribution by species over all sampling strat
    meanLD_bio_by_areaNOYEAR <- tapply.ID(df=meanLD_bio_by_area, datacols=DATACOLS, factorcols= c("FishLength_cm","SpeciesSciName"), mean,DATACOLS)  
    if(WRITE_LDs & (!BOOTSTRAP | (BOOTSTRAP & B==0))) write.csv(meanLD_bio_by_areaNOYEAR,file = paste(paste(FILENAM,sep='_'),"LD_tonnes_YEARave.csv",sep=".") )
    
    #third level biomass to regional sea (or survey extent i.e. coverage of 'L_REG' in year)
    for(SP in SPECIES){
      if(SP=="ALL") SP<-c("DEM","PEL")
      if(BY_LREG & BY_SREG){
        bio_spp_subdiv <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","L_REG","SpeciesSciName","Group"), sum,DATACOLS)  
        bio_spp_area <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","SpeciesSciName","Group"), sum,DATACOLS)  
        bio_by_subdiv <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","L_REG"), sum,DATACOLS)  
        bio_by_area <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year"), sum,DATACOLS)  
        bio_grp_subdiv <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","L_REG","Group"), sum,DATACOLS)  
        bio_grp_area <- tapply.ID(df=meanLD_bio_by_subdiv[meanLD_bio_by_subdiv$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","Group"), sum,DATACOLS)  
      } else {#one or other of BY_LREG | BY_SREG
        if(BY_LREG) bio_spp_subdiv <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","L_REG","SpeciesSciName","Group"), sum,DATACOLS)  
        if(BY_LREG) bio_by_subdiv <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","L_REG"), sum,DATACOLS)  
        if(BY_LREG) bio_grp_subdiv <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","L_REG","Group"), sum,DATACOLS)  
        bio_spp_area <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","SpeciesSciName","Group"), sum,DATACOLS)  
        bio_by_area  <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year"), sum,DATACOLS)  
        bio_grp_area <- tapply.ID(df=meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL %in% SP,], datacols=DATACOLS, factorcols=c("Year","Group"), sum,DATACOLS)  
      }
      if(length(SP)==2) SP<-"ALL"
      if(WRITE & (!BOOTSTRAP | (BOOTSTRAP & B==0)) ){
        #all
        if(BY_LREG) write.csv(bio_spp_subdiv,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_AllSpecies_subdivYear.csv",sep="") )
        if(BY_LREG) write.csv(bio_by_subdiv,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_subdivYear.csv",sep="") )
        if(BY_LREG & !BYGUILD) write.csv(bio_grp_subdiv,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_GroupsubdivYear.csv",sep="") )
        if(BY_LREG & BYGUILD) write.csv(bio_grp_subdiv,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_GuildsubdivYear.csv",sep="") )
        
        write.csv(bio_spp_area,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_AllSpecies_Year.csv",sep="") )
        if(BYGUILD) write.csv(bio_spp_area,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_AllSpecies_Guild_Year.csv",sep="") )
        write.csv(bio_by_area,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_Year.csv",sep="") )
        if(BYGUILD) write.csv(bio_grp_area,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_Guild_Year.csv",sep="") )
        if(!BYGUILD) write.csv(bio_grp_area,file = paste(paste(FILENAM,SP,sep=''), "_Surv_biotonnes_Group_Year.csv",sep="") )
        
        #plot biomass by area once
        if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD){
          windows(width=8, height=8); par(mfrow=c(2,1))
          plot(bio_by_area$Year, bio_by_area$CatCatchWgtSwept/1000,type="l",
               xlab="Year",ylab="Surveyed biomass with Q correction (kt)",main=paste(survey,SP) )
          points(bio_by_area$Year, bio_by_area$CatCatchWgtSwept/1000,pch=19)
          
          plot(bio_by_area$Year, bio_by_area$CatCatchWgtSwept_beforeQmult/1000,type="l",
               xlab="Year",ylab="Surveyed biomass (kt)",main=paste(survey,SP) )
          points(bio_by_area$Year, bio_by_area$CatCatchWgtSwept_beforeQmult/1000,pch=19)
          
        } else {
          windows(width=8, height=8)
          par(mfrow=c(2,1))
          plot(bio_by_area$Year, bio_by_area$CatCatchWgtSwept/1000,type="l",
               xlab="Year",ylab="Surveyed biomass (kt)",main=paste(survey,SP) )
          points(bio_by_area$Year, bio_by_area$CatCatchWgtSwept/1000,pch=19)
          
          plot(bio_by_area$Year, bio_by_area$Abund_N_Swept/1000,type="l",
               xlab="Year",ylab="Surveyed abundance",main=paste(survey,SP) )
          points(bio_by_area$Year, bio_by_area$Abund_N_Swept/1000,pch=19)
          
        }
        
        savePlot(filename= paste(FILENAM,SP,"BIO.bmp",sep='_'),type="bmp")
        dev.off()
      }#plot(bio_by_area[bio_by_area$SpeciesSciName=="Clupea harengus",2:1])
      
      #plot(bio_by_area[bio_by_area$SpeciesSciName=="Clupea harengus",2], bio_by_area[bio_by_area$SpeciesSciName=="Clupea harengus",1]/2000,main="000t, assume half spawning female",type='b')
      
      if(BY_LREG){
        bio_by_subdiv$L_REGName <- as.factor(bio_by_subdiv$L_REG)
        NL<- nlevels(bio_by_subdiv$L_REGName)
        PLOTCOLN<- ceiling(sqrt(NL))
        PLOTROWNbio <- ifelse(NL==2,1,PLOTCOLN)
        PLOTROWNbio <- ifelse(NL==5 | NL==6,2,PLOTCOLN)
        #windows(width=8*PLOTCOLN, height=4*PLOTROWNbio)
        
        #plot biomass by area once
        if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD){
          bmp(filename= paste(FILENAM,SP,"BIOstrata_beforeQmult.bmp",sep='_'))
          xy<- xyplot(data=bio_by_subdiv, CatCatchWgtSwept_beforeQmult~Year | L_REGName,type="b",
                      xlab="Year",ylab="Surveyed biomass (tonnes)",main=paste(survey,SP) )
          print(xy)
          dev.off()
          
          bmp(filename= paste(FILENAM,SP,"BIOstrata_withQmult.bmp",sep='_'))
          xy<- xyplot(data=bio_by_subdiv, CatCatchWgtSwept~Year | L_REGName,type="b",
                      xlab="Year",ylab="Surveyed biomass with Q correction (tonnes)",main=paste(survey,SP) )
          print(xy)
          dev.off()
          
        } else {
          bmp(filename= paste(FILENAM,SP,"BIOstrata.bmp",sep='_'))
          xy<- xyplot(data=bio_by_subdiv, CatCatchWgtSwept~Year | L_REGName,type="b",
                      xlab="Year",ylab="Surveyed biomass (tonnes)",main=paste(survey,SP) )
          print(xy)
          dev.off()
          
          bmp(filename= paste(FILENAM,SP,"ABUNDstrata.bmp",sep='_'))
          xy<- xyplot(data=bio_by_subdiv, Abund_N_Swept~Year | L_REGName,type="b",
                      xlab="Year",ylab="Surveyed abundance",main=paste(survey,SP) )
          print(xy)
          dev.off()
        }
        
        
        bio_by_subdiv <- bio_by_subdiv[,-which(names(bio_by_subdiv) == "L_REGName")]
      }
      
    } #end species loop
  } #END if(BY_LREG | BY_SREG)
  
  
  
  #indicators by species dem pel groups    h(species_bioL_by_area_DEMPEL[species_bioL_by_area_DEMPEL$DEMPEL=='PEL',])
  FILENAM_DEMPEL <- FILENAM # copy as overwrite later
  #meanLD_bio_by_area = meanLD_bio_by_area[!is.na(meanLD_bio_by_area$DEMPEL),]
  
  SPECIES_LOOP <- SPECIES #DEC 2021 TO AVOID HAVING TO RUN CODE TWICE
  if(SPECIES=="ALL") SPECIES_LOOP<-c("DEM","PEL")
  for(SP in SPECIES_LOOP){ # SP<-"ALL" 
    print(SP)
    FILENAM <-  paste(FILENAM_DEMPEL,SP,sep='') 
    #meanLD_bio_by_area# is species_bio_by_area_DEMPEL but raised by area to lowest S_REG
    if(BY_LREG | BY_SREG){
      if(SP == "ALL"){ species_bio_by_area <- meanLD_bio_by_area
      } else { species_bio_by_area <- meanLD_bio_by_area[meanLD_bio_by_area$DEMPEL==SP,]; }
    } else { #not raised by area above
      if(SP == "ALL"){ species_bio_by_area <- species_bio_by_area_DEMPEL
      } else { species_bio_by_area <- species_bio_by_area_DEMPEL[species_bio_by_area_DEMPEL$DEMPEL==SP,]; }
    }
    if(nrow(species_bio_by_area)==0){ print(paste("no species_bio_by_area data for ",SP,sep=''));  break}
    if( length( unique(species_bio_by_area$SpeciesSciName) )<5){ if( length( unique(species_bio_by_area$SpeciesSciName) )<5) print(paste("<5 species recorded in group ",SP,sep='')); print(paste("<5 species recorded in group ",SP,sep=''));  break}
    # include correction for area of strata here so have CPUE_estimates * area of S_REG (or L_REG if lowest level)
    if(BY_LREG & BY_SREG){ # merge in scaling factor
      species_bio_by_area <- merge(x=species_bio_by_area, y=area_by_subdiv, by.x = c("Year","L_REG"),by.y = c("Year","L_REG"),all.x=T)
      if(AREASCALE){
        species_bio_by_area$CatCatchWgtSwept <- species_bio_by_area$CatCatchWgtSwept*species_bio_by_area$scale
        if(CATCHABILITY_COR_WALKER | CATCHABILITY_COR_MOD) species_bio_by_area$CatCatchWgtSwept_beforeQmult <- species_bio_by_area$CatCatchWgtSwept_beforeQmult*species_bio_by_area$scale
      }
    }# corrected for any change in sampling between L_REG, but not scaled up to include missing L_REG
    
    #just elasmos etc
    if(!is.null(GROUP)){
      if(SP=="PEL" & GROUP!="Other" & !BYICESGROUP & !BYGUILD){ next; 
      } else { species_bio_by_area <- species_bio_by_area[species_bio_by_area$Group==GROUP,]; } 
      if(SP!= "ALL" & nrow( species_bio_by_area[species_bio_by_area$DEMPEL==SP,])==0 ) next
    }
    
    
    
    #Large Fish Indicator
    if(!BY_SREG) numS_REG_by_sea <- 0*numhaulsyr
    if(LFI & SP!="PEL"){   IND_LFI <- INDfn_LFI( species_bio_by_area=species_bio_by_area, SP=SP,
                                                 numhaulsyr=numhaulsyr,numS_REG_by_sea=numS_REG_by_sea, WRITE=WRITE, 
                                                 FILENAM=paste(FILENAM,GROUP,LFI_THRESHOLD,sep="_"),BY_LREG=BY_LREG,LFI_THRESHOLD = LFI_THRESHOLD)
    } else { IND_LFI <-NULL }
    #Mean Length cm by S_REGa and year 
    if(MeanL) IND_MeanL <- INDfn_MeanL(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=paste(FILENAM,GROUP,sep="_"),BY_SREG=BY_SREG,BY_LREG=BY_LREG,SP=SP)
    
    #MaxL by rectangle and year 
    if(MaxL) IND_MaxL <- INDfn_Mtrait(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=paste(FILENAM,GROUP,sep="_"),BY_SREG=BY_SREG,BY_LREG=BY_LREG,SP=SP, IND.NAM = "MaxL")
    if(Loo) IND_Loo <- INDfn_Mtrait(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=paste(FILENAM,GROUP,sep="_"),BY_SREG=BY_SREG,BY_LREG=BY_LREG,SP=SP, IND.NAM = "Loo")
    if(Lm) IND_Lm <- INDfn_Mtrait(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=paste(FILENAM,GROUP,sep="_"),BY_SREG=BY_SREG,BY_LREG=BY_LREG,SP=SP, IND.NAM = "Lm")
    
    #TL by rectangle and year 
    if(MEANTL) IND_MeanTL <- INDfn_MeanTL(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=paste(FILENAM,GROUP,sep="_"),BY_SREG=BY_SREG,BY_LREG=BY_LREG,SP=SP)
    
    #Geometric mean length (Typical Length cm)
    # Weight in kg raised to 60 min haul 
    # log[ length(cm) ]
    if(TyL_GeoM) IND_TyL_GeoM <-INDfn_TyL_GeoM(species_bio_by_area=species_bio_by_area, WRITE=WRITE, FILENAM=paste(FILENAM,GROUP,sep="_"),BY_SREG=BY_SREG,BY_LREG=BY_LREG,SP=SP,TyL_SPECIES=TyL_SPECIES)
    
    #noddy way to make a list with all indicators for pelagic and demersals
    if(SP=="ALL"){      
      if(TyL_GeoM) TyL.cm.sea_all <- IND_TyL_GeoM;
      if(MeanL) FishLength_cmsea_all<-IND_MeanL; 
      if(MaxL)  MaxLsea_all<-IND_MaxL; 
      if(Loo)  Loosea_all<-IND_Loo; 
      if(Lm)   Lmsea_all<-IND_Lm; 
      if(MEANTL) TLsea_all<-IND_MeanTL; 
      if(LFI & !is.null(IND_LFI[[1]]) ){
        LFIbind_all <-   IND_LFI[[1]]       #all years
        rownames(LFIbind_all) <- LFIbind_all$Year
        # add ncol(LFI_by_sub_all) cols
        if(!is.null(IND_LFI[[2]])){ 
          LFI_by_sub_all <- IND_LFI[[2]]; 
          for(i in 1:ncol(LFI_by_sub_all)) LFIbind_all <- cbind(LFIbind_all,NA)
          names(LFIbind_all)[(ncol(LFIbind_all)-ncol(LFI_by_sub_all)+1):ncol(LFIbind_all)]  <- colnames(LFI_by_sub_all) 
          
          LFIbind_all[ rownames(LFIbind_all) %in% rownames(LFI_by_sub_all),
                       (ncol(LFIbind_all)-ncol(LFI_by_sub_all)+1):ncol(LFIbind_all)] <- (LFI_by_sub_all)
        }
      }
    }
    if(SP=="DEM"){      
      if(TyL_GeoM) TyL.cm.sea_dem <- IND_TyL_GeoM;
      if(MeanL) FishLength_cmsea_dem<-IND_MeanL; 
      if(MaxL)  MaxLsea_dem<-IND_MaxL;  
      if(Loo)  Loosea_dem<-IND_Loo; 
      if(Lm)   Lmsea_dem<-IND_Lm; 
      if(MEANTL) TLsea_dem<-IND_MeanTL; 
      if(LFI & !is.null(IND_LFI[[1]]) ){
        LFIbind_dem <-   IND_LFI[[1]]       #all years
        rownames(LFIbind_dem) <- LFIbind_dem$Year
        # add ncol(LFI_by_sub_dem) cols
        if(!is.null(IND_LFI[[2]])){ 
          LFI_by_sub_dem <- IND_LFI[[2]]; 
          for(i in 1:ncol(LFI_by_sub_dem)) LFIbind_dem <- cbind(LFIbind_dem,NA)
          names(LFIbind_dem)[(ncol(LFIbind_dem)-ncol(LFI_by_sub_dem)+1):ncol(LFIbind_dem)]  <- colnames(LFI_by_sub_dem) 
          
          LFIbind_dem[ rownames(LFIbind_dem) %in% rownames(LFI_by_sub_dem),
                       (ncol(LFIbind_dem)-ncol(LFI_by_sub_dem)+1):ncol(LFIbind_dem)] <- (LFI_by_sub_dem)
        }
      }
    }
    if(SP=="PEL"){
      
      if(TyL_GeoM)    TyL.cm.sea_pel<-IND_TyL_GeoM
      if(MeanL)  FishLength_cmsea_pel<-IND_MeanL
      if(MaxL)   MaxLsea_pel <- IND_MaxL 
      if(Loo)  Loosea_pel<-IND_Loo; 
      if(Lm)   Lmsea_pel<-IND_Lm; 
      if(MEANTL) TLsea_pel   <- IND_MeanTL
      if(LFI & !is.null(IND_LFI[[1]]) ){
        LFIbind_pel <-  IND_LFI[[1]]     #all years
        rownames(LFIbind_pel) <- LFIbind_pel$Year
        # add ncol(LFI_by_sub_dem) cols
        if(!is.null(IND_LFI[[2]])){ 
          LFI_by_sub_pel<-IND_LFI[[2]]; 
          for(i in 1:ncol(LFI_by_sub_pel)) LFIbind_pel <- cbind(LFIbind_pel,NA)
          names(LFIbind_pel)[(ncol(LFIbind_pel)-ncol(LFI_by_sub_pel)+1):ncol(LFIbind_pel)]  <- colnames(LFI_by_sub_pel)              
          LFIbind_pel[rownames(LFIbind_pel) %in% rownames(LFI_by_sub_pel),(ncol(LFIbind_pel)-ncol(LFI_by_sub_pel)+1):ncol(LFIbind_pel)] <- (LFI_by_sub_pel)              
        }
      }
    }
  }#species set loop 
  #        IND_OUT<-
  if(!BOOTSTRAP) return(list(                                               
    LFI_by_sub_all = LFIbind_all, TyL.cm.sea_all = TyL.cm.sea_all, FishLength_cmsea_all =FishLength_cmsea_all, MaxLsea_all =MaxLsea_all, Loosea_all =Loosea_all, Lmsea_all =Lmsea_all, TLsea_all =TLsea_all 
    , LFI_by_sub_dem = LFIbind_dem, TyL.cm.sea_dem = TyL.cm.sea_dem, FishLength_cmsea_dem =FishLength_cmsea_dem, MaxLsea_dem =MaxLsea_dem, Loosea_dem =Loosea_dem, Lmsea_dem =Lmsea_dem, TLsea_dem =TLsea_dem 
    , LFI_by_sub_pel = LFIbind_pel, TyL.cm.sea_pel = TyL.cm.sea_pel, FishLength_cmsea_pel =FishLength_cmsea_pel, MaxLsea_pel =MaxLsea_pel, Loosea_pel =Loosea_pel, Lmsea_pel =Lmsea_pel, TLsea_pel =TLsea_pel
    , species_bio_by_area=species_bio_by_area_DEMPEL, 
    numhauls=numhauls, numhaulsyr=numhaulsyr, numhaulsBYS_REG=numhaulsBYS_REG, numhaulsBYsubdiv=numhaulsBYsubdiv) )
  #dont save everything in bootstrap
  if(BOOTSTRAP) return(list(                    
    LFI_by_sub_all = LFIbind_all, TyL.cm.sea_all = TyL.cm.sea_all, FishLength_cmsea_all =FishLength_cmsea_all, MaxLsea_all =MaxLsea_all, Loosea_all =Loosea_all, Lmsea_all =Lmsea_all, TLsea_all =TLsea_all 
    , LFI_by_sub_dem = LFIbind_dem, TyL.cm.sea_dem = TyL.cm.sea_dem, FishLength_cmsea_dem =FishLength_cmsea_dem, MaxLsea_dem =MaxLsea_dem, Loosea_dem =Loosea_dem, Lmsea_dem =Lmsea_dem, TLsea_dem =TLsea_dem
    , LFI_by_sub_pel = LFIbind_pel, TyL.cm.sea_pel = TyL.cm.sea_pel, FishLength_cmsea_pel =FishLength_cmsea_pel, MaxLsea_pel =MaxLsea_pel, Loosea_pel =Loosea_pel, Lmsea_pel =Lmsea_pel, TLsea_pel =TLsea_pel 
  ) )
}
