INDfn_TyL_GeoM <- function(species_bio_by_area, WRITE=F, FILENAM="",SAMP_STRAT=T,BYSUBDIV=F,SP="DEM",TyL_SPECIES=F){
  #may 2018 added if PLOTTED to avoid error if no plot#  
  #for mean calc: product catch*logL in each rect_stsq       
  species_bio_by_area$FishLength_cm.Swept <- species_bio_by_area$CatCatchWgtSwept * log(species_bio_by_area$FishLength_cm) #NAs cascade
  TyL.cm.subdiv_est<-NULL
  
  if(SAMP_STRAT){
    FACT <- c("Year","sampstrat")
    if(BYSUBDIV) FACT <- c(FACT,"subdiv","scale") #Jul2020 added scale
    
    #calc typ length by rect in 2 steps with stsqs and sub-regions
    TyL.cm_by_samp <- tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$FishLength_cm),],
                              datacols=c("CatCatchWgtSwept","FishLength_cm.Swept"),
                              factorcols=FACT,sum,c("CatCatchWgtSwept","FishLength_cm.Swept"))    #step 1
  
          #typ len-region (by rectangle) #[,1:4] CatCatchWgtSwept	Year	ICESStSq	subdiv
          TyL.cm_by_samp$FishLength_cm.SweptEst <- (TyL.cm_by_samp$FishLength_cm.Swept / TyL.cm_by_samp$CatCatchWgtSwept)             #step 2
  
          #EXP output cm #Jul2020 replaced c(1:4) by '-ncol(TyL.cm_by_samp)'
          TyL.cm_by_sampOUT <-cbind(TyL.cm_by_samp[,-ncol(TyL.cm_by_samp)], TyL.cm = exp(TyL.cm_by_samp[,"FishLength_cm.SweptEst"]) )
          if(WRITE) write.table(TyL.cm_by_sampOUT,
                        paste(FILENAM,'TyL.cm_swept_bysampstrat.csv',sep="_"),sep=',',append=F,col.names=T,row.names=F)         #output by rect
  }
  
  if(BYSUBDIV){
    #and for subdiv by year # here correct weighted average on log scale
    TyL.cm.subdiv <-tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$FishLength_cm),],
                         datacols=c("CatCatchWgtSwept","FishLength_cm.Swept"),
                         factorcols=c("Year","subdiv"),sum,c("CatCatchWgtSwept","FishLength_cm.Swept")) 
          #now back to natural length scale #i.e. the geometric mean which is the median of a lognormal distn
          TyL.cm.subdiv$TyL_subdiv_mean <- exp(TyL.cm.subdiv$FishLength_cm.Swept/TyL.cm.subdiv$CatCatchWgtSwept)
          #and reshape since have one value per year and subdiv combination
          TyL.cm.subdiv_est <- (tapply(TyL.cm.subdiv$TyL_subdiv_mean,list(TyL.cm.subdiv$Year, TyL.cm.subdiv$subdiv), FUN=mean, na.rm=T))
  }
  
  #and for sea by year
  # corrected for any change in sampling between subdiv but not for missing subdiv
  
  # here correct weighted average on log scale
  TyL.cm.sea <-tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$FishLength_cm),],
                         datacols=c("CatCatchWgtSwept","FishLength_cm.Swept"),
                         factorcols=c("Year"),sum,c("CatCatchWgtSwept","FishLength_cm.Swept")) 
  #now back to natural length scale #i.e. the geometric mean which is the median of a lognormal distn
  TyL.cm.sea$TyL_subdiv_mean <- exp(TyL.cm.sea$FishLength_cm.Swept/TyL.cm.sea$CatCatchWgtSwept)
  
  #and for sea by year and regions, average across rects then EXP
  if(BYSUBDIV) TyL.cm.sea <- cbind(as.numeric(as.character(TyL.cm.sea[,"Year"])),TyL.cm.subdiv_est, TyL.cm.sea[,"TyL_subdiv_mean"])
  colnames(TyL.cm.sea)[ncol(TyL.cm.sea)] <- "sea"
  if(WRITE) write.table(TyL.cm.sea,paste(FILENAM,'TyL.cm_swept_subdiv_sea.csv',sep="_"),sep=',',append=F,col.names=T,row.names=F)
  if(BYSUBDIV) TyL.cm.sea <- TyL.cm.sea[,-1]
  
  
  ##
  if(TyL_SPECIES){
    #and for sea by year as above but output by species
    
    # here correct weighted average on log scale
    TyL.cm.sea.SP <-tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$FishLength_cm),],
                           datacols=c("CatCatchWgtSwept","FishLength_cm.Swept"),
                           factorcols=c("Year","SpeciesSciName"),sum,c("CatCatchWgtSwept","FishLength_cm.Swept")) 
    #now back to natural length scale #i.e. the geometric mean which is the median of a lognormal distn
    TyL.cm.sea.SP$TyL_subdiv_mean <- exp(TyL.cm.sea.SP$FishLength_cm.Swept/TyL.cm.sea.SP$CatCatchWgtSwept)
    
    if(WRITE) write.table(TyL.cm.sea.SP,paste(FILENAM,'TyL.cm_SP.swept_subdiv_sea.csv',sep="_"),sep=',',append=F,col.names=T,row.names=F)
  }
  ##
  
  
  #plot code
  if(!BOOTSTRAP & BYSUBDIV & !is.null(TyL.cm.subdiv_est)){
    if(nrow(TyL.cm.subdiv_est)>3  & length(u(ac(TyL.cm.subdiv_est)))>1 ){ 
      
    YLAB <- "TyL (cm)"
    
    DATA2PLOT <- TyL.cm.subdiv_est 
    YLIM<-c(floor(min(DATA2PLOT,na.rm=T)*.95), ceiling(max(DATA2PLOT,na.rm=T)*1.05) )
    N<-ncol(DATA2PLOT)
    winwidth <- ceiling(sqrt(N))
    winhgt <- ifelse(N==2,1,winwidth)
    windows(width=winwidth*8, height=winhgt*8)
    par(mfrow=c(winhgt,winwidth),mar=c(2,4,2,2),oma=c(1,1,3,1))
    
    for(n in 1:N){
      if(length(u(ac(DATA2PLOT[!is.na(DATA2PLOT[,n]),n])))<4) next 
      PLOTTED<-T
      
      TITLE <- colnames(DATA2PLOT)[n]
      YRSPLOT <- an(names(DATA2PLOT[!is.na(DATA2PLOT[,n]),n]))
      #if(length(YRSPLOT)==0) next # blue2_lam missing from "CSScoOT4"!
      INDPLOTFN(DATA2PLOT=DATA2PLOT[!is.na(DATA2PLOT[,n]),n], BOOTDATA2PLOT=NULL,YRS=YRSPLOT,YLIM=YLIM,TITLE=TITLE, GAMMOD=NA, YLAB=YLAB
                ,BOOTSTRAP=F, ADDBOOTTREND=F, ADDBOOTTREND_CI=F, ADDBOOT_ERRBAR=F, ADDGAM=F, BEST_AND_BOOT=T)
    }
    
    if(SP=="DEM") TITLE <-paste(survey, "Demersal fish", sep=" ")
    if(SP=="PEL") TITLE <-paste(survey, "Pelagic fish", sep=" ")
    if(SP=="ALL") TITLE <-paste(survey, "All fish", sep=" ")
    if( PLOTTED & dev.cur()>1 ){ mtext(TITLE,line=0.5,outer=T)
    savePlot(filename= paste(FILENAM,"TyL.bmp",sep=''),type="bmp")
    }
    dev.off()  
    }
  }
  return(TyL.cm.sea)
}



