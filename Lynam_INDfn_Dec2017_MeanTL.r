INDfn_MeanTL <- function(species_bio_by_area, WRITE=F, FILENAM="",SAMP_STRAT=T,BYSUBDIV=F,SP=SP){
  #may 2018 added if PLOTTED to avoid error if no plot#
  species_bio_by_area$TL.Swept <- species_bio_by_area$CatCatchWgtSwept * species_bio_by_area$TL  #NAs cascade
  TLsubdiv<-NULL
  
  if(SAMP_STRAT){
    FACT <- c("Year","sampstrat")
    if(BYSUBDIV) FACT <- c(FACT,"subdiv")
    #sum by rectangles and subdivisions
    TLbio_by_samp <- tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$TL),], 
                               datacols=c("CatCatchWgtSwept", "TL.Swept", "TL"),
                               factorcols=FACT,sum,c("CatCatchWgtSwept","TL.Swept","TL"))
    #MTL by rectangle
    TLbio_by_samp$TL.Swept <- TLbio_by_samp$TL.Swept / TLbio_by_samp$CatCatchWgtSwept
    if(WRITE) write.table(TLbio_by_samp,paste(FILENAM,'TL_swept_sampstrat.csv',sep="_"),sep=',',append=F,col.names=T,row.names=F)
  }
  if(BYSUBDIV){
    #MTL for subdivision by year
    TLsubdiv <- tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$TL),], 
                          datacols=c("CatCatchWgtSwept", "TL.Swept", "TL"),
                          factorcols=c("Year","subdiv"),sum,c("CatCatchWgtSwept","TL.Swept","TL"))
    TLsubdiv$TL <- TLsubdiv$TL.Swept / TLsubdiv$CatCatchWgtSwept
    if(WRITE) write.table(TLsubdiv,paste(FILENAM,'TLbio_by_subdiv.csv',sep="_"),sep=',',append=F,col.names=T,row.names=F)
    #reshape
    TLsubdiv <- tapply(TLsubdiv$TL,list(TLsubdiv$Year, TLsubdiv$subdiv), FUN=mean, na.rm=T)
  }
  
  #MTL for sea by year
  # corrected for any change in sampling between subdiv
  
  suppressWarnings( #NAs introduced by coercion since some TL are NA
    TLsea <- tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$TL),],  
                     datacols=c("CatCatchWgtSwept","TL.Swept","TL"),
                     factorcols=c("Year"),sum,c("CatCatchWgtSwept","TL.Swept","TL"))
  )
  TLsea$TL <- TLsea$TL.Swept / TLsea$CatCatchWgtSwept
  TLsea <- data.frame(Year= TLsea$Year, TL=as.numeric(TLsea$TL))
  
  
  #MTL for sea by year and sub-divisions
  if(BYSUBDIV) TLsea <- cbind(as.numeric(as.character(TLsea[,1])), TLsubdiv, TLsea[,2])
  colnames(TLsea)[1] <- "Year"; colnames(TLsea)[ncol(TLsea)] <- "sea"
  if(WRITE) write.table(TLsea,paste(FILENAM,'TL_swept_sea.csv',sep="_"),sep=',',append=F,col.names=T,row.names=F)
  if(BYSUBDIV) TLsea <- TLsea[,-1]

  
  #plot code
  if(!BOOTSTRAP & BYSUBDIV & !is.null(TLsubdiv)){
    if(nrow(TLsubdiv)>3 & length(u(ac(TLsubdiv)))>1){ 
    
    YLAB <- "MTL"
    DATA2PLOT <- TLsubdiv
    
    YLIM<-c( (min(DATA2PLOT,na.rm=T)*.95), (max(DATA2PLOT,na.rm=T)*1.05) )
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
      INDPLOTFN(DATA2PLOT=DATA2PLOT[!is.na(DATA2PLOT[,n]),n], BOOTDATA2PLOT=NULL,YRS=YRSPLOT,YLIM=YLIM,TITLE=TITLE, GAMMOD=NA, YLAB=YLAB
                ,BOOTSTRAP=F, ADDBOOTTREND=F, ADDBOOTTREND_CI=F, ADDBOOT_ERRBAR=F, ADDGAM=F, BEST_AND_BOOT=F)
    }
    
    if(SP=="DEM") TITLE <-paste(survey, "Demersal fish", sep=" ")
    if(SP=="PEL") TITLE <-paste(survey, "Pelagic fish", sep=" ")
    if(SP=="ALL") TITLE <-paste(survey, "All fish", sep=" ")
    if(PLOTTED & dev.cur()>1 ){ mtext(TITLE,line=0.5,outer=T)
      savePlot(filename= paste(FILENAM,"_",SP,"_MTL.bmp",sep=''),type="bmp")
    }
    dev.off()  
    }
  }
  return(TLsea)
}