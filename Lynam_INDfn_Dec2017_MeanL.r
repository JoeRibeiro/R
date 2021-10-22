INDfn_MeanL <- function(species_bio_by_area, WRITE=F, FILENAM="",SAMP_STRAT=T,BYSUBDIV=F,SP=SP){
  #may 2018 added if PLOTTED to avoid error if no plot#  
  species_bio_by_area$FishLength_cm.Swept <- species_bio_by_area$CatCatchWgtSwept * species_bio_by_area$FishLength_cm  #NAs cascade
  FishLength_cmsubdiv <- NULL
  #mean len by rectangle in 2 steps (with sub-regions lose species etc)
  if(SAMP_STRAT){
    FACT <- c("Year","sampstrat")
    if(BYSUBDIV) FACT <- c(FACT,"subdiv")
    FishLength_cmbio_by_samp <- tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$FishLength_cm),],
                                      datacols=c("CatCatchWgtSwept","FishLength_cm.Swept","FishLength_cm"),
                                      factorcols=FACT,  
                                      sum,c("CatCatchWgtSwept","FishLength_cm.Swept","FishLength_cm")
                                      )  #step 1 #if BYSUBDIV is F then subdiv is just NA in fn
    #Mean Length by sampstrata (e.g. rectangle)
    FishLength_cmbio_by_samp$FishLength_cm <- FishLength_cmbio_by_samp$FishLength_cm.Swept / FishLength_cmbio_by_samp$CatCatchWgtSwept   #step 2
    if(WRITE) write.table(FishLength_cmbio_by_samp,paste(FILENAM,'Len.cm_swept_sampstrat.csv',sep="_"),sep=',',append=F,col.names=T,row.names=F)              #output by rect
  }
  
  #Mean Length by subdiv
  if(BYSUBDIV){
    FishLength_cmsubdiv <- tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$FishLength_cm),],
                                  datacols=c("CatCatchWgtSwept","FishLength_cm.Swept","FishLength_cm"),
                                  factorcols=c("Year","subdiv"),sum,c("CatCatchWgtSwept","FishLength_cm.Swept","FishLength_cm")
                                  )  #step 1
    FishLength_cmsubdiv$FishLength_cm <- FishLength_cmsubdiv$FishLength_cm.Swept / FishLength_cmsubdiv$CatCatchWgtSwept   #step 2
    #reshape
    FishLength_cmsubdiv <- tapply(FishLength_cmsubdiv$FishLength_cm,list(FishLength_cmsubdiv$Year, FishLength_cmsubdiv$subdiv), FUN=mean, na.rm=T)
  }
  
  #Mean Length for sea by year
  # corrected for any change in sampling between subdiv
  FishLength_cmsea <- tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area$FishLength_cm),],
                                datacols=c("CatCatchWgtSwept","FishLength_cm.Swept","FishLength_cm"),
                                factorcols=c("Year"),sum,c("CatCatchWgtSwept","FishLength_cm.Swept","FishLength_cm")
                                )
  FishLength_cmsea$FishLength_cm <- FishLength_cmsea$FishLength_cm.Swept / FishLength_cmsea$CatCatchWgtSwept   #step 2
  FishLength_cmsea <- data.frame(Year=FishLength_cmsea$Year, FishLength_cm= FishLength_cmsea$FishLength_cm)
  
  
  #Mean Length for sea by year and regions
  if(BYSUBDIV) FishLength_cmsea <- cbind(as.numeric(as.character(FishLength_cmsea[,1])),FishLength_cmsubdiv,FishLength_cmsea[,2])
  colnames(FishLength_cmsea)[1] <- "Year"; colnames(FishLength_cmsea)[ncol(FishLength_cmsea)] <- "sea"
  if(WRITE) write.table(FishLength_cmsea,paste(FILENAM,'Len.cm_swept_sea.csv',sep="_"),sep=',',append=F,col.names=T,row.names=F)
  if(BYSUBDIV) FishLength_cmsea <- FishLength_cmsea[,-1]
  
  #plot code
  if(!BOOTSTRAP & BYSUBDIV & !is.null(FishLength_cmsubdiv)){
    if(nrow(FishLength_cmsubdiv)>3){ 
    YLAB <- "Mean Length (cm)"
    DATA2PLOT <- FishLength_cmsubdiv
    
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
      INDPLOTFN(DATA2PLOT=DATA2PLOT[!is.na(DATA2PLOT[,n]),n], BOOTDATA2PLOT=NULL,YRS=YRSPLOT,YLIM=YLIM,TITLE=TITLE, GAMMOD=NA, YLAB=YLAB
                ,BOOTSTRAP=F, ADDBOOTTREND=F, ADDBOOTTREND_CI=F, ADDBOOT_ERRBAR=F, ADDGAM=F, BEST_AND_BOOT=F)
    }
    
    if(SP=="DEM") TITLE <-paste(survey, "Demersal fish", sep=" ")
    if(SP=="PEL") TITLE <-paste(survey, "Pelagic fish", sep=" ")
    if(SP=="ALL") TITLE <-paste(survey, "All fish", sep=" ")
      if(PLOTTED & dev.cur()>1 ){ mtext(TITLE,line=0.5,outer=T)
      savePlot(filename= paste(FILENAM,"_",SP,"_MeanL.bmp",sep=''),type="bmp")
      dev.off() 
      }
    }
  }
  return(FishLength_cmsea)
}