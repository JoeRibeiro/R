INDfn_Mtrait <- function(species_bio_by_area, WRITE=F, FILENAM="",SAMP_STRAT=T, BYSUBDIV=F,SP=SP, IND.NAM="MaxL"){#IND.NAM<-"MaxL"
  #Dec 2021 added as.numeric since have character for Loo and Lm
  IND.NAM.SWEPT<-paste(IND.NAM,"Swept",sep="")
  #intermediate calc for MML ## jul2020 added line 69 as had poorly sampled strata in SWBEAM
  species_bio_by_area[,IND.NAM] <- as.numeric(species_bio_by_area[,IND.NAM]) 
  
  species_bio_by_area[,IND.NAM.SWEPT] <- species_bio_by_area$CatCatchWgtSwept * species_bio_by_area[,IND.NAM]  #NAs cascade
  INDbio_by_subdiv<-NULL
  
  if(SAMP_STRAT){
    FACT <- c("Year","sampstrat")
    if(BYSUBDIV) FACT <- c(FACT,"subdiv","scale") #Jul2020 added scale
    
    #sum by sampstrat (rectangle) and by sub-divisions
    INDbio_by_samp <- tapply.ID( df = species_bio_by_area[!is.na(species_bio_by_area[,IND.NAM]),], 
                                datacols = c("CatCatchWgtSwept",IND.NAM.SWEPT,IND.NAM),
                                factorcols = FACT, sum, c("CatCatchWgtSwept",IND.NAM.SWEPT,IND.NAM))
  
    # MML by sampstrat (rectangle)
    INDbio_by_samp[,IND.NAM] <- INDbio_by_samp[,IND.NAM.SWEPT] / INDbio_by_samp$CatCatchWgtSwept
    if(WRITE) write.table(INDbio_by_samp,paste(FILENAM,IND.NAM,'swept_sampstrat.csv',sep="_"), sep=',', append=F, col.names=T, row.names=F)
  }
  
  # MML by subdiv
  if(BYSUBDIV){
    INDbio_by_subdiv <- tapply.ID( df = species_bio_by_area[!is.na(species_bio_by_area[,IND.NAM]),], 
                                    datacols = c("CatCatchWgtSwept",IND.NAM.SWEPT,IND.NAM),
                                    factorcols = c("Year","subdiv"),sum,c("CatCatchWgtSwept",IND.NAM.SWEPT,IND.NAM))
    
    INDbio_by_subdiv[,IND.NAM] <- INDbio_by_subdiv[,IND.NAM.SWEPT] / INDbio_by_subdiv$CatCatchWgtSwept
    
    if(WRITE) write.table(INDbio_by_subdiv,paste(FILENAM,IND.NAM,'bio_by_subdiv.csv',sep="_"), sep=',', append=F, col.names=T, row.names=F)
    #reshape
    INDbio_by_subdiv <- tapply(INDbio_by_subdiv[,IND.NAM],list(INDbio_by_subdiv$Year, INDbio_by_subdiv$subdiv), FUN=mean, na.rm=T)
  }
  
  # MML for sub-regional sea by year 
  # corrected for any change in sampling between subdiv, but not scaled up to include missing subdiv
  
  suppressWarnings( #NAs introduced by coercion since some MaxL are NA
    INDsea <- tapply.ID(df=species_bio_by_area[!is.na(species_bio_by_area[,IND.NAM]),],
                       datacols=c("CatCatchWgtSwept",IND.NAM.SWEPT,IND.NAM),
                       factorcols=c("Year"),sum,c("CatCatchWgtSwept",IND.NAM.SWEPT,IND.NAM))
  )
  INDsea[,IND.NAM] <- INDsea[,IND.NAM.SWEPT] / INDsea$CatCatchWgtSwept
  
  INDsea <- data.frame(Year=INDsea$Year, IND=INDsea[,IND.NAM])                
  
  
  #combine for output# mean MML for sub-regional sea by year and now with subdivisions (sub.reg)
  if(BYSUBDIV) INDsea <- cbind(as.numeric(as.character(INDsea[,1])),  INDbio_by_subdiv, INDsea[,2])
  colnames(INDsea)[1] <- "Year"; colnames(INDsea)[ncol(INDsea)] <- "sea"
  if(WRITE) write.table(INDsea,paste(FILENAM,IND.NAM,'_swept_sea.csv',sep="_"),sep=',',append=F,col.names=T,row.names=F)
  if(BYSUBDIV) INDsea <- INDsea[,-1]

  
  #plot code
  if(!BOOTSTRAP & BYSUBDIV & !is.null(INDbio_by_subdiv)){
    if(nrow(INDbio_by_subdiv)>3){ 
    YLAB <- "Length (cm)"
    DATA2PLOT <- INDbio_by_subdiv
    
    YLIM<-c(floor(min(DATA2PLOT,na.rm=T)*.95), ceiling(max(DATA2PLOT,na.rm=T)*1.05) )
    N<-ncol(DATA2PLOT)
    winwidth <- ceiling(sqrt(N))
    winhgt <- ifelse(N==2,1,winwidth)
    windows(width=winwidth*8, height=winhgt*8)
    par(mfrow=c(winhgt,winwidth),mar=c(2,4,2,2),oma=c(1,1,3,1))
    for(n in 1:N){
      if(length(DATA2PLOT[!is.na(DATA2PLOT[,n]),n])==1) next #added jul2020
      TITLE <- colnames(DATA2PLOT)[n]
      YRSPLOT <- an(names(DATA2PLOT[!is.na(DATA2PLOT[,n]),n]))
      INDPLOTFN(DATA2PLOT=DATA2PLOT[!is.na(DATA2PLOT[,n]),n], BOOTDATA2PLOT=NULL,YRS=YRSPLOT,YLIM=YLIM,TITLE=TITLE, 
                GAMMOD=NA, YLAB=YLAB
                ,BOOTSTRAP=F, ADDBOOTTREND=F, ADDBOOTTREND_CI=F, ADDBOOT_ERRBAR=F, ADDGAM=F, BEST_AND_BOOT=T)
    }
    
    if(SP=="DEM") TITLE <-paste(survey, "Demersal fish", sep=" ")
    if(SP=="PEL") TITLE <-paste(survey, "Pelagic fish", sep=" ")
    if(SP=="ALL") TITLE <-paste(survey, "All fish", sep=" ")
    mtext(TITLE,line=0.5,outer=T)
    
    savePlot(filename= paste(FILENAM,IND.NAM,".bmp",sep=''),type="bmp")
    dev.off()  
      }
  }
  return(INDsea)
  
}