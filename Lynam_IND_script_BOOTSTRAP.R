# Author: C Lynam, Cefas
# Contact: chris.lynam@cefas.co.uk
# Version: 1 
# Date: May 2020 
#if(BOOTSTRAP){ print("Bootstrap dataset")
  # for surveys with many samples per sampling unit only
  # from best estimate output above retain info on numhaul used and send to bootstrap function 
  NUMHAULS <- IND_OUT$numhauls
  NUMHAULSYR <- IND_OUT$numhaulsyr
  NUMHAULSbySAMPSTRAT <- IND_OUT$numhaulsBYsampstrat
  NUMHAULSbySUBDIV <- IND_OUT$numhaulsBYsubdiv
  #NUMHAULSbySAMPSTRAT[NUMHAULSbySAMPSTRAT$sampstrat=="23E6" & NUMHAULSbySAMPSTRAT$Year==2008,]
  #numhaulsBYsampstrat[numhaulsBYsampstrat$sampstrat=="23E6" & numhaulsBYsampstrat$Year==2008,]
  
  #and save output if WRITE_BOOT <- T
  BOOTOUTPATH <- paste(OUTPATH,"BOOT/",sep='')
  
  #now run the bootstrap B times 
  #ok if  SAMP_STRAT <- T   and   BYSUBDIV <- T    above
  
  for(B in 1:NBOOT){  # B<-1
    print(B)
    
    #bootstrap the dataset
    #SAMP_STRAT <- T and BYSUBDIV <- T   runs ok
    
    if(SAMP_STRAT) boot_hspp <- BOOTfn(DATA=dhspp, NUMHAULS=NUMHAULS, NUMHAULSbySAMPSTRAT=NUMHAULSbySAMPSTRAT, LON_RANGE=1,LAT_RANGE=0.5)
    if(!SAMP_STRAT) boot_hspp <- BOOTfnNoRect(DATA=dhspp, NUMHAULS=NUMHAULS, NUMHAULSbySUBDIV=NUMHAULSbySUBDIV)
    if(WRITE_BOOT){
      if(!dir.exists(BOOTOUTPATH)) dir.create(BOOTOUTPATH)
      save(boot_hspp,file=paste(BOOTOUTPATH,"BOOT_DAT_",B,"_",Sys.Date(),".RData",sep=""))
    }
    #run the Indicator function       #DATA<-boot_hspp$dat_boot; WRITE<-F; FILENAM<-"bootS"; BOOT<-T
    #simple geometric mean log length with averaging of catch rates by rectangles or other spatial strata (subdiv)
    #TyL_GeoM<-T; LFI    <-F; MeanL  <-F; MaxL   <-F; 
    #SPECIES <- c("DEM", "PEL")
    
    #boot_hspp$dat_boot <- boot_hspp$dat_boot[!is.na(boot_hspp$dat_boot$MaxL),]
    
    BOOT<-T #dont save everything in bootstrap
    try(
      BOOT_OUT <-  INDfn(DATA=boot_hspp$dat_boot, WRITE=WRITE_BOOT,SPECIES=SPECIES,FILENAM=FILENAM, BOOT=BOOT, SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV,
                         LFI=LFI,LFI_THRESHOLD=LFI_THRESHOLD, MEANTL=MEANTLs, MaxL=MaxL, Loo=Loo, Lm=Lm, MeanL=MeanL, TyL_GeoM=TyL_GeoM
                         ,TyL_SPECIES=F, BYGUILD=F)
      ,silent=TRUE)
    if(WRITE_BOOT) save(BOOT_OUT,file=paste(BOOTOUTPATH,"BOOT_IND_",B,"_",Sys.Date(),".RData",sep=""))   
    
    if(BYGROUP){
      BOOT_OUT_BYGROUP<-list()
      for(SPGROUP in levels(boot_hspp$Group)){ 
        FILENAM2 <- paste(OUTPATH,survey,"_",format(Sys.time(), "%d%b%Y"),sep="")
        
        if(nrow(boot_hspp[boot_hspp$Group==SPGROUP,])>3){
          try(
            BOOT_OUT_BYGROUP[[SPGROUP]] <- INDfn( DATA=boot_hspp, WRITE=T, SPECIES=SPECIES, GROUP=SPGROUP,
                                                  LFI_THRESHOLD=LFI_THRESHOLD, SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, FILENAM=FILENAM2,
                                                  LFI=LFI, MEANTL=MEANTLs, MaxL=MaxL, Loo=Loo, Lm=Lm, MeanL=MeanL, TyL_GeoM=TyL_GeoM,
                                                  TyL_SPECIES=F, BYGUILD=F)
            ,silent=TRUE)
        }
      }
    }
    
    if(BYICESGROUP){
      BOOT_OUT_BYICESGROUP<-list()
      for(SPGROUP in u(boot_hspp$habitat.guild)){
        FILENAM_GROUP <- paste(FILENAM,"_",format(Sys.time(), "%d%b%Y"),sep="")
        if(nrow(boot_hspp[boot_hspp$Group==SPGROUP,])>3){
          try(
            BOOT_OUT_BYICESGROUP[[SPGROUP]] <- INDfn( DATA=boot_hspp, WRITE=T, SPECIES=SPECIES, GROUP=SPGROUP,
                                                      LFI_THRESHOLD=LFI_THRESHOLD, SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, FILENAM=FILENAM_GROUP,
                                                      LFI=LFI, MEANTL=MEANTLs, MaxL=MaxL, Loo=Loo, Lm=Lm, MeanL=MeanL, TyL_GeoM=TyL_GeoM,
                                                      TyL_SPECIES=F, BYGUILD=F)
            ,silent=TRUE)
        }
      }
    }
    
    if(BYGUILD){#not run boot_hspp does not contain fguild
      next;
      BOOT_OUT_BYGUILD<-list()
      for(SPGROUP in u(guild_dat$fguild)){ # guilds 1:6
        FILENAM_GROUP <- paste(FILENAM,"_",format(Sys.time(), "%d%b%Y"),sep="")
        #select data for guild
        #to do
        if(nrow(boot_hspp[boot_hspp$Group==SPGROUP,])>3){
          try(
            BOOT_OUT_BYGUILD[[SPGROUP]] <- INDfn( DATA=boot_hspp, WRITE=T, SPECIES=SPECIES, GROUP=SPGROUP,
                                                  LFI_THRESHOLD=LFI_THRESHOLD, SAMP_STRAT=SAMP_STRAT, BYSUBDIV=BYSUBDIV, FILENAM=FILENAM_GROUP,
                                                  LFI=LFI, MEANTL=MEANTLs, MaxL=MaxL, Loo=Loo, Lm=Lm, MeanL=MeanL, TyL_GeoM=TyL_GeoM,
                                                  TyL_SPECIES=F, BYGUILD=T)
            ,silent=TRUE)
        }
      }
    }
    
    #x11(); plot(IND_OUT$TyL.cm.sea_dem[,'sea'],    IND_OUT[["TyL_meanln_dem"]][,"Region_mu"], pch=19); abline(0,1)
    #x11(); plot(BOOT_OUT$TyL.cm.sea_dem[,'sea'],   BOOT_OUT[["TyL_meanln_dem"]][,"Region_mu"]); abline(0,1)# why different now?
    #points(BOOT_OUT$TyL.cm.sea_dem[,'sea'],   BOOT_OUT[["TyL_meanln_dem"]][,"Region_mu"],col=B+1,pch=19); 
    #points(IND_OUT$TyL.cm.sea_dem[,'sea'],BOOT_OUT$TyL.cm.sea_dem[,'sea'],col=4,pch=4) # so problem with boot meanln_dem
    #points(IND_OUT$TyL.cm.sea_dem[,'sea'],BOOT_OUT[["TyL_meanln_dem"]][,"Region_mu"],col=6,pch=2)
    #BOOT_OUT$TyL.cm.sea_dem #seems fine 
    #BOOT_OUT[["TyL_meanln_dem"]][,"Region_mu"] #too low
    
    #collate each indicator output in a list
    #indicators from bootstrap by region and whole area  dem / pel ? #what is average of tyl mean by stsq; plot(YRS,apply(BOOT_OUT$TyL_meanln_dem[ ,1:(ncol(BOOT_OUT$TyL_meanln_dem)-4)],1,mean,na.rm=T),col=1,lwd=4)
    
    ### need to add SPGROUP here so not only BOOT_OUT saved
    for(GROUP in SPECIES){
      if(GROUP=="ALL") GROUP <- "all" # edit Jul2017
      if(GROUP=="PEL") GROUP <- "pel"
      if(GROUP=="DEM") GROUP <- "dem"
      print(GROUP) 
      #select and rename to simplify                                                                                  
      if( LFI ){  LFI_BOOT_OUT <- data.frame(BOOT_OUT[paste("LFI_by_sub_",GROUP,sep="")]); 
      if(nrow(LFI_BOOT_OUT)>0  & !is.null(nrow(LFI_BOOT_OUT)) ){
        ss<- unlist( strsplit( names(LFI_BOOT_OUT), paste(GROUP,".",sep="") ) );      
        names(LFI_BOOT_OUT) <- ss[-seq(1,length(ss),2)]
        LFI_regional[[GROUP]] <- rbind(LFI_regional[[GROUP]],  LFI_BOOT_OUT[["LFIregional"]] )
      }
      }
      if(TyL_GeoM){ 
        TyLrect_BOOT_OUT <- data.frame(BOOT_OUT[paste("TyL.cm.sea_",GROUP,sep="")]   );
        if(nrow(TyLrect_BOOT_OUT)>0  & !is.null(nrow(TyLrect_BOOT_OUT))  ){
          ss<- unlist( strsplit( names(TyLrect_BOOT_OUT), paste(GROUP,".",sep="") ) );  names(TyLrect_BOOT_OUT) <- ss[-seq(1,length(ss),2)]
          TyLrect_regional[[GROUP]] <- rbind(TyLrect_regional[[GROUP]],  TyLrect_BOOT_OUT$sea )
        }
      }
      if(MeanL ){
        Len_BOOT_OUT <- data.frame(BOOT_OUT[paste("FishLength_cmsea_",GROUP,sep="")] );
        if(nrow(Len_BOOT_OUT)>0  & !is.null(nrow(Len_BOOT_OUT))){
          ss<- unlist( strsplit( names(Len_BOOT_OUT), paste(GROUP,".",sep="") ) );      names(Len_BOOT_OUT) <- ss[-seq(1,length(ss),2)]
          Len_regional[[GROUP]] <- rbind(Len_regional[[GROUP]],  Len_BOOT_OUT$sea )
        }
      }
      if( MaxL ){ 
        MaxL_BOOT_OUT <- data.frame(BOOT_OUT[paste("MaxLsea_",GROUP,sep="")]);
        if(nrow(MaxL_BOOT_OUT)>0  & !is.null(nrow(MaxL_BOOT_OUT))){
          ss<- unlist( strsplit( names(MaxL_BOOT_OUT), paste(GROUP,".",sep="") ) );     names(MaxL_BOOT_OUT) <- ss[-seq(1,length(ss),2)]
          MaxL_regional[[GROUP]] <- rbind(MaxL_regional[[GROUP]],  MaxL_BOOT_OUT$sea )
        }
      }
      if( Loo ){ 
        Loo_BOOT_OUT <- data.frame(BOOT_OUT[paste("Loosea_",GROUP,sep="")]);
        if(nrow(Loo_BOOT_OUT)>0  & !is.null(nrow(Loo_BOOT_OUT))){
          ss<- unlist( strsplit( names(Loo_BOOT_OUT), paste(GROUP,".",sep="") ) );     names(Loo_BOOT_OUT) <- ss[-seq(1,length(ss),2)]
          Loo_regional[[GROUP]] <- rbind(Loo_regional[[GROUP]],  Loo_BOOT_OUT$sea )
        }
      }
      if( Lm ){ 
        Lm_BOOT_OUT <- data.frame(BOOT_OUT[paste("Lmsea_",GROUP,sep="")]);
        if(nrow(Lm_BOOT_OUT)>0  & !is.null(nrow(Lm_BOOT_OUT))){
          ss<- unlist( strsplit( names(Lm_BOOT_OUT), paste(GROUP,".",sep="") ) );     names(Lm_BOOT_OUT) <- ss[-seq(1,length(ss),2)]
          Lm_regional[[GROUP]] <- rbind(Lm_regional[[GROUP]],  Lm_BOOT_OUT$sea )
        }
      }
      if(MEANTLs){ 
        TL_BOOT_OUT <- data.frame(BOOT_OUT[paste("TLsea_",GROUP,sep="")])
        if(nrow(TL_BOOT_OUT)>0 & !is.null(nrow(TL_BOOT_OUT))){
          ss<- unlist( strsplit( names(TL_BOOT_OUT), paste(GROUP,".",sep="") ) );       names(TL_BOOT_OUT) <- ss[-seq(1,length(ss),2)]
          TL_regional[[GROUP]] <- rbind(TL_regional[[GROUP]],  TL_BOOT_OUT$sea )
        }
      }
      #by strata?
      #if(STRATA) BOOTSTRAP_STRATA_OUT(GROUP=GROUP, Len_BOOT_OUT,MaxL_BOOT_OUT,TL_BOOT_OUT,TyLrect_BOOT_OUT,TyL_BOOT_OUT,LFI_BOOT_OUT,STRATA_SET="subdiv") 
    }#finish group
  } 
  # finish B loop
  
  #now sort out column names
  if(LFI & !is.null(LFI_regional[[GROUP]]) )  colnames(LFI_regional[[GROUP]])<-YRS
  if(TyL_GeoM & !is.null(TyLrect_regional[[GROUP]]) ) colnames(TyLrect_regional[[GROUP]])<-YRS
  if(MeanL  & !is.null(Len_regional[[GROUP]]) ) colnames(Len_regional[[GROUP]])<-YRS
  if( MaxL  & !is.null(MaxL_regional[[GROUP]]) ) colnames(MaxL_regional[[GROUP]])<-YRS
  if( Loo  & !is.null(Loo_regional[[GROUP]]) ) colnames(Loo_regional[[GROUP]])<-YRS
  if( Lm  & !is.null(Lm_regional[[GROUP]]) ) colnames(Lm_regional[[GROUP]])<-YRS
  if(MEANTLs & !is.null(TL_regional[[GROUP]]) ) colnames(TL_regional[[GROUP]])<-YRS
  if(STRATA){##not run
    if(exists("LFI_LIST")){ for(i in 1:length(LFI_LIST)){ if(!is.null(dim(LFI_LIST[[i]]))) colnames(LFI_LIST[[i]])<-YRS} }
    if(exists("TL_LIST")){ for(i in 1:length(TL_LIST)){ if(!is.null(dim(TL_LIST[[i]]))) colnames(TL_LIST[[i]])<-YRS} }
    if(exists("Len_LIST")){ for(i in 1:length(Len_LIST)){ if(!is.null(dim(Len_LIST[[i]]))) colnames(Len_LIST[[i]])<-YRS} }
    if(exists("MaxL_LIST")){ for(i in 1:length(MaxL_LIST)){ if(!is.null(dim(MaxL_LIST[[i]]))) colnames(MaxL_LIST[[i]])<-YRS} }
    if(exists("Loo_LIST")){ for(i in 1:length(Loo_LIST)){ if(!is.null(dim(Loo_LIST[[i]]))) colnames(Loo_LIST[[i]])<-YRS} }
    if(exists("Lm_LIST")){ for(i in 1:length(Lm_LIST)){ if(!is.null(dim(Lm_LIST[[i]]))) colnames(Lm_LIST[[i]])<-YRS} }
    if(exists("TyL_LIST")){ for(i in 1:length(TyL_LIST)){ if(!is.null(dim(TyL_LIST[[i]]))) colnames(TyL_LIST[[i]])<-YRS} }
    if(exists("TyLrect_LIST")){ for(i in 1:length(TyLrect_LIST)){ if(!is.null(dim(TyLrect_LIST[[i]]))) colnames(TyLrect_LIST[[i]])<-YRS}}
  }
  #done  
  if(SAVE) save.image(paste(BOOTOUTPATH,"Indicators_",format(Sys.time(), "%d%b %H%M"),"_",B,"_BOOTS.RData",sep=''))
##end BOOTSTRAP