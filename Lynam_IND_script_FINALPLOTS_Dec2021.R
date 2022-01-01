# Author: C Lynam, Cefas
# Contact: chris.lynam@cefas.co.uk
# Version: 2 
# Date: DEC 2021 
#Lynam_IND_script_FINALPLOTS.R

#DEC 2021 HAVE ALTERED INDfn to AVOID HAVING TO RUN CODE TWICE and no longer creating ALL combined
if(SPECIES=="ALL") SPECIES<-c("DEM","PEL")
RATIOS<-F #tyl/mml and tyl/lm

PLOTROWN <- sum(TyL_GeoM,MaxL,Loo,Lm,MEANTL,LFI) + ifelse(sum(TyL_GeoM,MaxL,Lm)==3,1,0)
if(RATIOS & ( (Lm & TyL_GeoM) | (MaxL & TyL_GeoM & !Lm)) ) PLOTROWN<-PLOTROWN+1
PLOTCOLN <- length(SPECIES)
windows(width=PLOTCOLN*8, height=4*PLOTROWN)
par(mfrow=c(PLOTROWN,PLOTCOLN),mar=c(2,4,2,2),oma=c(1,1,3,1))

YRS<- YRS
if(BOOTSTRAP){ ADDBOOTTREND<-F; ADDBOOTTREND_CI<-F; ADDBOOT_ERRBAR<-T; } else { ADDBOOT_ERRBAR<-F; ADDBOOTTREND<-F; ADDBOOTTREND_CI<-F}
ADDGAM<-F
BEST_AND_BOOT<- T # add crosses to plot
TITAdd<-T; ADDLOESS=T; ADDLAST6LM=T


if(LFI){
  #### plot LFI
  YLAB <- "Large Fish Indicator"
  for(plotgroups in 1:length(SPECIES)){
    
    if(SPECIES[plotgroups]=="ALL"){ DATA2PLOT <- IND_OUT[["LFI_by_sub_all"]][,"LFIregional"]; YLIM<- c(0,0.05);  BOOTDATA2PLOT <- LFI_regional$all;if(TITAdd){ TITLE<-"All fish"}  }
    if(SPECIES[plotgroups]=="PEL"){ plot(1,1,col="white",axes=F,ylab="",main="Pelagic fish"); next; 
                                    DATA2PLOT <- IND_OUT[["LFI_by_sub_pel"]][,"LFIregional"]; YLIM<- c(0,0.05);  BOOTDATA2PLOT <- LFI_regional$pel;  }
    if(SPECIES[plotgroups]=="DEM"){ DATA2PLOT <- IND_OUT[["LFI_by_sub_dem"]][,"LFIregional"]; YLIM<- c(0.1,0.45);BOOTDATA2PLOT <- LFI_regional$dem;if(TITAdd){ TITLE<-"Demersal fish"} }
    #print(DATA2PLOT)
    if(!is.null(DATA2PLOT)){
      summary(lmseaLFI<-gam(DATA2PLOT ~ (YRS))); GAMMOD<-predict(lmseaLFI,se=T) 
      if(ADDGAM){ summary(gseaLFI<-gam(DATA2PLOT ~ s(YRS,k=6))); GAMMOD<-predict(gseaLFI,se=T) }
      YLIM<-c( (min(DATA2PLOT,na.rm=T)*.95),  (max(DATA2PLOT,na.rm=T)*1.05) )
      if(exists("BOOT_OUT") & !is.null(BOOTDATA2PLOT) ){ 
        INDPLOTFN(DATA2PLOT[is.finite(DATA2PLOT)], BOOTDATA2PLOT[is.finite(DATA2PLOT)],YRS[is.finite(DATA2PLOT)],YLIM=YLIM,TITLE=TITLE,BEST_AND_BOOT=T,GAMMOD=GAMMOD, YLAB=YLAB
                  ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM)
      } else {#no boot
        INDPLOTFN(DATA2PLOT[is.finite(DATA2PLOT)], BOOTDATA2PLOT[is.finite(DATA2PLOT)],YRS[is.finite(DATA2PLOT)],YLIM=YLIM,TITLE=TITLE,GAMMOD=GAMMOD, YLAB=YLAB
                  ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM, BEST_AND_BOOT=BEST_AND_BOOT)
      }
    }
  }
  TITAdd<-F; TITLE<-""
}

if(TyL_GeoM){
  #### plot tyl geometric
  YLAB <- "TyL geomean (cm)"
  for(plotgroups in 1:length(SPECIES)){
    if(SPECIES[plotgroups]=="ALL"){ DATA2PLOT <- IND_OUT[["TyL.cm.sea_all"]][,"sea"]; BOOTDATA2PLOT <- TyLrect_regional$all; if(TITAdd){TITLE<-"All fish"} }
    if(SPECIES[plotgroups]=="PEL"){ DATA2PLOT <- IND_OUT[["TyL.cm.sea_pel"]][,"sea"]; BOOTDATA2PLOT <- TyLrect_regional$pel; if(TITAdd){TITLE<-"Pelagic fish"} }
    if(SPECIES[plotgroups]=="DEM"){ DATA2PLOT <- IND_OUT[["TyL.cm.sea_dem"]][,"sea"]; BOOTDATA2PLOT <- TyLrect_regional$dem; if(TITAdd){TITLE<-"Demersal fish"} }
    summary(lmTyL<-gam(DATA2PLOT ~ (YRS))); GAMMOD<-predict(lmTyL,se=T) 
    if(ADDGAM){ summary(gTyL<-gam(DATA2PLOT ~ s(YRS,k=6))); GAMMOD<-predict(gTyL,se=T) }
    YLIM<-c(floor(min(DATA2PLOT,na.rm=T)*.95), ceiling(max(DATA2PLOT,na.rm=T)*1.05) )
    if(exists("BOOT_OUT") & !is.null(BOOTDATA2PLOT) ){ 
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,BEST_AND_BOOT=T,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM)
    } else {#no boot
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE, GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM, BEST_AND_BOOT=BEST_AND_BOOT)
    }
  }
  TITAdd<-F; TITLE<-""
}

if(MaxL){
  #### plot MML
  YLAB <- "Mean Max Length (cm)"
  for(plotgroups in 1:length(SPECIES)){
    if(SPECIES[plotgroups]=="ALL"){ DATA2PLOT <- IND_OUT[["MaxLsea_all"]][,"sea"]; YLIM<- c(35,50); BOOTDATA2PLOT <- MaxL_regional$all; if(TITAdd){TITLE<-"All fish"} }
    if(SPECIES[plotgroups]=="PEL"){ DATA2PLOT <- IND_OUT[["MaxLsea_pel"]][,"sea"]; YLIM<- c(35,50); BOOTDATA2PLOT <- MaxL_regional$pel; if(TITAdd){TITLE<-"Pelagic fish"}}
    if(SPECIES[plotgroups]=="DEM"){ DATA2PLOT <- IND_OUT[["MaxLsea_dem"]][,"sea"]; YLIM<- c(45,125); BOOTDATA2PLOT <- MaxL_regional$dem; if(TITAdd){TITLE<-"Demersal fish"} }
    summary(lmseaMML<-gam(DATA2PLOT ~ (YRS))); GAMMOD<-predict(lmseaMML,se=T) 
    if(ADDGAM){ summary(gseaMML<-gam(DATA2PLOT ~ s(YRS,k=6))); GAMMOD<-predict(gseaMML,se=T) }
    YLIM<-c(floor(min(DATA2PLOT,na.rm=T)*.95), ceiling(max(DATA2PLOT,na.rm=T)*1.05) )
    if(exists("BOOT_OUT") & !is.null(BOOTDATA2PLOT) ){
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,BEST_AND_BOOT=T,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM)
    } else {#no boot
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM, BEST_AND_BOOT=BEST_AND_BOOT)
    }
  }
  TITAdd<-F; TITLE<-""
}
if(Loo){
  YLAB <- "Asymptotic Length (cm)"
  for(plotgroups in 1:length(SPECIES)){
    if(SPECIES[plotgroups]=="ALL"){ DATA2PLOT <- IND_OUT[["Loosea_all"]][,"sea"]; YLIM<- c(35,50); BOOTDATA2PLOT <- Loo_regional$all; if(TITAdd){TITLE<-"All fish"} }
    if(SPECIES[plotgroups]=="PEL"){ DATA2PLOT <- IND_OUT[["Loosea_pel"]][,"sea"]; YLIM<- c(35,50); BOOTDATA2PLOT <- Loo_regional$pel; if(TITAdd){TITLE<-"Pelagic fish"}}
    if(SPECIES[plotgroups]=="DEM"){ DATA2PLOT <- IND_OUT[["Loosea_dem"]][,"sea"]; YLIM<- c(45,125); BOOTDATA2PLOT <- Loo_regional$dem; if(TITAdd){TITLE<-"Demersal fish"} }
    summary(lmseaMML<-gam(DATA2PLOT ~ (YRS))); GAMMOD<-predict(lmseaMML,se=T) 
    if(ADDGAM){ summary(gseaMML<-gam(DATA2PLOT ~ s(YRS,k=6))); GAMMOD<-predict(gseaMML,se=T) }
    YLIM<-c(floor(min(DATA2PLOT,na.rm=T)*.95), ceiling(max(DATA2PLOT,na.rm=T)*1.05) )
    if(exists("BOOT_OUT") & !is.null(BOOTDATA2PLOT) ){
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,BEST_AND_BOOT=T,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM)
    } else {#no boot
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM, BEST_AND_BOOT=BEST_AND_BOOT)
    }
  }
  TITAdd<-F; TITLE<-""
}
if(Lm){
  YLAB <- "Length at Maturity (cm)"
  for(plotgroups in 1:length(SPECIES)){
    if(SPECIES[plotgroups]=="ALL"){ DATA2PLOT <- IND_OUT[["Lmsea_all"]][,"sea"]; YLIM<- c(35,50); BOOTDATA2PLOT <- Lm_regional$all; if(TITAdd){TITLE<-"All fish"} }
    if(SPECIES[plotgroups]=="PEL"){ DATA2PLOT <- IND_OUT[["Lmsea_pel"]][,"sea"]; YLIM<- c(35,50); BOOTDATA2PLOT <- Lm_regional$pel; if(TITAdd){TITLE<-"Pelagic fish"}}
    if(SPECIES[plotgroups]=="DEM"){ DATA2PLOT <- IND_OUT[["Lmsea_dem"]][,"sea"]; YLIM<- c(45,125); BOOTDATA2PLOT <- Lm_regional$dem; if(TITAdd){TITLE<-"Demersal fish"} }
    summary(lmseaMML<-gam(DATA2PLOT ~ (YRS))); GAMMOD<-predict(lmseaMML,se=T) 
    if(ADDGAM){ summary(gseaMML<-gam(DATA2PLOT ~ s(YRS,k=6))); GAMMOD<-predict(gseaMML,se=T) }
    YLIM<-c(floor(min(DATA2PLOT,na.rm=T)*.95), ceiling(max(DATA2PLOT,na.rm=T)*1.05) )
    if(exists("BOOT_OUT") & !is.null(BOOTDATA2PLOT) ){
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,BEST_AND_BOOT=T,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM)
    } else {#no boot
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM, BEST_AND_BOOT=BEST_AND_BOOT)
    }
  }
  TITAdd<-F; TITLE<-""
}

if(MEANTL){
  #### plot tl 
  YLAB <- "MTL"
  for(plotgroups in 1:length(SPECIES)){
    if(SPECIES[plotgroups]=="ALL"){ DATA2PLOT <- IND_OUT[["TLsea_all"]][,"sea"]; BOOTDATA2PLOT <- TL_regional$all; if(TITAdd){TITLE<-"All fish"}  }
    if(SPECIES[plotgroups]=="PEL"){ DATA2PLOT <- IND_OUT[["TLsea_pel"]][,"sea"]; BOOTDATA2PLOT <- TL_regional$pel; if(TITAdd){TITLE<-"Pelagic fish"}  }
    if(SPECIES[plotgroups]=="DEM"){ DATA2PLOT <- IND_OUT[["TLsea_dem"]][,"sea"]; BOOTDATA2PLOT <- TL_regional$dem; if(TITAdd){TITLE<-"Demersal fish"}  }
    summary(lmTL<-gam(DATA2PLOT ~ (YRS))); GAMMOD<-predict(lmTL,se=T) 
    if(ADDGAM){ summary(gTL<-gam(DATA2PLOT ~ s(YRS,k=6))); GAMMOD<-predict(gTL,se=T) }
    YLIM<-c((min(DATA2PLOT,na.rm=T)*.975), (max(DATA2PLOT,na.rm=T)*1.025) )
    if(exists("BOOT_OUT") & !is.null(BOOTDATA2PLOT) ){ 
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,BEST_AND_BOOT=T,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM)
    } else {#no boot
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE, GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=ADDBOOTTREND_CI, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, ADDGAM=ADDGAM, BEST_AND_BOOT=BEST_AND_BOOT)
    }
  }
  TITAdd<-F; TITLE<-""
}

if(RATIOS & MaxL & TyL_GeoM & !Lm){
  ####plot TyL / MML
  YLAB <- "TyL/MML"
  for(plotgroups in 1:length(SPECIES)){
    if(SPECIES[plotgroups]=="ALL"){ DATA2PLOT <- IND_OUT[["TyL.cm.sea_all"]][,"sea"]/IND_OUT[["MaxLsea_all"]][,"sea"]; BOOTDATA2PLOT <- TyLrect_regional$all/MaxL_regional$all }
    if(SPECIES[plotgroups]=="PEL"){ DATA2PLOT <- IND_OUT[["TyL.cm.sea_pel"]][,"sea"]/IND_OUT[["MaxLsea_pel"]][,"sea"]; BOOTDATA2PLOT <- TyLrect_regional$pel/MaxL_regional$pel }
    if(SPECIES[plotgroups]=="DEM"){ DATA2PLOT <- IND_OUT[["TyL.cm.sea_dem"]][,"sea"]/IND_OUT[["MaxLsea_dem"]][,"sea"]; BOOTDATA2PLOT <- TyLrect_regional$dem/MaxL_regional$dem}
    summary(lmTyL<-gam(DATA2PLOT ~ (YRS))); GAMMOD<-predict(lmTyL,se=T) 
    if(ADDGAM){ summary(gTyL<-gam(DATA2PLOT ~ s(YRS,k=6))); GAMMOD<-predict(gTyL,se=T) }
    YLIM<-c((min(DATA2PLOT,na.rm=T)*.95), (max(DATA2PLOT,na.rm=T)*1.05) )
    if(exists("BOOT_OUT") & !is.null(BOOTDATA2PLOT) ){ 
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,BEST_AND_BOOT=T,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,
                BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=F, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, 
                ADDGAM=ADDGAM)
    } else {#no boot
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE, BEST_AND_BOOT=BEST_AND_BOOT, GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM, BOOTSTRAP=F,ADDGAM=ADDGAM)
    }
  }
}
if(RATIOS & Lm & TyL_GeoM){
  ####plot TyL / Lm
  YLAB <- "TyL/Lm"
  for(plotgroups in 1:length(SPECIES)){
    if(SPECIES[plotgroups]=="ALL"){ DATA2PLOT <- IND_OUT[["TyL.cm.sea_all"]][,"sea"]/IND_OUT[["Lmsea_all"]][,"sea"]; BOOTDATA2PLOT <- TyLrect_regional$all/Lm_regional$all }
    if(SPECIES[plotgroups]=="PEL"){ DATA2PLOT <- IND_OUT[["TyL.cm.sea_pel"]][,"sea"]/IND_OUT[["Lmsea_pel"]][,"sea"]; BOOTDATA2PLOT <- TyLrect_regional$pel/Lm_regional$pel }
    if(SPECIES[plotgroups]=="DEM"){ DATA2PLOT <- IND_OUT[["TyL.cm.sea_dem"]][,"sea"]/IND_OUT[["Lmsea_dem"]][,"sea"]; BOOTDATA2PLOT <- TyLrect_regional$dem/Lm_regional$dem}
    summary(lmTyL<-gam(DATA2PLOT ~ (YRS))); GAMMOD<-predict(lmTyL,se=T) 
    if(ADDGAM){ summary(gTyL<-gam(DATA2PLOT ~ s(YRS,k=6))); GAMMOD<-predict(gTyL,se=T) }
    YLIM<-c((min(DATA2PLOT,na.rm=T)*.95), (max(DATA2PLOT,na.rm=T)*1.05) )
    if(exists("BOOT_OUT") & !is.null(BOOTDATA2PLOT) ){ 
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE,BEST_AND_BOOT=T,GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM,
                BOOTSTRAP=BOOTSTRAP, ADDBOOTTREND=ADDBOOTTREND, ADDBOOTTREND_CI=F, ADDBOOT_ERRBAR=ADDBOOT_ERRBAR, 
                ADDGAM=ADDGAM)
    } else {#no boot
      INDPLOTFN(DATA2PLOT, BOOTDATA2PLOT,YRS,YLIM=YLIM,TITLE=TITLE, BEST_AND_BOOT=BEST_AND_BOOT, GAMMOD=GAMMOD, YLAB=YLAB
                ,ADDLOESS=ADDLOESS,ADDLAST6LM=ADDLAST6LM, BOOTSTRAP=F,ADDGAM=ADDGAM)
    }
    abline(h=1,col=2)
  }
}
mtext(paste(survey,sep=' '),line=0,outer=T)
WHICHIND<- "_"
if(LFI) WHICHIND<- paste(WHICHIND, "LFI.",sep='')
if(MaxL) WHICHIND<- paste(WHICHIND, "MaxL.",sep='')
if(Loo) WHICHIND<- paste(WHICHIND, "Loo.",sep='')
if(Lm) WHICHIND<- paste(WHICHIND, "Lm.",sep='')
if(TyL_GeoM) WHICHIND<- paste(WHICHIND, "TyL.",sep='')
if(MEANTL) WHICHIND<- paste(WHICHIND, "MTL.",sep='')
if(CATCHABILITY_COR_WALKER){qcor<-"Qwalk"} else {qcor<-NULL}

if(length(SPECIES)>1) SPECIES<-c("DEMPEL")
savePlot(filename= paste(FILENAM,SPECIES,qcor,WHICHIND,"bmp",sep=''),type="bmp")
