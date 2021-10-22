INDPLOTFN<-function(DATA2PLOT=DATA2PLOT, BOOTDATA2PLOT=BOOTDATA2PLOT,YRS=YRS,YLIM=NULL,YLAB=YLAB,TITLE=TITLE,GAMMOD=GAMMOD,
         BOOTSTRAP=F, ADDBOOTTREND=F, ADDBOOTTREND_CI=F, ADDBOOT_ERRBAR=F, ADDGAM=F, BEST_AND_BOOT=F,
         ADDLOESS=T,ADDLAST6LM=T){
  #BOOTSTRAP<-T; ADDBOOTTREND<-T; ADDBOOTTREND_CI<-F; ADDBOOT_ERRBAR<-T; ADDGAM<-T; BEST_AND_BOOT<-T
  #BOOTSTRAP<-F; ADDBOOTTREND<-F; ADDBOOTTREND_CI<-F; ADDBOOT_ERRBAR<-F; ADDGAM<-F; BEST_AND_BOOT<-F;  ADDLOESS=T; ADDLAST6LM=T
  #DATA2PLOT<-DATA2PLOT[,n]; YRS<-YRSPLOT#an(names(DATA2PLOT))
  plot(YRS,DATA2PLOT,col='white',ylim=YLIM,ylab=YLAB)
  title(TITLE)
  if( (ADDGAM & (!ADDBOOTTREND_CI | !BOOTSTRAP)) | BEST_AND_BOOT)  points(YRS,   DATA2PLOT,col=4, pch=1); 
  #jul2020 changed from green triangles to blue circles
  
  #add gam fit
  if(ADDGAM){
    lines(YRS, GAMMOD$fit,lwd=4)
    lines(YRS, GAMMOD$fit-GAMMOD$se.fit,col='dark grey',lwd=4)
    lines(YRS, GAMMOD$fit+GAMMOD$se.fit,col='dark grey',lwd=4)
  }
  if(ADDLOESS & length(DATA2PLOT)>=5){#changed to = and >5 not just >5 in Jul2020
    l<-loess(DATA2PLOT~YRS, span= 10/length(DATA2PLOT), degree=2, control = loess.control(surface = "direct"))#      S<-c(S,l$s)  # the residual standard error
    #LOESS with se in line
    PRED <-  predict(l,newdata=data.frame(YRS=YRS), se = TRUE)
    graphics::polygon( c(YRS, YRS[length(YRS):1]), c(PRED$fit+PRED$se.fit, (PRED$fit-PRED$se.fit)[length(YRS):1] ),col='light blue')
    lines(YRS, PRED$fit,col=4,lwd=1)
    points(YRS, DATA2PLOT,pch=19,col=4,lwd=1)
    #abline(h=min(DATA2PLOT[-((length(DATA2PLOT)-5):length(DATA2PLOT))],na.rm=T),lwd=2,col=4)#excluding last 6
  }
  if(ADDLAST6LM & length(DATA2PLOT)>7){
    LMFIT <- lm(DATA2PLOT[((length(DATA2PLOT)-5):length(DATA2PLOT))]~YRS[((length(DATA2PLOT)-5):length(DATA2PLOT))])
    if(! is.na(summary(LMFIT)$coefficients[2,4] )& summary(LMFIT)$coefficients[2,4] < 0.05) lines(YRS[((length(DATA2PLOT)-5):length(DATA2PLOT))],LMFIT$fit,lwd=1,lty=2,col=2)
  }
  #add all trendlines from bootstrap 
  if(BOOTSTRAP & !is.null(BOOTDATA2PLOT)){ #add !is.null(BOOTDATA2PLOT) since may select to bootstrap some but not all indicators!
    modout<-NULL
    for(i in 1:nrow(BOOTDATA2PLOT)){ 
      mod <- lm(BOOTDATA2PLOT[i,] ~ YRS); 
      if(ADDBOOTTREND) abline(mod,col=i) 
      modout <-  rbind(modout, t(data.frame(coef(mod))) )
    }
    #add 95% in trendline from bootstrap 
    if(ADDBOOTTREND_CI){
      modoutperc<-NULL
      for(YR in YRS) modout <- cbind(modout, modout[,1] + (YR*modout[,2]) )
      colnames(modout)[3:ncol(modout)] <- YRS
      if(nrow(modout)>1){
        modoutperc <- apply(modout[,3:ncol(modout)],2,quantile,probs = c(0.025, 0.5, 0.975))
        lines(YRS, modoutperc[1,],lwd=2,col=4)
        lines(YRS, modoutperc[2,],lwd=2,col=4)
        lines(YRS, modoutperc[3,],lwd=2,col=4)
      }
    }
    #add 95% CI estimates for each year from bootstrap 
    if(ADDBOOT_ERRBAR & nrow(BOOTDATA2PLOT)>=1){ 
      perc <- apply(BOOTDATA2PLOT,2,quantile,probs = c(0.025, 0.5, 0.975))
      errbar(x=YRS,y=perc[2,],yplus=perc[3,],yminus=perc[1,],add=T)
    }
  }
  box()
}
