
MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR_JR/"
RDIR<- paste0(MAINDIR,"R/")
QADIR_pre_update<- paste0(MAINDIR,"app_QA_downloaded_from_zapp/preupdated/")
QADIR_update<- paste0(MAINDIR,"app_QA_downloaded_from_zapp/updated/")
outdir<- paste0(MAINDIR,"app_QA_downloaded_from_zapp/")

checklist_pre = list.files(QADIR_pre_update,include.dirs=T,full.names = T)
checklist_post = list.files(QADIR_update,include.dirs=T,full.names = T)
outdf = data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c("speciesarea","Survey","EEZgroup","EEZJointRegime","Year","ZAPercentBeforeUpdate","ZAPercentAfterUpdate"))))

for(combrow in 1:length(checklist_pre)){
  prefile=read.csv(checklist_pre[combrow])
  postfile=read.csv(checklist_post[combrow])
  
  prefile$assessmentyear = 2021
  postfile$assessmentyear = 2022

  prefile$link=paste0(prefile$Year,prefile$EEZ.or.joint.regime,prefile$Survey,prefile$EEZ.group)
  postfile$link=paste0(postfile$Year,postfile$EEZ.or.joint.regime,postfile$Survey,postfile$EEZ.group)
  
  combfiles=merge(prefile,postfile,by='link',suffixes = c("2021","2022"))
  arealst = strsplit(checklist_post[combrow],"/");   arealst = arealst[[1]][length(arealst[[1]])]
  areai = substr(arealst,27,nchar(arealst)-4)
  combfiles$speciesarea = areai

  combfiles = combfiles[,c("speciesarea","Survey2022","Year2022","EEZ.group2022","EEZ.or.joint.regime2022","ZA..2021","ZA..2022")]
  colnames(combfiles) = c("speciesarea","Survey","Year","EEZgroup","EEZJointRegime","ZAPercentBeforeUpdate","ZAPercentAfterUpdate")

  outdf=rbind(outdf,combfiles)
}

outdf <- outdf[order(outdf$speciesarea,outdf$Survey,outdf$EEZgroup,outdf$EEZJointRegime,outdf$Year),]

outdf$percent_difference = 100*(outdf$ZAPercentAfterUpdate/outdf$ZAPercentBeforeUpdate)-100
write.csv(outdf,paste0(outdir,"differences_caused_by_update.csv"),row.names = F)
