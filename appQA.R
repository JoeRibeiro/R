
MAINDIR<- "C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR_JR/"
RDIR<- paste0(MAINDIR,"R/")
QADIR_pre_update<- paste0(MAINDIR,"app_QA_downloaded_from_zapp/pre_update/")
QADIR_update<- paste0(MAINDIR,"app_QA_downloaded_from_zapp/updated/")

checklist_pre = list.files(QADIR_pre_update,include.dirs=T,full.names = T)
checklist_post = list.files(QADIR_update,include.dirs=T,full.names = T)
outdf = data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Survey","EEZgroup","EEZJointRegime","ZAPercentBeforeUpdate","ZAPercentAfterUpdate"))))

for(combrow in 1:length(checklist_pre)){
  prefile=read.csv(checklist_pre[combrow])
  postfile=read.csv(checklist_post[combrow])
  
  prefile$assessmentyear = 2021
  postfile$assessmentyear = 2022

  prefile$link=paste0(prefile$Year,prefile$EEZ.or.joint.regime,prefile$Survey,prefile$EEZ.group)
  postfile$link=paste0(postfile$Year,postfile$EEZ.or.joint.regime,postfile$Survey,postfile$EEZ.group)
  
  combfiles=merge(prefile,postfile,by='link',suffixes = c("2021","2022"))
  combfiles = combfiles[,c("Survey2022","EEZ.group2022","EEZ.or.joint.regime2022","ZA..2021","ZA..2022")]
  colnames(combfiles) = c("Survey","EEZgroup","EEZJointRegime","ZAPercentBeforeUpdate","ZAPercentAfterUpdate")
  outdf=rbind(outdf,combfiles)
}
  
  