# Author: C Lynam, Cefas
# Contact: chris.lynam@cefas.co.uk
# Version: 1 
# Date: May 2020 
#Lynam_IND_script_MAKEGUILDIND.R
#guild_dat$taxa_PRED
dhspp$pred.size.class <- "M" #if >Lmat and < half way between Lmat and Lmax
dhspp$pred.size.class[ ( dhspp$FishLength_cm <  dhspp$Lm )]  <- "Sm" #if <Lmat  
dhspp$pred.size.class[ ( dhspp$FishLength_cm < (0.5*dhspp$Lm))]  <- "S" #if <.5*Lmat 
dhspp$pred.size.class[ ( dhspp$FishLength_cm <  3 )]  <- "Lv" #larv 
dhspp$pred.size.class[ ( dhspp$FishLength_cm >  (dhspp$Lm + (dhspp$MaxL-dhspp$Lm)/2) )] <- "L" #if >Lmat and >= half way between Lmat and Lmax
#
# make taxonomic info between data compatible
# sandeels are in by family in trawl data'
guild_dat$taxa_PRED <- ac(guild_dat$taxa_PRED)
#guild_dat$family_PRED <- ac(guild_dat$family_PRED)
guild_dat$taxa_PRED[guild_dat$taxa_PRED=="Ammodytes"]<- "Ammodytidae"
#dhsppf<- merge(dhspp, guild_dat, by.x=c("SpeciesSciName","pred.size.class"), by.y=c("taxa_PRED","pred.size.class"),all=F)

guild_dat$taxa_PRED[guild_dat$taxa_PRED == 'Diplecogaster bimaculata'] <- "Diplecogaster bimaculata bimaculata"
guild_dat$taxa_PRED[guild_dat$taxa_PRED=='Liparis liparis']<- "Liparis liparis liparis"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == 'Mullus barbatus'] <- "Mullus barbatus barbatus"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == 'Salmo trutta'] <- "Salmo trutta trutta"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == 'Scomberesox saurus'] <- "Scomberesox saurus saurus"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == 'Gasterosteus aculeatus'] <- "Gasterosteus aculeatus aculeatus"# Gasterosteus aculeatus williamsoni

guild_dat$taxa_PRED[guild_dat$taxa_PRED == 'Psetta maxima'] <- "Scophthalmus maximus"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == 'Raja fullonica'] <- "Leucoraja fullonica"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == "Callionymidae" ]<- "Callionymus reticulatus"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == "Argentinidae" ]<- "Argentina silus"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == "Syngnathidae" ]<- "Argentinidae"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == "Sebastes marinus" ]<- "Sebastes norvegicus"
guild_dat$taxa_PRED[guild_dat$taxa_PRED == "Sebastes" ]<- "Sebastes viviparus"
#guild_dat$taxa_PRED[guild_dat$taxa_PRED == "Raja microocellata" ]<- #small eyed ray
#not matched: "Dipturus linteus"#white skate   #and #     "Leucoraja circularis" #Sandy Skate 
dhspp<- merge(dhspp, guild_dat, by.x=c("SpeciesSciName","pred.size.class"), 
              by.y=c("taxa_PRED","pred.size.class"), all.x=T)#taxa not SpeciesSciName if using trawldf

#sometimes missing guild for groups as no DAPSTOM Data 
# sort(u(guild_dat$taxa_PRED[!guild_dat$taxa_PRED %in% intersect(guild_dat$taxa_PRED,dhspp$SpeciesSciName)]))
#sort(u(dhspp$SpeciesSciName[!dhspp$SpeciesSciName %in% intersect(guild_dat$taxa_PRED,dhspp$SpeciesSciName)]))
#which are missing
if(nrow( dhspp[is.na(dhspp$fguild),])>0){ 
  missing <- dhspp[is.na(dhspp$fguild),] #93988 # 17180 from nsea sites
  NCOL_DELETE_MISSING<- ncol(missing)-ncol(guild_dat)+2
  NROW <- nrow(dhspp)
  print(paste(survey, signif(100*nrow( dhspp[is.na(dhspp$fguild),])/NROW,3), "% missing records from guild analyses"),sep=',')#dhspp_raw
  #fill in blanks...
  if(!FILL_GUILD) dhspp[is.na(dhspp$fguild),"fguild"] <- "NG"
  if(FILL_GUILD){
    dhspp <- dhspp[!is.na(dhspp$fguild),]#remove missing
    missing_L<-NULL
    if(nrow( missing[missing$pred.size.class=="L",]) >0){ 
      missing_L <- missing[missing$pred.size.class=="L",1:NCOL_DELETE_MISSING]
      missing_L$pred.size.class <- "M" 
      missing_L<- merge(missing_L, guild_dat, by.x=c("SpeciesSciName","pred.size.class"), by.y=c("SpeciesSciName","pred.size.class"), all.x=T)
      
      filled_LwithM<- missing_L[!is.na(missing_L$fguild),]
      missing_L <- missing_L[is.na(missing_L$fguild),]
      #if still missing try S? good idea or not? could be ok for spome species and not others! no data 
      if(nrow(missing_L)>0){
        missing_L$pred.size.class <- "S" 
        missing_L<- merge(missing_L[,1:NCOL_DELETE_MISSING], guild_dat, by.x=c("SpeciesSciName","pred.size.class"), by.y=c("SpeciesSciName","pred.size.class"), all.x=T)
        filled_LwithS<- missing_L[!is.na(missing_L$fguild),]
        missing_L <- missing_L[is.na(missing_L$fguild),]
      }
      #if still missing have to lose
    }
    
    missing_S<-NULL
    if(nrow( missing[missing$pred.size.class=="S",])>0){ 
      missing_S <- missing[missing$pred.size.class=="S",1:NCOL_DELETE_MISSING]
      missing_S$pred.size.class <- "M" 
      missing_S<- merge(missing_S, guild_dat, by.x=c("SpeciesSciName","pred.size.class"), by.y=c("SpeciesSciName","pred.size.class"), all.x=T)
      #
      filled_SwithM<- missing_S[!is.na(missing_S$fguild),]
      missing_S <- missing_S[is.na(missing_S$fguild),]
      #if still missing try L ? good idea or not? could be ok for spome species and not others! no data 
      if(nrow(missing_S)>0){
        missing_S$pred.size.class <- "L" 
        missing_S<- merge(missing_S[,1:NCOL_DELETE_MISSING], guild_dat, by.x=c("SpeciesSciName","pred.size.class"), by.y=c("SpeciesSciName","pred.size.class"), all.x=T)
        filled_SwithL<- missing_S[!is.na(missing_S$fguild),]
        missing_S <- missing_S[is.na(missing_S$fguild),]
      }
      #if still missing have to lose
    }
    
    missing_M<-NULL
    if(nrow( missing[missing$pred.size.class=="M",])>0){ 
      missing_M <- missing[missing$pred.size.class=="M",1:NCOL_DELETE_MISSING]
      #assusme M fish more similar to other mature L fish
      missing_M$pred.size.class <- "L" 
      missing_M<- merge(missing_M, guild_dat, by.x=c("SpeciesSciName","pred.size.class"), by.y=c("SpeciesSciName","pred.size.class"), all.x=T)
      filled_MwithL<- missing_M[!is.na(missing_M$fguild),]
      missing_M <- missing_M[is.na(missing_M$fguild),]
      #if still missing try S
      if(nrow( missing_M)>0){ 
        missing_M$pred.size.class <- "S" 
        missing_M<- merge(missing_M[,1:NCOL_DELETE_MISSING], guild_dat, by.x=c("SpeciesSciName","pred.size.class"), by.y=c("SpeciesSciName","pred.size.class"), all.x=T)
        filled_MwithS<- missing_M[!is.na(missing_M$fguild),]
        missing_M <- missing_M[is.na(missing_M$fguild),]
      }
    }
    missing<- rbind(missing_L,missing_M,missing_S)
    filled<- rbind(filled_LwithM,filled_LwithS, filled_MwithL,filled_MwithS, filled_SwithL,filled_SwithM)
    print(paste(survey, signif(100*nrow( missing )/NROW,3), "% missing records after filling lost from guild analyses "),sep=',')
    print(paste(survey, signif(100*nrow( filled  )/NROW,3), "% missing records filled with other size classes for guild analyses"),sep=',')
    
    dhspp <- rbind(dhspp, filled)
    write.csv(file=paste(OUTPATH,survey,"_filled guilds.csv",sep=""),x=filled)
  }#FILL
  dhspp$Group<- dhspp$fguild
  #note the missing records
  write.csv(file=paste(OUTPATH,survey,"_missing guilds.csv",sep=""),x=missing)
}#IF is.nas in GUILD


