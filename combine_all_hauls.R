

if(WRITE_LDs | IEO_FW4){
  # append all the data that are required for (a) Chibuzor's modelling work an (b) an output to go into the BX5 app
  fileslisted=list.files(OUTPATHstem,'haul_by_spp',full.names=T,recursive = T);
  # In case of reruns done on previous date
  fileslisted = fileslisted[grepl(time_of_run, fileslisted)]
  file1=read.csv(fileslisted[1]) 
  # Also get survey name from filename, crude approach
  file1$survey_name = strsplit(fileslisted[1],'/')[[1]][length(strsplit(fileslisted[1],'/')[[1]])-1]
  for(file in 2:length(fileslisted)){file2=read.csv(fileslisted[file]); 
    file2$survey_name = strsplit(fileslisted[file],'/')[[1]][length(strsplit(fileslisted[file],'/')[[1]])-1]
    # SurvStratum is required but not present for some surveys, these must be GNS surveys
    if(!"SurvStratum" %in% colnames(file2)){file2$SurvStratum = file2$ICESStSq}
    #file2_names = colnames(file2); file1_names = colnames(file1); common_names = intersect(file2_names, file1_names); file1 = rbind(file2[common_names], file1[common_names])
    file1 = rbind(file2, file1)
    }
  
  # Add species ID from lookup table originating from database table bx005.public.species
  specieslookup=read.csv("C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/public_species_table.csv")
  specieslookup$SpeciesSciName = specieslookup$latin_name
  specieslookup$latin_name <- NULL
  specieslookup = specieslookup[,c("SpeciesSciName","species_id")]
  file1=merge(file1,specieslookup,all.x=T)
  
  # File is too big - aggregate and remove some columns. And rename some to be more consistent with the expected file
  desiredcols=colnames(cod_final)
  desiredcols[!desiredcols %in% colnames(file1)]
  file1$Survey_Acronym = file1$survey_name
  file1$survey_name <- NULL
  file1$GearType = file1$Gear
  file1$Gear <- NULL
  file1$YearShot = file1$Year
  file1$Year <- NULL
  file1$DayShot = file1$Day
  file1$Day <- NULL
  file1=file1[,!colnames(file1) %in% c("KM2_LAM","checkLen","Order","Group","Loo","Lm","MaxL","Max.L..cm.")]
  file1$numhauls=1
  file1=aggregate(cbind(DensBiom_kg_Sqkm,DensAbund_N_Sqkm,DensAbund_N_perhr,DensBiom_kg_perhr,numhauls) ~ SpeciesSciName + HaulID + SensFC1 + DEMPEL + HaulDur_min + ICESStSq + WingSwpArea_sqkm + WingSwpVol_CorF + NetOpen_m + Ship + MonthShot + TimeShot + ShootLong_degdec + ShootLat_degdec + habitat.guild + ScientificName_WoRMS + SurvStratum + species_id + Survey_Acronym + GearType + YearShot + DayShot , data = file1, FUN = sum, na.rm = TRUE)
  file1$sumDensBiom_kg_Sqkm = file1$DensBiom_kg_Sqkm
  file1$sumDensBiom_kg_perhr = file1$DensBiom_kg_perhr
  file1$sumDensAbund_N_Sqkm = file1$DensAbund_N_Sqkm
  file1$sumDensAbund_N_perhr = file1$DensAbund_N_perhr
  file1$DensBiom_kg_Sqkm <- NULL
  file1$DensBiom_kg_perhr <- NULL
  file1$DensAbund_N_Sqkm <- NULL
  file1$DensAbund_N_perhr <- NULL
    
  # Write file for Chibuzor's modelling work with lats and longs
  write.csv(file1,paste0(OUTPATHstem,"hauls_by_spp_all.csv"))
}
