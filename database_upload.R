
if(WRITE_to_DB){
  # upload to bx5 postgres database
  library(DBI)

  # append all the data that are required for an output to go into the BX5 app
  fileslisted=list.files(OUTPATHstem,'LD_tonnes_Year_W.by.subdiv',full.names=T,recursive = T);
  # In case of reruns done on previous date
  fileslisted = fileslisted[grepl(time_of_run, fileslisted)]
  file1=read.csv(fileslisted[1]) 
  # Also get survey name from filename, crude approach
  file1$survey_name = strsplit(fileslisted[1],'/')[[1]][length(strsplit(fileslisted[1],'/')[[1]])-1]
  for(file in 2:length(fileslisted)){file2=read.csv(fileslisted[file]); 
    file2$survey_name = strsplit(fileslisted[file],'/')[[1]][length(strsplit(fileslisted[file],'/')[[1]])-1]
    file2_names = colnames(file2); file1_names = colnames(file1); common_names = intersect(file2_names, file1_names); file1 = rbind(file2[common_names], file1[common_names])
  }
  
  # Add species ID from lookup table originating from database table bx005.public.species
  specieslookup=read.csv("C:/Users/JR13/OneDrive - CEFAS/Fish_dataproduct_QSR/SweptArea_29Oct2021/public_species_table.csv")
  specieslookup$SpeciesSciName = specieslookup$latin_name
  specieslookup$latin_name <- NULL
  specieslookup = specieslookup[,c("SpeciesSciName","species_id")]
  file1=merge(file1,specieslookup,all.x=T)

  # chosen column names to lower, this next file will be going into a postgresql database for the app:
  colnames(file1) <-tolower(colnames(file1))
  fileout = file1[,c("survstratum","year","catcatchwgtswept","fishlength_cm","dempel","order","group","survey_name","species_id")]
#   fileout=fileout[c(!is.na(fileout$species_id) & !is.na(fileout$catcatchwgtswept) & !is.na(fileout$year) & !is.na(fileout$survey_name) & !is.na(fileout$wingswparea_sqkm)) ,]
  fileout$fishlength_cm=as.integer(fileout$fishlength_cm)

  connect_t0_db <- function(connection_info_text="W:/WPx_Application/appBX005IK/ZAPP beta/bxappserver.txt"){
    connection_info <- read.csv(connection_info_text, header=T, stringsAsFactors=FALSE, sep="\t")
    #establish connection to the database we want to make the table in
    drv <- RPostgres::Postgres()
    con <- RPostgres::dbConnect( #before there was DBI
      drv,
      dbname = connection_info$res[connection_info$var == "db"],
      host = connection_info$res[connection_info$var == "host"],
      port = connection_info$res[connection_info$var == "port"],
      user = connection_info$res[connection_info$var == "user"],
      password = connection_info$res[connection_info$var == "pword"]
    )
    return(c(con, drv,connection_info))
  }
  conF<-connect_t0_db()
  con1 <-conF[[1]]
  drv <- conF[[2]]
  con_info <- conF[[3]]

  # There seem to be surveys that don't exist in the ICES data product. I don't want the app to lose functionality by losing these surveys, so download them and merge before reuploading
  cpua_original_2021 <- DBI::dbReadTable(con1, Id(schema = "wp2", table = "cpua_original_2021"))
  cpua_original_2021 = cpua_original_2021[cpua_original_2021$survey_name %in% c("heras","peltic","Q1SWOTTER","Q1SWBEAM"),]
  fileout$cpua_id = rownames(fileout)
  fileout=rbind(fileout,cpua_original_2021)

  # Also make a cpua_id
  fileout$cpua_id = rownames(fileout)

  # Write to table so field types are consistent with other table called cpua
  DBI::dbWriteTable(con1, Id(schema = "wp2", table = "cpua"), value = fileout, field.types = c(cpua_id="INTEGER PRIMARY KEY",survstratum="varchar",year="int",catcatchwgtswept="real",fishlength_cm="int",dempel="varchar",order="varchar",group="varchar",survey_name="varchar",species_id="int"), row.names=FALSE)
}

print("script complete")


# ## Code for some comparisons vs the old cpua table. This is a dataproduct from quite a different derivation so will be different due to different correction factors and SSA. But there should at least be a correlation
# cpua_original_2021 <- DBI::dbReadTable(con1, Id(schema = "wp2", table = "cpua_original_2021"))
# fileout = fileout[complete.cases(fileout), ]; cpua_original_2021 = cpua_original_2021[complete.cases(cpua_original_2021), ]
# cpua_original_2021$mergeon = paste0(cpua_original_2021$survstratum,cpua_original_2021$year,cpua_original_2021$fishlength_cm,cpua_original_2021$dempel,cpua_original_2021$order,cpua_original_2021$group,cpua_original_2021$survey_name,cpua_original_2021$species_id)
# fileout$mergeon = paste0(fileout$survstratum,fileout$year,fileout$fishlength_cm,fileout$dempel,fileout$order,fileout$group,fileout$survey_name,fileout$species_id)
# fileout$newcpua <- fileout$catcatchwgtswept; fileout$catcatchwgtswept <- NULL
# plotme=merge(fileout,cpua_original_2021,by="mergeon")
# plotme=plotme[!duplicated(plotme$mergeon),] # unique joins only, otherwise something has gone wrong in the setup of the join
# 
# # Reasonable correlation with previous
# plot(plotme$catcatchwgtswept[plotme$survey_name.x!='GNSIntOT1'],plotme$newcpua[plotme$survey_name.x!='GNSIntOT1'],ylim=c(0,50000), xlim=c(0,50000))
# 
# # GNSINTOT1 survey has a poorer agreement but still seems consistent with a corrected dataproduct
# # Drill into gnsintot1 - problem seems to be specific to 2018 and 2019 in this survey. Original cpua estimates for 2018/19 were too low, they must have been erroneous. The replacement is more consistent with other years.
# plot(plotme$catcatchwgtswept[plotme$survey_name.x=='GNSIntOT1'],plotme$newcpua[plotme$survey_name.x=='GNSIntOT1'],ylim=c(0,50000), xlim=c(0,50000))
# for(yr in unique(plotme$year.x)){plot(plotme$catcatchwgtswept[plotme$survey_name.x=='GNSIntOT1' & plotme$year.x==yr],plotme$newcpua[plotme$survey_name.x=='GNSIntOT1' & plotme$year.x==yr],ylim=c(0,10000), xlim=c(0,10000))}