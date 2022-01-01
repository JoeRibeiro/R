library(rfishbase)
RDIR<- "C:/Users/JR13/Desktop/Fish_dataproduct_QSR/SweptArea_29Oct2021/R/"
SPPFILE<-paste(RDIR,"SpeciesInfoSG.csv",sep="")
filein=read.csv(SPPFILE)

fromfishbase=rfishbase::length_weight(species_list = filein$CorrSciName, server='fishbase')

fishbaseLW=rfishbase::length_weight()
write.table(fishbaseLW,paste(RDIR,"rfishbase.txt",sep=""))
fishbaseLW=fishbaseLW[is.finite(fishbaseLW$a),]


# rfishbase is surprisingly much empty, trying FishLife.
# devtools::install_github("james-thorson/FishLife")
# library(FishLife)
Predict = Plot_taxa( Search_species(Genus="Acantholabrus",Species = "palloni")$match_taxonomy, mfrow=c(2,2) )
# Fishlife doesn't seem to return a and b factors in any of its functions.



  