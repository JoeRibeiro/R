#######################################################
###Script code et TL
#######################################################
library(tidyr)
library(dplyr)
library(ggplot2)
library(tRophicPosition)
library(nycflights13)
library(forcats)
library(labelled)
library(questionr)
library(skimr)
library(data.table)

setwd("C:/PRÁCTICAS/TFM 2021 Univ Brest/Datos/")


##Data
diet.dem <- read.table("diet_DEMERSALES.csv", sep = ";", dec = ".",header = T)
diet.ind <- read.table("diet_INDEMARES.csv", sep = ";", dec = ".",header = T)
diet.ind$area[diet.ind$area %in% "Avilés"]<-"Aviles"
diet.pel <- read.table("diet_PELACUS.csv", sep = ";", dec = ".",header = T)
diet.pel <- subset(diet.pel, select = -c(X, X.1))
tab1 <- read.table("Tabla IE2014-0600.csv", sep = ";", dec = ".",header = T)
tab2 <- read.table("Tabla IE2014-0601.csv", sep = ";", dec = ".",header = T)
tab.tl <- read.table("TL_complete.csv", sep = ";", dec = ".",header = T)
demersal <- read.table("isotopes_DEMERSALES.csv", sep = ";", dec = ".",header = T)
#demersal <- subset(demersal, select = -c(X))
indemares <- read.table("isotopes_INDEMARES.csv", sep = ";", dec = ".",header = T)
pelagic <- read.table("isotopes_PELACUS.csv", sep = ";", dec = ".",header = T)
pelagic$Delta.15N = as.numeric(pelagic$Delta.15N)


#####
###Pelagic :
str(pelagic)
names(pelagic)
Merlu.p <- subset(pelagic, pelagic$species == "Merluccius merluccius" & pelagic$Total.length < 17, select = survey : Delta.13C)
Merlu.p$species[Merlu.p$species %in% "Merluccius merluccius"]<-"Merluccius merluccius (< 17 cm)"
Merlu.m <- subset(pelagic, pelagic$species == "Merluccius merluccius" & pelagic$Total.length > 16 & pelagic$species == "Merluccius merluccius" & pelagic$Total.length < 35 , select = survey : Delta.13C)
Merlu.m$species[Merlu.m$species %in% "Merluccius merluccius"]<-"Merluccius merluccius (18 - 34cm)"
Merlu.g <- subset(pelagic, pelagic$species == "Merluccius merluccius" & pelagic$Total.length > 34, select = survey : Delta.13C)
Merlu.g$species[Merlu.g$species %in% "Merluccius merluccius"]<-"Merluccius merluccius ( > 35 cm)"

#Méthode supression plus courte avec substr
pelagic2.2 <- pelagic[substr(pelagic$species, 1, 126) != "Merluccius merluccius", ]

##Création tableau final avec le changement de noms
pelagic <- rbind(pelagic2.2, Merlu.p, Merlu.m, Merlu.g)


###Demersal
str(diet.dem)
names(diet.dem)
diet.merlu <- subset(diet.dem, diet.dem$consumer_sp == "Merluccius merluccius", select = haul : consumer_weight)
Merlu.p <- subset(diet.dem, diet.dem$consumer_sp == "Merluccius merluccius" & diet.dem$length < 17, select = haul : consumer_weight)
Merlu.p$consumer_sp[Merlu.p$consumer_sp %in% "Merluccius merluccius"]<-"Merluccius merluccius (< 17 cm)"

Merlu.m <- subset(diet.dem, diet.dem$consumer_sp == "Merluccius merluccius" & diet.dem$length > 16 & diet.dem$consumer_sp == "Merluccius merluccius" & diet.dem$length < 35 , select = haul : consumer_weight)
Merlu.m$consumer_sp[Merlu.m$consumer_sp %in% "Merluccius merluccius"]<-"Merluccius merluccius (18 - 34cm)"

Merlu.g <- subset(diet.dem, diet.dem$consumer_sp == "Merluccius merluccius" & diet.dem$length > 34, select = haul : consumer_weight)
Merlu.g$consumer_sp[Merlu.g$consumer_sp %in% "Merluccius merluccius"]<-"Merluccius merluccius ( > 35 cm)"

diet.dem2 <- diet.dem[substr(diet.dem$consumer_sp, 1, 33391) != "Merluccius merluccius", ]
diet.dem <- rbind(diet.dem2, Merlu.p, Merlu.m, Merlu.g)

rm(Merlu.g,Merlu.m, Merlu.p, diet.merlu, diet.dem2)

###Indemares
str(diet.ind)
names(diet.ind)
diet.merlu <- subset(diet.ind, diet.ind$consumer_sp == "Merluccius merluccius", select = area : consumer_weight)
Merlu.p <- subset(diet.ind, diet.ind$consumer_sp == "Merluccius merluccius" & diet.ind$length < 17, select = area : consumer_weight)
Merlu.p$consumer_sp[Merlu.p$consumer_sp %in% "Merluccius merluccius"]<-"Merluccius merluccius (< 17 cm)"
Merlu.g <- subset(diet.ind, diet.ind$consumer_sp == "Merluccius merluccius" & diet.ind$length > 34, select = area : consumer_weight)
Merlu.g$consumer_sp[Merlu.g$consumer_sp %in% "Merluccius merluccius"]<-"Merluccius merluccius ( > 35 cm)"

diet.ind2 <- diet.ind[substr(diet.ind$consumer_sp, 1, 1930) != "Merluccius merluccius", ]
diet.ind <- rbind(diet.ind2, Merlu.p, Merlu.g)

rm(Merlu.g, Merlu.p, diet.merlu, diet.ind2)

###Pelagic
str(diet.pel)
names(diet.pel)
diet.merlu <- subset(diet.pel, diet.pel$consumer_sp == "Merluccius merluccius", select = haul : consumer_weight)
Merlu.p <- subset(diet.pel, diet.pel$consumer_sp == "Merluccius merluccius" & diet.pel$length < 17, select = haul : consumer_weight)
Merlu.p$consumer_sp[Merlu.p$consumer_sp %in% "Merluccius merluccius"]<-"Merluccius merluccius (< 17 cm)"

Merlu.m <- subset(diet.pel, diet.pel$consumer_sp == "Merluccius merluccius" & diet.pel$length > 16 & diet.pel$consumer_sp == "Merluccius merluccius" & diet.pel$length < 35 , select = haul : consumer_weight)
Merlu.m$consumer_sp[Merlu.m$consumer_sp %in% "Merluccius merluccius"]<-"Merluccius merluccius (18 - 34cm)"

Merlu.g <- subset(diet.pel, diet.pel$consumer_sp == "Merluccius merluccius" & diet.pel$length > 34, select = haul : consumer_weight)
Merlu.g$consumer_sp[Merlu.g$consumer_sp %in% "Merluccius merluccius"]<-"Merluccius merluccius ( > 35 cm)"

diet.pel2 <- diet.pel[substr(diet.pel$consumer_sp, 1, 4323) != "Merluccius merluccius", ]

diet.pel <- rbind(diet.pel2, Merlu.p, Merlu.m, Merlu.g)

rm(Merlu.g, Merlu.m, Merlu.p, diet.merlu, diet.pel2)

#########################
######Etape 1 :
#########################
names(tab1)
names(tab2)

##Manipulating data
#Permet de separer une colonne en deux
tab2 <- separate(tab2,Muestra, c("Especie","Year"), sep = "_")

#moyenne POM et SOM
mean.15N <- tapply(tab1$Delta.15N,tab1$Especie,mean,na.rm = T)
mean.15N
mean.13C <- tapply(tab2$Delta.13C,tab2$Especie,mean,na.rm = T)
mean.13C
mean.iso <- rbind(mean.13C,mean.15N)
mean.iso <- t(mean.iso)
mean.iso = as.data.frame(mean.iso)


##Scatterplot only mean
#Traitement des data
mean.dem.C <- tapply(demersal$Delta.13C,demersal$species,mean,na.rm = T)
sd.dem.C <- tapply(demersal$Delta.13C,demersal$species,sd,na.rm = T)
mean.dem.N <- tapply(demersal$Delta.15N,demersal$species,mean,na.rm = T)
sd.dem.N <- tapply(demersal$Delta.15N,demersal$species,sd,na.rm = T)
demersal.mean.sd <- cbind(mean.dem.C,sd.dem.C,mean.dem.N,sd.dem.N)

demersal.mean.sd = as.data.frame(demersal.mean.sd)
demersal.mean.sd <- add_rownames(demersal.mean.sd, "Species")
names(demersal.mean.sd)


mean.pel.C <- tapply(pelagic$Delta.13C, pelagic$species,mean,na.rm = T)
sd.pel.C <- tapply(pelagic$Delta.13C, pelagic$species,sd,na.rm = T)
mean.pel.N <- tapply(pelagic$Delta.15N, pelagic$species,mean,na.rm = T)
sd.pel.N <- tapply(pelagic$Delta.15N,pelagic$species,sd,na.rm = T)
pelagic.mean.sd <- cbind(mean.pel.C,sd.pel.C,mean.pel.N,sd.pel.N)

pelagic.mean.sd = as.data.frame(pelagic.mean.sd)
pelagic.mean.sd <- add_rownames(pelagic.mean.sd, "Species")
names(pelagic.mean.sd)


mean.ind.C <- tapply(indemares$Delta.13C, indemares$species,mean,na.rm = T)
sd.ind.C <- tapply(indemares$Delta.13C, indemares$species,sd,na.rm = T)
mean.ind.N <- tapply(indemares$Delta.15N, indemares$species,mean,na.rm = T)
sd.ind.N <- tapply(indemares$Delta.15N, indemares$species,sd,na.rm = T)
indemares.mean.sd <- cbind(mean.ind.C,sd.ind.C,mean.ind.N,sd.ind.N)

indemares.mean.sd = as.data.frame(indemares.mean.sd)
indemares.mean.sd <- add_rownames(indemares.mean.sd, "species")
names(indemares.mean.sd)


########################
###Création d'une fonction permettant de calculer les Level trophic
TL <- function(x, POM, lamb = 1, delta = 3.4)
{
  p1 <- x - POM
  p2 <- p1/delta
  fin <- lamb + p2
  fin
}
#########################
##Calcul des TL
dem.TL <- TL(x = demersal$Delta.15N, POM = 4.535714 )
dem.TL
dem.TL = as.data.frame(dem.TL)
colnames(dem.TL)<- c("TL") 
demersal.2 <- cbind(demersal, dem.TL)


pel.TL <- TL(x = pelagic$Delta.15N, POM = 4.096667 )
pel.TL
pel.TL = as.data.frame(pel.TL)
colnames(pel.TL)<- c("TL") 
pelagic.2 <- cbind(pelagic, pel.TL)


ind.TL <- TL(x = indemares$Delta.15N, POM = 4.535714 )
ind.TL
ind.TL = as.data.frame(ind.TL)
colnames(ind.TL)<- c("TL") 
indemares.2 <- cbind(indemares, ind.TL)


#Supression des données inutiles
rm(dem.TL, demersal.mean.sd, ind.TL,indemares.mean.sd, mean.iso, pel.TL, pelagic.mean.sd, mean.13C,
   mean.15N, mean.dem.C, mean.dem.N, mean.ind.C, mean.ind.N, mean.pel.C, mean.pel.N, sd.dem.C, 
   sd.dem.N, sd.ind.C, sd.ind.N, sd.pel.C, sd.pel.N, demersal, pelagic, indemares, tab1, tab2)

############################################################
##Obtention tableau code et TL pour les isotopes

####Renommer des variables
indemares.2$zone[indemares.2$zone %in% "Avilés"]<-"Aviles"

###Retirer une lettre dans une colonne du jeu de données
indemares.2.2 <- separate(indemares.2, haul, c("r", "haul"), sep = "G")
indemares.2 <- subset(indemares.2.2, select = -c(r))

#Retirer une colonne
demersal.2 <- subset(demersal.2, select = -c(survey))

#Ajout de colonnes
demersal.2["survey"]= "DEMERSALES"
indemares.2["survey"]= "INDEMARES"

#####Fusion des colonnes pour créer une colonne code
code <- paste(demersal.2$haul, demersal.2$year, demersal.2$survey,sep = "_")
hauls.d3 <- cbind(code, demersal.2)

code <- paste(pelagic.2$haul, pelagic.2$year, pelagic.2$survey, sep = "_")
hauls.p3 <- cbind(code, pelagic.2)

code <- paste(indemares.2$haul, indemares.2$year, indemares.2$zone, indemares.2$survey, sep = "_")
hauls.ind3<- cbind(code, indemares.2)

#Fusionner en un seul tableau
names(hauls.d3)
names(hauls.p3)
names(hauls.ind3)
hauls.ind3 <- subset(hauls.ind3, select = -c(zone, Sex, Total.length, Preanal.length, Cephalotorax.length, depth, X))
hauls.p3 <- subset(hauls.p3, select = -c(Sex, Total.length))

tab.fin <- rbind(hauls.ind3, hauls.d3, hauls.p3)

write.table(tab.fin, "Data_iso_code_TL.csv", sep=";",row.names=FALSE)

rm(hauls.d3,hauls.ind3, hauls.p3, indemares.2.2, tab.fin)
#############################################################
######################
#######Etape 2 :
#####################

###################
####### Demersal dataset :
##################
##Manipulating data
names(diet.dem)

dietd <- diet.dem[substr(diet.dem$prey_sp, 1, 33391) != "", ]   #Supprimer les cases vides de la colonne prey_sp
dietd <- subset(dietd, select = c("haul", "year", "consumer", "prey_n", "consumer_sp", "volume_total", "prey_sp", "prey_percent"))
str(dietd)

#permet de rajouter les TL provenant d'un autre jeu de données
test <- left_join(dietd, tab.tl, by = c("prey_sp" = "species"))

#observation des NA et création d'un nouveau tableau sans NA
row.has.na <- apply(test, 1, function(x){any(is.na(x))}) 
sum(row.has.na)
test2 <- na.omit(test)

#Fusion de 3 colonnes en 1 :
test$ID <- paste(test$haul, test$year, test$consumer, sep ="_")

if(any(is.na(test$TL))) test$Presence_NA <-ifelse(is.na(test$TL), "TRUE", "FALSE")

#Retirer les individus avec des données NA
test$ID[test$Presence_NA == "TRUE"]
test3 <- test[!test$ID %in% test$ID[test$Presence_NA == "TRUE"],]

#Ajout d'une colonne avec un changement d'unité au niveau des pourcentages
Prey_mod_percent <- test3$prey_percent * 0.01
test3<- cbind(test3,Prey_mod_percent)
str(test3)

##Calcul du level trophic avec la formule
Mult <- test3$Prey_mod_percent * test3$TL
test.mult<- cbind(test3,Mult)
TL_calc <- tapply(test.mult$Mult, test.mult$ID, sum, na.rm =T)
TL_calc = as.data.frame(TL_calc)
TL_calc <- add_rownames(TL_calc, "ID")

TL.deme <- left_join(test3, TL_calc, by = "ID")
TL.deme <- left_join(TL_calc, test3, by = "ID")

#Identifier et retirer les lignes en double du tableau
TL.deme$Dup <- duplicated(TL.deme$ID)
TL.demersal <- TL.deme[!TL.deme$Dup, ]

#Supression des data inutiles
rm(Mult, test3, test, test2, test.mult, TL_calc, dietd, Prey_mod_percent, TL.deme)


############Creation tableau demersal TL estomac
#Ajout de colonnes
TL.demersal["Survey"]= "DEMERSALES"

#####Fusion des colonnes pour créer une colonne code
code <- paste(TL.demersal$haul, TL.demersal$year, TL.demersal$Survey,sep = "_")
hauls.d3 <- cbind(code, TL.demersal)

#Supprimer les colonnes inutiles
names(hauls.d3)
hauls.d3 <- subset(hauls.d3, select = -c(ID,TL, SE, Dup, Presence_NA))


##################
####Etape 3 :
##################
##Fusion des données des deux TL
#Data de l'étape 1 :
mean.TL.dem <- tapply(demersal.2$TL, demersal.2$species,mean,na.rm = T)
sd.TL.dem <- tapply(demersal.2$TL, demersal.2$species,sd,na.rm = T)
dem.TL.mean.sd <- cbind(mean.TL.dem,sd.TL.dem)
dem.TL.mean.sd = as.data.frame(dem.TL.mean.sd)
dem.TL.mean.sd <- add_rownames(dem.TL.mean.sd, "species")
dem.TL.mean.sd["TL"] = "TL_type_1"
colnames(dem.TL.mean.sd) <- c("species","TL","TL_sd", "TL_type")

#Data de l'étape 2 :
mean.TL2.dem <- tapply(TL.demersal$TL_calc, TL.demersal$consumer_sp,mean,na.rm = T)
sd.TL2.dem <- tapply(TL.demersal$TL_calc, TL.demersal$consumer_sp,sd,na.rm = T)
dem.TL2.mean.sd <- cbind(mean.TL2.dem,sd.TL2.dem)
dem.TL2.mean.sd = as.data.frame(dem.TL2.mean.sd)
dem.TL2.mean.sd <- add_rownames(dem.TL2.mean.sd, "species")
dem.TL2.mean.sd["TL"] = "TL_type_2"
colnames(dem.TL2.mean.sd) <- c("species","TL","TL_sd", "TL_type")


#permet de rajouter les TL provenant d'un autre jeu de données
tot <- rbind(dem.TL.mean.sd, dem.TL2.mean.sd)
intersect(dem.TL.mean.sd$species, dem.TL2.mean.sd$species)
#Com.sp.dem <- subset(tot, tot$Species == "Chelidonichthys cuculus" | tot$Species == "Conger conger" | tot$Species == "Helicolenus dactylopterus" | tot$Species == "Lepidorhombus boscii"|  tot$Species == "Lepidorhombus whiffiagonis"|  tot$Species == "Lepidotrigla dieuzeidei"|  tot$Species == "Micromesistius poutassou", select = Species : TL_type)

#Extrait espèces en commun de TL2 puis de TL
tet <- subset(dem.TL2.mean.sd, species %in% dem.TL.mean.sd$species)
tet2 <- subset(dem.TL.mean.sd, species %in% dem.TL2.mean.sd$species)
Com.sp <- rbind(tet2, tet)


#seulement les espèces en commun avec colonnes séparées
Com.sep <- merge(dem.TL.mean.sd, dem.TL2.mean.sd, by = "species")
colnames(Com.sep) <- c("species","TL_type1","TL_type1_sd", "TL_cat1", "TL_type2", "TL_type2_sd", "TL_cat2")

rm(diet.dem, Com.sp, dem.TL.mean.sd, dem.TL2.mean.sd, demersal.2, diet.dem, tet, tet2, TL.demersal, tot, Com.sep)
rm(mean.TL.dem, mean.TL2.dem, sd.TL.dem, sd.TL2.dem, row.has.na, TL, code)

######################################################################################################
###################
####### Pelagic dataset :
##################
##Manipulating data
names(diet.pel)

dietp <- diet.pel[substr(diet.pel$prey_sp, 1, 4323) != "", ]   #Supprimer les cases vides de la colonne prey_sp
dietp <- subset(dietp, select = c("haul", "year", "consumer", "prey_n", "consumer_sp", "volume_total", "prey_sp", "prey_percent"))
str(dietp)

#permet de rajouter les TL provenant d'un autre jeu de données
test <- left_join(dietp, tab.tl, by = c("prey_sp" = "species"))
names(test)
test <- subset(test, select = -c(X, X.1,X.2))

#observation des NA et création d'un nouveau tableau sans NA
row.has.na <- apply(test, 1, function(x){any(is.na(x))}) 
sum(row.has.na)
test2 <- na.omit(test)

#Fusion de 3 colonnes en 1 :
test$ID <- paste(test$haul, test$year, test$consumer, sep ="_")

if(any(is.na(test$TL))) test$Presence_NA <-ifelse(is.na(test$TL), "TRUE", "FALSE")

#Retirer les individus avec des données NA
test$ID[test$Presence_NA == "TRUE"]
test3 <- test[!test$ID %in% test$ID[test$Presence_NA == "TRUE"],]

#Ajout d'une colonne avec un changement d'unité au niveau des pourcentages
Prey_mod_percent <- test3$prey_percent * 0.01
test3<- cbind(test3,Prey_mod_percent)
str(test3)

##Calcul du level trophic avec la formule
Mult <- test3$Prey_mod_percent * test3$TL
test.mult<- cbind(test3,Mult)
TL_calc <- tapply(test.mult$Mult, test.mult$ID, sum, na.rm =T)
TL_calc = as.data.frame(TL_calc)
TL_calc <- add_rownames(TL_calc, "ID")

TL.pela <- left_join(test3, TL_calc, by = "ID")
TL.pela <- left_join(TL_calc, test3, by = "ID")

#Identifier et retirer les lignes en double du tableau
TL.pela$Dup <- duplicated(TL.pela$ID)
TL.pelagic <- TL.pela[!TL.pela$Dup, ]

#Supression des data inutiles
rm(Mult, test3, test, test2, test.mult, TL_calc, dietp, Prey_mod_percent, TL.pela)

############Creation tableau pelagique TL estomac
#Ajout de colonnes
TL.pelagic["Survey"]= "PELACUS"

#####Fusion des colonnes pour créer une colonne code
code <- paste(TL.pelagic$haul, TL.pelagic$year, TL.pelagic$Survey,sep = "_")
hauls.p3 <- cbind(code, TL.pelagic)

#Supprimer les colonnes inutiles
names(hauls.p3)
hauls.p3 <- subset(hauls.p3, select = -c(ID, TL, SE, Presence_NA, Dup))


##################
####Etape 3 :
##################
##Fusion des données des deux TL
#Data de l'étape 1 :
mean.TL.pel <- tapply(pelagic.2$TL, pelagic.2$species,mean,na.rm = T)
sd.TL.pel <- tapply(pelagic.2$TL, pelagic.2$species,sd,na.rm = T)
pel.TL.mean.sd <- cbind(mean.TL.pel,sd.TL.pel)
pel.TL.mean.sd = as.data.frame(pel.TL.mean.sd)
pel.TL.mean.sd <- add_rownames(pel.TL.mean.sd, "Species")
pel.TL.mean.sd["TL"] = "TL_type_1"
colnames(pel.TL.mean.sd) <- c("Species","TL","TL_sd", "TL_type")

#Data de l'étape 2 :
mean.TL2.pel <- tapply(TL.pelagic$TL_calc, TL.pelagic$consumer_sp,mean,na.rm = T)
sd.TL2.pel <- tapply(TL.pelagic$TL_calc, TL.pelagic$consumer_sp,sd,na.rm = T)
pel.TL2.mean.sd <- cbind(mean.TL2.pel,sd.TL2.pel)
pel.TL2.mean.sd = as.data.frame(pel.TL2.mean.sd)
pel.TL2.mean.sd <- add_rownames(pel.TL2.mean.sd, "Species")
pel.TL2.mean.sd["TL"] = "TL_type_2"
colnames(pel.TL2.mean.sd) <- c("Species","TL","TL_sd", "TL_type")


#permet de rajouter les TL provenant d'un autre jeu de données
tot <- rbind(pel.TL.mean.sd, pel.TL2.mean.sd)
intersect(pel.TL.mean.sd$Species, pel.TL2.mean.sd$Species)

#Extrait espèces en commun de TL2 puis de TL
tet <- subset(pel.TL2.mean.sd, Species %in% pel.TL.mean.sd$Species)
tet2 <- subset(pel.TL.mean.sd, Species %in% pel.TL2.mean.sd$Species)
Com.sp <- rbind(tet2, tet)

#seulement les espèces en commun avec colonnes séparées
Com.sep <- merge(pel.TL.mean.sd, pel.TL2.mean.sd, by = "Species")
colnames(Com.sep) <- c("Species","TL_type1","TL_type1_sd", "TL_cat1", "TL_type2", "TL_type2_sd", "TL_cat2")

rm(diet.pel, Com.sp, pel.TL.mean.sd, pel.TL2.mean.sd, pelagic.2, tet, tet2, TL.pelagic, tot, Com.sep)
rm(mean.TL.pel, mean.TL2.pel, sd.TL.pel, sd.TL2.pel, row.has.na, code)

######################################################################################################
###################
####### Indemares dataset :
##################
##Manipulating data
names(diet.ind)

dieti <- diet.ind[substr(diet.ind$prey_sp, 1, 1930) != "", ]   #Supprimer les cases vides de la colonne prey_sp
dieti <- subset(dieti, select = c("haul", "area", "year", "consumer", "prey_n", "consumer_sp", "volume_total", "prey_sp", "prey_percent"))
str(dieti)

#permet de rajouter les TL provenant d'un autre jeu de données
test <- left_join(dieti, tab.tl, by = c("prey_sp" = "species"))
names(test)
test <- subset(test, select = -c(X, X.1,X.2))

#observation des NA et création d'un nouveau tableau sans NA
row.has.na <- apply(test, 1, function(x){any(is.na(x))}) 
sum(row.has.na)
test2 <- na.omit(test)

#Fusion de 3 colonnes en 1 :
test$ID <- paste(test$haul, test$year, test$consumer, sep ="_")

if(any(is.na(test$TL))) test$Presence_NA <-ifelse(is.na(test$TL), "TRUE", "FALSE")

#Retirer les individus avec des données NA
test$ID[test$Presence_NA == "TRUE"]
test3 <- test[!test$ID %in% test$ID[test$Presence_NA == "TRUE"],]

#Ajout d'une colonne avec un changement d'unité au niveau des pourcentages
Prey_mod_percent <- test3$prey_percent * 0.01
test3<- cbind(test3,Prey_mod_percent)
str(test3)

##Calcul du level trophic avec la formule
Mult <- test3$Prey_mod_percent * test3$TL
test.mult<- cbind(test3,Mult)
TL_calc <- tapply(test.mult$Mult, test.mult$ID, sum, na.rm =T)
TL_calc = as.data.frame(TL_calc)
TL_calc <- add_rownames(TL_calc, "ID")

TL.inde <- left_join(test3, TL_calc, by = "ID")
TL.inde <- left_join(TL_calc, test3, by = "ID")

#Identifier et retirer les lignes en double du tableau
TL.inde$Dup <- duplicated(TL.inde$ID)
TL.indemares <- TL.inde[!TL.inde$Dup, ]

#Supression des data inutiles
rm(Mult, test3, test, test2, test.mult, TL_calc, dietp, Prey_mod_percent)

############Creation tableau indemares TL estomac
#Ajout de colonnes
TL.indemares["Survey"]= "INDEMARES"

#####Fusion des colonnes pour créer une colonne code
code <- paste(TL.indemares$haul, TL.indemares$year, TL.indemares$area, TL.indemares$Survey, sep = "_")
hauls.ind3<- cbind(code, TL.indemares)

#Supprimer les colonnes inutiles
names(hauls.ind3)
hauls.ind3 <- subset(hauls.ind3, select = -c(ID, TL, SE, Dup, Presence_NA, area))


##################
####Etape 3 :
##################
##Fusion des données des deux TL
#Data de l'étape 1 :
mean.TL.ind <- tapply(indemares.2$TL, indemares.2$species,mean,na.rm = T)
sd.TL.ind <- tapply(indemares.2$TL, indemares.2$species,sd,na.rm = T)
ind.TL.mean.sd <- cbind(mean.TL.ind,sd.TL.ind)
ind.TL.mean.sd = as.data.frame(ind.TL.mean.sd)
ind.TL.mean.sd <- add_rownames(ind.TL.mean.sd, "Species")
ind.TL.mean.sd["TL"] = "TL_type_1"
colnames(ind.TL.mean.sd) <- c("Species","TL","TL_sd", "TL_type")

#Data de l'étape 2 :
mean.TL2.ind <- tapply(TL.indemares$TL_calc, TL.indemares$consumer_sp,mean,na.rm = T)
sd.TL2.ind <- tapply(TL.indemares$TL_calc, TL.indemares$consumer_sp,sd,na.rm = T)
ind.TL2.mean.sd <- cbind(mean.TL2.ind,sd.TL2.ind)
ind.TL2.mean.sd = as.data.frame(ind.TL2.mean.sd)
ind.TL2.mean.sd <- add_rownames(ind.TL2.mean.sd, "Species")
ind.TL2.mean.sd["TL"] = "TL_type_2"
colnames(ind.TL2.mean.sd) <- c("Species","TL","TL_sd", "TL_type")


#permet de rajouter les TL provenant d'un autre jeu de données
tot <- rbind(ind.TL.mean.sd, ind.TL2.mean.sd)
intersect(ind.TL.mean.sd$Species, ind.TL2.mean.sd$Species)

#Extrait espèces en commun de TL2 puis de TL
tet <- subset(ind.TL2.mean.sd, Species %in% ind.TL.mean.sd$Species)
tet2 <- subset(ind.TL.mean.sd, Species %in% ind.TL2.mean.sd$Species)
Com.sp <- rbind(tet2, tet)

#seulement les espèces en commun avec colonnes séparées
Com.sep <- merge(ind.TL.mean.sd, ind.TL2.mean.sd, by = "Species")
colnames(Com.sep) <- c("Species","TL_type1","TL_type1_sd", "TL_cat1", "TL_type2", "TL_type2_sd", "TL_cat2")


#####################################
#########Fusion de tout en 1 tableau
############################

tab.conclu <- rbind(hauls.ind3, hauls.d3, hauls.p3)
write.table(tab.conclu, "Data_code_TL_stomac.csv", sep=";",row.names=FALSE)





