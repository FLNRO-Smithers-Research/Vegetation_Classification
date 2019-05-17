####Script for Importing BEC vegetatation data. 
###for rule based classification,testing noise clustering, and creating edatopic grids
##Kiri Daust, July 2018
##MacKenzie, August 2018 extensive updates

.libPaths("E:/R packages351")
#install.packages("Hmisc")
require(reshape)
require(reshape2)
require(vegan)
require(caret)
require(tcltk)
require(randomForest)
require(Matrix)
require(labdsv)
require(gdata)
require(MASS)
require(openxlsx)
require (C50)
require(tidyr)
require(stringr)
require(rpart)
require(tree)
require(rattle)
require(rpart.plot)
require(partykit)
require(vegclust)
require(standardize)
require(dplyr)
require(tictoc)
require(plyr)
require(Hmisc)

rm(list=ls())
wd=tk_choose.dir(); setwd(wd)

####trialing Goldstream#############
goldstream <- function(x, type){ ##type either goldstream or prominence
  if(grepl("prom", type)){
    return(x[1]*sqrt(x[2]))
  }else{
    return(x[2]*sqrt(x[1]))
  }
}
#####IMPORT AND CLEAN DATA##############################
#################################################

#####################importing veg data########################################
############### Uses 4 column R export FORMAT FROM Vpro with Lifeform option selected (Plot - Species - Cover - Lifeform)
vegData <- read.table("BECMasterVeg_Feb7_2019.txt", header = TRUE) 
vegData <- separate(vegData, Species, c("Species","Type"), "-", remove = TRUE)
vegData <- mutate_all(vegData, funs(toupper)) ### converts lower case characters to upper
vegData$Cover <- as.numeric(vegData$Cover)
###return on significant decimal places
vegData$Cover <- round(vegData$Cover, digits = 3)
## remove species with zero cover
vegData <- vegData[vegData$Cover > 0,]
save(vegData, file = "VegDat_Raw.RData")##includes type field for lifeform
lifeform <- unique(vegData[,c(2:3)])
save(lifeform, file= "SppLifeForm.RData")###library list of species lifeform codes
load("VegDat_Raw.RData")


##########update old codes###################
masterList <- read.csv("USysAllSpecs.csv", stringsAsFactors = FALSE)
noMatch <- masterList[masterList$OldCode != masterList$Code,3:4]
temp <- merge(vegData,noMatch,by.x = "Species", by.y = "OldCode")
temp$Species <- temp$Code
temp <- temp[,-5]
vegData <- rbind(vegData,temp) ###Add section with new names
vegData <- vegData[!vegData$Species %in% noMatch$OldCode,] ##remove old codes

###remove codes not in master list
notIn <- vegData[!vegData$Species %in% masterList$Code,]
vegData <- vegData[vegData$Species %in% masterList$Code,]

if(length(notIn$Species) > 0){
  notIn <- dcast(notIn, PlotNumber ~ Species, value.var = "Species", fun.aggregate = length)
  write.csv(notIn, file = "CodesNotInMasterList.csv")
}
vegRemove <- vegData[vegData$Cover <= 0,] 
write.csv(vegRemove, file = "Plots_0Cover.csv") ## output of species with zero cover for review
##############Some stats on imported data
####Count number of columns (species) required
NTaxa = as.data.frame (length(unique(vegData$Species)))
NSpp <- as.data.frame (length(str_sub(vegData$Species,1,7)))
Taxa <- as.data.frame (unique(vegData$Species))
Spp <- as.data.frame (unique(str_sub(vegData$Species,1,7)))
Plots <- as.data.frame (unique(vegData$PlotNumber))
####Counts the number of instances of each unique species
CountsTaxa <- ddply(vegData,~Species,summarise,sppcount=length(Plots))
CountsSpp <- ddply(vegData,~str_sub(vegData$Species,1,7),summarise,sppcount=length(unique(PlotNumber)))
CountsSpp <- arrange(CountsSpp, -sppcount)
colnames (CountsSpp) [1] <- "Species"
rareSpp <- as.data.frame (CountsSpp[CountsSpp$sppcount <= 3,]) # set level for rare species removal
CountsSppReduced <- CountsSpp[!CountsSpp$Species %in% rareSpp$Species,]

#############Optional -- Remove subtaxa coding
  #vegData$Species <-str_sub(vegData$Species,1,7) ### adds field with only species codes (no spp or varieties)
  #vegData$Species <- as.factor(vegData$Species)

###Plot number of species in density graph
ggplot(data=CountsSpp,aes(x=reorder(Species, -sppcount), y=sppcount))+
 geom_point()
vegDatareduced <- vegData[!vegData$Species %in% rareSpp$Species,]
ggplot(vegDatareduced,aes(Species))+
  geom_density(adjust = 5)+
  theme(axis.text.x=element_text(angle = 90, hjust=1))
#####reduce vegData by eliminating rareSpp
vegData <- vegData[!vegData$Species %in% rareSpp$SppOnly,]

#####reduce vegData to only lifeforms listede (lifeform 1 and 2 = trees; 6 = grasses) only
#vegData <- vegData[vegData$Type %in% c(3,4,6),]
#vegData <- vegData[!is.na(vegData$Species),]
#treeSpp <- as.character(unique(vegData$Species))



#####################################################

#vegData$Species <- unlist(lapply(vegData$Species, toupper))
vegData3c <- vegData[,-3]##removes type field and setsback to 3-column format
save(vegData3c, file = "VegDat_Raw_3column.RData")
load("VegDat_Raw3column.RData")


vegData <- vegData[vegData$Cover > 0,]## remove records with zero cover
save(vegData, file = "VegDat_Clean.RData")
load("VegDat_Clean.RData")

###Optional application of lump species
lump <- read.csv("NewSppLump27Oct2018_Lump.csv", stringsAsFactors = FALSE)
lump$Lump <- unlist(lapply(lump$Lump, tolower))
lump <- lump[,1:2]
colnames(lump)[1:2] <- c("Lump","Species")
vegData <- merge(vegData, lump, by.x = "Species", all.x = TRUE) ##lump data
vegData$Species <- as.character(vegData$Species)
vegData$Species <- ifelse(!is.na(vegData$Lump), vegData$Lump, vegData$Species)
vegData <- vegData[-5]
save(vegData, file = "VegDat_Lumped.RData")
load("VegDat_Lumped.RData")


