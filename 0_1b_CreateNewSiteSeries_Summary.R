####Script for producing veg unit summary datasets from plot data and hierarchy units. Then various methods
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
require(ggplot2)
require(ggdendro)
require(pvclust)
require(dendextend)
require(ape)
require(doParallel)

rm(list=ls())
wd=tk_choose.dir(); setwd(wd)

###Import BECMaster vegetation plot data produced from 0_BEC_data_import_clean.R script
load("VegDat_Lumped.RData")
vegDat2 <- vegData
##Import SU table for new Site Series of Interest
NewSS <- read.csv("NewSiteSeries_SU.csv", stringsAsFactors = FALSE)

vegData2 <- merge(vegData, NewSS, by = "PlotNumber", all.x = TRUE)###Must select SiteSeries as level for now
#colnames(vegData)[5] <- "SiteSeries" ##depending on previous optional runs the column number will need to be changed
if(any(is.na(vegData2$SiteUnit))){
  warning("Data contains Plots not in hierachy table. These will be removed.")
}
vegData2 <- vegData2[!is.na(vegData2$SiteUnit),]
vegData2$PlotNumber <- as.character(vegData2$PlotNumber)
vegData2 <- vegData2[vegData2$SiteUnit != "",]

##roll up into site series summary data
set.seed(123321)
coreNo <- makeCluster(detectCores() - 1)
registerDoParallel(coreNo, cores = detectCores() - 1)
Cores <- as.numeric(detectCores()-1)
clusterEvalQ(coreNo, .libPaths("E:/R packages351"))

constCut <- 0 ##remove species less than cutoff
tic()
temp <- foreach(SS = unique(vegData2$SiteUnit), .combine = rbind, .packages = "foreach") %dopar% {
  sub <- vegData2[vegData2$SiteUnit == SS,]
  num <- length(unique(sub$PlotNumber))
  foreach(Spp = unique(sub$Species), .combine = rbind) %do% {
    sub2 <- sub[sub$Species == Spp,]
    numSpp <- dim(unique(sub2[,1:2]))[1]
    covsum <- sum(sub2$Cover)
    mean <- covsum/num
    const <- numSpp/num
    if(const >= constCut){
      out <- data.frame(SiteUnit = SS, Species = Spp, MeanCov = mean, Constancy = const*100, NoPlots = num)
    }
    
  }
}
toc()
NewSSSum <- temp
colnames(NewSSSum) [1] <- "SiteSeries"
save(NewSSSum, file = "NewSS_SummaryData.RData")








#####Statistics and Evaluations of new site series
NewSSUnits <- NewSSSum[,c(1,5)]
UnitPlots <- unique(NewSSUnits)# plots per unit
numspp <- ddply(NewSSSum,~SiteUnit,summarise,sppcount=length(unique(Species))) # number of species per unit
NewSSUnits2 <-merge(NewSSSum, numspp, by = "SiteUnit")
numspp2<-ddply(NewSSSum[NewSSSum$Constancy > 20,],~ SiteUnit,summarise,sppcount=length(unique(Species))) # number of non-accidental species per unit
NewSSUnits2 <-merge(NewSSUnits2, numspp2, by = "SiteUnit" )
constspp <-ddply(NewSSSum[NewSSSum$Constancy > 59,],~SiteUnit,summarise,sppcount=length(unique(Species))) # number of constant species per unit
NewSSUnits2 <-merge(NewSSUnits2, constspp, by= "SiteUnit" )
  
##### create potential differential valueslist
  #CovConst <- melt(temp)
  selectUnits <- as.data.frame(unique(NewSSSum[,1]))
len <- length(unique(NewSSSum[,1]))
load("SppLifeForm.RData")
differential <- foreach(rowNum = 1:len, .combine = rbind, .packages = c("foreach","reshape2")) %dopar% {
  select <- as.factor (selectUnits[rowNum,])
  CovTemp <- NewSSSum[NewSSSum$SiteSeries %in% select,]
  Cov1 <- CovTemp
  Cov1[is.na(Cov1)] <- 0
  
  # Cov1 <- Cov1[Cov1[,4] >= 60 & !is.na(Cov1[,4]),] ## remove non constant species
  #Potential differntial value by constancy
  Cov1$PotDiff <- ifelse(Cov1$Constancy >= 80,3,
                         ifelse(Cov1$Constancy >= 60,2,0))
  #add differential value by cover                                   
  Cov1$PotDiff <- ifelse(Cov1$Constancy < 60,0,
                         ifelse((Cov1$MeanCov <= 1 & Cov1$Constancy >= 60),(Cov1$PotDiff - 1),
                                ifelse((Cov1$MeanCov <= 10 & Cov1$Constancy >= 60),(Cov1$PotDiff),(Cov1$PotDiff + 7))))
  #reduce differential value for non-vascular spp by 1/2
  Cov1 <- merge(Cov1, lifeform, by = "Species", all.x = TRUE)
  Cov1$Type <- as.integer(Cov1$Type)
  Cov1$PotDiff <- ifelse((Cov1$Type >8 & Cov1$Type <12) ,(Cov1$PotDiff*0.5),(Cov1$PotDiff))
  #adjust differential value by constancy
  Cov1$PotDiff <- (Cov1$Constancy * Cov1$PotDiff)/100
  Cov1
  
}
differential[is.na(differential)] <- 0 ##without this trap some units end with N/A sum in summaryBy function
###Add Differential potential
diffSum <-summaryBy(PotDiff + NoPlots ~ SiteUnit, data = differential, FUN = c(length, sum))
diffSum[,c(5)] <- diffSum[,c(5)]/diffSum[,c(3)]
diffSum <- diffSum[,c(1,4)]
######
VegUnits2 <-merge(VegUnits2, diffSum, by= "SiteUnit" )
VegUnits2 <- VegUnits2[,c(1,5:9)]
colnames(VegUnits2)[3:5] <- c("TotalSpp","NonrareSpp","ConstSpp")
VegUnits2 <- unique(VegUnits2)
VegUnits2$Issues <- ifelse(VegUnits2$NoPlots < 6, "Low Plots", 
                           ifelse(VegUnits2$NoPlots <10, "Few Plots",
                                  ifelse(VegUnits2$ConstSpp <6, "Few Constants", "OK")))
write.csv(VegUnits2, paste(level,"_Unit_stats.csv", sep = ""))
save(VegUnits2, file = paste(level,"_SU_Stats.RData", sep = ""))
#load("Order_SummaryData.RData")

#### Optional remove units will less than X plots
#    SUsumData <- SUsumData[SUsumData$NoPlots > 3,]
##Retain only species with a constancy >= 60%
SUsumData <- NewSSSum
maxConst <- aggregate(Constancy ~ Species, SUsumData, FUN = max)
SpKeep <- maxConst[maxConst$Constancy >= 60,] ## identify species with >60% constancy in at least one unit
SpKeep <- as.character(SpKeep$Species)
#SpKeep <- SpKeep[SpKeep != "luzu+p"] ##this code causes issues further down - possibly the + sign
SUsumData$Species <- as.character(SUsumData$Species)
SUsumData <- SUsumData[SUsumData$Species %in% SpKeep,]
SUsumData$Species <- as.factor(SUsumData$Species)
#keep only left 4 columns of data
SUsumData<- SUsumData[,1:4]
SUcovData <- SUsumData[,1:3]
SUcovData <- SUcovData[SUcovData$MeanCov > 0.05,]
diffspp <- as.matrix(unique(SUsumData$Species)) ##list of species possibly differential in at least one unit
#factor (SUsumData$Group)
droplevels(SUsumData$SiteUnit, SUsumData$Species)
droplevels(SUcovData$SiteUnit, SUcovData$Species)
save(SUsumData, file = "SS_Differential_data.RData")
save(SUcovData, file = "SS_Differential_coverdata.RData") ###to be used for site series cluster analysis

 ###################################END of Summarizing Data#################################################################







########################################
#####SELECT LEVEL TO GROUP BY AND CONVERT TO MATRIX###############
SUsumData2 <- melt(SSsumData, id.vars = c("Group","Species"))
SUsumData3 <- dcast(SUsumData2, Group ~ Species + variable)###Convert to site by species matrix
##SUsumData <- dcast(SUsumData, PlotNumber ~ Species, value.var = "Cover", fun.aggregate = mean)
SUhier <- read.csv("GrasslandHierReduced.csv")
SUhier <- SUhierALL
level <- select.list(choices = colnames(SUhier), graphics = TRUE)###Select which level to test classification
SUhier <- SUhier[,c("11SiteSeries", level)]
SUhier <- unique(SUhier)
colnames(vegData)[1] <- "11SiteSeries"
vegData <- merge(vegData, SUhier, by = "11SiteSeries", all.x = TRUE)###Merge classification with matrix
vegData <- unique(vegData)
colnames(vegData)[length(vegData)] <- "Class"
if(any(is.na(vegData$Class))){
  warning("Data contains Plots not in hierachy table. These will be removed.")
}
vegData <- vegData[!is.na(vegData$Class),]###Remove missing groups
vegData <- vegData[,c(length(vegData),1:(length(vegData)-1))]
unique(vegData$Class)
##vegData <- vegData[!grepl("SUB|ALLIANCE",vegData$Class),] ##removes units from other classes
vegData$Class <- as.factor(vegData$Class)
vegData[is.na(vegData)] <- 0 ###set NAs to 0