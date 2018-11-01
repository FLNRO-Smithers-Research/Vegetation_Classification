###This script imports vegetation data, and calculates pairwise diagnostics either between
###user specified units or all units to find bad associations.
###Kiri Daust, July 2018

.libPaths("E:/R packages351")
require(reshape2)
require(plyr)
require(dplyr)
require(tidyr)
require(ggplot2)
require(magrittr)
require(foreach)
require(tcltk)
require(openxlsx)
require(doParallel)
require(doBy)
require(doParallel)

#install.packages("tidyr")
rm(list=ls())
wd <- tk_choose.dir(); setwd(wd)


set.seed(123321)
coreNum <- as.numeric(detectCores()-1)
coreNo <- makeCluster(coreNum)
registerDoParallel(coreNo, cores = coreNum)
clusterEvalQ(coreNo, .libPaths("E:/R packages351"))

############################################################################################
#load Summary data from 01 summary script
load("Association_SummaryData.RData")
load("SppLifeForm.RData")
maxConst <- aggregate(Constancy ~ Species, SUsumData, FUN = max)
SpKeep <- maxConst[maxConst$Constancy >= 60,] ## identify species with >60% constancy in at least one unit
SpKeep <- as.character(SpKeep$Species)
#SpKeep <- SpKeep[SpKeep != "luzu+p"] ##this code causes issues further down - possibly the + sign
SUsumData$Species <- as.character(SUsumData$Species)
SUsumData <- SUsumData[SUsumData$Species %in% SpKeep,]
SUsumData$Species <- as.factor(SUsumData$Species)
#keep only left 4 columns of data
SUsumData<- SUsumData[,1:4]
diffspp <- as.matrix(unique(SUsumData$Species)) ##list of species possibly differential in at least one unit
#factor (SUsumData$SiteUnit)
droplevels(SUsumData$SiteUnit, SUsumData$Species)
save(SUsumData, file = paste(level,"_Differential_data.RData", sep = ""))
###Lookup tables
domDiffCls <- data.frame(SigDiff = c(9,8,7,6,5,4,3,2), DomDiff = c("dd1","dd1","dd1","dd1","dd1","dd2","dd3","dd4"))
scoreVals <- data.frame(Code = c("d1","d2","d3","dd1","dd2","dd3","dd4","c","cd","cm"),
                        Value = c(3,2,1,5,4,3,2,1,2,-1))
differential <- SUsumData[,c("SiteUnit", "Species", "MeanCov", "Constancy" )]
CovConst <- melt(differential)
selectUnits <- unique(as.character(CovConst$SiteUnit))
len <- length(selectUnits)
j = 23
k = 24

 ##Loop to calculate pairwise diagnostics for every possible combination (returns score for each pair)
out <- foreach(j = (1:(len-1)), .combine = rbind, .packages = c("foreach","reshape2")) %:% 
  foreach(k = ((j+1):len), .combine = rbind, .packages = c("foreach","reshape2")) %dopar% {
  select <- selectUnits[c(j,k)]
  CovTemp <- CovConst[CovConst$SiteUnit %in% select,] ##subset

  CovMatrix <- dcast(CovTemp, Species ~ SiteUnit + variable, value.var = "value", fun.aggregate = mean)
  CovMatrix[is.na(CovMatrix)] <- 0
 
  # Unit j vs k
  Cov1 <- CovMatrix[CovMatrix[,3] >= 60 & !is.na(CovMatrix[,3]),]###remove < 60% constancy for j unit
  Cov1$ConstDiff <- abs(Cov1[,3] - Cov1[,5])
  Cov1$Differential <- ifelse(Cov1$ConstDiff >= 80,"d1",
                              ifelse(Cov1$ConstDiff >= 60,"d2",
                                     ifelse(Cov1$ConstDiff >= 40,"d3",NA)))
  Cov1$SigA <- ifelse(Cov1[,2] <= 0.1,0,
                      ifelse(Cov1[,2] <= 0.3,1,
                             ifelse(Cov1[,2] <= 1,2,
                                    ifelse(Cov1[,2] <= 2.2,3,
                                           ifelse(Cov1[,2] <= 5,4,
                                                  ifelse(Cov1[,2] <= 10, 5,
                                                         ifelse(Cov1[,2] <= 20,6,
                                                                ifelse(Cov1[,2] <= 33,7,
                                                                       ifelse(Cov1[,2] <= 50,8,
                                                                              ifelse(Cov1[,2] <= 75,9,10))))))))))
  Cov1$SigB <- ifelse(Cov1[,4] <= 0.1,0,
                      ifelse(Cov1[,4] <= 0.3,1,
                             ifelse(Cov1[,4] <= 1,2,
                                    ifelse(Cov1[,4] <= 2.2,3,
                                           ifelse(Cov1[,4] <= 5,4,
                                                  ifelse(Cov1[,4] <= 10, 5,
                                                         ifelse(Cov1[,4] <= 20,6,
                                                                ifelse(Cov1[,4] <= 33,7,
                                                                       ifelse(Cov1[,4] <= 50,8,
                                                                              ifelse(Cov1[,4] <= 75,9,10))))))))))
  Cov1$SigDiff <- Cov1$SigA - Cov1$SigB
  Cov1 <- merge(Cov1, domDiffCls, by = "SigDiff", all.x = TRUE)
  #Cov1$DomDiff <- ifelse(Cov1$SigDiff < 6, NA, Cov1$SigDiff)
  Cov1$Const <- ifelse(Cov1$SigA >= 6,"cd",
                       ifelse(Cov1$SigA >= 3, "c","cm"))
  
  ###sum values
  Cov1$Value <- apply(Cov1[,c(8,11,12)],1,FUN = function(x){sum(scoreVals$Value[scoreVals$Code %in% x])})
  Cov1$Value <- ifelse(Cov1$Value>0 & is.na(Cov1$Differential), Cov1$Value -2, Cov1$Value)
  Cov1$Value <- ifelse(Cov1$Value<0, 0, Cov1$Value)
  Cov1$Value <- Cov1$Value*(Cov1[,4]/100)
  Cov1 <- merge(Cov1, lifeform, by = "Species", all.x = TRUE)
  Cov1$Value <- ifelse(Cov1$Type %in% c(9,10,11,13), Cov1$Value/2, Cov1$Value)
  Cov1$Value <- ifelse(is.na(Cov1$Differential ) & is.na(Cov1$DomDiff), 0, Cov1$Value )
  sumA <- sum(Cov1$Value)
 
   # Unit k vs j 
  Cov2 <- CovMatrix[CovMatrix[,5] >= 60 & !is.na(CovMatrix[,5]),]###remove < 60% constancy for k unit
  Cov2$ConstDiff <- abs(Cov2[,5] - Cov2[,3])
  Cov2$Differential <- ifelse(Cov2$ConstDiff >= 80,"d1",
                              ifelse(Cov2$ConstDiff >= 60,"d2",
                                     ifelse(Cov2$ConstDiff >= 40,"d3",NA)))
  Cov2$SigA <- ifelse(Cov2[,2] <= 0.1,0,
                      ifelse(Cov2[,2] <= 0.3,1,
                             ifelse(Cov2[,2] <= 1,2,
                                    ifelse(Cov2[,2] <= 2.2,3,
                                           ifelse(Cov2[,2] <= 5,4,
                                                  ifelse(Cov2[,2] <= 10, 5,
                                                         ifelse(Cov2[,2] <= 20,6,
                                                                ifelse(Cov2[,2] <= 33,7,
                                                                       ifelse(Cov2[,2] <= 50,8,
                                                                              ifelse(Cov2[,2] <= 75,9,10))))))))))
  Cov2$SigB <- ifelse(Cov2[,4] <= 0.1,0,
                      ifelse(Cov2[,4] <= 0.3,1,
                             ifelse(Cov2[,4] <= 1,2,
                                    ifelse(Cov2[,4] <= 2.2,3,
                                           ifelse(Cov2[,4] <= 5,4,
                                                  ifelse(Cov2[,4] <= 10, 5,
                                                         ifelse(Cov2[,4] <= 20,6,
                                                                ifelse(Cov2[,4] <= 33,7,
                                                                       ifelse(Cov2[,4] <= 50,8,
                                                                              ifelse(Cov2[,4] <= 75,9,10))))))))))
  Cov2$SigDiff <- Cov2$SigB - Cov2$SigA
  Cov2 <- merge(Cov2, domDiffCls, by = "SigDiff", all.x = TRUE)
  #Cov2$DomDiff <- ifelse(Cov2$SigDiff < 6, NA, Cov2$SigDiff)
  Cov2$Const <- ifelse(Cov2$SigB >= 6,"cd",
                       ifelse(Cov2$SigB >= 3, "c","cm"))
  
  ###sum values
  Cov2$Value <- apply(Cov2[,c(8,11,12)],1,FUN = function(x){sum(scoreVals$Value[scoreVals$Code %in% x])})
  Cov2$Value <- ifelse(Cov2$Value>0 & is.na(Cov2$Differential), Cov2$Value -2, Cov2$Value)
  Cov2$Value <- ifelse(Cov2$Value<0, 0, Cov2$Value)
  Cov2$Value <- Cov2$Value*(Cov2[,6]/100)
  Cov2 <- merge(Cov2, lifeform, by = "Species", all.x = TRUE)
  Cov2$Value <- ifelse(Cov2$Type %in% c(9,10,11,13), Cov2$Value/2, Cov2$Value)
  Cov2$Value <- ifelse(is.na(Cov2$Differential ) & is.na(Cov2$DomDiff), 0, Cov2$Value )
  sumB <- sum(Cov2$Value)
  
  totDiff <- sumA+sumB
  
  outTemp <- data.frame(SiteUnits = paste(select[1],"|",select[2]), Score = totDiff)
  outTemp
  
  }  
   
 

badAssoc <- out[out$Score < 15,]

badAssoc <- separate(badAssoc, SiteUnits, c("G1","G2"), " \\| ", remove = TRUE)
out <- separate(out, SiteUnits, c("G1","G2"), " \\| ", remove = TRUE)
save(out, file = "AllForestSiteSeriesScores.RData")
write.csv(badAssoc, "AllForestSiteSeriesDifferentialSum_Low.csv", row.names = FALSE)
write.csv(out, "AllForestSiteSeriesDifferentialSum_ALL.csv", row.names = FALSE)
          ###identify those subassociations that are similar but in different associations
##Import hierachy table
SUhier <- read.csv("AllForestHier.csv", stringsAsFactors = FALSE)
colnames(SUhier)[1:12]=c("PlotNumber", "Region", "Class", "Order", "SubOrder", "Alliance", "SubAlliance", "Association", "SubAssociation", "Facies", "Working", "SiteSeries")
hierunits <-unique (SUhier[,c(2:12)])

colnames(out)[1] <-c("Association")
out2 <- merge(out, hierunits[7:8], by = "Association", all = FALSE)
out2 <- out2[out2$Score<10,] #retain only closely related units
colnames(out2)[1:4] <-c("SubAssociation1","SubAssociation" , "Score", "Association1")
out3 <- merge(out2, hierunits[7:8], by = "Association", all = FALSE)
out3$movetosubass <-ifelse(out3$Association1 != out3$Association, 1, 0 )
out3 <-out3[(out3$movetosubass == 1),]
write.csv(out3, "GrasslandsPairedUnits_tomoveAssociations.csv", row.names = FALSE)

          ###identify those siteseries that are similar but in different associations
colnames(out)[1] <-c("SiteSeries")
out2 <- merge(out, hierunits[7:11], by = "SiteSeries", all = FALSE)
out2 <- out2[out2$Score<10,] #retain only closely related units
colnames(out2)[1:4] <-c("SiteSeries1","SiteSeries" , "Score", "Association1")
out3 <- merge(out2, hierunits[7:11], by = "SiteSeries", all = FALSE)
out3$movetosubass <-ifelse(out3$Association1 != out3$Association, 1, 0 )
out3 <-out3[(out3$movetosubass == 1),]
out3 <- out3[,c("SiteSeries", "SiteSeries1", "Score", "Association", "Association1", "SubAssociation.x", "SubAssociation.y")]
write.csv(out3, "ForestSiteSeries_tomoveAssociations.csv", row.names = FALSE)


##How many bad associations does each group have?
lenG1 <- aggregate(Score ~ G1, badAssoc, FUN = length)
lenG2 <- aggregate(Score ~ G2, badAssoc, FUN = length)
len <- merge(lenG1, lenG2, by.x = "G1", by.y = "G2", all = TRUE)
len$Total <- apply(len[,2:3],1,FUN = sum, na.rm = TRUE)
len <- len[,-(2:3)]
write.csv(len, "GrasslandnumBad1.csv", row.names = FALSE)
bad <- unique(badAssoc[,1:2])

###Function to determine closest parent branch between two bad associations
hierUnit <- function(x,hierClass){
  if(hierClass$Association[hierClass$SubAssoc == x[1]] == hierClass$Association[hierClass$SubAssoc == x[2]]){
    return(hierClass$Association[hierClass$SubAssoc == x])
  }else if(hierClass$Alliance[hierClass$SubAssoc == x[1]] == hierClass$Alliance[hierClass$SubAssoc == x[2]]){
    return(hierClass$Alliance[hierClass$SubAssoc == x])
  }else if(hierClass$Order[hierClass$SubAssoc == x[1]] == hierClass$Order[hierClass$SubAssoc == x[2]]){
    return(hierClass$Order[hierClass$SubAssoc == x])
  }else if(hierClass$Class[hierClass$SubAssoc == x[1]] == hierClass$Class[hierClass$SubAssoc == x[2]]){
    return(hierClass$Class[hierClass$SubAssoc == x])
  }else{
    return("Not")
  }
} 

badAssoc$HierSiteUnit <- apply(badAssoc[,1:2],1,FUN = hierUnit, hierClass = hierClass)
write.csv(badAssoc, "NearestHierarchyJoin.csv", row.names = FALSE)
####Display full data for each bad association
###for checking specific non-bad groups, create dataframe called "bad" with columns
####G1 and G2 which have the names of the units to compare. 
fullData <- foreach(rowNum = 1:length(bad$G1), .combine = rbind, .packages = c("foreach","reshape2")) %dopar% {
    select <- as.character(bad[rowNum,])
    CovTemp <- CovConst[CovConst$SiteUnit %in% select,]
    
    Cov1 <- dcast(CovTemp, Species ~ SiteUnit + variable, value.var = "value")
    Cov1[is.na(Cov1)] <- 0
    
    Cov1 <- Cov1[Cov1[,3] >= 60 & !is.na(Cov1[,3]),]
    Cov1$ConstDiff <- abs(Cov1[,3] - Cov1[,6])
    Cov1$Differential <- ifelse(Cov1$ConstDiff >= 80,"d1",
                                ifelse(Cov1$ConstDiff >= 60,"d2",
                                       ifelse(Cov1$ConstDiff >= 40,"d3",NA)))
    Cov1$SigA <- ifelse(Cov1[,2] <= 0.1,0,
                        ifelse(Cov1[,2] <= 0.3,1,
                               ifelse(Cov1[,2] <= 1,2,
                                      ifelse(Cov1[,2] <= 2.2,3,
                                             ifelse(Cov1[,2] <= 5,4,
                                                    ifelse(Cov1[,2] <= 10, 5,
                                                           ifelse(Cov1[,2] <= 20,6,
                                                                  ifelse(Cov1[,2] <= 33,7,
                                                                         ifelse(Cov1[,2] <= 50,8,
                                                                                ifelse(Cov1[,2] <= 75,9,10))))))))))
    Cov1$SigB <- ifelse(Cov1[,5] <= 0.1,0,
                        ifelse(Cov1[,5] <= 0.3,1,
                               ifelse(Cov1[,5] <= 1,2,
                                      ifelse(Cov1[,5] <= 2.2,3,
                                             ifelse(Cov1[,5] <= 5,4,
                                                    ifelse(Cov1[,5] <= 10, 5,
                                                           ifelse(Cov1[,5] <= 20,6,
                                                                  ifelse(Cov1[,5] <= 33,7,
                                                                         ifelse(Cov1[,5] <= 50,8,
                                                                                ifelse(Cov1[,5] <= 75,9,10))))))))))
    Cov1$SigDiff <- Cov1$SigA - Cov1$SigB
    Cov1 <- merge(Cov1, domDiffCls, by = "SigDiff", all.x = TRUE)
    Cov1$DomDiff <- ifelse(Cov1$SigDiff < 6, NA, Cov1$SigDiff)
    Cov1$Const <- ifelse(Cov1$SigA >= 6,"cd",
                         ifelse(Cov1$SigA >= 3, "c","cm"))
    
    Cov1$Value <- apply(Cov1[,c(10,13,14)],1,FUN = function(x){sum(scoreVals$Value[scoreVals$Code %in% x])})
    Cov1$Value <- Cov1$Value*(Cov1[,4]/100)
    Cov1 <- merge(Cov1, typeCodes, by = "Species", all.x = TRUE)
    Cov1$Value <- ifelse(Cov1$Type %in% c(9,10,11,13), Cov1$Value/2, Cov1$Value)
    colnames(Cov1)[3:4] <- c("G1Cov","G1Const")
    colnames(Cov1)[6:7] <- c("G2Cov","G2Const")
    Cov1$G1 <- select[1]
    Cov1$G2 <- select[2]
    Cov1
    Cov1 <- Cov1[c("Species", "G1","G1Cov","G1Const", "G2","G2Cov","G2Const")]   

}

write.csv(fullData, "CompareUnitsTooSimilar.csv", row.names = FALSE)

##############Not Required#################

##### create potential differential list

differential <- foreach(rowNum = 1:length(bad$G1), .combine = rbind, .packages = c("foreach","reshape2")) %dopar% {
  select <- as.character(bad[rowNum,])
  CovTemp <- CovConst[CovConst$SiteUnit %in% select,]
  
  Cov1 <- dcast(CovTemp, Species ~ SiteUnit + variable, value.var = "value")
  Cov1[is.na(Cov1)] <- 0
  
  Cov1 <- Cov1[Cov1[,3] >= 60 & !is.na(Cov1[,3]),]
  Cov1$ConstDiff <- abs(Cov1[,3] - Cov1[,5])
  Cov1$Differential <- ifelse(Cov1$ConstDiff >= 80,"d1",
                              ifelse(Cov1$ConstDiff >= 60,"d2",
                                     ifelse(Cov1$ConstDiff >= 40,"d3",NA)))
  Cov1$SigA <- ifelse(Cov1[,2] <= 0.1,0,
                      ifelse(Cov1[,2] <= 0.3,1,
                             ifelse(Cov1[,2] <= 1,2,
                                    ifelse(Cov1[,2] <= 2.2,3,
                                           ifelse(Cov1[,2] <= 5,4,
                                                  ifelse(Cov1[,2] <= 10, 5,
                                                         ifelse(Cov1[,2] <= 20,6,
                                                                ifelse(Cov1[,2] <= 33,7,
                                                                       ifelse(Cov1[,2] <= 50,8,
                                                                              ifelse(Cov1[,2] <= 75,9,10))))))))))
  Cov1$SigB <- ifelse(Cov1[,4] <= 0.1,0,
                      ifelse(Cov1[,4] <= 0.3,1,
                             ifelse(Cov1[,4] <= 1,2,
                                    ifelse(Cov1[,4] <= 2.2,3,
                                           ifelse(Cov1[,4] <= 5,4,
                                                  ifelse(Cov1[,4] <= 10, 5,
                                                         ifelse(Cov1[,4] <= 20,6,
                                                                ifelse(Cov1[,4] <= 33,7,
                                                                       ifelse(Cov1[,4] <= 50,8,
                                                                              ifelse(Cov1[,4] <= 75,9,10))))))))))
  Cov1$SigDiff <- Cov1$SigA - Cov1$SigB
  Cov1 <- merge(Cov1, domDiffCls, by = "SigDiff", all.x = TRUE)
  Cov1$DomDiff <- ifelse(Cov1$SigDiff < 6, NA, Cov1$SigDiff)
  Cov1$Const <- ifelse(Cov1$SigA >= 6,"cd",
                       ifelse(Cov1$SigA >= 3, "c","cm"))
  
  Cov1$Value <- apply(Cov1[,c(8,11,12)],1,FUN = function(x){sum(scoreVals$Value[scoreVals$Code %in% x])})
  Cov1$Value <- Cov1$Value*(Cov1[,4]/100)
  Cov1 <- merge(Cov1, typeCodes, by = "Species", all.x = TRUE)
  Cov1$Value <- ifelse(Cov1$Type %in% c(9,10,11,13), Cov1$Value/2, Cov1$Value)
  colnames(Cov1)[3:6] <- c("G1Cov","G1Const","G2Cov","G2Const")
  Cov1$G1 <- select[1]
  Cov1$G2 <- select[2]
  Cov1
  
}

differential <- differential[,c(15,16,1,3:6,2,8:14)]
