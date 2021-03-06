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
require(ape)
require(dendextend)
require (cba)

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
load("SS_Differential_data.RData")
colnames (SUsumData) [1] <- "SiteUnit"
load("SppLifeForm.RData")
SUsumData <- SUsumData[!grepl("[?]" ,SUsumData$SiteUnit),] ##removes units flagged as questionable
SUsumData <- SUsumData[!grepl("%" ,SUsumData$SiteUnit),] # removes units flagged as provisional
SUsumData <- SUsumData[!grepl("[$]" ,SUsumData$SiteUnit),] #removes units flagged as seral
#SpKeep <- SpKeep[SpKeep != "luzu+p"]
maxConst <- aggregate(Constancy ~ Species, SUsumData, FUN = max)
SpKeep <- maxConst[maxConst$Constancy >= 60,] ## identify species with >60% constancy in at least one unit
SpKeep <- as.character(SpKeep$Species)
#SpKeep <- SpKeep[SpKeep != "luzu+p"] ##this code causes issues further down - possibly the + sign
SUsumData$Species <- as.character(SUsumData$Species)
SUsumData <- SUsumData[SUsumData$Species %in% SpKeep,]
SUsumData$Species <- as.factor(SUsumData$Species)
#keep only left 4 columns of data
SUsumData<- SUsumData[,1:4]
SUsumData[3:4] <-round(SUsumData[3:4],2)
diffspp <- as.matrix(unique(SUsumData$Species)) ##list of species possibly differential in at least one unit
#factor (SUsumData$SiteUnit)
droplevels(SUsumData$SiteUnit, SUsumData$Species)
save(SUsumData, file = paste(level,"_Differential_data.RData", sep = ""))
###Lookup tables
domDiffCls <- data.frame(SigDiff = c(9,8,7,6,5,4,3,2), DomDiff = c("dd1","dd1","dd1","dd1","dd1","dd2","dd3","dd4"))
scoreVals <- data.frame(Code = c("d1","d2","d3","dd1","dd2","dd3","dd4","cd","c","cm"),
                        Value = c(4,3,1,0,0,0,0,0,0,0))
differential <- SUsumData[,c("SiteUnit", "Species", "MeanCov", "Constancy" )]
CovConst <- melt(differential)
selectUnits <- unique(as.character(CovConst$SiteUnit))
selectUnits[sort.list(selectUnits)]
len <- length(selectUnits)
jname  <- select.list(choices = selectUnits[sort.list(selectUnits)], graphics = TRUE)
j <- which(selectUnits == jname)[1]
kname = select.list(choices = selectUnits[sort.list(selectUnits)], graphics = TRUE)
k <- which(selectUnits == kname)[1]
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
 ####Points for differntial based on constancy differernces
   Cov1$Differential <- ifelse(Cov1$ConstDiff >= 80, "d1",
                              ifelse(Cov1$ConstDiff >= 60, "d2",
                                     ifelse(Cov1$ConstDiff >= 40, "d3" ,NA)))

  #####Points for dominant differential 
  Cov1$CovDiff <- ifelse (Cov1[,2] >10, (Cov1[,2] / Cov1[,4])/10, 
                          ifelse(Cov1[,2] >5 , ((Cov1[,2] / (Cov1[,4]*1.5) / 10)), 0))
  Cov1$CovDiff <- ifelse(Cov1$CovDiff >1, 1, 
                         ifelse(Cov1$CovDiff <0.2, 0, Cov1$CovDiff ))
  Cov1$DDpoints <- ifelse(Cov1[,2] <10, (Cov1$CovDiff * 4) *(Cov1[,2]/10),Cov1$CovDiff * 4)## adjust for cover
  Cov1$DDclass <- ifelse(Cov1$DDpoints >=4, "dd1",
                       ifelse(Cov1$DDpoints >=2, "dd2",
                              ifelse(Cov1$DDpoints >=1.2, "dd3",
                                     ifelse(Cov1$DDpoints > 0, "dd4", NA))))
  Cov1$Const <- ifelse(Cov1[,2] >= 10,"cd",
  ifelse(Cov1[,2] >= 0.3, "c","cm"))
  ###sum values
  Cov1$Value <- apply(Cov1[,c("Differential","Const","DDpoints")],1,FUN = function(x){sum(scoreVals$Value[scoreVals$Code %in% x])})
 Cov1$Value <- Cov1$Value + Cov1$DDpoints
 Cov1$Value <- ifelse (((!is.na(Cov1$Differential)  & (!is.na(Cov1$DDclass)))), Cov1$Value + 2, Cov1$Value)
 Cov1$Value <- Cov1$Value * (Cov1[,3]/100) ### adjust for constancy
 #  Cov1$Value <- ifelse(Cov1$Value>0 & is.na(Cov1$Differential), Cov1$Value -2, Cov1$Value)
 # Cov1$Value <- ifelse(Cov1$Value<0, 0, Cov1$Value)
 # Cov1$Value <- Cov1$Value*(Cov1[,3]/100)
 # Cov1$Value <- Cov1$Value + Cov1$DDpoints
  Cov1 <- merge(Cov1, lifeform, by = "Species", all.x = TRUE)
  Cov1$Value <- ifelse(Cov1$Type %in% c(9,10,11,13), Cov1$Value/2, Cov1$Value)
  Cov1$Value <- ifelse(is.na(Cov1$Differential ) & is.na(Cov1$DDclass), 0, Cov1$Value )

   sumA <- sum(Cov1$Value)

   # Unit k vs j 
   Cov2 <- CovMatrix[CovMatrix[,5] >= 60 & !is.na(CovMatrix[,5]),]###remove < 60% constancy for j unit
   Cov2$ConstDiff <- abs(Cov2[,5] - Cov2[,3])
   ####Points for differntial based on constancy differernces
   Cov2$Differential <- ifelse(Cov2$ConstDiff >= 80, "d1",
                               ifelse(Cov2$ConstDiff >= 60, "d2",
                                      ifelse(Cov2$ConstDiff >= 40, "d3" ,NA)))
   
   #####Points for dominant differential 
   Cov2$CovDiff <- ifelse (Cov2[,4] >10, (Cov2[,4] / Cov2[,2])/10, 
                           ifelse(Cov2[,4] >5 , ((Cov2[,4] / (Cov2[,2]*1.5) / 10)), 0))
   Cov2$CovDiff <- ifelse(Cov2$CovDiff >1, 1, 
                          ifelse(Cov2$CovDiff <0.2, 0, Cov2$CovDiff ))
   Cov2$DDpoints <- ifelse(Cov2[,4] <10, (Cov2$CovDiff * 4) *(Cov2[,4]/10),Cov2$CovDiff * 4)
   Cov2$DDclass <- ifelse(Cov2$DDpoints >=4, "dd1",
                          ifelse(Cov2$DDpoints >=2, "dd2",
                                 ifelse(Cov2$DDpoints >=1.2, "dd3",
                                        ifelse(Cov2$DDpoints > 0, "dd4", NA))))
   Cov2$Const <- ifelse(Cov2[,4] >= 10,"cd",
                        ifelse(Cov2[,4] >= 0.3, "c","cm"))
   ###sum values
   Cov2$Value <- apply(Cov2[,c("Differential","Const","DDpoints")],1,FUN = function(x){sum(scoreVals$Value[scoreVals$Code %in% x])})
   Cov2$Value <- Cov2$Value + Cov2$DDpoints
   Cov2$Value <- ifelse (((!is.na(Cov2$Differential)  & (!is.na(Cov2$DDclass)))), Cov2$Value + 2, Cov2$Value)
   Cov2$Value <- Cov2$Value * (Cov2[,5]/100) ### adjust for constancy
   #  Cov2$Value <- ifelse(Cov2$Value>0 & is.na(Cov2$Differential), Cov2$Value -2, Cov2$Value)
   # Cov2$Value <- ifelse(Cov2$Value<0, 0, Cov2$Value)
   # Cov2$Value <- Cov2$Value*(Cov2[,5]/100)
   # Cov2$Value <- Cov2$Value + Cov2$DDpoints
   Cov2 <- merge(Cov2, lifeform, by = "Species", all.x = TRUE)
   Cov2$Value <- ifelse(Cov2$Type %in% c(9,10,11,13), Cov2$Value/2, Cov2$Value)
   Cov2$Value <- ifelse(is.na(Cov2$Differential ) & is.na(Cov2$DDclass), 0, Cov2$Value )
   
   sumB <- sum(Cov2$Value)
  
  totDiff <- sumA+sumB ####total diagnostic sum
  
  #########Ratio of constants to differentials
  allconst <- rbind(Cov1, Cov2)
  countconst <- length (unique(allconst$Species))
  countdiag <- allconst[allconst$Value > 0,]
  countdiag <- length (countdiag$Value)
  const_ratio <- 1-(countdiag/countconst)
  
  #################Ratio of potential diagnostic to actual diagnostic

  Cov1$PotDiffPoints <- ifelse(Cov1[,3] >= 80, 4,
                              ifelse(Cov1[,3] >= 60, 3 ,
                                     ifelse(Cov1[,3] >= 40, 1 ,0)))
  
  Cov1$PotDomPoints <- ifelse (Cov1[,2] >10, 6, 
                          ifelse(Cov1[,2] >5 , ((6 * (Cov1[,2])/ 10)), 0))

  Cov1$PotPoints <- Cov1$PotDiffPoints + Cov1$PotDomPoints
  Cov1$PotPoints <- ifelse(Cov1$Type %in% c(9,10,11,13), Cov1$PotPoints/2, Cov1$PotPoints)

  
  Cov2$PotDiffPoints <- ifelse(Cov2[,5] >= 80, 4,
                               ifelse(Cov2[,5] >= 60, 3 ,
                                      ifelse(Cov2[,5] >= 40, 1 ,0)))
  
  Cov2$PotDomPoints <- ifelse (Cov2[,4] >10, 6, 
                               ifelse(Cov2[,4] >5 , ((6 * (Cov2[,4])/ 10)), 0))
  
  Cov2$PotPoints <- Cov2$PotDiffPoints + Cov2$PotDomPoints
  Cov2$PotPoints <- ifelse(Cov2$Type %in% c(9,10,11,13), Cov2$PotPoints/2, Cov2$PotPoints)
 
  PotCompareA <- Cov1[,c(1,16)]
  PotCompareB <- Cov2[,c(1,16)]
  PotCompare <- merge(PotCompareA, PotCompareB, by = "Species", all=TRUE)
  PotCompare [is.na(PotCompare)] <-0
  PotCompare$PotMax <- pmax(PotCompare$PotPoints.x, PotCompare$PotPoints.y)
  PotCompare$PotMean <- rowMeans(PotCompare[2:3])
  PotCompare$PotDiff <- abs(PotCompare$PotPoints.x - PotCompare$PotPoints.y)
  potsummax <- sum(PotCompare$PotMax)  
  potsummean <- sum (PotCompare$PotMean)
  ConstmaxDiff_ratio <- 1-(totDiff/potsummax)
  ConstmeanDiff_ratio <- 1-(totDiff/potsummean)
  potdiffdiff <- sum (PotCompare$PotDiff)
  actpotratio <- totDiff/potdiffdiff
  outTemp <- data.frame(SiteUnits = paste(select[1],"|",select[2]), Score = totDiff, DifferenceinPotential = potdiffdiff, ActualPotential = actpotratio,  ConstRatio = const_ratio,   MeanPotPoints = potsummean, ConstDiffmax = ConstmaxDiff_ratio, ConstDiffmean = ConstmeanDiff_ratio)
  outTemp
  
  }  
   
 

badAssoc <- out[out$Score < 15,]

badAssoc <- separate(badAssoc, SiteUnits, c("G1","G2"), " \\| ", remove = TRUE)
out <- separate(out, SiteUnits, c("G1","G2"), " \\| ", remove = TRUE)
save(out, file = "AllForestSiteSeriesScores.RData")
write.csv(out, "Abies_ForestSiteSeriesScores.csv")
#####convert to matrix
outdup <- out[,c(1,2,5)]
colnames(outdup) [1:2] <- c("G2", "G1")
out2 <- rbind(out[,c(1,2,5)], outdup)
out2 <- melt(out2, id.vars = c("G1","G2"))
PairScoreMatrix <- dcast(out2, G1 ~ G2 + variable)###Convert to site unit by species matrix
PairScoreMatrix [is.na(PairScoreMatrix)] <- 0
rownames(PairScoreMatrix) <- PairScoreMatrix [,1]
PairScoreMatrix <- PairScoreMatrix [-1]
write.csv(PairScoreMatrix, "Abies_ActualPotential_SimilarityMatric.csv")
save(PairScoreMatrix, file = "PairScoreMatrixScores.RData")
write.csv(badAssoc, "NewMethodAllForestSiteSeriesDifferentialSum_Low.csv", row.names = FALSE)
write.csv(out, "NewMethodAllForestSiteSeriesDifferentialSum_ALL.csv", row.names = FALSE)
          ###identify those subassociations that are similar but in different associations
#PairScoreMatrix <- as.matrix(PairScoreMatrix, labels=TRUE)
DissPair <- dist(PairScoreMatrix)
DissPairhc <- hclust(DissPair)

DissPairdend <-as.dendrogram(DissPairhc)
plot(DissPairdend,  hang = -1, cex = .3)
heights_per_k.dendrogram(DissPairdend) #### examine heights at various cut levels
op = par(mfrow = c(2, 1))
plot(cut(DissPairdend, h = 2)$upper, main = "Upper tree of cut at h=75")
plot(cut(DissPairdend, h = 3)$lower[[2]], main = "Second branch of lower tree with cut at h=75") ###this will show the site series of a lower branch

#DissDendo <- as.dendrogram(DissPairhc)
#plot(DissDendo, type = "rectangle", xlab = "Height", horiz = TRUE, cex = 0.3)
####plots from ape package
plot(as.phylo(DissPairhc), type = "phylo", cex = .35 ,no.margin = FALSE, font = 3, lab4ut = "axial", 
     align.tip.label=TRUE, label.offset = .1, edge.width =.2)#, rotate.tree = 80
title ("Abies Site Series")

###Output cluster of cluster cuts
clustgroups <- cutree(DissPairhc, h=3)
######options to work in
# plot dendrogram with some cuts
op = par(mfrow = c(2, 1))
plot(cut(hcd, h = 75)$upper, main = "Upper tree of cut at h=75")
plot(cut(hcd, h = 75)$lower[[2]], main = "Second branch of lower tree with cut at h=75")

hc <- hclust(dist(USArrests[1:4,]), "ave")
dend <- as.dendrogram(hc)
heights_per_k.dendrogram(dend)
    
dend1 <- color_branches(dend, k = 3)
    dend2 <- color_labels(dend, k = 3)

    par(mfrow = c(1,2))
    plot(dend1, main = "Colored branches")
    plot(dend2, main = "Colored labels")
   
     # horiz mirror version
    par(mar = c(3,7,1,1))
    plot_horiz.dendrogram(d, side = TRUE)

write.csv(clustgroups, "ClusterAnalysisGroups_AbiesSiteSeries.csv")

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


##############removed from loops - required for significance class calculations

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
                                  ifelse(Cov1[,4] <= 5,4,
                                         ifelse(Cov1[,4] <= 10, 5,
                                                ifelse(Cov1[,4] <= 20,6,
                                                       ifelse(Cov1[,4] <= 33,7,
                                                              ifelse(Cov1[,4] <= 50,8,
                                                                     ifelse(Cov1[,4] <= 75,9,10)))))))))
Cov1$SigDiff <- Cov1$SigA - Cov1$SigB

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


###sum values
Cov1$Value <- apply(Cov1[,c("Differential","Const","DDpoints")],1,FUN = function(x){sum(scoreVals$Value[scoreVals$Code %in% x])})
Cov1$Value <- ifelse(Cov1$Value>0 & is.na(Cov1$Differential), Cov1$Value -2, Cov1$Value)
Cov1$Value <- ifelse(Cov1$Value<0, 0, Cov1$Value)
Cov1$Value <- Cov1$Value*(Cov1[,3]/100)
Cov1$Value <- Cov1$Value + Cov1$DDpoints
Cov1 <- merge(Cov1, lifeform, by = "Species", all.x = TRUE)
Cov1$Value <- ifelse(Cov1$Type %in% c(9,10,11,13), Cov1$Value/2, Cov1$Value)
Cov1$Value <- ifelse(is.na(Cov1$Differential ) & is.na(Cov1$DDclass), 0, Cov1$Value )
sumA <- sum(Cov1$Value)

# Unit k vs j 
Cov2 <- CovMatrix[CovMatrix[,5] >= 60 & !is.na(CovMatrix[,5]),]###remove < 60% constancy for k unit
Cov2$ConstDiff <- abs(Cov2[,5] - Cov2[,3])
Cov2$Differential <- ifelse(Cov2$ConstDiff >= 80,"d1",
                            ifelse(Cov2$ConstDiff >= 60,"d2",
                                   ifelse(Cov2$ConstDiff >= 40,"d3",NA)))

Cov2$CovDiff <- ifelse (Cov2[,4] >10, (Cov2[,4] / Cov2[,2])/10, 
                        ifelse(Cov2[,4] >5 , ((Cov2[,4] / (Cov2[,2]*1.5) / 10)), 0))
Cov2$CovDiff <- ifelse(Cov2$CovDiff >1, 1, 
                       ifelse(Cov2$CovDiff <0.2, 0, Cov2$CovDiff ))
Cov2$DDpoints <- ifelse(Cov2[,4] <10, (Cov2$CovDiff * 4) *(Cov2[,4]/10),Cov2$CovDiff * 4)
Cov2$DDclass <- ifelse(Cov2$DDpoints >=4, "dd1",
                       ifelse(Cov2$DDpoints >=2, "dd2",
                              ifelse(Cov2$DDpoints >=1.2, "dd3",
                                     ifelse(Cov2$DDpoints > 0, "dd4", NA))))

Cov2$Const <- ifelse(Cov2[,4] >= 10,"cd",
                     ifelse(Cov2[,4] >= 0.3, "c","cm"))
###sum values
Cov2$Value <- apply(Cov2[,c("Differential","Const","DDpoints")],1,FUN = function(x){sum(scoreVals$Value[scoreVals$Code %in% x])})
Cov2$Value <- ifelse(Cov2$Value>0 & is.na(Cov2$Differential), Cov2$Value -2, Cov2$Value)
Cov2$Value <- ifelse(Cov2$Value<0, 0, Cov2$Value)
Cov2$Value <- Cov2$Value*(Cov2[,5]/100)
Cov2$Value <- Cov2$Value + Cov2$DDpoints
Cov2 <- merge(Cov2, lifeform, by = "Species", all.x = TRUE)
Cov2$Value <- ifelse(Cov2$Type %in% c(9,10,11,13), Cov2$Value/2, Cov2$Value)
Cov2$Value <- ifelse(is.na(Cov2$Differential ) & is.na(Cov2$DDclass), 0, Cov2$Value )
sumB <- sum(Cov2$Value)


###########some clustering options
DissPairOrder <- order.optimal (DissPair,DissPairhc$merge )
md <- as.dist(crossprod(as.matrix(DissPair, diag = 0)))   # Murtagh's distances
hg <- order.greedy(md)
go <- order.optimal(md, hg$merge)
d <- DissPair
hc <- DissPairhc
co <- DissPairOrder
### compare images
op <- par(mfrow=c(2,2), pty="s")
implot(d[[hc$order]], main="hclust")
implot(d[[co$order]], main="hlcust + optimal")
implot(d[[hg$order]], main="greedy")
implot(d[[go$order]], main="greedy + optimal")
par(op)
# compare lengths
order.length(d, hc$order)
order.length(d, co$order)
order.length(d, hg$order)
order.length(d, go$order)






#########################old code

# Unit k vs j 
Cov2 <- CovMatrix[CovMatrix[,5] >= 60 & !is.na(CovMatrix[,5]),]###remove < 60% constancy for k unit
Cov2$ConstDiff <- abs(Cov2[,5] - Cov2[,3])
Cov2$Differential <- ifelse(Cov2$ConstDiff >= 80,"d1",
                            ifelse(Cov2$ConstDiff >= 60,"d2",
                                   ifelse(Cov2$ConstDiff >= 40,"d3",NA)))

Cov2$CovDiff <- ifelse (Cov2[,4] >10, (Cov2[,4] / Cov2[,2])/10, 
                        ifelse(Cov2[,4] >5 , ((Cov2[,4] / (Cov2[,2]*1.5) / 10)), 0))
Cov2$CovDiff <- ifelse(Cov2$CovDiff >1, 1, 
                       ifelse(Cov2$CovDiff <0.2, 0, Cov2$CovDiff ))
Cov2$DDpoints <- ifelse(Cov2[,4] <10, (Cov2$CovDiff * 4) *(Cov2[,4]/10),Cov2$CovDiff * 4)
Cov2$DDclass <- ifelse(Cov2$DDpoints >=4, "dd1",
                       ifelse(Cov2$DDpoints >=2, "dd2",
                              ifelse(Cov2$DDpoints >=1.2, "dd3",
                                     ifelse(Cov2$DDpoints > 0, "dd4", NA))))

Cov2$Const <- ifelse(Cov2[,4] >= 10,"cd",
                     ifelse(Cov2[,4] >= 0.3, "c","cm"))
###sum values
Cov2$Value <- apply(Cov2[,c("Differential","Const","DDpoints")],1,FUN = function(x){sum(scoreVals$Value[scoreVals$Code %in% x])})
Cov2$Value <- ifelse(Cov2$Value>0 & is.na(Cov2$Differential), Cov2$Value -2, Cov2$Value)
Cov2$Value <- ifelse(Cov2$Value<0, 0, Cov2$Value)
Cov2$Value <- Cov2$Value*(Cov2[,5]/100)
Cov2$Value <- Cov2$Value + Cov2$DDpoints
Cov2 <- merge(Cov2, lifeform, by = "Species", all.x = TRUE)
Cov2$Value <- ifelse(Cov2$Type %in% c(9,10,11,13), Cov2$Value/2, Cov2$Value)
Cov2$Value <- ifelse(is.na(Cov2$Differential ) & is.na(Cov2$DDclass), 0, Cov2$Value )
sumB <- sum(Cov2$Value)