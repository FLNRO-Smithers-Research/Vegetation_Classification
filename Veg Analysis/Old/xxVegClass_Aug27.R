####Script for various veg classification analyses and test. First part cleans data, then various methods
###for rule based classification,testing noise clustering, and creating edatopic grids
###Kiri Daust, July 2018

.libPaths("E:/R packages")
library(reshape)
library(reshape2)
library(vegan)
library(caret)
library(tcltk)
library(ddply)
library(randomForest)
library(Matrix)
require(labdsv)
require(vegan)
require(gdata)
require(MASS)
require(xlsx)
require (C50)
require(tidyr)
require(stringr)
require(rpart)
library(tree)
library(rattle)
library(rpart.plot)
library(partykit)
library(vegclust)
library(standardize)
library(dplyr)
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
vegData <- read.table("BECMasterALLVegData.txt", header = TRUE) #3 column

vegData <- separate(vegData, Species, c("Species","Type"), "-", remove = TRUE)
vegData <- vegData[vegData$Type %in% c(1,2),]
vegData <- vegData[!is.na(vegData$Species),]

treeSpp <- as.character(unique(vegData$Species))
#####################################################

vegData$Species <- unlist(lapply(vegData$Species, toupper))
vegData <- vegData[,-3]
save(vegData, file = "VegDat_Raw.RData")
load("VegDat_Raw.RData")

###update old codes
masterList <- read.csv("MasterSpeciesCodeList.csv", stringsAsFactors = FALSE)
noMatch <- masterList[masterList$OldCode != masterList$Code,4:5]
temp <- merge(vegData,noMatch,by.x = "Species", by.y = "OldCode")
temp$Species <- temp$Code
temp <- temp[,-4]
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
write.csv(vegRemove, file = "Plots_0Cover.csv")
vegData <- vegData[vegData$Cover > 0,]
save(vegData, file = "VegDat_Clean.RData")
load("VegDat_Clean.RData")

###lump species
lump <- read.csv("CombineCodeCorr5Jan2018_Lump.csv", stringsAsFactors = FALSE)
lump$Lump <- unlist(lapply(lump$Lump, tolower))

vegData <- merge(vegData, lump, by.x = "Species", by.y = "Code", all.x = TRUE) ##lump data
vegData$Species <- ifelse(!is.na(vegData$Lump), vegData$Lump, vegData$Species)

save(vegData, file = "VegDat_Lumped.RData")
load("VegDat_Lumped.RData")

##Import hierachy table
SUhier <- read.xlsx("HierachySUTable.xlsx", sheet = 1)
level <- select.list(choices = colnames(SUhier), graphics = TRUE)
SUhier <- SUhier[,c("PlotNumber", level)]
vegData <- merge(vegData, SUhier, by = "PlotNumber", all.x = TRUE)###Must select SiteSeries as level for now
colnames(vegData)[5] <- "Group"
if(any(is.na(vegData$Group))){
  warning("Data contains Plots not in hierachy table. These will be removed.")
}
vegData <- vegData[!is.na(vegData$Group),]
vegData$PlotNumber <- as.character(vegData$PlotNumber)
vegData <- vegData[vegData$Group != "",]
constCut <- 0 ##remove species less than cutoff

##roll up into groups
require(doParallel)
set.seed(123321)
coreNo <- makeCluster(detectCores() - 1)
registerDoParallel(coreNo, cores = detectCores() - 1)
Cores <- as.numeric(detectCores()-1)
clusterEvalQ(coreNo, .libPaths("E:/R packages"))

temp <- foreach(SS = unique(vegData$Group), .combine = rbind, .packages = "foreach") %dopar% {
  sub <- vegData[vegData$Group == SS,]
  num <- length(unique(sub$PlotNumber))
  foreach(Spp = unique(sub$Species), .combine = rbind) %do% {
    sub2 <- sub[sub$Species == Spp,]
    numSpp <- dim(unique(sub2[,1:2]))[1]
    mean <- mean(sub2$Cover)
    const <- numSpp/num
    if(const >= constCut){
      out <- data.frame(Group = SS, Species = Spp, MeanCov = mean, Constancy = const*100)
    }
    
  }
}

vegData <- temp

##Remove species with a maximum constancy < 60%
maxConst <- aggregate(Constancy ~ Species, vegData, FUN = max)
SpKeep <- maxConst[maxConst$Constancy >= 60,]
SpKeep <- as.character(SpKeep$Species)
SpKeep <- SpKeep[SpKeep != "luzu+p"] ##this code causes issues further down - possibly the + sign
vegData$Species <- as.character(vegData$Species)
vegData <- vegData[vegData$Species %in% SpKeep,]

save(vegData, file = "VegSU_DifferentialSpecies.RData")
###################################END of CLeaning Data#################################################################

####################Can usually start here####################################

####This part creates a list of all tree species for subsetting##################
temp <- read.table("BECMasterALLVegData.txt", header = TRUE)
temp <- separate(temp, Species, c("Species","Type"), "-", remove = TRUE)
temp <- temp[temp$Type %in% c(1,2),]
temp <- temp[!is.na(temp$Species),]
treeSpp <- as.character(unique(temp$Species))
##############################################################################

load("VegSU_DifferentialSpecies.RData")

######Subset for just tree species (optional)#################
vegData <- vegData[vegData$Species %in% treeSpp,]

#####SELECT LEVEL TO GROUP BY AND CONVERT TO MATRIX###############
library(openxlsx)
vegData <- melt(vegData, id.vars = c("Group","Species"))
vegData <- dcast(vegData, Group ~ Species + variable)###Convert to site by species matrix
##vegData <- dcast(vegData, PlotNumber ~ Species, value.var = "Cover", fun.aggregate = mean)
SUhier <- read.xlsx("HierachySUTable.xlsx", sheet = 1)
level <- select.list(choices = colnames(SUhier), graphics = TRUE)###Select which level to test classification
SUhier <- SUhier[,c("SiteSeries", level)]

SUhier <- unique(SUhier)
colnames(vegData)[1] <- "SiteSeries"
vegData <- merge(vegData, SUhier, by = "SiteSeries", all.x = TRUE)###Merge classification with matrix
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

###FOR NOISE CLUSTERING USE vegData AND GO TO SCRIPT BELOW RULE BASED CLASSIFICATION

###RULE BASED CLASSIFICATION#####################
####RPart##############
RPcontrol <- rpart.control(maxcompete = 1, maxsurrogate = 2, usesurrogate = 1, surrogatestyle = 1)
fit.rp <- rpart(Class ~ ., data = vegData[,-2], method = "class")
write(capture.output(summary(fit.rp)), file = "RPartClassLevel.txt")
rpart.plot(fit.rp, type = 0, cex = 0.7, extra = 0)

##Test quality
vegData$Pred <- predict(fit.rp, newdata = vegData[,-c(1,2)], type = "class")
compare <- vegData[,c(1,length(vegData),2)]
compare$Same <- ifelse(compare$Class == compare$Pred,1,0)
rp.qual <- sum(compare$Same)/length(compare$Same)

#####Random Forest######################
vegData <- vegData[,-length(vegData)]
vegData[is.na(vegData)] <- 0
vegData$Class <- as.factor(as.character(vegData$Class))

tic()
rf.model <- randomForest(Class ~ .,data=vegData[,-2], nodesize = 2, 
                         do.trace = 10, ntree=501, na.action=na.fail, importance=TRUE, proximity=TRUE)
toc()
MDSplot(rf.model, vegData$Class)
summary(rf.model)

###test quality
vegData$Pred <- predict(rf.model, newdata = vegData[,-c(1,2)])
compareRF <- vegData[,c(1,length(vegData),2)]
compareRF$Same <- ifelse(compareRF$Class == compareRF$Pred,1,0)
rf.qual <- sum(compareRF$Same)/length(compareRF$Same)

varImpPlot(rf.model, n.var = 30)
varImpMat <- importance(rf.model)

#####C50############################
#vegData[is.na(vegData)] <- 0
#vegData <- vegData[,-length(vegData)]

c50.fit <- C5.0(Class ~ ., data=vegData[,-2], rules = FALSE)
plot(c50.fit)
summary(c50.fit)

vegData$Pred <- predict(c50.fit, newdata = vegData[,-c(1,2)])
compareC5 <- vegData[,c(1,length(vegData),2)]
compareC5$Same <- ifelse(compareC5$Class == compareC5$Pred,1,0)
temp <- compareC5[compareC5$Same == 0,]
write.csv(temp, "C50WrongSS_Order.csv")
c5.qual <- sum(compareC5$Same)/length(compareC5$Same)
write(capture.output(summary(c50.fit)), "c50OrderMod.txt")

#####End of C50################################

###########################################################################
#####Noise Clustering with predefined classes#######################
###Create training and testing data

vegDat.chord <- decostand(vegData[,-(1:2)], "normalize") ##standardise data
vegMat <- cbind(vegData[,1:2],vegDat.chord)

####create training and testing data sets
vegNew <- vegMat[rownames(vegMat) %in% sample(rownames(vegMat), 100, replace = FALSE),] 
vegOld <- vegMat[!(rownames(vegMat) %in% rownames(vegNew)),]
actualClass <- vegNew[,1:2] ###siteseries to grouping lookup
rownames(vegNew)  <- vegNew$SiteSeries ###set rownames to siteseries

###grouping of training set
grouping <- vegOld$Class
vegOld <- vegOld[,-(1:2)]
vegNew <- vegNew[,-(1:2)]

grouping <- as.vector.factor(grouping)
vegOld.clust <- as.vegclust(x = vegOld, y = grouping, method = "HNC", dnoise = 70)###create noise clustering with grouping as classes
vegComb <- vegclass(vegOld.clust, vegNew) ##classify vegNew
vegComb.memb <- vegComb[["memb"]]
newGroup <- dematrify(vegComb.memb) ##extract classification
newGroup <- newGroup[,1:2]
colnames(newGroup) <- c("SiteSeries","Class")
newGroup <- cbind(newGroup, actualClass$Class)####merge actual classification for comparison
colnames(newGroup)[3] <- "Actual"

###MDS for visualisation
MDS <- metaMDS(vegNew, distance = "bray", k = 2, trymax = 200)
MDS.df <- as.data.frame(scores(MDS, display = "sites"))###extract mds scores
MDS.df$SiteSeries <- rownames(MDS.df)
MDS.df <- merge(MDS.df,newGroup, by = "SiteSeries") ##merge predicted and actual classification
MDS.df <- MDS.df[,-1]

colnames(MDS.df)[3:4] <- c("Actual","Predict")

ggplot(MDS.df)+
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Predict), size = 2.5, shape = 17)+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none")


#####NOISE CLUSTERING################################
vegDat.chord <- decostand(vegMat, "normalize") ##standardise data
vegData <- cbind(covMat,constMat)
veg.nc <- vegclust(vegOld, mobileCenters = 6, m = 1.5, dnoise = 1, method = "NC", nstart = 20)
round(t(veg.nc$memb), digits = 2)
groups <- as.data.frame(defuzzify(veg.nc)$cluster)
groups
memb <- veg.nc$memb
memb <- round(t(memb), digits = 2)
memb$BGC <- rownames(memb)
memb <- merge(classTable, memb, by = "BGC")
groups$BGC <- rownames(groups)
groups <- merge(classTable, groups, by = "BGC")
groups$Class <- as.factor(groups$Class)
table(groups[,2:3])
x <- sample(seq(1,1123,1), 150, replace = FALSE)
newDat <- vegDat.chord[x,]
oldDat <- vegDat.chord[-x,]


######################################################################
#####CREATE EDATOPIC GRIDS WITH VEG STATS IN EACH CELL################
###############################################################

vegAll <- read.table("BECMasterALLVegData.txt", header = TRUE)
codeCross <- read.csv("CodeCrosswalk.csv", stringsAsFactors = FALSE)

vegAll <- separate(vegAll, Species, c("Species","Type"), "-", remove = TRUE)
vegAll <- vegAll[vegAll$Type %in% c(1,2),]

BGCLookup <- read.csv("BGCLookup.csv", stringsAsFactors = FALSE)
BGCLookup <- BGCLookup[,c(3,12)]
BGCLookup <- BGCLookup[BGCLookup$BGC_LABEL != "",]
colnames(BGCLookup)[1] <- "PlotNumber"

##import edatopic data
plotEnv <- read.csv("KiriEnvDat.csv", stringsAsFactors = FALSE)
plotEnv <- plotEnv[plotEnv$NutrientRegime %in% c("A","B","C","D","E"),]
plotEnv <- plotEnv[plotEnv$MoistureRegime %in% c(0,1,2,3,4,5,6,7,8),]
plotEnv <- plotEnv[,-2]
plotEnv <- merge(plotEnv,BGCLookup, by = "PlotNumber", all.x = TRUE)
plotEnv$BGC_LABEL <- gsub("[[:space:]]","",plotEnv$BGC_LABEL)
plotEnv <- plotEnv[plotEnv$BGC_LABEL != "",]
colnames(plotEnv)[4] <- "Unit"

modBGC <- read.csv("ModelledBGC.csv", stringsAsFactors = FALSE)

for(i in 1:length(modBGC$BGC)){
  Unit <- modBGC$BGC[i]
  envSub <- plotEnv[plotEnv$Unit == Unit,]
  vegData <- merge(vegAll, envSub, by = "PlotNumber", all.y = TRUE)
  
  if(length(vegData$PlotNumber) > 2){
  vegData$Group <- paste(vegData$MoistureRegime,"-",vegData$NutrientRegime, sep = "")
  vegData <- vegData[!is.na(vegData$Species),]
  vegData <- vegData[,-c(3,5:7)]
  constCut <- 0.2
  ###roll up
    temp <- foreach(SS = unique(vegData$Group), .combine = rbind, .packages = "foreach") %dopar% {
      sub <- vegData[vegData$Group == SS,]
      num <- length(unique(sub$PlotNumber))
      foreach(Spp = unique(sub$Species), .combine = rbind) %do% {
        sub2 <- sub[sub$Species == Spp,]
        numSpp <- dim(unique(sub2[,1:2]))[1]
        mean <- mean(sub2$Cover)
        const <- numSpp/num
        if(const >= constCut){
          out <- data.frame(Group = SS, Species = Spp, MeanCov = mean, Constancy = const*100, Num = num)
        }
        
      }
    }
    
    vegGrid <- temp
    ###classify as 1,2 or 3
    vegGrid$Pres <- ifelse(vegGrid$MeanCov > 20 & vegGrid$Constancy > 50, 1,
                           ifelse(vegGrid$MeanCov > 10 & vegGrid$Constancy > 25,2,3))
    vegGrid$Order <- vegGrid$MeanCov*vegGrid$Constancy
    vegGrid <- merge(vegGrid,codeCross, by.x = "Species", by.y = "Code", all.x = TRUE)
    vegGrid <- vegGrid[,c(2,5:8)]
    colnames(vegGrid)[5] <- "Species"
    
    Lab <- ave(vegGrid[,c("Group","Pres","Species","Order")], vegGrid$Group, FUN = combineSpp)
    vegGrid <- separate(vegGrid, Group, c("Numeric","Alph"), "-", remove = TRUE)
    vegGrid$Lab <- Lab$Species
    vegGrid <- unique(vegGrid[,c(1:3,7)])
    
    ##plot
    pdf(file = paste("EdaGrid_",Unit,".pdf",sep = ""), height = 10.5, paper = "letter")
    print(ggplot(data = vegGrid)+
      geom_tile(aes(x= Alph, y = Numeric), color = "black", fill = "white")+
      geom_text(aes(x = Alph, y = Numeric, label = Lab), size = 3)+
      geom_text(aes(x = Alph, y = Numeric, label = Num,hjust = -4, vjust = -6), size = 2, color = "red")+
      scale_y_discrete(limits = c("8","7","6","5","4","3","2","1","0"))+
      scale_x_discrete(limits = c("A","B","C","D","E"))+
      labs(x = "Relative Soil Nutrient Regime", y = "Relative Soil Moisture Regime", cex = 1.5, title = paste(Unit," (",sum(vegGrid$Num)," plots)", sep = ""))+
      theme_bw(base_size = 10)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      coord_fixed())
    dev.off()
    
  }
}

###Function used above to combine data into formatted string
combineSpp <- function(x){
  x <- x[order(-x$Order),]
  if(any(x$Pres == 1)){
    dom <- paste(x$Species[x$Pres == 1], collapse = ",")
    dom <- paste("*",dom,"*", sep = "")
  }else{
    dom <- ""
  }
  
  if(any(x$Pres == 2)){
    sub <- x$Species[x$Pres == 2]
    if(length(sub) > 5){
      sec <- paste("(", paste(sub[1:4], collapse = ","),"\n",paste(sub[5:length(sub)], collapse = ","),")",sep = "")
    }else{
    sec <- paste(sub, collapse = ",")
    }
  }else{
    sec <- ""
  }
  
  if(any(x$Pres == 3)){
    sub <- x$Species[x$Pres == 3]
    if(length(sub) > 5){
      un <- paste("(", paste(sub[1:4], collapse = ","),"\n",paste(sub[5:length(sub)], collapse = ","),")",sep = "")
    }else{
      un <- paste(sub, collapse = ",")
      un <- paste("(",un,")", sep = "")
    }
  }else{
    un <- ""
  }
  
  return(paste(dom,"\n", sec, "\n", un, sep = ""))
}


########################################
