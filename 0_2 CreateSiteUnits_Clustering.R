####Script for producing veg units from density based clusering
### from http://www.sthda.com/english/wiki/wiki.php?id_contents=7940summary 
# Welcome! Related Books: `Practical Guide To Cluster Analysis in R` at https://goo.gl/13EFCZ
##MacKenzie, march 2019

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
require (dbscan)
require(factoextra)
require(fpc)
rm(list=ls())
wd=tk_choose.dir(); setwd(wd)


###Import vegetation plot data produced from 0_BEC_data_import_clean.R script
load("VegDat_Lumped.RData") ###variable = vegData

##Import SU table of plots desired for cluster analysis
SUhier <- read.csv("FebHierarchyLevelSU.csv", stringsAsFactors = FALSE)
#### Reduce vegdata to those listed in SUHier

########Function to define optimal eps
dbscan::kNNdistplot(iris, k =  4)
abline(h = 0.4, lty = 2)
###compute DBSAB using two different packages
set.seed(123)
# fpc package
res.fpc <- fpc::dbscan(iris, eps = 0.4, MinPts = 4)
# dbscan package
res.db <- dbscan::dbscan(iris, 0.4, 4)
#Make sure that both version produce the same results:
  
  all(res.fpc$cluster == res.db)

 # The result can be visualized as follow:
    
    fviz_cluster(res.fpc, iris, geom = "point")
    
    
  dbscan(data, eps, MinPts = 5, scale = FALSE, 
       method = c("hybrid", "raw", "dist"))
# Compute DBSCAN using fpc package
set.seed(123)
db <- fpc::dbscan(df, eps = 0.15, MinPts = 5)
# Plot DBSCAN results
plot(db, df, main = "DBSCAN", frame = FALSE)
# Print DBSCAN
print(db)
## or plot with factoextra
fviz_cluster(db, df, stand = FALSE, frame = FALSE, geom = "point")
km.res <- kmeans(df, 5, nstart = 25)
fviz_cluster(km.res, df, frame = FALSE, geom = "point")




###########################################################################
#####Noise Clustering with predefined classes#######################
###Create training and testing data

vegDat.chord <- decostand(SUsumMatrix[-1], "normalize") ##standardise data
vegMat <- vegDat.chord
vegMat <- cbind(SUsumMatrix[,1],vegDat.chord)

####create training and testing data sets
vegNew <- vegMat[rownames(vegMat) %in% sample(rownames(vegMat), 100, replace = FALSE),] 
vegOld <- vegMat[!(rownames(vegMat) %in% rownames(vegNew)),]
#actualClass <- vegNew[,1:2] ###siteseries to grouping lookup
#rownames(vegNew)  <- vegNew[1] ###set rownames to siteseries
k <- ncol(vegOld)
n <- nrow(vegOld)
###grouping of training set
grouping <- vegOld[1]
vegOld <- vegOld[,-(1)]
vegNew <- vegNew[,-(1)]
k <- ncol(vegOld)
n <- nrow(vegOld)

grouping <- as.vector.factor(grouping)
vegOld.clust <- as.vegclust(vegOld, grouping)###create noise clustering with grouping as classes
###Kmeans Clustering
vegOld.kmclst <- vegclust(x = vegOld[,-1], mobileCenters=5, method = "KM", nstart=20)###create noise clustering with grouping as classes
t(vegOld.kmclst$memb)
###NC Clustering
vegOld.kmclst <- vegclust(x = vegOld[,-1], mobileCenters=5, method = "NC", m=1.2, dnoise=0.8, nstart=20)###create noise clustering with grouping as classes
round(t(vegOld.kmclst$memb), dig=2)
groups = defuzzify(vegOld.kmclst, method="cut", alpha=0.8)$cluster
table(groups)


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