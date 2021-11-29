setwd("~/R/dataAnalysis_UFRN/BooleanNetworkBioME/TBN_IBN/EGEOD18494_KM2/ZhouGenes/AvgBin_Heatmap")

library(ellipsis)
library(gplots)
library(RColorBrewer)
library(BoolNet)
library(reshape2)
library(ggplot2)
library(igraph)
library(dplyr)

########################################################################
# determine the effect of stochasticity in KM3 binarization on node states
# by taking 1000 KM3 binarizations of the same input sample, averaging node
# states and constructing a heatmap of the data. This file is produced 
# by the python script "KM_binarization.py" in the subfolder "BinarizationScripts"
########################################################################

# ~/R/dataAnalysis_UFRN/BooleanNetworkBioME/TBN_IBN/EGEOD18494_KM2/ZhouGenes/binarized_colab/proxy/KM2/group8/proxy_grupo8_ 
# C:\Users\AlexandreSarmento\Documents\R\dataAnalysis_UFRN\BooleanNetworkBioME\TBN_IBN\EGEOD18494_KM2\ZhouGenes\binarized_colab\breast
filename = "~/R/dataAnalysis_UFRN/BooleanNetworkBioME/TBN_IBN/EGEOD18494_KM2/ZhouGenes/binarized_colab/breast/breast_raw_data"
#allTimePoints = list(-4,4,8,12)
allTimePoints = list(0,1,2,3,4,5,6,7,8,9,10,11)
mergedSampleList = list()
for(i in 1:500){
  #print(i)
  sampleName = paste("~/R/dataAnalysis_UFRN/BooleanNetworkBioME/TBN_IBN/EGEOD18494_KM2/ZhouGenes/binarized_colab/breast/breast_raw_data", toString(i), ".csv",sep="")
  mergedSample <- read.csv(sampleName, comment.char="#", header=FALSE, check.names=FALSE, row.names=1)
  mergedSampleList[[i]] = mergedSample
}
s1 = mergedSampleList[[1]]
noRows = length(s1[,1])
noColumns = length(s1[1,])
rowList = c(1:noRows)
colList = c(1:noColumns)
popAverage = data.frame(row.names = rownames(data.frame(s1)))
popCumulative = data.frame(row.names = rownames(data.frame(s1)))
rownames(popAverage) = rownames(data.frame(s1))
for(rowVal in rowList){
  for(colVal in colList){
    cumulativeVal = 0
    for(sample in mergedSampleList){
      tmpVal = sample[rowVal,colVal]
      cumulativeVal = tmpVal + cumulativeVal 
    }
    popCumulative[rowVal,colVal]= cumulativeVal
    popAverage[rowVal,colVal]= round(cumulativeVal/length(mergedSampleList),1)
  }
}


########################################################################
# round the sample averages to their closest integer value 
# value < .5 = 0 value > 0.5 = 1
########################################################################

avgSampleList = list('S1' = popAverage)
binAvgSampleList = avgSampleList
name = 'S1'
for(name in names(avgSampleList)){
  
  tmpSample = avgSampleList[[name]]
  noRows = length(tmpSample[,1])
  noColumns = length(tmpSample[1,])
  rowList = c(1:noRows)
  colList = c(1:noColumns)
  for(rowVal in rowList){
    for(colVal in colList){
      tmpVal = tmpSample[rowVal,colVal]
      if(tmpVal > .5){
        tmpSample[rowVal,colVal] = 1
      }
      if(tmpVal < .5){
        tmpSample[rowVal,colVal] = 0
      }
      # the binarized version of avg sample list
      binAvgSampleList[[name]] = tmpSample
    }
  }
}

# "HIF1A","CDKN1A","HIF1AN","ATR","TP53","MDM2","AKT1","BBC3","MIR17HG","PTEN",
# "BCL2L11","FOXO3","CASP3","CYCS","EP300"

# PUMA CASP3 FOXO HIF1A AKT PTEN miR1792 BIM ATR CYTOC P21 MDM2 PHD P53

myorder <-c("HIF1A","P21","ATR","TP53","MDM2","AKT","PUMA","miR1792",
            "PTEN","BIM","FOXO","CASP3","CYTOC","EP300")


binAvgSampleList$S1 <- binAvgSampleList$S1[myorder,]

rownames(binAvgSampleList$S1) <- c("HIF1A","P21","ATR","TP53","MDM2","AKT","PUMA","miR1792",
                                   "PTEN","BIM","FOXO","CASP3","CYTOC","EP300")


for(name in names(binAvgSampleList)){
  fileName =paste("~/R/dataAnalysis_UFRN/BooleanNetworkBioME/TBN_IBN/EGEOD18494_KM2/ZhouGenes/binarized_colab/breast/binAvg",name,".csv",sep="")
  write.table(binAvgSampleList[[name]], file =fileName, sep= ','
              ,quote=FALSE, row.names = rownames(binAvgSampleList[[name]])) 
}

# creates a own color palette from red to green for the heatmaps
my_palette <- colorRampPalette(c("blue", "gray", "yellow"))(n = 299)
col_breaks = c(seq(0,.1,length=100),  # for green
               seq(.11,0.9,length=100),              # for yellow
               seq(0.91,1,length=100))



plotHeatMap <- function(sample,titleDescripter,timePoints){  
  #i, paste("test#",counter,sep=""),equiDistantTimePoints
  #sample = avgSample1
  #titleDescripter = "Healthy-1"
  #timePoints = equiDistantTimePoints2
  rnames <- rownames(sample)          # assign labels in column 1 to "rnames"                           
  mat_data <- data.matrix(sample[,1:ncol(sample)])  # transform column 2- # of timepts into a matrix
  rownames(mat_data) <- rnames                  # assign row names 
  colnames(mat_data) <- timePoints
  # creates a 5 x 5 inch image
  outputName =paste(titleDescripter,".png",sep="") 
  png(outputName,   # create PNG for the heat map        
      width = 4*500,        # 5 x 300 pixels
      height = 3*500,
      res = 500,            # 300 pixels per inch
      pointsize = 6)      # smaller font size
  par(cex.main=0.8,cex.lab = 1)
  heatmap.2(mat_data,
            #cellnote = mat_data,  # same data set for cell labels
            #main = titleDescripter, # heat map title
            offsetCol = 2,
            offsetRow = 0,
            srtCol=270,
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(9,12),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier 
            breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="none",    # don't draw dendrograms
            Colv="NA",            # turn off column clustering
            Rowv="NA",            # turn off row clustering
            #srtCol=0,
            adjCol= c(.5,.5),
            adjRow= c(0,.5),
            sepwidth=c(0.001,0.001),
            sepcolor="black",
            colsep=1:ncol(mat_data),
            rowsep=1:nrow(mat_data),
            cexRow = 1.0,
            cexCol=.8,
            # color key + density info
            #key = TRUE,
            keysize = 1,
            #density.info=c("histogram","density","none"),
            key.title = NA,
            key.ylab = NA,
            key.xlab = "Binarized Abundance",
            key.ytickfun = NA
  )
  dev.off()               # close the PNG device  
}

plotHeatMap(binAvgSampleList$S1,"breast_KM2_Zhou",allTimePoints)

########################################################################################
# infer most probably Boolean rules using BoolNet package
# see text for how multiple rules of same error were integrated into consensus rules
########################################################################################
geneNames = rownames(binAvgSampleList$S1)
ibn <- reconstructNetwork(binAvgSampleList$S1, 
                          method="reveal",
                          maxK = length(geneNames),
                          returnPBN=FALSE,
                          readableFunctions=TRUE)

print(ibn)