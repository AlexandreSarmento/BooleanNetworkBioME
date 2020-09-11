#install.packages("BoolNet")
require(BoolNet)
#devtools::install_github("clsdavid/FBNNet2_public")
require(FBNNet)
setwd("~/d4/modelling/FBNNet")
#generate BoolNet type of cell cycle network
cellcyclenetwork<-loadNetwork("hif1a.bn")
trainingseries<-generateTimeSeries(cellcyclenetwork,1000,42)
#reduce the duplicate samples
trainingseries<-FBNDataReduction(trainingseries)
#generate initial states, i.e., get the first state from all samples
startstates<-lapply(trainingseries,function(x)x[,1])
#construct FBM cube, maxK=3, temporal=1
cube<-constructFBNCube(cellcyclenetwork$genes,cellcyclenetwork$genes,trainingseries,3,1,TRUE)
#mine FBN type of network
FBNcellcyclenetwork<-mineFBNNetwork(cube,cellcyclenetwork$genes)
#plot the static FBN type of Network
FBNNetwork.Graph(FBNcellcyclenetwork)
#print out the Network
print(FBNcellcyclenetwork)
#reconstruct the timeseries with initial states with the original time series generated with BoolNet
resultfile<-FBNBenchmark(FBNcellcyclenetwork,NULL,trainingseries,trainingseries,interval=1,temporal=1,FALSE,cube)
#genereate similarity report, i.e., compare the reconstructed time series with the original one
similar<-checkSimilarity(trainingseries,resultfile$FBNResult$reconstructed)
similarreport<-generateSimilarReport(similar)
print("Similar Report")
print(similarreport)
print("Benchmark result")
print(paste("ErrorRate=",resultfile$ErrorRate,sep="",collapse = ""))
print(paste("AccurateRate=",resultfile$AccurateRate,sep="",collapse = ""))
print(paste("MissMatchedRate=",resultfile$MissMatchedRate,sep="",collapse = ""))
print(paste("PerfectMatchedRate=",resultfile$PrfectMatchedRate,sep="",collapse = ""))
#get attractors
genes<-rownames(trainingseries [[1]])
attractor<-searchForAttractors(FBNcellcyclenetwork,startstates,genes)
#display the dynamic trajectory of the attactor 2
FBNNetwork.Graph.DrawAttractor(FBNcellcyclenetwork,attractor,2)
