Fundamental Boolean Network
================



# Installing and Loading Libraries

``` r
# packages_bioconductor = c("limma")
# 
# #use this function to check if each package is on the local machine
# #if a package is installed, it will be loaded
# #if any are not, the missing package(s) will be installed from Bioconductor and loaded
# package.check <- lapply(packages_bioconductor, FUN = function(x) {
#   if (!require(x, character.only = TRUE)) {
#     BiocManager::install(x, dependencies = TRUE)
#     library(x, character.only = TRUE)
#   }
# })

library(knitr)
library(BoolNet)
library(utils)
library(FBNNet)
library(visNetwork)


packages_cran = c("BoolNet","utils", "visNetwork", "FBNNet","tictoc")
  
#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from CRAN and loaded
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

#devtools::install_github("clsdavid/FBNNet2_public")
#devtools::install_github("collectivemedia/tictoc")

rm(package.check, packages_cran)
```

## Introduction

Fundamental Boolean modelling(FBM) has been proposed to draw insights
into gene activation, inhibition, and protein decay. This novel Boolean
model provides an intuitive definition of activation and inhibition
pathways and includes mechanisms to handle protein decay issues. To
prove the concept of the novel model, we implemented a platform using R
language, called FBNNet. Our experimental results show that the proposed
FBM could explicitly display the internal connections of the mammalian
cell cycle between genes separated into the connection types of
activation, inhibition and protein decay.

<https://clsdavid.github.io/FBNNet2_public/>

## Experiment

The experiments conducted and described here intend to prove the concept
of the new Boolean Model, i.e., the FBM. To verify the results, we apply
the general processes described in Figure 1 as a benchmark to compare
the results generated via BoolNet, with these consequently reconstructed
from the new R package, FBNNet.

## Data

# Theoretical Netwoork

``` r
hif1a.net <- loadNetwork("../data/ATOTS.bn")
hif1a.net
#> Boolean network with 6 genes
#> 
#> Involved genes:
#> EP300 HIF1A MDM2 TP53 VHL APOP
#> 
#> Transition functions:
#> EP300 = ((TP53 & HIF1A) & !VHL) | (!(TP53 & HIF1A) & VHL)
#> HIF1A = !VHL | (!TP53 & MDM2)
#> MDM2 = TP53 & !VHL
#> TP53 = !MDM2 | VHL
#> VHL = HIF1A & !TP53
#> APOP = TP53 & EP300 & (!HIF1A | !VHL | !MDM2)
```

## Extract the Fundamental Boolean Network

``` r

tic("Total - Extract FBN")
tic("generateAllCombinationBinary")  
   initialStates <- generateAllCombinationBinary(hif1a.net$genes)
toc()
#> generateAllCombinationBinary: 0.004 sec elapsed
tic("genereateBoolNetTimeseries") 
   trainingseries <- genereateBoolNetTimeseries(hif1a.net,
                                             initialStates,43,
                                             type = "synchronous")
toc()
#> genereateBoolNetTimeseries: 0.325 sec elapsed
tic("generateFBMNetwork") 
FBNcellcyclenetwork <- generateFBMNetwork(timeseries_data = trainingseries,
                                 maxK = 4,
                                 max_deep_temporal = 1,
                                 useParallel = FALSE,
                                 verbose = TRUE)
#> INFO [2021-02-05 17:29:09] Enter generateFBMNetwork zone: 
#>           method=kmeans,
#>           maxK = 4, 
#>           useParallel = FALSE, 
#>           max_deep_temporal = 1,
#>           threshold_confidence = 1,
#>           threshold_error = 0,
#>           threshold_support = 1e-05,
#>           maxFBNRules = 5,
#>           network_only = TRUE,
#>           verbose = TRUE
#> INFO [2021-02-05 17:29:09] Run generateFBMNetwork with a single cube
#> INFO [2021-02-05 17:29:09] Enter constructFBNCube zone: 
#>         target_genes=6 genes and they are EP300, HIF1A, MDM2, TP53, VHL, APOP,
#>         conditional_genes=6 genes and they are EP300, HIF1A, MDM2, TP53, VHL, APOP,
#>         data_length=64,
#>         maxK=4, 
#>         temporal=1,
#>         useParallel=FALSE
#> INFO [2021-02-05 17:29:10] Leave constructFBNCube zone.
#> INFO [2021-02-05 17:29:10] Enter mineFBNNetwork zone: genes=, useParallel=FALSE
#> INFO [2021-02-05 17:29:10] Enter search_FBN_core zone: useParallel=FALSE
#> INFO [2021-02-05 17:29:11] Leave search_FBN_core zone
#> INFO [2021-02-05 17:29:11] Enter mineFBNNetworkWithCores zone
#> INFO [2021-02-05 17:29:11] Enter mineFBNNetworkStage2 zone
#> INFO [2021-02-05 17:29:11] Leave mineFBNNetworkStage2 zone
#> INFO [2021-02-05 17:29:11] Enter convertMinedResultToFBNNetwork zone
#> INFO [2021-02-05 17:29:11] Leave convertMinedResultToFBNNetwork zone
#> INFO [2021-02-05 17:29:11] Leave mineFBNNetworkWithCores zone
#> INFO [2021-02-05 17:29:11] Leave mineFBNNetwork zone
toc()
#> generateFBMNetwork: 2.081 sec elapsed
toc()
#> Total - Extract FBN: 2.413 sec elapsed

print(FBNcellcyclenetwork)
#> Fundamental Boolean Network with  6 genes
#> Genes involved:
#> EP300, HIF1A, MDM2, TP53, VHL, APOP
#> 
#> Networks:
#> Multiple Transition Functions for EP300 with decay value = 1:
#> EP300_1_Activator: EP300 = !TP53&VHL (Confidence: 1, TimeStep: 1)
#> EP300_2_Activator: EP300 = !HIF1A&VHL (Confidence: 1, TimeStep: 1)
#> EP300_1_Inhibitor: EP300 = !HIF1A&!VHL (Confidence: 1, TimeStep: 1)
#> EP300_2_Inhibitor: EP300 = !TP53&!VHL (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for HIF1A with decay value = 1:
#> HIF1A_1_Activator: HIF1A = !VHL (Confidence: 1, TimeStep: 1)
#> HIF1A_2_Activator: HIF1A = MDM2&!TP53 (Confidence: 1, TimeStep: 1)
#> HIF1A_1_Inhibitor: HIF1A = !MDM2&VHL (Confidence: 1, TimeStep: 1)
#> HIF1A_2_Inhibitor: HIF1A = TP53&VHL (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for MDM2 with decay value = 1:
#> MDM2_1_Activator: MDM2 = TP53&!VHL (Confidence: 1, TimeStep: 1)
#> MDM2_1_Inhibitor: MDM2 = !TP53 (Confidence: 1, TimeStep: 1)
#> MDM2_2_Inhibitor: MDM2 = VHL (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for TP53 with decay value = 1:
#> TP53_1_Activator: TP53 = !MDM2 (Confidence: 1, TimeStep: 1)
#> TP53_2_Activator: TP53 = VHL (Confidence: 1, TimeStep: 1)
#> TP53_1_Inhibitor: TP53 = MDM2&!VHL (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for VHL with decay value = 1:
#> VHL_1_Activator: VHL = HIF1A&!TP53 (Confidence: 1, TimeStep: 1)
#> VHL_1_Inhibitor: VHL = TP53 (Confidence: 1, TimeStep: 1)
#> VHL_2_Inhibitor: VHL = !HIF1A (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for APOP with decay value = 1:
#> APOP_1_Activator: APOP = EP300&!MDM2&TP53 (Confidence: 1, TimeStep: 1)
#> APOP_2_Activator: APOP = EP300&!HIF1A&TP53 (Confidence: 1, TimeStep: 1)
#> APOP_1_Inhibitor: APOP = !EP300 (Confidence: 1, TimeStep: 1)
#> APOP_2_Inhibitor: APOP = !TP53 (Confidence: 1, TimeStep: 1)
```

``` r
   resultfile <- reconstructTimeseries(FBNcellcyclenetwork,
                                    initialStates,
                                    type = "synchronous",
                                    maxTimepoints = 43,
                                    useParallel = FALSE)
   similarreport <- generateSimilaryReport(trainingseries,resultfile)
   print(paste("ErrorRate=",similarreport$ErrorRate,sep = "",collapse = ""))
#> [1] "ErrorRate=0.0333091085271322"
   print(paste("AccurateRate=",similarreport$AccurateRate,sep = "",collapse = ""))
#> [1] "AccurateRate=0.966690891472868"
   print(paste("MissMatchedRate=",similarreport$MissMatchedRate,sep = "",collapse = ""))
#> [1] "MissMatchedRate=1"
   print(paste("PerfectMatchedRate=",similarreport$PerfectMatchedRate,sep = "",collapse = ""))
#> [1] "PerfectMatchedRate=0"
   #get attractors
   genes <- rownames(trainingseries[[1]])
   attractor <- searchForAttractors(FBNcellcyclenetwork,initialStates,genes)
   print(attractor)
#> Discovered Attractors via Fundamental Boolean Model :
#> Genes are encoded in the following order::
#> EP300 HIF1A MDM2 TP53 VHL APOP:
#> 
#> Attractor 1 is a complex attractor consisting of 5 state(s):
#> 
#> | --< - - - - - - |
#> v                 ^
#> 0 1 1 0 0 0       |
#> |                 |
#> 0 1 0 0 1 0       |
#> |                 |
#> 1 0 0 1 1 0       |
#> |                 |
#> 1 0 0 1 0 1       |
#> |                 |
#> 0 1 1 1 0 1       |
#> |                 |
#> v                 ^
#> | - - - - - - >-- |
   #display the dynamic trajectory of the attactor 2
   FBNNetwork.Graph.DrawAttractor(FBNcellcyclenetwork,attractor,1)
```

<!--html_preserve-->

<div id="htmlwidget-a251dbe4279bfbbd000a" class="visNetwork html-widget" style="width:100%;height:768px;">

</div>

<script type="application/json" data-for="htmlwidget-a251dbe4279bfbbd000a">{"x":{"nodes":{"id":["APOP_1","EP300_2","HIF1A_2","MDM2_2","TP53_2","VHL_2","APOP_2","EP300_1","HIF1A_1","MDM2_1","TP53_1","VHL_1","EP300_2_Inhibitor_2","HIF1A_1_Activator_2","HIF1A_2_Activator_2","MDM2_1_Inhibitor_2","TP53_1_Inhibitor_2","VHL_1_Activator_2","APOP_1_Inhibitor_2","APOP_2_Inhibitor_2","EP300_3","HIF1A_3","MDM2_3","TP53_3","VHL_3","APOP_3","EP300_1_Activator_3","HIF1A_1_Inhibitor_3","MDM2_1_Inhibitor_3","MDM2_2_Inhibitor_3","TP53_1_Activator_3","TP53_2_Activator_3","VHL_1_Activator_3","APOP_1_Inhibitor_3","APOP_2_Inhibitor_3","EP300_4","HIF1A_4","MDM2_4","TP53_4","VHL_4","APOP_4","EP300_2_Activator_4","HIF1A_1_Inhibitor_4","HIF1A_2_Inhibitor_4","MDM2_2_Inhibitor_4","TP53_1_Activator_4","TP53_2_Activator_4","VHL_1_Inhibitor_4","VHL_2_Inhibitor_4","APOP_1_Activator_4","APOP_2_Activator_4","EP300_5","HIF1A_5","MDM2_5","TP53_5","VHL_5","APOP_5","EP300_1_Inhibitor_5","HIF1A_1_Activator_5","MDM2_1_Activator_5","TP53_1_Activator_5","VHL_1_Inhibitor_5","VHL_2_Inhibitor_5","APOP_1_Activator_5","APOP_2_Activator_5","EP300_6","TP53_6","APOP_6","HIF1A_6","MDM2_6","VHL_6","HIF1A_1_Activator_6","MDM2_1_Activator_6","VHL_1_Inhibitor_6"],"shape":["ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","box","box","box","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","box","box","box","box","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","box","box","box","box","box","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","box","box","box","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box"],"color":["pink","pink","lightblue","pink","pink","lightblue","pink","pink","lightblue","lightblue","pink","pink","orange","LightGreen","LightGreen","orange","orange","LightGreen","orange","orange","lightblue","pink","pink","lightblue","lightblue","pink","LightGreen","orange","orange","orange","LightGreen","LightGreen","LightGreen","orange","orange","lightblue","pink","pink","lightblue","pink","lightblue","LightGreen","orange","orange","orange","LightGreen","LightGreen","orange","orange","LightGreen","LightGreen","pink","lightblue","lightblue","lightblue","pink","lightblue","orange","LightGreen","LightGreen","LightGreen","orange","orange","LightGreen","LightGreen","pink","pink","pink","lightblue","lightblue","pink","LightGreen","LightGreen","orange"],"type":["gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","TF","TF","TF","gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","TF","TF","TF","TF","gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","TF","TF","TF","gene","gene","gene","gene","gene","gene","TF","TF","TF"],"value":[4,4,4,4,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,4,4,4,4,4,4,2,2,2,2,2,2,2,2,2,4,4,4,4,4,4,2,2,2,2,2,2,2,2,2,2,4,4,4,4,4,4,2,2,2,2,2,2,2,2,4,4,4,4,4,4,2,2,2],"label":["APOP_1","EP300_2","HIF1A_2","MDM2_2","TP53_2","VHL_2","APOP_2","EP300_1","HIF1A_1","MDM2_1","TP53_1","VHL_1","-, 1","+, 1","+, 1","-, 1","-, 1","+, 1","-, 1","-, 1","EP300_3","HIF1A_3","MDM2_3","TP53_3","VHL_3","APOP_3","+, 1","-, 1","-, 1","-, 1","+, 1","+, 1","+, 1","-, 1","-, 1","EP300_4","HIF1A_4","MDM2_4","TP53_4","VHL_4","APOP_4","+, 1","-, 1","-, 1","-, 1","+, 1","+, 1","-, 1","-, 1","+, 1","+, 1","EP300_5","HIF1A_5","MDM2_5","TP53_5","VHL_5","APOP_5","-, 1","+, 1","+, 1","+, 1","-, 1","-, 1","+, 1","+, 1","EP300_6","TP53_6","APOP_6","HIF1A_6","MDM2_6","VHL_6","+, 1","+, 1","-, 1"],"shadow":[true,true,true,true,true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,true,true,true,true,true,true,false,false,false,false,false,false,false,false,true,true,true,true,true,true,false,false,false],"group":["Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Inhibit Function","Activate Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Inhibit Function","Inhibit Function","Gene","Gene","Gene","Gene","Gene","Gene","Activate Function","Inhibit Function","Inhibit Function","Inhibit Function","Activate Function","Activate Function","Activate Function","Inhibit Function","Inhibit Function","Gene","Gene","Gene","Gene","Gene","Gene","Activate Function","Inhibit Function","Inhibit Function","Inhibit Function","Activate Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Activate Function","Gene","Gene","Gene","Gene","Gene","Gene","Inhibit Function","Activate Function","Activate Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Activate Function","Gene","Gene","Gene","Gene","Gene","Gene","Activate Function","Activate Function","Inhibit Function"],"title":["APOP","EP300","HIF1A","MDM2","TP53","VHL","APOP","EP300","HIF1A","MDM2","TP53","VHL","!TP53&!VHL","!VHL","MDM2&!TP53","!TP53","MDM2&!VHL","HIF1A&!TP53","!EP300","!TP53","EP300","HIF1A","MDM2","TP53","VHL","APOP","!TP53&VHL","!MDM2&VHL","!TP53","VHL","!MDM2","VHL","HIF1A&!TP53","!EP300","!TP53","EP300","HIF1A","MDM2","TP53","VHL","APOP","!HIF1A&VHL","!MDM2&VHL","TP53&VHL","VHL","!MDM2","VHL","TP53","!HIF1A","EP300&!MDM2&TP53","EP300&!HIF1A&TP53","EP300","HIF1A","MDM2","TP53","VHL","APOP","!HIF1A&!VHL","!VHL","TP53&!VHL","!MDM2","TP53","!HIF1A","EP300&!MDM2&TP53","EP300&!HIF1A&TP53","EP300","TP53","APOP","HIF1A","MDM2","VHL","!VHL","TP53&!VHL","TP53"],"level":[1,3,3,3,3,3,3,1,1,1,1,1,2,2,2,2,2,2,2,2,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,9,9,9,9,9,9,8,8,8,8,8,8,8,8,11,11,11,11,11,11,10,10,10]},"edges":{"from":["EP300_2_Inhibitor_2","TP53_1","VHL_1","HIF1A_1_Activator_2","VHL_1","HIF1A_2_Activator_2","MDM2_1","TP53_1","MDM2_1_Inhibitor_2","TP53_1","TP53_1_Inhibitor_2","MDM2_1","VHL_1","VHL_1_Activator_2","HIF1A_1","TP53_1","APOP_1_Inhibitor_2","EP300_1","APOP_2_Inhibitor_2","TP53_1","EP300_1_Activator_3","TP53_2","VHL_2","HIF1A_1_Inhibitor_3","MDM2_2","VHL_2","MDM2_1_Inhibitor_3","TP53_2","MDM2_2_Inhibitor_3","VHL_2","TP53_1_Activator_3","MDM2_2","TP53_2_Activator_3","VHL_2","VHL_1_Activator_3","HIF1A_2","TP53_2","APOP_1_Inhibitor_3","EP300_2","APOP_2_Inhibitor_3","TP53_2","EP300_2_Activator_4","HIF1A_3","VHL_3","HIF1A_1_Inhibitor_4","MDM2_3","VHL_3","HIF1A_2_Inhibitor_4","TP53_3","VHL_3","MDM2_2_Inhibitor_4","VHL_3","TP53_1_Activator_4","MDM2_3","TP53_2_Activator_4","VHL_3","VHL_1_Inhibitor_4","TP53_3","VHL_2_Inhibitor_4","HIF1A_3","APOP_1_Activator_4","EP300_3","MDM2_3","TP53_3","APOP_2_Activator_4","EP300_3","HIF1A_3","TP53_3","EP300_1_Inhibitor_5","HIF1A_4","VHL_4","HIF1A_1_Activator_5","VHL_4","MDM2_1_Activator_5","TP53_4","VHL_4","TP53_1_Activator_5","MDM2_4","VHL_1_Inhibitor_5","TP53_4","VHL_2_Inhibitor_5","HIF1A_4","APOP_1_Activator_5","EP300_4","MDM2_4","TP53_4","APOP_2_Activator_5","EP300_4","HIF1A_4","TP53_4","HIF1A_1_Activator_6","VHL_5","MDM2_1_Activator_6","TP53_5","VHL_5","VHL_1_Inhibitor_6","TP53_5"],"to":["EP300_2","EP300_2_Inhibitor_2","EP300_2_Inhibitor_2","HIF1A_2","HIF1A_1_Activator_2","HIF1A_2","HIF1A_2_Activator_2","HIF1A_2_Activator_2","MDM2_2","MDM2_1_Inhibitor_2","TP53_2","TP53_1_Inhibitor_2","TP53_1_Inhibitor_2","VHL_2","VHL_1_Activator_2","VHL_1_Activator_2","APOP_2","APOP_1_Inhibitor_2","APOP_2","APOP_2_Inhibitor_2","EP300_3","EP300_1_Activator_3","EP300_1_Activator_3","HIF1A_3","HIF1A_1_Inhibitor_3","HIF1A_1_Inhibitor_3","MDM2_3","MDM2_1_Inhibitor_3","MDM2_3","MDM2_2_Inhibitor_3","TP53_3","TP53_1_Activator_3","TP53_3","TP53_2_Activator_3","VHL_3","VHL_1_Activator_3","VHL_1_Activator_3","APOP_3","APOP_1_Inhibitor_3","APOP_3","APOP_2_Inhibitor_3","EP300_4","EP300_2_Activator_4","EP300_2_Activator_4","HIF1A_4","HIF1A_1_Inhibitor_4","HIF1A_1_Inhibitor_4","HIF1A_4","HIF1A_2_Inhibitor_4","HIF1A_2_Inhibitor_4","MDM2_4","MDM2_2_Inhibitor_4","TP53_4","TP53_1_Activator_4","TP53_4","TP53_2_Activator_4","VHL_4","VHL_1_Inhibitor_4","VHL_4","VHL_2_Inhibitor_4","APOP_4","APOP_1_Activator_4","APOP_1_Activator_4","APOP_1_Activator_4","APOP_4","APOP_2_Activator_4","APOP_2_Activator_4","APOP_2_Activator_4","EP300_5","EP300_1_Inhibitor_5","EP300_1_Inhibitor_5","HIF1A_5","HIF1A_1_Activator_5","MDM2_5","MDM2_1_Activator_5","MDM2_1_Activator_5","TP53_5","TP53_1_Activator_5","VHL_5","VHL_1_Inhibitor_5","VHL_5","VHL_2_Inhibitor_5","APOP_5","APOP_1_Activator_5","APOP_1_Activator_5","APOP_1_Activator_5","APOP_5","APOP_2_Activator_5","APOP_2_Activator_5","APOP_2_Activator_5","HIF1A_6","HIF1A_1_Activator_6","MDM2_6","MDM2_1_Activator_6","MDM2_1_Activator_6","VHL_6","VHL_1_Inhibitor_6"],"support":[0.20089,0.20089,0.20089,0.60119,0.60119,0.1994,0.1994,0.1994,0.39881,0.39881,0.39286,0.39286,0.39286,0.39286,0.39286,0.39286,0.40774,0.40774,0.39881,0.39881,0.19792,0.19792,0.19792,0.39286,0.39286,0.39286,0.39881,0.39881,0.39881,0.39881,0.60119,0.60119,0.39881,0.39881,0.39286,0.39286,0.39286,0.40774,0.40774,0.39881,0.39881,0.19643,0.19643,0.19643,0.39286,0.39286,0.39286,0.20089,0.20089,0.20089,0.39881,0.39881,0.60119,0.60119,0.39881,0.39881,0.60119,0.60119,0.4003,0.4003,0.3869,0.3869,0.3869,0.3869,0.38393,0.38393,0.38393,0.38393,0.20387,0.20387,0.20387,0.60119,0.60119,0.4003,0.4003,0.4003,0.60119,0.60119,0.60119,0.60119,0.4003,0.4003,0.3869,0.3869,0.3869,0.3869,0.38393,0.38393,0.38393,0.38393,0.60119,0.60119,0.4003,0.4003,0.4003,0.60119,0.60119],"targetNode":["EP300","EP300","EP300","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","MDM2","MDM2","TP53","TP53","TP53","VHL","VHL","VHL","APOP","APOP","APOP","APOP","EP300","EP300","EP300","HIF1A","HIF1A","HIF1A","MDM2","MDM2","MDM2","MDM2","TP53","TP53","TP53","TP53","VHL","VHL","VHL","APOP","APOP","APOP","APOP","EP300","EP300","EP300","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","MDM2","MDM2","TP53","TP53","TP53","TP53","VHL","VHL","VHL","VHL","APOP","APOP","APOP","APOP","APOP","APOP","APOP","APOP","EP300","EP300","EP300","HIF1A","HIF1A","MDM2","MDM2","MDM2","TP53","TP53","VHL","VHL","VHL","VHL","APOP","APOP","APOP","APOP","APOP","APOP","APOP","APOP","HIF1A","HIF1A","MDM2","MDM2","MDM2","VHL","VHL"],"title":["EP300_2_Inhibitor","EP300_2_Inhibitor","EP300_2_Inhibitor","HIF1A_1_Activator","HIF1A_1_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","MDM2_1_Inhibitor","MDM2_1_Inhibitor","TP53_1_Inhibitor","TP53_1_Inhibitor","TP53_1_Inhibitor","VHL_1_Activator","VHL_1_Activator","VHL_1_Activator","APOP_1_Inhibitor","APOP_1_Inhibitor","APOP_2_Inhibitor","APOP_2_Inhibitor","EP300_1_Activator","EP300_1_Activator","EP300_1_Activator","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","MDM2_1_Inhibitor","MDM2_1_Inhibitor","MDM2_2_Inhibitor","MDM2_2_Inhibitor","TP53_1_Activator","TP53_1_Activator","TP53_2_Activator","TP53_2_Activator","VHL_1_Activator","VHL_1_Activator","VHL_1_Activator","APOP_1_Inhibitor","APOP_1_Inhibitor","APOP_2_Inhibitor","APOP_2_Inhibitor","EP300_2_Activator","EP300_2_Activator","EP300_2_Activator","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","HIF1A_2_Inhibitor","HIF1A_2_Inhibitor","HIF1A_2_Inhibitor","MDM2_2_Inhibitor","MDM2_2_Inhibitor","TP53_1_Activator","TP53_1_Activator","TP53_2_Activator","TP53_2_Activator","VHL_1_Inhibitor","VHL_1_Inhibitor","VHL_2_Inhibitor","VHL_2_Inhibitor","APOP_1_Activator","APOP_1_Activator","APOP_1_Activator","APOP_1_Activator","APOP_2_Activator","APOP_2_Activator","APOP_2_Activator","APOP_2_Activator","EP300_1_Inhibitor","EP300_1_Inhibitor","EP300_1_Inhibitor","HIF1A_1_Activator","HIF1A_1_Activator","MDM2_1_Activator","MDM2_1_Activator","MDM2_1_Activator","TP53_1_Activator","TP53_1_Activator","VHL_1_Inhibitor","VHL_1_Inhibitor","VHL_2_Inhibitor","VHL_2_Inhibitor","APOP_1_Activator","APOP_1_Activator","APOP_1_Activator","APOP_1_Activator","APOP_2_Activator","APOP_2_Activator","APOP_2_Activator","APOP_2_Activator","HIF1A_1_Activator","HIF1A_1_Activator","MDM2_1_Activator","MDM2_1_Activator","MDM2_1_Activator","VHL_1_Inhibitor","VHL_1_Inhibitor"],"arrowtail":["none","tee","tee","none","tee","none","none","tee","none","tee","none","none","tee","none","none","tee","none","tee","none","tee","none","tee","none","none","tee","none","none","tee","none","none","none","tee","none","none","none","none","tee","none","tee","none","tee","none","tee","none","none","tee","none","none","none","none","none","none","none","tee","none","none","none","none","none","tee","none","none","tee","none","none","none","tee","none","none","tee","tee","none","tee","none","none","tee","none","tee","none","none","none","tee","none","none","tee","none","none","none","tee","none","none","tee","none","none","tee","none","none"],"arrowhead":["odot","none","none","open","none","open","none","none","odot","none","odot","none","none","open","none","none","odot","none","odot","none","open","none","none","odot","none","none","odot","none","odot","none","open","none","open","none","open","none","none","odot","none","odot","none","open","none","none","odot","none","none","odot","none","none","odot","none","open","none","open","none","odot","none","odot","none","open","none","none","none","open","none","none","none","odot","none","none","open","none","open","none","none","open","none","odot","none","odot","none","open","none","none","none","open","none","none","none","open","none","open","none","none","odot","none"],"color":["darkred","red","red","darkblue","red","darkblue","green","red","darkred","red","darkred","green","red","darkblue","green","red","darkred","red","darkred","red","darkblue","red","green","darkred","red","green","darkred","red","darkred","green","darkblue","red","darkblue","green","darkblue","green","red","darkred","red","darkred","red","darkblue","red","green","darkred","red","green","darkred","green","green","darkred","green","darkblue","red","darkblue","green","darkred","green","darkred","red","darkblue","green","red","green","darkblue","green","red","green","darkred","red","red","darkblue","red","darkblue","green","red","darkblue","red","darkred","green","darkred","red","darkblue","green","red","green","darkblue","green","red","green","darkblue","red","darkblue","green","red","darkred","green"],"lty":["longdash","longdash","longdash","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","longdash","longdash","longdash","longdash","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","solid","longdash","longdash","longdash","longdash","solid","solid","solid","solid","solid","solid","solid","solid","longdash","longdash","longdash","solid","solid","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash","solid","solid","solid","solid","solid","solid","solid","solid","solid","solid","solid","solid","solid","longdash","longdash"],"arrow.mode":[2,0,0,2,0,2,0,0,2,0,2,0,0,2,0,0,2,0,2,0,2,0,0,2,0,0,2,0,2,0,2,0,2,0,2,0,0,2,0,2,0,2,0,0,2,0,0,2,0,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,2,0,0,0,2,0,0,2,0,2,0,0,2,0,2,0,2,0,2,0,0,0,2,0,0,0,2,0,2,0,0,2,0],"arrow.size":[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2],"type":["TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF"],"arrows":["to","none","none","to","none","to","none","none","to","none","to","none","none","to","none","none","to","none","to","none","to","none","none","to","none","none","to","none","to","none","to","none","to","none","to","none","none","to","none","to","none","to","none","none","to","none","none","to","none","none","to","none","to","none","to","none","to","none","to","none","to","none","none","none","to","none","none","none","to","none","none","to","none","to","none","none","to","none","to","none","to","none","to","none","none","none","to","none","none","none","to","none","to","none","none","to","none"],"label":["","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"dashes":[true,true,true,false,true,false,false,true,true,true,true,false,true,false,false,true,true,true,true,true,false,true,false,true,true,false,true,true,true,false,false,true,false,false,false,false,true,true,true,true,true,false,true,false,true,true,false,true,false,false,true,false,false,true,false,false,true,false,true,true,false,false,true,false,false,false,true,false,true,true,true,false,true,false,false,true,false,true,true,false,true,true,false,false,true,false,false,false,true,false,false,true,false,false,true,true,false],"shadow":[true,false,false,true,false,true,false,false,true,false,true,false,false,true,false,false,true,false,true,false,true,false,false,true,false,false,true,false,true,false,true,false,true,false,true,false,false,true,false,true,false,true,false,false,true,false,false,true,false,false,true,false,true,false,true,false,true,false,true,false,true,false,false,false,true,false,false,false,true,false,false,true,false,true,false,false,true,false,true,false,true,false,true,false,false,false,true,false,false,false,true,false,true,false,false,true,false]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot"},"manipulation":{"enabled":false},"interaction":{"dragNodes":true,"dragView":true,"zoomView":true},"edges":{"arrows":"from"},"layout":{"hierarchical":{"enabled":true,"levelSeparation":200,"direction":"LR"}}},"groups":["Gene","Inhibit Function","Activate Function"],"width":null,"height":null,"idselection":{"enabled":false},"byselection":{"enabled":false},"main":{"text":"Dynamic Fundamental Boolean Networks from the time point of 1 to 6","style":"font-family:Georgia, Times New Roman, Times, serif;font-weight:bold;font-size:20px;text-align:center;"},"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","legend":{"width":0.1,"useGroups":false,"position":"left","ncol":1,"stepX":100,"stepY":100,"zoom":true,"edges":{"color":["darkblue","darkred","grey","green","red"],"label":["activate","inhibit","decay","activated input","deactivated input"],"arrows":["to","to","to","none","none"],"dashes":[false,true,true,false,true],"shadow":[true,true,true,false,false]},"edgesToDataframe":true,"nodes":{"label":["Activated Gene","Inhibited Gene","Activate Function (+)","Inhibit Function (-)"],"shape":["ellipse","ellipse","box","box"],"color":["lightblue","pink","lightgreen","orange"],"shadow":[true,true,false,false]},"nodesToDataframe":true},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);"},"evals":[],"jsHooks":[]}</script>

<!--/html_preserve-->

``` r

FBNNet::FBNNetwork.Graph(FBNcellcyclenetwork)
```

<!--html_preserve-->

<div id="htmlwidget-28231faff821b6c640b9" class="visNetwork html-widget" style="width:100%;height:768px;">

</div>

<script type="application/json" data-for="htmlwidget-28231faff821b6c640b9">{"x":{"nodes":{"id":["EP300","HIF1A","MDM2","TP53","VHL","APOP","EP300_1_Activator","EP300_2_Activator","EP300_1_Inhibitor","EP300_2_Inhibitor","HIF1A_1_Activator","HIF1A_2_Activator","HIF1A_1_Inhibitor","HIF1A_2_Inhibitor","MDM2_1_Activator","MDM2_1_Inhibitor","MDM2_2_Inhibitor","TP53_1_Activator","TP53_2_Activator","TP53_1_Inhibitor","VHL_1_Activator","VHL_1_Inhibitor","VHL_2_Inhibitor","APOP_1_Activator","APOP_2_Activator","APOP_1_Inhibitor","APOP_2_Inhibitor"],"shape":["ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box"],"color":["lightblue","lightblue","lightblue","lightblue","lightblue","lightblue","LightGreen","LightGreen","orange","orange","LightGreen","LightGreen","orange","orange","LightGreen","orange","orange","LightGreen","LightGreen","orange","LightGreen","orange","orange","LightGreen","LightGreen","orange","orange"],"type":["gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF"],"value":[4,4,4,4,4,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],"label":["EP300","HIF1A","MDM2","TP53","VHL","APOP","+, 1","+, 1","-, 1","-, 1","+, 1","+, 1","-, 1","-, 1","+, 1","-, 1","-, 1","+, 1","+, 1","-, 1","+, 1","-, 1","-, 1","+, 1","+, 1","-, 1","-, 1"],"shadow":[true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false],"group":["Gene","Gene","Gene","Gene","Gene","Gene","Activate Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Activate Function","Inhibit Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Activate Function","Inhibit Function","Inhibit Function"],"title":["E1A binding protein p300(EP300)","hypoxia inducible factor 1 alpha subunit(HIF1A)","MDM2 proto-oncogene(MDM2)","tumor protein p53(TP53)","von Hippel-Lindau tumor suppressor(VHL)",null,"!TP53&VHL","!HIF1A&VHL","!HIF1A&!VHL","!TP53&!VHL","!VHL","MDM2&!TP53","!MDM2&VHL","TP53&VHL","TP53&!VHL","!TP53","VHL","!MDM2","VHL","MDM2&!VHL","HIF1A&!TP53","TP53","!HIF1A","EP300&!MDM2&TP53","EP300&!HIF1A&TP53","!EP300","!TP53"]},"edges":{"from":["EP300_1_Activator","TP53","VHL","EP300_2_Activator","HIF1A","VHL","EP300_1_Inhibitor","HIF1A","VHL","EP300_2_Inhibitor","TP53","VHL","HIF1A_1_Activator","VHL","HIF1A_2_Activator","MDM2","TP53","HIF1A_1_Inhibitor","MDM2","VHL","HIF1A_2_Inhibitor","TP53","VHL","MDM2_1_Activator","TP53","VHL","MDM2_1_Inhibitor","TP53","MDM2_2_Inhibitor","VHL","TP53_1_Activator","MDM2","TP53_2_Activator","VHL","TP53_1_Inhibitor","MDM2","VHL","VHL_1_Activator","HIF1A","TP53","VHL_1_Inhibitor","TP53","VHL_2_Inhibitor","HIF1A","APOP_1_Activator","EP300","MDM2","TP53","APOP_2_Activator","EP300","HIF1A","TP53","APOP_1_Inhibitor","EP300","APOP_2_Inhibitor","TP53"],"to":["EP300","EP300_1_Activator","EP300_1_Activator","EP300","EP300_2_Activator","EP300_2_Activator","EP300","EP300_1_Inhibitor","EP300_1_Inhibitor","EP300","EP300_2_Inhibitor","EP300_2_Inhibitor","HIF1A","HIF1A_1_Activator","HIF1A","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","HIF1A","HIF1A_2_Inhibitor","HIF1A_2_Inhibitor","MDM2","MDM2_1_Activator","MDM2_1_Activator","MDM2","MDM2_1_Inhibitor","MDM2","MDM2_2_Inhibitor","TP53","TP53_1_Activator","TP53","TP53_2_Activator","TP53","TP53_1_Inhibitor","TP53_1_Inhibitor","VHL","VHL_1_Activator","VHL_1_Activator","VHL","VHL_1_Inhibitor","VHL","VHL_2_Inhibitor","APOP","APOP_1_Activator","APOP_1_Activator","APOP_1_Activator","APOP","APOP_2_Activator","APOP_2_Activator","APOP_2_Activator","APOP","APOP_1_Inhibitor","APOP","APOP_2_Inhibitor"],"support":[0.19792,0.19792,0.19792,0.19643,0.19643,0.19643,0.20387,0.20387,0.20387,0.20089,0.20089,0.20089,0.60119,0.60119,0.1994,0.1994,0.1994,0.39286,0.39286,0.39286,0.20089,0.20089,0.20089,0.4003,0.4003,0.4003,0.39881,0.39881,0.39881,0.39881,0.60119,0.60119,0.39881,0.39881,0.39286,0.39286,0.39286,0.39286,0.39286,0.39286,0.60119,0.60119,0.4003,0.4003,0.3869,0.3869,0.3869,0.3869,0.38393,0.38393,0.38393,0.38393,0.40774,0.40774,0.39881,0.39881],"targetNode":["EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","MDM2","MDM2","MDM2","MDM2","MDM2","MDM2","MDM2","TP53","TP53","TP53","TP53","TP53","TP53","TP53","VHL","VHL","VHL","VHL","VHL","VHL","VHL","APOP","APOP","APOP","APOP","APOP","APOP","APOP","APOP","APOP","APOP","APOP","APOP"],"title":["EP300_1_Activator","EP300_1_Activator","EP300_1_Activator","EP300_2_Activator","EP300_2_Activator","EP300_2_Activator","EP300_1_Inhibitor","EP300_1_Inhibitor","EP300_1_Inhibitor","EP300_2_Inhibitor","EP300_2_Inhibitor","EP300_2_Inhibitor","HIF1A_1_Activator","HIF1A_1_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","HIF1A_2_Inhibitor","HIF1A_2_Inhibitor","HIF1A_2_Inhibitor","MDM2_1_Activator","MDM2_1_Activator","MDM2_1_Activator","MDM2_1_Inhibitor","MDM2_1_Inhibitor","MDM2_2_Inhibitor","MDM2_2_Inhibitor","TP53_1_Activator","TP53_1_Activator","TP53_2_Activator","TP53_2_Activator","TP53_1_Inhibitor","TP53_1_Inhibitor","TP53_1_Inhibitor","VHL_1_Activator","VHL_1_Activator","VHL_1_Activator","VHL_1_Inhibitor","VHL_1_Inhibitor","VHL_2_Inhibitor","VHL_2_Inhibitor","APOP_1_Activator","APOP_1_Activator","APOP_1_Activator","APOP_1_Activator","APOP_2_Activator","APOP_2_Activator","APOP_2_Activator","APOP_2_Activator","APOP_1_Inhibitor","APOP_1_Inhibitor","APOP_2_Inhibitor","APOP_2_Inhibitor"],"arrowtail":["none","tee","none","none","tee","none","none","tee","tee","none","tee","tee","none","tee","none","none","tee","none","tee","none","none","none","none","none","none","tee","none","tee","none","none","none","tee","none","none","none","none","tee","none","none","tee","none","none","none","tee","none","none","tee","none","none","none","tee","none","none","tee","none","tee"],"arrowhead":["open","none","none","open","none","none","odot","none","none","odot","none","none","open","none","open","none","none","odot","none","none","odot","none","none","open","none","none","odot","none","odot","none","open","none","open","none","odot","none","none","open","none","none","odot","none","odot","none","open","none","none","none","open","none","none","none","odot","none","odot","none"],"color":["darkblue","red","green","darkblue","red","green","darkred","red","red","darkred","red","red","darkblue","red","darkblue","green","red","darkred","red","green","darkred","green","green","darkblue","green","red","darkred","red","darkred","green","darkblue","red","darkblue","green","darkred","green","red","darkblue","green","red","darkred","green","darkred","red","darkblue","green","red","green","darkblue","green","red","green","darkred","red","darkred","red"],"lty":["solid","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","longdash","longdash","longdash","longdash","solid","solid","solid","solid","longdash","longdash","longdash","solid","solid","solid","longdash","longdash","longdash","longdash","solid","solid","solid","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash"],"arrow.mode":[2,0,0,2,0,0,2,0,0,2,0,0,2,0,2,0,0,2,0,0,2,0,0,2,0,0,2,0,2,0,2,0,2,0,2,0,0,2,0,0,2,0,2,0,2,0,0,0,2,0,0,0,2,0,2,0],"arrow.size":[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2],"type":["TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF"],"arrows":["to","none","none","to","none","none","to","none","none","to","none","none","to","none","to","none","none","to","none","none","to","none","none","to","none","none","to","none","to","none","to","none","to","none","to","none","none","to","none","none","to","none","to","none","to","none","none","none","to","none","none","none","to","none","to","none"],"label":["","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"dashes":[false,true,false,false,true,false,true,true,true,true,true,true,false,true,false,false,true,true,true,false,true,false,false,false,false,true,true,true,true,false,false,true,false,false,true,false,true,false,false,true,true,false,true,true,false,false,true,false,false,false,true,false,true,true,true,true],"shadow":[true,false,false,true,false,false,true,false,false,true,false,false,true,false,true,false,false,true,false,false,true,false,false,true,false,false,true,false,true,false,true,false,true,false,true,false,false,true,false,false,true,false,true,false,true,false,false,false,true,false,false,false,true,false,true,false]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot"},"manipulation":{"enabled":false},"interaction":{"dragNodes":true,"dragView":true,"zoomView":true}},"groups":["Gene","Activate Function","Inhibit Function"],"width":null,"height":null,"idselection":{"enabled":false},"byselection":{"enabled":false},"main":{"text":"Fundamental Boolean Networks","style":"font-family:Georgia, Times New Roman, Times, serif;font-weight:bold;font-size:20px;text-align:center;"},"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","legend":{"width":0.2,"useGroups":false,"position":"left","ncol":1,"stepX":100,"stepY":100,"zoom":true,"edges":{"color":["darkblue","darkred","green","red"],"label":["activate","inhibit","activated input","deactivated input"],"arrows":["to","to","none","none"],"dashes":[false,true,false,true],"shadow":[true,true,false,false]},"edgesToDataframe":true,"nodes":{"label":["Gene","Activate Function (+, Timestep)","Inhibit Function (-, Timestep)"],"shape":["ellipse","box","box"],"color":["lightblue","lightgreen","orange"],"shadow":[true,false,false]},"nodesToDataframe":true},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);"},"evals":[],"jsHooks":[]}</script>

<!--/html_preserve-->

# Inferred Netwoork

``` r
hif1a.net <- loadNetwork("../data/ATOTS_inferred.bn")
hif1a.net
#> Boolean network with 5 genes
#> 
#> Involved genes:
#> EP300 HIF1A MDM2 TP53 VHL
#> 
#> Transition functions:
#> EP300 = (!HIF1A & TP53) | (HIF1A & !TP53)
#> HIF1A = EP300
#> MDM2 = EP300
#> TP53 = HIF1A
#> VHL = !TP53
```

## Extract the Fundamental Boolean Network

``` r

tic("Total - Extract FBN")
tic("generateAllCombinationBinary")  
   initialStates <- generateAllCombinationBinary(hif1a.net$genes)
toc()
#> generateAllCombinationBinary: 0.001 sec elapsed
tic("genereateBoolNetTimeseries") 
   trainingseries <- genereateBoolNetTimeseries(hif1a.net,
                                             initialStates,43,
                                             type = "synchronous")
toc()
#> genereateBoolNetTimeseries: 0.148 sec elapsed
tic("generateFBMNetwork") 
FBNcellcyclenetwork <- generateFBMNetwork(timeseries_data = trainingseries,
                                 maxK = 4,
                                 max_deep_temporal = 1,
                                 useParallel = FALSE,
                                 verbose = TRUE)
#> INFO [2021-02-05 17:29:15] Enter generateFBMNetwork zone: 
#>           method=kmeans,
#>           maxK = 4, 
#>           useParallel = FALSE, 
#>           max_deep_temporal = 1,
#>           threshold_confidence = 1,
#>           threshold_error = 0,
#>           threshold_support = 1e-05,
#>           maxFBNRules = 5,
#>           network_only = TRUE,
#>           verbose = TRUE
#> INFO [2021-02-05 17:29:15] Run generateFBMNetwork with a single cube
#> INFO [2021-02-05 17:29:15] Enter constructFBNCube zone: 
#>         target_genes=5 genes and they are EP300, HIF1A, MDM2, TP53, VHL,
#>         conditional_genes=5 genes and they are EP300, HIF1A, MDM2, TP53, VHL,
#>         data_length=32,
#>         maxK=4, 
#>         temporal=1,
#>         useParallel=FALSE
#> INFO [2021-02-05 17:29:15] Leave constructFBNCube zone.
#> INFO [2021-02-05 17:29:15] Enter mineFBNNetwork zone: genes=, useParallel=FALSE
#> INFO [2021-02-05 17:29:15] Enter search_FBN_core zone: useParallel=FALSE
#> INFO [2021-02-05 17:29:15] Leave search_FBN_core zone
#> INFO [2021-02-05 17:29:15] Enter mineFBNNetworkWithCores zone
#> INFO [2021-02-05 17:29:15] Enter mineFBNNetworkStage2 zone
#> INFO [2021-02-05 17:29:15] Leave mineFBNNetworkStage2 zone
#> INFO [2021-02-05 17:29:15] Enter convertMinedResultToFBNNetwork zone
#> INFO [2021-02-05 17:29:15] Leave convertMinedResultToFBNNetwork zone
#> INFO [2021-02-05 17:29:15] Leave mineFBNNetworkWithCores zone
#> INFO [2021-02-05 17:29:15] Leave mineFBNNetwork zone
toc()
#> generateFBMNetwork: 0.094 sec elapsed
toc()
#> Total - Extract FBN: 0.246 sec elapsed

print(FBNcellcyclenetwork)
#> Fundamental Boolean Network with  5 genes
#> Genes involved:
#> HIF1A, MDM2, TP53, VHL, EP300
#> 
#> Networks:
#> Multiple Transition Functions for EP300 with decay value = 1:
#> 
#> Multiple Transition Functions for HIF1A with decay value = 1:
#> HIF1A_1_Activator: HIF1A = EP300 (Confidence: 1, TimeStep: 1)
#> HIF1A_1_Inhibitor: HIF1A = !EP300 (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for MDM2 with decay value = 1:
#> MDM2_1_Activator: MDM2 = EP300 (Confidence: 1, TimeStep: 1)
#> MDM2_1_Inhibitor: MDM2 = !EP300 (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for TP53 with decay value = 1:
#> TP53_1_Activator: TP53 = HIF1A (Confidence: 1, TimeStep: 1)
#> TP53_1_Inhibitor: TP53 = !HIF1A (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for VHL with decay value = 1:
#> VHL_1_Activator: VHL = !TP53 (Confidence: 1, TimeStep: 1)
#> VHL_1_Inhibitor: VHL = TP53 (Confidence: 1, TimeStep: 1)
```

``` r
   resultfile <- reconstructTimeseries(FBNcellcyclenetwork,
                                    initialStates,
                                    type = "synchronous",
                                    maxTimepoints = 43,
                                    useParallel = FALSE)
   similarreport <- generateSimilaryReport(trainingseries,resultfile)
   print(paste("ErrorRate=",similarreport$ErrorRate,sep = "",collapse = ""))
#> [1] "ErrorRate=0.488372093023256"
   print(paste("AccurateRate=",similarreport$AccurateRate,sep = "",collapse = ""))
#> [1] "AccurateRate=0.511627906976744"
   print(paste("MissMatchedRate=",similarreport$MissMatchedRate,sep = "",collapse = ""))
#> [1] "MissMatchedRate=1"
   print(paste("PerfectMatchedRate=",similarreport$PerfectMatchedRate,sep = "",collapse = ""))
#> [1] "PerfectMatchedRate=0"
   #get attractors
   #genes <- rownames(trainingseries[[1]])
   #attractor <- searchForAttractors(FBNcellcyclenetwork,initialStates,genes)
   print(attractor)
#> Discovered Attractors via Fundamental Boolean Model :
#> Genes are encoded in the following order::
#> EP300 HIF1A MDM2 TP53 VHL APOP:
#> 
#> Attractor 1 is a complex attractor consisting of 5 state(s):
#> 
#> | --< - - - - - - |
#> v                 ^
#> 0 1 1 0 0 0       |
#> |                 |
#> 0 1 0 0 1 0       |
#> |                 |
#> 1 0 0 1 1 0       |
#> |                 |
#> 1 0 0 1 0 1       |
#> |                 |
#> 0 1 1 1 0 1       |
#> |                 |
#> v                 ^
#> | - - - - - - >-- |
   #display the dynamic trajectory of the attactor 2
   #FBNNetwork.Graph.DrawAttractor(FBNcellcyclenetwork,attractor,1)
```

``` r

FBNNet::FBNNetwork.Graph(FBNcellcyclenetwork)
```

<!--html_preserve-->

<div id="htmlwidget-d563ad4ef4f35bc31a80" class="visNetwork html-widget" style="width:100%;height:768px;">

</div>

<script type="application/json" data-for="htmlwidget-d563ad4ef4f35bc31a80">{"x":{"nodes":{"id":["HIF1A","MDM2","TP53","VHL","EP300","HIF1A_1_Activator","HIF1A_1_Inhibitor","MDM2_1_Activator","MDM2_1_Inhibitor","TP53_1_Activator","TP53_1_Inhibitor","VHL_1_Activator","VHL_1_Inhibitor"],"shape":["ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","box","box","box"],"color":["lightblue","lightblue","lightblue","lightblue","lightblue","LightGreen","orange","LightGreen","orange","LightGreen","orange","LightGreen","orange"],"type":["gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","TF","TF","TF"],"value":[4,4,4,4,4,2,2,2,2,2,2,2,2],"label":["HIF1A","MDM2","TP53","VHL","EP300","+, 1","-, 1","+, 1","-, 1","+, 1","-, 1","+, 1","-, 1"],"shadow":[true,true,true,true,true,false,false,false,false,false,false,false,false],"group":["Gene","Gene","Gene","Gene","Gene","Activate Function","Inhibit Function","Activate Function","Inhibit Function","Activate Function","Inhibit Function","Activate Function","Inhibit Function"],"title":["hypoxia inducible factor 1 alpha subunit(HIF1A)","MDM2 proto-oncogene(MDM2)","tumor protein p53(TP53)","von Hippel-Lindau tumor suppressor(VHL)","E1A binding protein p300(EP300)","EP300","!EP300","EP300","!EP300","HIF1A","!HIF1A","!TP53","TP53"]},"edges":{"from":["HIF1A_1_Activator","EP300","HIF1A_1_Inhibitor","EP300","MDM2_1_Activator","EP300","MDM2_1_Inhibitor","EP300","TP53_1_Activator","HIF1A","TP53_1_Inhibitor","HIF1A","VHL_1_Activator","TP53","VHL_1_Inhibitor","TP53"],"to":["HIF1A","HIF1A_1_Activator","HIF1A","HIF1A_1_Inhibitor","MDM2","MDM2_1_Activator","MDM2","MDM2_1_Inhibitor","TP53","TP53_1_Activator","TP53","TP53_1_Inhibitor","VHL","VHL_1_Activator","VHL","VHL_1_Inhibitor"],"support":[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],"targetNode":["HIF1A","HIF1A","HIF1A","HIF1A","MDM2","MDM2","MDM2","MDM2","TP53","TP53","TP53","TP53","VHL","VHL","VHL","VHL"],"title":["HIF1A_1_Activator","HIF1A_1_Activator","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","MDM2_1_Activator","MDM2_1_Activator","MDM2_1_Inhibitor","MDM2_1_Inhibitor","TP53_1_Activator","TP53_1_Activator","TP53_1_Inhibitor","TP53_1_Inhibitor","VHL_1_Activator","VHL_1_Activator","VHL_1_Inhibitor","VHL_1_Inhibitor"],"arrowtail":["none","none","none","tee","none","none","none","tee","none","none","none","tee","none","tee","none","none"],"arrowhead":["open","none","odot","none","open","none","odot","none","open","none","odot","none","open","none","odot","none"],"color":["darkblue","green","darkred","red","darkblue","green","darkred","red","darkblue","green","darkred","red","darkblue","red","darkred","green"],"lty":["solid","solid","longdash","longdash","solid","solid","longdash","longdash","solid","solid","longdash","longdash","solid","solid","longdash","longdash"],"arrow.mode":[2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0],"arrow.size":[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2],"type":["TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF"],"arrows":["to","none","to","none","to","none","to","none","to","none","to","none","to","none","to","none"],"label":["","","","","","","","","","","","","","","",""],"dashes":[false,false,true,true,false,false,true,true,false,false,true,true,false,true,true,false],"shadow":[true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot"},"manipulation":{"enabled":false},"interaction":{"dragNodes":true,"dragView":true,"zoomView":true}},"groups":["Gene","Activate Function","Inhibit Function"],"width":null,"height":null,"idselection":{"enabled":false},"byselection":{"enabled":false},"main":{"text":"Fundamental Boolean Networks","style":"font-family:Georgia, Times New Roman, Times, serif;font-weight:bold;font-size:20px;text-align:center;"},"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","legend":{"width":0.2,"useGroups":false,"position":"left","ncol":1,"stepX":100,"stepY":100,"zoom":true,"edges":{"color":["darkblue","darkred","green","red"],"label":["activate","inhibit","activated input","deactivated input"],"arrows":["to","to","none","none"],"dashes":[false,true,false,true],"shadow":[true,true,false,false]},"edgesToDataframe":true,"nodes":{"label":["Gene","Activate Function (+, Timestep)","Inhibit Function (-, Timestep)"],"shape":["ellipse","box","box"],"color":["lightblue","lightgreen","orange"],"shadow":[true,false,false]},"nodesToDataframe":true},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);"},"evals":[],"jsHooks":[]}</script>

<!--/html_preserve-->
