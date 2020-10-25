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

``` r
hif1a.net<-loadNetwork("../data/reduced_HIFaxis.bn")
hif1a.net
#> Boolean network with 6 genes
#> 
#> Involved genes:
#> EP300 HIF1A MDM2 TP53 VHL O2
#> 
#> Transition functions:
#> EP300 = (!VHL & TP53) | (VHL & !TP53)
#> HIF1A = !VHL & ((!O2 & EP300 & !TP53) | (!O2 & !MDM2))
#> MDM2 = TP53 & !VHL
#> TP53 = !MDM2 | (!O2 & EP300 & VHL)
#> VHL = HIF1A & !TP53
#> O2 = 0
#> 
#> Knocked-out and over-expressed genes:
#> O2 = 0
```

## Extract the Fundamental Boolean Network

``` r

tic("Total - Extract FBN")
tic("generateAllCombinationBinary")  
   initialStates <- generateAllCombinationBinary(hif1a.net$genes)
toc()
#> generateAllCombinationBinary: 0.009 sec elapsed
tic("genereateBoolNetTimeseries") 
   trainingseries <- genereateBoolNetTimeseries(hif1a.net,
                                             initialStates,43,
                                             type = "synchronous")
toc()
#> genereateBoolNetTimeseries: 0.302 sec elapsed
tic("generateFBMNetwork") 
FBNcellcyclenetwork <- generateFBMNetwork(timeseries_data = trainingseries,
                                 maxK = 4,
                                 max_deep_temporal = 1,
                                 useParallel = FALSE,
                                 verbose = TRUE)
#> INFO [2020-10-25 19:15:42] Enter generateFBMNetwork zone: 
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
#> INFO [2020-10-25 19:15:42] Run generateFBMNetwork with a single cube
#> INFO [2020-10-25 19:15:42] Enter constructFBNCube zone: 
#>         target_genes=6 genes and they are EP300, HIF1A, MDM2, TP53, VHL, O2,
#>         conditional_genes=6 genes and they are EP300, HIF1A, MDM2, TP53, VHL, O2,
#>         data_length=64,
#>         maxK=4, 
#>         temporal=1,
#>         useParallel=FALSE
#> INFO [2020-10-25 19:15:43] Leave constructFBNCube zone.
#> INFO [2020-10-25 19:15:43] Enter mineFBNNetwork zone: genes=, useParallel=FALSE
#> INFO [2020-10-25 19:15:43] Enter search_FBN_core zone: useParallel=FALSE
#> INFO [2020-10-25 19:15:43] Leave search_FBN_core zone
#> INFO [2020-10-25 19:15:43] Enter mineFBNNetworkWithCores zone
#> INFO [2020-10-25 19:15:43] Enter mineFBNNetworkStage2 zone
#> INFO [2020-10-25 19:15:43] Leave mineFBNNetworkStage2 zone
#> INFO [2020-10-25 19:15:43] Enter convertMinedResultToFBNNetwork zone
#> INFO [2020-10-25 19:15:43] Leave convertMinedResultToFBNNetwork zone
#> INFO [2020-10-25 19:15:43] Leave mineFBNNetworkWithCores zone
#> INFO [2020-10-25 19:15:43] Leave mineFBNNetwork zone
toc()
#> generateFBMNetwork: 1.105 sec elapsed
toc()
#> Total - Extract FBN: 1.419 sec elapsed

print(FBNcellcyclenetwork)
#> Fundamental Boolean Network with  6 genes
#> Genes involved:
#> EP300, HIF1A, MDM2, TP53, VHL, O2
#> 
#> Networks:
#> Multiple Transition Functions for EP300 with decay value = 1:
#> EP300_1_Activator: EP300 = TP53&!VHL (Confidence: 1, TimeStep: 1)
#> EP300_2_Activator: EP300 = !TP53&VHL (Confidence: 1, TimeStep: 1)
#> EP300_1_Inhibitor: EP300 = !TP53&!VHL (Confidence: 1, TimeStep: 1)
#> EP300_2_Inhibitor: EP300 = TP53&VHL (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for HIF1A with decay value = 1:
#> HIF1A_1_Activator: HIF1A = !MDM2&!VHL&!O2 (Confidence: 1, TimeStep: 1)
#> HIF1A_2_Activator: HIF1A = EP300&!TP53&!VHL&!O2 (Confidence: 1, TimeStep: 1)
#> HIF1A_1_Inhibitor: HIF1A = VHL (Confidence: 1, TimeStep: 1)
#> HIF1A_2_Inhibitor: HIF1A = O2 (Confidence: 1, TimeStep: 1)
#> HIF1A_3_Inhibitor: HIF1A = MDM2&TP53 (Confidence: 1, TimeStep: 1)
#> HIF1A_4_Inhibitor: HIF1A = !EP300&MDM2 (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for MDM2 with decay value = 1:
#> MDM2_1_Activator: MDM2 = TP53&!VHL (Confidence: 1, TimeStep: 1)
#> MDM2_1_Inhibitor: MDM2 = !TP53 (Confidence: 1, TimeStep: 1)
#> MDM2_2_Inhibitor: MDM2 = VHL (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for TP53 with decay value = 1:
#> TP53_1_Activator: TP53 = !MDM2 (Confidence: 1, TimeStep: 1)
#> TP53_1_Inhibitor: TP53 = MDM2&!VHL (Confidence: 1, TimeStep: 1)
#> TP53_2_Inhibitor: TP53 = !EP300&MDM2 (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for VHL with decay value = 1:
#> VHL_1_Activator: VHL = HIF1A&!TP53 (Confidence: 1, TimeStep: 1)
#> VHL_1_Inhibitor: VHL = TP53 (Confidence: 1, TimeStep: 1)
#> VHL_2_Inhibitor: VHL = !HIF1A (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for O2 with decay value = 1:
```

``` r
   resultfile <- reconstructTimeseries(FBNcellcyclenetwork,
                                    initialStates,
                                    type = "synchronous",
                                    maxTimepoints = 43,
                                    useParallel = FALSE)
   similarreport <- generateSimilaryReport(trainingseries,resultfile)
   print(paste("ErrorRate=",similarreport$ErrorRate,sep = "",collapse = ""))
#> [1] "ErrorRate=0.0181080426356589"
   print(paste("AccurateRate=",similarreport$AccurateRate,sep = "",collapse = ""))
#> [1] "AccurateRate=0.981891957364341"
   print(paste("MissMatchedRate=",similarreport$MissMatchedRate,sep = "",collapse = ""))
#> [1] "MissMatchedRate=0.0625"
   print(paste("PerfectMatchedRate=",similarreport$PerfectMatchedRate,sep = "",collapse = ""))
#> [1] "PerfectMatchedRate=0.9375"
   #get attractors
   genes <- rownames(trainingseries[[1]])
   attractor <- searchForAttractors(FBNcellcyclenetwork,initialStates,genes)
   print(attractor)
#> Discovered Attractors via Fundamental Boolean Model :
#> Genes are encoded in the following order::
#> EP300 HIF1A MDM2 TP53 VHL O2:
#> 
#> Attractor 1 is a complex attractor consisting of 5 state(s):
#> 
#> | --< - - - - - - |
#> v                 ^
#> 1 1 1 1 0 0       |
#> |                 |
#> 1 0 1 0 0 0       |
#> |                 |
#> 0 1 0 0 0 0       |
#> |                 |
#> 0 1 0 1 1 0       |
#> |                 |
#> 0 0 0 1 0 0       |
#> |                 |
#> v                 ^
#> | - - - - - - >-- |
   #display the dynamic trajectory of the attactor 2
   FBNNetwork.Graph.DrawAttractor(FBNcellcyclenetwork,attractor,1)
```

<!--html_preserve-->

<div id="htmlwidget-42c35c1c61ed2bc2a91c" class="visNetwork html-widget" style="width:100%;height:768px;">

</div>

<script type="application/json" data-for="htmlwidget-42c35c1c61ed2bc2a91c">{"x":{"nodes":{"id":["EP300_1","HIF1A_1","O2_1","O2_2","EP300_2","HIF1A_2","MDM2_2","TP53_2","VHL_2","MDM2_1","TP53_1","VHL_1","EP300_1_Activator_2","HIF1A_3_Inhibitor_2","MDM2_1_Activator_2","TP53_1_Inhibitor_2","VHL_1_Inhibitor_2","O2_3","EP300_3","HIF1A_3","MDM2_3","TP53_3","VHL_3","EP300_1_Inhibitor_3","HIF1A_2_Activator_3","MDM2_1_Inhibitor_3","TP53_1_Inhibitor_3","VHL_2_Inhibitor_3","O2_4","EP300_4","HIF1A_4","MDM2_4","TP53_4","VHL_4","EP300_1_Inhibitor_4","HIF1A_1_Activator_4","MDM2_1_Inhibitor_4","TP53_1_Activator_4","VHL_1_Activator_4","O2_5","EP300_5","HIF1A_5","MDM2_5","TP53_5","VHL_5","EP300_2_Inhibitor_5","HIF1A_1_Inhibitor_5","MDM2_2_Inhibitor_5","TP53_1_Activator_5","VHL_1_Inhibitor_5","EP300_6","HIF1A_6","MDM2_6","O2_6","TP53_6","VHL_6","TP53_1_Activator_6","VHL_1_Inhibitor_6"],"shape":["ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box"],"color":["lightblue","lightblue","pink","pink","lightblue","pink","lightblue","pink","pink","lightblue","lightblue","pink","LightGreen","orange","LightGreen","orange","orange","pink","pink","lightblue","pink","pink","pink","orange","LightGreen","orange","orange","orange","pink","pink","lightblue","pink","lightblue","lightblue","orange","LightGreen","orange","LightGreen","LightGreen","pink","pink","pink","pink","lightblue","pink","orange","orange","orange","LightGreen","orange","lightblue","lightblue","lightblue","pink","lightblue","pink","LightGreen","orange"],"type":["gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","gene","gene","gene","gene","gene","gene","TF","TF"],"value":[4,4,4,4,4,4,4,4,4,4,4,4,2,2,2,2,2,4,4,4,4,4,4,2,2,2,2,2,4,4,4,4,4,4,2,2,2,2,2,4,4,4,4,4,4,2,2,2,2,2,4,4,4,4,4,4,2,2],"label":["EP300_1","HIF1A_1","O2_1","O2_2","EP300_2","HIF1A_2","MDM2_2","TP53_2","VHL_2","MDM2_1","TP53_1","VHL_1","+, 1","-, 1","+, 1","-, 1","-, 1","O2_3","EP300_3","HIF1A_3","MDM2_3","TP53_3","VHL_3","-, 1","+, 1","-, 1","-, 1","-, 1","O2_4","EP300_4","HIF1A_4","MDM2_4","TP53_4","VHL_4","-, 1","+, 1","-, 1","+, 1","+, 1","O2_5","EP300_5","HIF1A_5","MDM2_5","TP53_5","VHL_5","-, 1","-, 1","-, 1","+, 1","-, 1","EP300_6","HIF1A_6","MDM2_6","O2_6","TP53_6","VHL_6","+, 1","-, 1"],"shadow":[true,true,true,true,true,true,true,true,true,true,true,true,false,false,false,false,false,true,true,true,true,true,true,false,false,false,false,false,true,true,true,true,true,true,false,false,false,false,false,true,true,true,true,true,true,false,false,false,false,false,true,true,true,true,true,true,false,false],"group":["Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Activate Function","Inhibit Function","Activate Function","Inhibit Function","Inhibit Function","Gene","Gene","Gene","Gene","Gene","Gene","Inhibit Function","Activate Function","Inhibit Function","Inhibit Function","Inhibit Function","Gene","Gene","Gene","Gene","Gene","Gene","Inhibit Function","Activate Function","Inhibit Function","Activate Function","Activate Function","Gene","Gene","Gene","Gene","Gene","Gene","Inhibit Function","Inhibit Function","Inhibit Function","Activate Function","Inhibit Function","Gene","Gene","Gene","Gene","Gene","Gene","Activate Function","Inhibit Function"],"title":["EP300","HIF1A","O2","O2","EP300","HIF1A","MDM2","TP53","VHL","MDM2","TP53","VHL","TP53&!VHL","MDM2&TP53","TP53&!VHL","MDM2&!VHL","TP53","O2","EP300","HIF1A","MDM2","TP53","VHL","!TP53&!VHL","EP300&!TP53&!VHL&!O2","!TP53","MDM2&!VHL","!HIF1A","O2","EP300","HIF1A","MDM2","TP53","VHL","!TP53&!VHL","!MDM2&!VHL&!O2","!TP53","!MDM2","HIF1A&!TP53","O2","EP300","HIF1A","MDM2","TP53","VHL","TP53&VHL","VHL","VHL","!MDM2","TP53","EP300","HIF1A","MDM2","O2","TP53","VHL","!MDM2","TP53"],"level":[1,1,1,3,3,3,3,3,3,1,1,1,2,2,2,2,2,5,5,5,5,5,5,4,4,4,4,4,7,7,7,7,7,7,6,6,6,6,6,9,9,9,9,9,9,8,8,8,8,8,11,11,11,11,11,11,10,10]},"edges":{"from":["EP300_1_Activator_2","TP53_1","VHL_1","HIF1A_3_Inhibitor_2","MDM2_1","TP53_1","MDM2_1_Activator_2","TP53_1","VHL_1","TP53_1_Inhibitor_2","MDM2_1","VHL_1","VHL_1_Inhibitor_2","TP53_1","EP300_1_Inhibitor_3","TP53_2","VHL_2","HIF1A_2_Activator_3","EP300_2","TP53_2","VHL_2","O2_2","MDM2_1_Inhibitor_3","TP53_2","TP53_1_Inhibitor_3","MDM2_2","VHL_2","VHL_2_Inhibitor_3","HIF1A_2","EP300_1_Inhibitor_4","TP53_3","VHL_3","HIF1A_1_Activator_4","MDM2_3","VHL_3","O2_3","MDM2_1_Inhibitor_4","TP53_3","TP53_1_Activator_4","MDM2_3","VHL_1_Activator_4","HIF1A_3","TP53_3","EP300_2_Inhibitor_5","TP53_4","VHL_4","HIF1A_1_Inhibitor_5","VHL_4","MDM2_2_Inhibitor_5","VHL_4","TP53_1_Activator_5","MDM2_4","VHL_1_Inhibitor_5","TP53_4","TP53_1_Activator_6","MDM2_5","VHL_1_Inhibitor_6","TP53_5"],"to":["EP300_2","EP300_1_Activator_2","EP300_1_Activator_2","HIF1A_2","HIF1A_3_Inhibitor_2","HIF1A_3_Inhibitor_2","MDM2_2","MDM2_1_Activator_2","MDM2_1_Activator_2","TP53_2","TP53_1_Inhibitor_2","TP53_1_Inhibitor_2","VHL_2","VHL_1_Inhibitor_2","EP300_3","EP300_1_Inhibitor_3","EP300_1_Inhibitor_3","HIF1A_3","HIF1A_2_Activator_3","HIF1A_2_Activator_3","HIF1A_2_Activator_3","HIF1A_2_Activator_3","MDM2_3","MDM2_1_Inhibitor_3","TP53_3","TP53_1_Inhibitor_3","TP53_1_Inhibitor_3","VHL_3","VHL_2_Inhibitor_3","EP300_4","EP300_1_Inhibitor_4","EP300_1_Inhibitor_4","HIF1A_4","HIF1A_1_Activator_4","HIF1A_1_Activator_4","HIF1A_1_Activator_4","MDM2_4","MDM2_1_Inhibitor_4","TP53_4","TP53_1_Activator_4","VHL_4","VHL_1_Activator_4","VHL_1_Activator_4","EP300_5","EP300_2_Inhibitor_5","EP300_2_Inhibitor_5","HIF1A_5","HIF1A_1_Inhibitor_5","MDM2_5","MDM2_2_Inhibitor_5","TP53_5","TP53_1_Activator_5","VHL_5","VHL_1_Inhibitor_5","TP53_6","TP53_1_Activator_6","VHL_6","VHL_1_Inhibitor_6"],"support":[0.39695,0.39695,0.39695,0.1994,0.1994,0.1994,0.39695,0.39695,0.39695,0.39286,0.39286,0.39286,0.59673,0.59673,0.39472,0.39472,0.39472,0.19606,0.19606,0.19606,0.19606,0.19606,0.40327,0.40327,0.39286,0.39286,0.39286,0.40885,0.40885,0.39472,0.39472,0.39472,0.39583,0.39583,0.39583,0.39583,0.40327,0.40327,0.60119,0.60119,0.19717,0.19717,0.19717,0.19978,0.19978,0.19978,0.20833,0.20833,0.20833,0.20833,0.60119,0.60119,0.59673,0.59673,0.60119,0.60119,0.59673,0.59673],"targetNode":["EP300","EP300","EP300","HIF1A","HIF1A","HIF1A","MDM2","MDM2","MDM2","TP53","TP53","TP53","VHL","VHL","EP300","EP300","EP300","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","MDM2","MDM2","TP53","TP53","TP53","VHL","VHL","EP300","EP300","EP300","HIF1A","HIF1A","HIF1A","HIF1A","MDM2","MDM2","TP53","TP53","VHL","VHL","VHL","EP300","EP300","EP300","HIF1A","HIF1A","MDM2","MDM2","TP53","TP53","VHL","VHL","TP53","TP53","VHL","VHL"],"title":["EP300_1_Activator","EP300_1_Activator","EP300_1_Activator","HIF1A_3_Inhibitor","HIF1A_3_Inhibitor","HIF1A_3_Inhibitor","MDM2_1_Activator","MDM2_1_Activator","MDM2_1_Activator","TP53_1_Inhibitor","TP53_1_Inhibitor","TP53_1_Inhibitor","VHL_1_Inhibitor","VHL_1_Inhibitor","EP300_1_Inhibitor","EP300_1_Inhibitor","EP300_1_Inhibitor","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","MDM2_1_Inhibitor","MDM2_1_Inhibitor","TP53_1_Inhibitor","TP53_1_Inhibitor","TP53_1_Inhibitor","VHL_2_Inhibitor","VHL_2_Inhibitor","EP300_1_Inhibitor","EP300_1_Inhibitor","EP300_1_Inhibitor","HIF1A_1_Activator","HIF1A_1_Activator","HIF1A_1_Activator","HIF1A_1_Activator","MDM2_1_Inhibitor","MDM2_1_Inhibitor","TP53_1_Activator","TP53_1_Activator","VHL_1_Activator","VHL_1_Activator","VHL_1_Activator","EP300_2_Inhibitor","EP300_2_Inhibitor","EP300_2_Inhibitor","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","MDM2_2_Inhibitor","MDM2_2_Inhibitor","TP53_1_Activator","TP53_1_Activator","VHL_1_Inhibitor","VHL_1_Inhibitor","TP53_1_Activator","TP53_1_Activator","VHL_1_Inhibitor","VHL_1_Inhibitor"],"arrowtail":["none","none","tee","none","none","none","none","none","tee","none","none","tee","none","none","none","tee","tee","none","none","tee","tee","tee","none","tee","none","none","tee","none","tee","none","tee","tee","none","tee","tee","tee","none","tee","none","tee","none","none","tee","none","none","none","none","none","none","none","none","tee","none","none","none","tee","none","none"],"arrowhead":["open","none","none","odot","none","none","open","none","none","odot","none","none","odot","none","odot","none","none","open","none","none","none","none","odot","none","odot","none","none","odot","none","odot","none","none","open","none","none","none","odot","none","open","none","open","none","none","odot","none","none","odot","none","odot","none","open","none","odot","none","open","none","odot","none"],"color":["darkblue","green","red","darkred","green","green","darkblue","green","red","darkred","green","red","darkred","green","darkred","red","red","darkblue","green","red","red","red","darkred","red","darkred","green","red","darkred","red","darkred","red","red","darkblue","red","red","red","darkred","red","darkblue","red","darkblue","green","red","darkred","green","green","darkred","green","darkred","green","darkblue","red","darkred","green","darkblue","red","darkred","green"],"lty":["solid","solid","solid","longdash","longdash","longdash","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","solid","longdash","longdash","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","longdash","longdash","solid","solid","longdash","longdash"],"arrow.mode":[2,0,0,2,0,0,2,0,0,2,0,0,2,0,2,0,0,2,0,0,0,0,2,0,2,0,0,2,0,2,0,0,2,0,0,0,2,0,2,0,2,0,0,2,0,0,2,0,2,0,2,0,2,0,2,0,2,0],"arrow.size":[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2],"type":["TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF"],"arrows":["to","none","none","to","none","none","to","none","none","to","none","none","to","none","to","none","none","to","none","none","none","none","to","none","to","none","none","to","none","to","none","none","to","none","none","none","to","none","to","none","to","none","none","to","none","none","to","none","to","none","to","none","to","none","to","none","to","none"],"label":["","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"dashes":[false,false,true,true,false,false,false,false,true,true,false,true,true,false,true,true,true,false,false,true,true,true,true,true,true,false,true,true,true,true,true,true,false,true,true,true,true,true,false,true,false,false,true,true,false,false,true,false,true,false,false,true,true,false,false,true,true,false],"shadow":[true,false,false,true,false,false,true,false,false,true,false,false,true,false,true,false,false,true,false,false,false,false,true,false,true,false,false,true,false,true,false,false,true,false,false,false,true,false,true,false,true,false,false,true,false,false,true,false,true,false,true,false,true,false,true,false,true,false]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot"},"manipulation":{"enabled":false},"interaction":{"dragNodes":true,"dragView":true,"zoomView":true},"edges":{"arrows":"from"},"layout":{"hierarchical":{"enabled":true,"levelSeparation":200,"direction":"LR"}}},"groups":["Gene","Activate Function","Inhibit Function"],"width":null,"height":null,"idselection":{"enabled":false},"byselection":{"enabled":false},"main":{"text":"Dynamic Fundamental Boolean Networks from the time point of 1 to 6","style":"font-family:Georgia, Times New Roman, Times, serif;font-weight:bold;font-size:20px;text-align:center;"},"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","legend":{"width":0.1,"useGroups":false,"position":"left","ncol":1,"stepX":100,"stepY":100,"zoom":true,"edges":{"color":["darkblue","darkred","grey","green","red"],"label":["activate","inhibit","decay","activated input","deactivated input"],"arrows":["to","to","to","none","none"],"dashes":[false,true,true,false,true],"shadow":[true,true,true,false,false]},"edgesToDataframe":true,"nodes":{"label":["Activated Gene","Inhibited Gene","Activate Function (+)","Inhibit Function (-)"],"shape":["ellipse","ellipse","box","box"],"color":["lightblue","pink","lightgreen","orange"],"shadow":[true,true,false,false]},"nodesToDataframe":true},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);"},"evals":[],"jsHooks":[]}</script>

<!--/html_preserve-->

``` r

FBNNet::FBNNetwork.Graph(FBNcellcyclenetwork)
```

<!--html_preserve-->

<div id="htmlwidget-d186e7ce427e59e502fc" class="visNetwork html-widget" style="width:100%;height:768px;">

</div>

<script type="application/json" data-for="htmlwidget-d186e7ce427e59e502fc">{"x":{"nodes":{"id":["EP300","HIF1A","MDM2","TP53","VHL","O2","EP300_1_Activator","EP300_2_Activator","EP300_1_Inhibitor","EP300_2_Inhibitor","HIF1A_1_Activator","HIF1A_2_Activator","HIF1A_1_Inhibitor","HIF1A_2_Inhibitor","HIF1A_3_Inhibitor","HIF1A_4_Inhibitor","MDM2_1_Activator","MDM2_1_Inhibitor","MDM2_2_Inhibitor","TP53_1_Activator","TP53_1_Inhibitor","TP53_2_Inhibitor","VHL_1_Activator","VHL_1_Inhibitor","VHL_2_Inhibitor"],"shape":["ellipse","ellipse","ellipse","ellipse","ellipse","ellipse","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box","box"],"color":["lightblue","lightblue","lightblue","lightblue","lightblue","lightblue","LightGreen","LightGreen","orange","orange","LightGreen","LightGreen","orange","orange","orange","orange","LightGreen","orange","orange","LightGreen","orange","orange","LightGreen","orange","orange"],"type":["gene","gene","gene","gene","gene","gene","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF","TF"],"value":[4,4,4,4,4,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],"label":["EP300","HIF1A","MDM2","TP53","VHL","O2","+, 1","+, 1","-, 1","-, 1","+, 1","+, 1","-, 1","-, 1","-, 1","-, 1","+, 1","-, 1","-, 1","+, 1","-, 1","-, 1","+, 1","-, 1","-, 1"],"shadow":[true,true,true,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false],"group":["Gene","Gene","Gene","Gene","Gene","Gene","Activate Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Activate Function","Inhibit Function","Inhibit Function","Inhibit Function","Inhibit Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Inhibit Function","Inhibit Function","Activate Function","Inhibit Function","Inhibit Function"],"title":["E1A binding protein p300(EP300)","hypoxia inducible factor 1 alpha subunit(HIF1A)","MDM2 proto-oncogene(MDM2)","tumor protein p53(TP53)","von Hippel-Lindau tumor suppressor(VHL)",null,"TP53&!VHL","!TP53&VHL","!TP53&!VHL","TP53&VHL","!MDM2&!VHL&!O2","EP300&!TP53&!VHL&!O2","VHL","O2","MDM2&TP53","!EP300&MDM2","TP53&!VHL","!TP53","VHL","!MDM2","MDM2&!VHL","!EP300&MDM2","HIF1A&!TP53","TP53","!HIF1A"]},"edges":{"from":["EP300_1_Activator","TP53","VHL","EP300_2_Activator","TP53","VHL","EP300_1_Inhibitor","TP53","VHL","EP300_2_Inhibitor","TP53","VHL","HIF1A_1_Activator","MDM2","VHL","O2","HIF1A_2_Activator","EP300","TP53","VHL","O2","HIF1A_1_Inhibitor","VHL","HIF1A_2_Inhibitor","O2","HIF1A_3_Inhibitor","MDM2","TP53","HIF1A_4_Inhibitor","EP300","MDM2","MDM2_1_Activator","TP53","VHL","MDM2_1_Inhibitor","TP53","MDM2_2_Inhibitor","VHL","TP53_1_Activator","MDM2","TP53_1_Inhibitor","MDM2","VHL","TP53_2_Inhibitor","EP300","MDM2","VHL_1_Activator","HIF1A","TP53","VHL_1_Inhibitor","TP53","VHL_2_Inhibitor","HIF1A"],"to":["EP300","EP300_1_Activator","EP300_1_Activator","EP300","EP300_2_Activator","EP300_2_Activator","EP300","EP300_1_Inhibitor","EP300_1_Inhibitor","EP300","EP300_2_Inhibitor","EP300_2_Inhibitor","HIF1A","HIF1A_1_Activator","HIF1A_1_Activator","HIF1A_1_Activator","HIF1A","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A","HIF1A_1_Inhibitor","HIF1A","HIF1A_2_Inhibitor","HIF1A","HIF1A_3_Inhibitor","HIF1A_3_Inhibitor","HIF1A","HIF1A_4_Inhibitor","HIF1A_4_Inhibitor","MDM2","MDM2_1_Activator","MDM2_1_Activator","MDM2","MDM2_1_Inhibitor","MDM2","MDM2_2_Inhibitor","TP53","TP53_1_Activator","TP53","TP53_1_Inhibitor","TP53_1_Inhibitor","TP53","TP53_2_Inhibitor","TP53_2_Inhibitor","VHL","VHL_1_Activator","VHL_1_Activator","VHL","VHL_1_Inhibitor","VHL","VHL_2_Inhibitor"],"support":[0.39695,0.39695,0.39695,0.00856,0.00856,0.00856,0.39472,0.39472,0.39472,0.19978,0.19978,0.19978,0.39583,0.39583,0.39583,0.39583,0.19606,0.19606,0.19606,0.19606,0.19606,0.20833,0.20833,0.0119,0.0119,0.1994,0.1994,0.1994,0.00595,0.00595,0.00595,0.39695,0.39695,0.39695,0.40327,0.40327,0.20833,0.20833,0.60119,0.60119,0.39286,0.39286,0.39286,0.00595,0.00595,0.00595,0.19717,0.19717,0.19717,0.59673,0.59673,0.40885,0.40885],"targetNode":["EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","EP300","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","HIF1A","MDM2","MDM2","MDM2","MDM2","MDM2","MDM2","MDM2","TP53","TP53","TP53","TP53","TP53","TP53","TP53","TP53","VHL","VHL","VHL","VHL","VHL","VHL","VHL"],"title":["EP300_1_Activator","EP300_1_Activator","EP300_1_Activator","EP300_2_Activator","EP300_2_Activator","EP300_2_Activator","EP300_1_Inhibitor","EP300_1_Inhibitor","EP300_1_Inhibitor","EP300_2_Inhibitor","EP300_2_Inhibitor","EP300_2_Inhibitor","HIF1A_1_Activator","HIF1A_1_Activator","HIF1A_1_Activator","HIF1A_1_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_2_Activator","HIF1A_1_Inhibitor","HIF1A_1_Inhibitor","HIF1A_2_Inhibitor","HIF1A_2_Inhibitor","HIF1A_3_Inhibitor","HIF1A_3_Inhibitor","HIF1A_3_Inhibitor","HIF1A_4_Inhibitor","HIF1A_4_Inhibitor","HIF1A_4_Inhibitor","MDM2_1_Activator","MDM2_1_Activator","MDM2_1_Activator","MDM2_1_Inhibitor","MDM2_1_Inhibitor","MDM2_2_Inhibitor","MDM2_2_Inhibitor","TP53_1_Activator","TP53_1_Activator","TP53_1_Inhibitor","TP53_1_Inhibitor","TP53_1_Inhibitor","TP53_2_Inhibitor","TP53_2_Inhibitor","TP53_2_Inhibitor","VHL_1_Activator","VHL_1_Activator","VHL_1_Activator","VHL_1_Inhibitor","VHL_1_Inhibitor","VHL_2_Inhibitor","VHL_2_Inhibitor"],"arrowtail":["none","none","tee","none","tee","none","none","tee","tee","none","none","none","none","tee","tee","tee","none","none","tee","tee","tee","none","none","none","none","none","none","none","none","tee","none","none","none","tee","none","tee","none","none","none","tee","none","none","tee","none","tee","none","none","none","tee","none","none","none","tee"],"arrowhead":["open","none","none","open","none","none","odot","none","none","odot","none","none","open","none","none","none","open","none","none","none","none","odot","none","odot","none","odot","none","none","odot","none","none","open","none","none","odot","none","odot","none","open","none","odot","none","none","odot","none","none","open","none","none","odot","none","odot","none"],"color":["darkblue","green","red","darkblue","red","green","darkred","red","red","darkred","green","green","darkblue","red","red","red","darkblue","green","red","red","red","darkred","green","darkred","green","darkred","green","green","darkred","red","green","darkblue","green","red","darkred","red","darkred","green","darkblue","red","darkred","green","red","darkred","red","green","darkblue","green","red","darkred","green","darkred","red"],"lty":["solid","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","solid","solid","solid","solid","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","longdash","longdash","longdash","longdash","solid","solid","longdash","longdash","longdash","longdash","longdash","longdash","solid","solid","solid","longdash","longdash","longdash","longdash"],"arrow.mode":[2,0,0,2,0,0,2,0,0,2,0,0,2,0,0,0,2,0,0,0,0,2,0,2,0,2,0,0,2,0,0,2,0,0,2,0,2,0,2,0,2,0,0,2,0,0,2,0,0,2,0,2,0],"arrow.size":[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2],"type":["TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","Gene_to_TF","TF_to_Gene","Gene_to_TF","TF_to_Gene","Gene_to_TF"],"arrows":["to","none","none","to","none","none","to","none","none","to","none","none","to","none","none","none","to","none","none","none","none","to","none","to","none","to","none","none","to","none","none","to","none","none","to","none","to","none","to","none","to","none","none","to","none","none","to","none","none","to","none","to","none"],"label":["","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""],"dashes":[false,false,true,false,true,false,true,true,true,true,false,false,false,true,true,true,false,false,true,true,true,true,false,true,false,true,false,false,true,true,false,false,false,true,true,true,true,false,false,true,true,false,true,true,true,false,false,false,true,true,false,true,true],"shadow":[true,false,false,true,false,false,true,false,false,true,false,false,true,false,false,false,true,false,false,false,false,true,false,true,false,true,false,false,true,false,false,true,false,false,true,false,true,false,true,false,true,false,false,true,false,false,true,false,false,true,false,true,false]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot"},"manipulation":{"enabled":false},"interaction":{"dragNodes":true,"dragView":true,"zoomView":true}},"groups":["Gene","Activate Function","Inhibit Function"],"width":null,"height":null,"idselection":{"enabled":false},"byselection":{"enabled":false},"main":{"text":"Fundamental Boolean Networks","style":"font-family:Georgia, Times New Roman, Times, serif;font-weight:bold;font-size:20px;text-align:center;"},"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","legend":{"width":0.2,"useGroups":false,"position":"left","ncol":1,"stepX":100,"stepY":100,"zoom":true,"edges":{"color":["darkblue","darkred","green","red"],"label":["activate","inhibit","activated input","deactivated input"],"arrows":["to","to","none","none"],"dashes":[false,true,false,true],"shadow":[true,true,false,false]},"edgesToDataframe":true,"nodes":{"label":["Gene","Activate Function (+, Timestep)","Inhibit Function (-, Timestep)"],"shape":["ellipse","box","box"],"color":["lightblue","lightgreen","orange"],"shadow":[true,false,false]},"nodesToDataframe":true},"tooltipStay":300,"tooltipStyle":"position: fixed;visibility:hidden;padding: 5px;white-space: nowrap;font-family: verdana;font-size:14px;font-color:#000000;background-color: #f5f4ed;-moz-border-radius: 3px;-webkit-border-radius: 3px;border-radius: 3px;border: 1px solid #808074;box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);"},"evals":[],"jsHooks":[]}</script>

<!--/html_preserve-->
