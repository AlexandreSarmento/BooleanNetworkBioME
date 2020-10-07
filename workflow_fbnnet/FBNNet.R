#' ---
#' title: "Fundamental Boolean Network"
#' pagetitle: "Fundamental Boolean Network"
#' author: "Leshi Chen"
#' affiliation: Lincoln University, New Zealand
#' date: "`r Sys.Date()`"
#' url: https://www.frontiersin.org/articles/10.3389/fphys.2018.01328/full
#' output: 
#'   pdf_document: default
#'   html_document: 
#'     default
#'   github_document: 
#'     df_print: paged
#'     html_preview: FALSE
#'     keep_html: TRUE
#' knit: (function(inputFile, encoding) {
#'   rmarkdown::render(inputFile, encoding = encoding, output_format = "all") })    
#' ---
#' 
#' 
#' 

#' 
#' 
#' # Installing and Loading Libraries            
#' 
## ----message=FALSE, warning=FALSE, paged.print=FALSE--------------------------
packages_bioconductor = c("limma")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from Bioconductor and loaded
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

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

rm(package.check, packages_bioconductor, packages_cran)

#' 
#' 
#' ## Introduction
#' Fundamental Boolean modelling(FBM) has been proposed to draw insights into gene activation, inhibition, and protein decay. This novel Boolean model provides an intuitive definition of activation and inhibition pathways and includes mechanisms to handle protein decay issues. To prove the concept of the novel model, we implemented a platform using R language, called FBNNet. Our experimental results show that the proposed FBM could explicitly display the internal
#' connections of the mammalian cell cycle between genes separated into the connection
#' types of activation, inhibition and protein decay.
#' 
#' https://clsdavid.github.io/FBNNet2_public/
#' 
#' 
#' ## Experiment
#' The experiments conducted and described here intend to prove the concept of the new Boolean Model, i.e., the FBM. To verify the results, we apply the general processes described in Figure 1 as a benchmark to compare the results generated via BoolNet, with these consequently reconstructed from the new R package, FBNNet.
#' 
#' ## Data
#' 
## ---- out.width = "100%"------------------------------------------------------

hif1a.net<-loadNetwork("boolean_network_HIFaxis.bn")
hif1a.net

#' 
#' ## Extract the Fundamental Boolean Network
## ---- out.width = "100%"------------------------------------------------------

tic("Total - Extract FBN")
tic("generateAllCombinationBinary")  
   initialStates <- generateAllCombinationBinary(hif1a.net$genes)
toc()
tic("genereateBoolNetTimeseries") 
   trainingseries <- genereateBoolNetTimeseries(hif1a.net,
                                             initialStates,43,
                                             type = "synchronous")
toc()
tic("generateFBMNetwork") 
FBNcellcyclenetwork <- generateFBMNetwork(timeseries_data = trainingseries,
                                 maxK = 4,
                                 max_deep_temporal = 1,
                                 useParallel = FALSE,
                                 verbose = TRUE)
toc()
toc()

   print(FBNcellcyclenetwork)

#' 
#' 
## ---- fig.width = 7.1, fig.height = 8-----------------------------------------
   FBNNet::FBNNetwork.Graph(FBNcellcyclenetwork)

#' 
#' 
## ---- out.width = "100%"------------------------------------------------------
   resultfile <- reconstructTimeseries(FBNcellcyclenetwork,
                                    initialStates,
                                    type = "synchronous",
                                    maxTimepoints = 10,
                                    useParallel = FALSE)
   similarreport <- generateSimilaryReport(trainingseries,resultfile)
   print(paste("ErrorRate=",similarreport$ErrorRate,sep = "",collapse = ""))
   print(paste("AccurateRate=",similarreport$AccurateRate,sep = "",collapse = ""))
   print(paste("MissMatchedRate=",similarreport$MissMatchedRate,sep = "",collapse = ""))
   print(paste("PerfectMatchedRate=",similarreport$PerfectMatchedRate,sep = "",collapse = ""))
   #get attractors
   genes <- rownames(trainingseries[[1]])
   attractor <- searchForAttractors(FBNcellcyclenetwork,initialStates,genes)
   print(attractor)
   #display the dynamic trajectory of the attactor 2
   FBNNetwork.Graph.DrawAttractor(FBNcellcyclenetwork,attractor,2)

