#' ---
#' title: "BoolNet Inference (E-GEOD-18494, GSE47533 and GSE41491)"
#' output:
#'   pdf_document: default
#'   github_document: 
#'     df_print: paged
#'     html_preview: FALSE
#'     keep_html: TRUE
#' knit: (function(inputFile, encoding) {
#'   rmarkdown::render(inputFile, encoding = encoding, output_format = "all") })    
#' ---

#' 
#' 
## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
#packages_cran = c("igraph", "BoolNet", "BiocManager", "tidyverse", "fs", "ff","RSQLite")
packages_cran = c("igraph", "BoolNet", "BiocManager", "tidyverse", "ggraph")

# Install and load packages
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

if(!require('cowplot')) {install.packages('cowplot')} # provides various features that help with creating publication-quality figures
if(!require('grid')){install.packages('grid')} # 
if(!require('gridExtra')) {install.packages('gridExtra')} #  Provides a number of user-level functions to work with "grid" graphics, notably to arrange multiple grid-based plots on a page, and draw tables.


# For oligo and ArrayExpress First install:
#install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz',repos=NULL)
# 
# if(!require("hgu133plus2hsentrezg.db")){ 
#   install.packages('http://mbni.org/customcdf/13.0.0/entrezg.download/hgu133plus2hsentrezg.db_13.0.0.tar.gz',repos=NULL) 
# }
# 
#packages_bioconductor = c("Biobase", "GEOquery", "affyPLM", "ArrayExpress", "illuminaHumanv3.db", "hgu133plus2.db")
packages_bioconductor = c()

# Install and load packages
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(package.check, packages_bioconductor, packages_cran)


#' 
#' # Integrated analysis of microRNA and mRNA expression and association with HIF binding in MCF-7 cells under hypoxia (GSE47533)
#' 
#' Camps C, Saini HK, Mole DR, Choudhry H et al. Integrated analysis of microRNA and mRNA expression and association with HIF binding reveals the complexity of microRNA expression regulation under hypoxia. Mol Cancer 2014 Feb 11;13:28. PMID: 24517586
#' 
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47533
#' 
#' 
## ----include=FALSE------------------------------------------------------------
load("../data/data.GSE47533.Rdata")

#' 
## ----include=FALSE------------------------------------------------------------
# Selected genes from HIF Axis
hif.symbols <- c("TP53", "HIF1A", "EP300", "MDM2", "VHL")

hif.probes <- anno.GSE47533$probes[anno.GSE47533$symbol %in% hif.symbols]

# Select the probes and genes
expr.GSE47533.hif <- data.frame(expr.GSE47533) %>% 
  rownames_to_column('probes') %>% 
  filter(probes %in% hif.probes) %>% 
  merge(anno.GSE47533[anno.GSE47533$symbol %in% hif.symbols, c("probes","symbol")], by = "probes") %>% 
  #distinct(symbol, .keep_all = TRUE) %>% # Take the first one
  dplyr::select(!(probes)) %>% 
  arrange(symbol)
  

#' 
## ----include=FALSE------------------------------------------------------------
# Function to binarize according an consensus mean of probes, add the O2 state and rename columns 

binNet <- function(b){
  
  cols <- data.GSE47533$codes %in% names(b)
  
  binarizeTimeSeries(b[,-1], method="kmeans")$binarizedMeasurements  %>% 
  as.data.frame(.)  %>% 
  aggregate(., list(symbol = b$symbol), mean) %>%  # mean of binarized probes
  mutate_at(vars(-symbol), funs(ifelse(. >= 0.5, 1, 0))) %>%  # consensus with a bies to 1 (>= 0.5)
  rbind(., c("O2", 1,0,0,0)) %>% 
    rename_at(vars(data.GSE47533$codes[cols] ),
            ~paste0(substr(data.GSE47533$condition[cols],1,2),".",
                    data.GSE47533$time[cols],".",
                    substr(data.GSE47533$cell_line[cols],1,2), ".",
                    data.GSE47533$rep[cols])) %>% 
  column_to_rownames("symbol")
  
}

breast_MCF7.1 <- 
expr.GSE47533.hif %>% 
  dplyr::select(c("symbol", data.GSE47533$codes[data.GSE47533$rep == 1]))  %>% 
  binNet(.) 

breast_MCF7.2 <- 
expr.GSE47533.hif %>% 
  dplyr::select(c("symbol", data.GSE47533$codes[data.GSE47533$rep == 2]))  %>% 
  binNet(.) 

breast_MCF7.3 <- 
expr.GSE47533.hif %>% 
  dplyr::select(c("symbol", data.GSE47533$codes[data.GSE47533$rep == 3]))  %>% 
  binNet(.) 


#' 
#' # Expression profiling of hypoxic HepG2 hepatoma, U87 glioma, and MDA-MB231 breast cancer cells: time course (E-GEOD-18494)
#' 
#' https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-18494/   
#' 
#' 
## ----include=FALSE------------------------------------------------------------
load("../data/data.EGEOD18494.Rdata")

#' 
#' 
## ----include=FALSE------------------------------------------------------------

# Selecting genes from HIF Axis
hif.symbols <- c("TP53", "HIF1A", "EP300", "MDM2", "VHL")

hif.probes <- anno.EGEOD18494$probes[anno.EGEOD18494$symbol %in% hif.symbols]

# Select the probes and genes
expr.EGEOD18494.hif <- as.data.frame(expr.EGEOD18494) %>% 
  rownames_to_column('probes') %>% 
  filter(probes %in% hif.probes) %>% 
  merge(anno.EGEOD18494[anno.EGEOD18494$symbol %in% hif.symbols, c("probes","symbol")], by = "probes") %>% 
  #distinct(symbol, .keep_all = TRUE) %>% # Take the first one
  dplyr::select(!(probes)) 
  

#' 
#' 
## ----include=FALSE------------------------------------------------------------
# Function to binarize according an consensus mean of probes, add the O2 state and rename columns 

binNet <- function(b){
  
  cols <- data.EGEOD18494$codes %in% names(b)
  
  binarizeTimeSeries(b[,-1], method="kmeans")$binarizedMeasurements  %>% 
  as.data.frame(.)  %>% 
  aggregate(., list(symbol = b$symbol), mean) %>%  # mean of binarized probes
  mutate_at(vars(-symbol), funs(ifelse(. >= 0.5, 1, 0))) %>%  # consensus with a bies to 1 (>= 0.5)
  rbind(., c("O2", 1,0,0,0)) %>% 
    rename_at(vars(data.EGEOD18494$codes[cols] ),
            ~paste0(substr(data.EGEOD18494$condition[cols],1,2),".",
                    data.EGEOD18494$time[cols],".",
                    substr(data.EGEOD18494$cell_line[cols],1,2), ".",
                    data.EGEOD18494$rep[cols])) %>% 
  column_to_rownames("symbol")
  
}

cellline.rep1 <- (data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &  data.EGEOD18494$rep == 1)
cellline.rep2 <- (data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &  data.EGEOD18494$rep == 2)
cellline.rep3 <- (data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &  data.EGEOD18494$rep == 3)

breast_MDA.1 <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep1])) %>% 
  binNet(.) 

breast_MDA.2 <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep2])) %>% 
  binNet(.) 

breast_MDA.3 <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep3])) %>% 
  binNet(.) 

#' 
#' 
#' # Hypoxia transcriptomic time-series data in three different cancer cell lines (GSE41491)
#' 
#' Tumour hypoxia exhibits a highly dynamic spatial and temporal distribution and is associated with increased malignancy and poor prognosis.
#' 
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41491
#' 
#' Exponentially growing prostate (DU145), colon (HT29) and breast (MCF7) carcinoma cells were seeded on glass dishes in McCoy, DMEM or RPMI media, respectively with 10% FCS.
#' 
#' 
## ----include=FALSE------------------------------------------------------------
load("../data/data.GSE41491.Rdata")

#' 
## ----include=FALSE------------------------------------------------------------
# Selected genes from HIF Axis
hif.symbols <- c("TP53", "HIF1A", "EP300", "MDM2", "VBP1") 

hif.probes <- unique(anno.GSE41491$probes[anno.GSE41491$symbol %in% hif.symbols])

# Select the probes and genes
expr.GSE41491.hif <- data.frame(expr.GSE41491) %>% 
  rownames_to_column('probes') %>% 
  filter(probes %in% hif.probes) %>% 
  merge(anno.GSE41491[anno.GSE41491$symbol %in% hif.symbols, c("probes","symbol")], by = "probes") %>% 
  #distinct(symbol, .keep_all = TRUE) %>% # Take the first one
  dplyr::select(!(probes)) %>% 
  arrange(symbol)

expr.GSE41491.hif$symbol[expr.GSE41491.hif$symbol == "VBP1"] <- "VHL"
  

#' 
## ----include=FALSE------------------------------------------------------------
# Function to binarize according an consensus mean of probes, add the O2 state and rename columns 

binNet <- function(b){
  
  cols <- data.GSE41491$codes %in% names(b)
  
  binarizeTimeSeries(b[,-1], method="kmeans")$binarizedMeasurements  %>% 
  as.data.frame(.)  %>% 
  aggregate(., list(symbol = b$symbol), mean) %>%  # mean of binarized probes
  mutate_at(vars(-symbol), funs(ifelse(. >= 0.5, 1, 0))) %>%  # consensus with a bies to 1 (>= 0.5)
  rbind(., c("O2", 1,0,0,0,0,0,0,0)) %>% 
    rename_at(vars(data.GSE41491$codes[cols] ),
            ~paste0(substr(data.GSE41491$condition[cols],1,2),".",
                    data.GSE41491$time[cols],".",
                    substr(data.GSE41491$cell_line[cols],1,2))) %>% 
  column_to_rownames("symbol")
  
}

breast_MCF7 <- 
expr.GSE41491.hif %>% 
  dplyr::select(c("symbol", data.GSE41491$codes[data.GSE41491$cell_line == "MCF7"]))  %>% 
  binNet(.) 


#' 
#' 
#' 
## ----include=FALSE------------------------------------------------------------

# MDA-MB231 breast cancer- 4 time-points
breast_MDA.1.net <- reconstructNetwork(breast_MDA.1, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)
breast_MDA.2.net <- reconstructNetwork(breast_MDA.2, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)
breast_MDA.3.net <- reconstructNetwork(breast_MDA.3, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)

breast_MDA.1.p <- plotNetworkWiring(breast_MDA.1.net, plotIt=F)
breast_MDA.2.p <- plotNetworkWiring(breast_MDA.2.net, plotIt=F)
breast_MDA.3.p <- plotNetworkWiring(breast_MDA.3.net, plotIt=F)

# MCF7 breast cancer - 4 time-points
breast_MCF7.1.net <- reconstructNetwork(breast_MCF7.1, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
breast_MCF7.2.net <- reconstructNetwork(breast_MCF7.2, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
breast_MCF7.3.net <- reconstructNetwork(breast_MCF7.3, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)

breast_MCF7.1.p <- plotNetworkWiring(breast_MCF7.1.net, plotIt=F)
breast_MCF7.2.p <- plotNetworkWiring(breast_MCF7.2.net, plotIt=F)
breast_MCF7.3.p <- plotNetworkWiring(breast_MCF7.3.net, plotIt=F)

# MCF7 breast cancer - 8 time-points
breast_MCF7.net <- reconstructNetwork(breast_MCF7, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
breast_MCF7.p <- plotNetworkWiring(breast_MCF7.net, plotIt=F)

# All breast cancer nets merged:
all.nets <- reconstructNetwork(list(breast_MCF7.1, breast_MCF7.2, breast_MCF7.3, breast_MCF7, breast_MDA.1, breast_MDA.2, breast_MDA.3),
                               method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)

all.p <- plotNetworkWiring(all.nets, plotIt=F)



#' 
#' # Network inference:
#' 
## ----fig.width=10-------------------------------------------------------------

# MDA-MB231 breast cancer - 4 time-points
par(mfrow = c(1,3))
plot(breast_MDA.1.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MDA-MB231 breast\n 4 steps, replicate 1")
plot(breast_MDA.2.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MDA-MB231 breast\n 4 steps, replicate 2")
plot(breast_MDA.3.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MDA-MB231 breast\n 4 steps, replicate 3")

# MCF7 breast cancer - 4 time-points 
par(mfrow = c(1,3))
plot(breast_MCF7.1.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MCF7 breast\n 4 steps, replicate 1")
plot(breast_MCF7.2.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MCF7 breast\n 4 steps, replicate 2")
plot(breast_MCF7.3.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MCF7 breast\n 4 steps, replicate 3")

# MCF7 breast cancer - 8 time-points
par(mfrow = c(1,3))
plot(breast_MCF7.p, vertex.label.color="#440154ff",  vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MCF7 breast\n 8 steps")

# all samples breast cancer merged
par(mfrow = c(1,1))
plot(all.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="All breast cell-lines")

#' 
## -----------------------------------------------------------------------------
HIFaxis.net <- loadNetwork("boolean_network_HIFaxis.bn")
HIFaxis.p <- plotNetworkWiring(HIFaxis.net, plotIt=F)
HIFaxis.net

#' 
## -----------------------------------------------------------------------------
plot(HIFaxis.p, vertex.label.color="#440154ff",  vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="HIF axis Network\n Theoretical")

#' 
#' 
