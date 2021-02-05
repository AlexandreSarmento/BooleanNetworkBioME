#' ---
#' title: "BoolNet Inference HepG2 hepatoma, U87 glioma, and MDA-MB231 breast cancer (E-GEOD-18494)"
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
#' Expression profiling of hypoxic HepG2 hepatoma, U87 glioma, and MDA-MB231 breast cancer cells: time course (E-GEOD-18494)
#' 
#' Analysis of expression changes of cultured HepG2 hepatoma, U87 glioma, and MDA-MB231 breast cancer cells subjected to hypoxia (0.5% O2) for 0, 4, 8, 12 hours . Results provide insight to cell type-specific response to hypoxia. HepG2 hepatoma, U87 glioma, and MDA-MB231 breast cancer cells were collected under normoxic conditions (~19% O2, 0 hours) and after 4, 8 and 12 hours of hypoxia treatment (0.5% O2). For each cell line, three replicates of total RNA at each time point were prepared using Trizol and submitted to the DFCI Microarray Core for labeling, hybridization to Affymetrix HG-U133Plus2 oligonucleotide arrays and image scanning.
#'     
#' https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-18494/   
#' 
#' 
## ----message=FALSE, warning=FALSE---------------------------------------------
packages_cran = c("igraph", "BoolNet", "BiocManager", "tidyverse", "fs")
# Install and load packages
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
packages_bioconductor = c("Biobase", "GEOquery", "vsn", "hgu133plus2.db")
# Install and load packages
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(package.check, packages_bioconductor, packages_cran)


#' 
#' <!-- ```{r message=FALSE, warning=FALSE } -->
#' <!-- download_dir <- fs::path(".data_tmp") -->
#' <!-- if (!dir_exists(download_dir)) { -->
#' <!--     dir_create(download_dir) -->
#' <!--     EGEOD18494 <- ArrayExpress( "E-GEOD-18494", save=TRUE, path=download_dir) -->
#' <!-- } else { -->
#' <!--     EGEOD18494 <- ArrayExpress( "E-GEOD-18494", save=TRUE, path=download_dir) -->
#' <!-- } -->
#' 
#' <!-- data.EGEOD18494 <- Biobase::pData(EGEOD18494) -->
#' 
#' <!-- data.EGEOD18494 <- data.frame( -->
#' <!--                   codes = substr(data.EGEOD18494$Source.Name,1,9), -->
#' <!--                   cell_line = data.EGEOD18494$Characteristics..cell.line., -->
#' <!--                   time = data.EGEOD18494$Characteristics..time, -->
#' <!--                   condition = data.EGEOD18494$Characteristics..stress. -->
#' <!--                   ) -->
#' <!-- data.EGEOD18494 <- data.EGEOD18494[order(data.EGEOD18494$codes),] -->
#' <!-- data.EGEOD18494$rep <- rep(1:3, n= length(data.EGEOD18494$codes)) -->
#' 
#' <!-- # Normalisation -->
#' <!-- eset.EGEOD18494 <- oligo::rma(EGEOD18494,  normalize = TRUE) -->
#' 
#' <!-- expr.EGEOD18494 <- exprs(eset.EGEOD18494) -->
#' 
#' <!-- # Convert to a data.frame -->
#' <!-- expr.EGEOD18494 <- as.data.frame(as.ffdf(expr.EGEOD18494)) -->
#' 
#' <!-- colnames(expr.EGEOD18494) <- substr(colnames(expr.EGEOD18494),1,9) -->
#' 
#' <!-- # Convert the probes to Symbol names -->
#' <!-- anno.EGEOD18494 <- AnnotationDbi::select(hgu133plus2.db,  -->
#' <!--                                          keys=rownames(expr.EGEOD18494),  -->
#' <!--                                          columns=c("ENSEMBL", "SYMBOL", "GENENAME"),  -->
#' <!--                                          keytype="PROBEID") -->
#' 
#' <!-- colnames(anno.EGEOD18494) <- c("probes", "ensgene", "symbol", "description") -->
#' 
#' <!-- rm(download_dir, EGEOD18494, eset.EGEOD18494) -->
#' 
#' <!-- save.image("../data/data.EGEOD18494.Rdata") -->
#' <!-- ``` -->
#' 
#' # Load the pre-processed
#' 
## -----------------------------------------------------------------------------
load("../data/data.EGEOD18494.Rdata")
eset <- ExpressionSet(assayData = as.matrix(expr.EGEOD18494), 
                      probeNames = row.names(expr.EGEOD18494))
expr.EGEOD18494 <- exprs(justvsn(eset))

#' 
#' # Selecting the HIF Genes 
#'     
## -----------------------------------------------------------------------------

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
## -----------------------------------------------------------------------------
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


#' 
#' # Exemplifying the Binarization 
#' 
## -----------------------------------------------------------------------------

cols <- (data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" & data.EGEOD18494$rep == 1)

breast1x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cols])) %>% arrange(symbol) %>% 
  arrange(symbol) %>% 
  rename_at(vars(data.EGEOD18494$codes[cols]),
            ~paste0(substr(data.EGEOD18494$condition[cols],1,2),".",
                    data.EGEOD18494$time[cols],".",
                    substr(data.EGEOD18494$cell_line[cols],1,2)))

breast1x %>% 
  knitr::kable(.)

binarizeTimeSeries(breast1x[,-1], method="kmeans")$binarizedMeasurements  %>% 
  data.frame(.)  %>% 
  add_column(symbol = breast1x$symbol, .before=0) %>% 
  knitr::kable(.)
  
binarizeTimeSeries(breast1x[,-1], method="kmeans")$binarizedMeasurements  %>% 
  data.frame(.)  %>% 
  aggregate(., list(symbol = breast1x$symbol), mean) %>% 
  mutate_at(vars(-symbol), funs(ifelse(. >= 0.5, 1, 0))) %>% 
  rbind(., c("O2", 1,0,0,0)) %>% 
  knitr::kable(.)


#' 
#' # MDA-MB231 breast cancer
#' 
## -----------------------------------------------------------------------------
cellline.rep1 <- (data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &  data.EGEOD18494$rep == 1)
cellline.rep2 <- (data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &  data.EGEOD18494$rep == 2)
cellline.rep3 <- (data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &  data.EGEOD18494$rep == 3)

breast1x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep1])) %>% binNet(.) 

breast1x %>% knitr::kable(.)

breast2x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep2])) %>% binNet(.) 

breast2x  %>% knitr::kable(.)

breast3x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep3])) %>% binNet(.) 

breast3x %>% knitr::kable(.)

#' 
#' # HepG2 hepatoma
#' 
## -----------------------------------------------------------------------------

cellline.rep1 <- (data.EGEOD18494$cell_line == "HepG2 hepatoma" &  data.EGEOD18494$rep == 1)
cellline.rep2 <- (data.EGEOD18494$cell_line == "HepG2 hepatoma" &  data.EGEOD18494$rep == 2)
cellline.rep3 <- (data.EGEOD18494$cell_line == "HepG2 hepatoma" &  data.EGEOD18494$rep == 3)

hepatoma1x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep1]))  %>% 
  binNet(.) 

hepatoma1x %>% 
  knitr::kable(.)

hepatoma2x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep2]))  %>% 
  binNet(.) 

hepatoma2x %>% 
  knitr::kable(.)

hepatoma3x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep3]))  %>% 
  binNet(.) 

hepatoma3x %>% 
  knitr::kable(.)


#' 
#' # U87 glioma
#' 
## -----------------------------------------------------------------------------
cellline.rep1 <- (data.EGEOD18494$cell_line == "U87 glioma" &  data.EGEOD18494$rep == 1)
cellline.rep2 <- (data.EGEOD18494$cell_line == "U87 glioma" &  data.EGEOD18494$rep == 2)
cellline.rep3 <- (data.EGEOD18494$cell_line == "U87 glioma" &  data.EGEOD18494$rep == 3)

glioma1x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep1]))  %>% 
  binNet(.) 

glioma1x %>% 
  knitr::kable(.)

glioma2x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep2]))  %>% 
  binNet(.) 

glioma2x %>% 
  knitr::kable(.)

glioma3x <- 
expr.EGEOD18494.hif %>% 
  dplyr::select(c("symbol", data.EGEOD18494$codes[cellline.rep3]))  %>% 
  binNet(.) 

glioma3x %>% 
  knitr::kable(.)


#' 
#' 
## ----include=FALSE------------------------------------------------------------
breast1x.net <- reconstructNetwork(breast1x, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)
breast2x.net <- reconstructNetwork(breast2x, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)
breast3x.net <- reconstructNetwork(breast3x, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)

hepatoma1x.net <- reconstructNetwork(hepatoma1x, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)
hepatoma2x.net <- reconstructNetwork(hepatoma2x, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)
hepatoma3x.net <- reconstructNetwork(hepatoma3x, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)

glioma1x.net <- reconstructNetwork(glioma1x, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)
glioma2x.net <- reconstructNetwork(glioma2x, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)
glioma3x.net <- reconstructNetwork(glioma3x, method="bestfit", returnPBN=TRUE, readableFunctions=TRUE)

breast1x.p <- plotNetworkWiring(breast1x.net, plotIt=F)
breast2x.p <- plotNetworkWiring(breast2x.net, plotIt=F)
breast3x.p <- plotNetworkWiring(breast3x.net, plotIt=F)

hepatoma1x.p <- plotNetworkWiring(hepatoma1x.net, plotIt=F)
hepatoma2x.p <- plotNetworkWiring(hepatoma2x.net, plotIt=F)
hepatoma3x.p <- plotNetworkWiring(hepatoma3x.net, plotIt=F)

glioma1x.p <- plotNetworkWiring(glioma1x.net, plotIt=F)
glioma2x.p <- plotNetworkWiring(glioma2x.net, plotIt=F)
glioma3x.p <- plotNetworkWiring(glioma3x.net, plotIt=F)

# All breast cancer nets merged:
breast.all <- reconstructNetwork(list(breast1x, breast2x, breast3x),
                               method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
breast.all.p <- plotNetworkWiring(breast.all, plotIt=F)

# All hepatoma cancer nets merged:
hepatoma.all <- reconstructNetwork(list(hepatoma1x, hepatoma2x, hepatoma3x),
                               method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
hepatoma.all.p <- plotNetworkWiring(hepatoma.all, plotIt=F)

# All hepatoma cancer nets merged:
glioma.all <- reconstructNetwork(list(glioma1x, glioma2x, glioma3x),
                               method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
glioma.all.p <- plotNetworkWiring(glioma.all, plotIt=F)

#' 
#' # Network inference:
#' 
## -----------------------------------------------------------------------------

# MDA-MB231 breast cancer - 4 time-points
par(mfrow = c(1,3))
plot(breast1x.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MDA-MB231 breast\n 4 steps, replicate 1")
plot(breast2x.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MDA-MB231 breast\n 4 steps, replicate 2")
plot(breast3x.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MDA-MB231 breast\n 4 steps, replicate 3")

par(mfrow = c(1,1))
plot(breast.all.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="MDA-MB231 breast\n 4 steps, replicate 3")

# HepG2 hepatoma
par(mfrow = c(1,3))
plot(hepatoma1x.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="HepG2 hepatoma\n 4 steps, replicate 1")
plot(hepatoma2x.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="HepG2 hepatoma\n 4 steps, replicate 2")
plot(hepatoma3x.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="HepG2 hepatoma\n 4 steps, replicate 3")

par(mfrow = c(1,1))
plot(hepatoma.all.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="HepG2 hepatoma\n 4 steps, replicate 3")

# U87 glioma
par(mfrow = c(1,3))
plot(glioma1x.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="U87 glioma\n 4 steps, replicate 1")
plot(glioma2x.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="U87 glioma\n 4 steps, replicate 2")
plot(glioma3x.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="U87 glioma\n 4 steps, replicate 3")

par(mfrow = c(1,1))
plot(glioma.all.p, vertex.label.color="#440154ff", vertex.color="lightblue", vertex.frame.color="white", layout=layout_in_circle, edge.curved=.3,
     main="U87 glioma\n 4 steps, replicate 3")

