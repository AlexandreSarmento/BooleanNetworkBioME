#' ---
#' title: "BoolNet Inference  (E-GEOD-18494)"
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

library(RSQLite, lib.loc = "/usr/local/lib/R/site-library")

# For oligo and ArrayExpress First install:
#install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz',repos=NULL)

packages_bioconductor = c("Biobase", "GEOquery", "ArrayExpress", "hgu133plus2.db", "preprocessCore")

# Install and load packages
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(package.check, packages_bioconductor, packages_cran)


#' 
## ----message=FALSE, warning=FALSE---------------------------------------------
download_dir <- fs::path(".data_tmp")
if (!dir_exists(download_dir)) {
    dir_create(download_dir)
    EGEOD18494 <- ArrayExpress( "E-GEOD-18494", save=TRUE, path=download_dir)
} else {
    EGEOD18494 <- ArrayExpress( "E-GEOD-18494", save=TRUE, path=download_dir)
}

data.EGEOD18494 <- Biobase::pData(EGEOD18494)

data.EGEOD18494 <- data.frame(
                  codes = substr(data.EGEOD18494$Source.Name,1,9),
                  cell_line = data.EGEOD18494$Characteristics..cell.line.,
                  time = data.EGEOD18494$Characteristics..time,
                  condition = data.EGEOD18494$Characteristics..stress.
                  )
data.EGEOD18494 <- data.EGEOD18494[order(data.EGEOD18494$codes),]
data.EGEOD18494$rep <- rep(1:3, n= length(data.EGEOD18494$codes))

# Normalisation
eset.EGEOD18494 <- oligo::rma(EGEOD18494,  normalize = TRUE)

exp.EGEOD18494 <- exprs(eset.EGEOD18494)

colnames(exp.EGEOD18494) <- substr(colnames(exp.EGEOD18494),1,9)

EGEOD18494@annotation

rm(download_dir)

#' 
#' # Convert the probes to Symbol names
#' 
## -----------------------------------------------------------------------------
anno.EGEOD18494 <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(exp.EGEOD18494), columns=c("ENSEMBL", "SYMBOL", "GENENAME"), keytype="PROBEID")

colnames(anno.EGEOD18494) <- c("probes", "ensgene", "symbol", "description")


#' 
#' # Selecting the HIF Genes 
#'     
## -----------------------------------------------------------------------------

# Selecting genes from HIF Axis
hif.symbols <- c("TP53", "HIF1A", "EP300", "MDM2", "VHL")

hif.probes <- anno.EGEOD18494$probes[anno.EGEOD18494$symbol %in% hif.symbols]

# Select the probes and genes
exp.EGEOD18494.hif <- as.data.frame(exp.EGEOD18494) %>% 
  rownames_to_column('probes') %>% 
  filter(probes %in% hif.probes) %>% 
  merge(anno.EGEOD18494[anno.EGEOD18494$symbol %in% hif.symbols, c("probes","symbol")], by = "probes") %>% 
  #distinct(symbol, .keep_all = TRUE) %>% # Take the first one
  dplyr::select(!(probes)) 
  

# Function to binarize according an consensus mean of probes, add the O2 state and rename columns 
binNet <- function(b){
  binarizeTimeSeries(b[,-5], method="kmeans")$binarizedMeasurements  %>% 
  as.data.frame(.)  %>% 
  aggregate(., list(symbol = b$symbol), mean) %>% 
  mutate_at(vars(-symbol), funs(ifelse(. > 0.4, 1, 0))) %>% 
  rbind(., c("O2", 1,0,0,0)) %>% 
    rename_at(vars(data.EGEOD18494$codes[data.EGEOD18494$codes %in% names(b)] ),
            ~paste0(substr(data.EGEOD18494$condition[data.EGEOD18494$codes %in% names(b)],1,4),".",
                    data.EGEOD18494$time[data.EGEOD18494$codes %in% names(b)],".",
                    substr(data.EGEOD18494$cell_line[data.EGEOD18494$codes %in% names(b)],1,1), ".",
                    data.EGEOD18494$rep[data.EGEOD18494$codes %in% names(b)])) %>% 
  column_to_rownames("symbol")
}


#' 
#' # Exemplifying the Binarization 
#' 
## -----------------------------------------------------------------------------

breast1x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &
                  data.EGEOD18494$rep == 1], "symbol")) %>% arrange(symbol)

names(breast1x) <- c("norm.control.M.1",  "hypo.4h.M.1", "hypo.8h.M.1", "hypo.12h.M.1", "symbol")

breast1x[, c("symbol","norm.control.M.1",  "hypo.4h.M.1", "hypo.8h.M.1", "hypo.12h.M.1")] %>% 
  knitr::kable(.)

binarizeTimeSeries(breast1x[,-5], method="kmeans")$binarizedMeasurements  %>% 
  as.data.frame(.)  %>% 
  add_column(symbol = breast1x$symbol) %>%   dplyr::select( c("symbol","norm.control.M.1",  "hypo.4h.M.1", "hypo.8h.M.1", "hypo.12h.M.1"))  %>% 
  knitr::kable(.)
  
binarizeTimeSeries(breast1x[,-5], method="kmeans")$binarizedMeasurements  %>% 
  as.data.frame(.)  %>% 
  aggregate(., list(symbol = breast1x$symbol), mean) %>% 
  mutate_at(vars(-symbol), funs(ifelse(. > 0.4, 1, 0))) %>% 
  rbind(., c("O2", 1,0,0,0)) %>% 
  knitr::kable(.)
  

#' 
#' # MDA-MB231 breast cancer
#' 
## -----------------------------------------------------------------------------

breast1x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &
                  data.EGEOD18494$rep == 1], "symbol"))  %>% 
  binNet(.) 

breast1x %>% 
  knitr::kable(.)

breast2x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &
                  data.EGEOD18494$rep == 2], "symbol"))  %>% 
  binNet(.) 

breast2x  %>% 
  knitr::kable(.)

breast3x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "MDA-MB231 breast cancer" &
                  data.EGEOD18494$rep == 3], "symbol"))  %>% 
  binNet(.) 

breast3x %>% 
  knitr::kable(.)

# All breast cancer nets merged:

net <- reconstructNetwork(list(breast1x, breast2x, breast3x), method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

# Individual nets of each replica:

net <- reconstructNetwork(breast1x, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

net <- reconstructNetwork(breast2x, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

net <- reconstructNetwork(breast3x, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

#' 
#' 
#' # HepG2 hepatoma
#' 
## -----------------------------------------------------------------------------

hepatoma1x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "HepG2 hepatoma" &
                  data.EGEOD18494$rep == 1], "symbol"))  %>% 
  binNet(.) 

hepatoma1x %>% 
  knitr::kable(.)

hepatoma2x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "HepG2 hepatoma" &
                  data.EGEOD18494$rep == 2], "symbol"))  %>% 
  binNet(.) 

hepatoma2x %>% 
  knitr::kable(.)

hepatoma3x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "HepG2 hepatoma" &
                  data.EGEOD18494$rep == 3], "symbol"))  %>% 
  binNet(.) 

hepatoma3x %>% 
  knitr::kable(.)

# All nets hepatoma merged:

net <- reconstructNetwork(list(hepatoma1x, hepatoma2x, hepatoma3x), method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

# Individual nets of each replica:

net <- reconstructNetwork(hepatoma1x, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

net <- reconstructNetwork(hepatoma2x, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

net <- reconstructNetwork(hepatoma3x, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

#' 
#' # U87 glioma
#' 
## -----------------------------------------------------------------------------

glioma1x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "U87 glioma" &
                  data.EGEOD18494$rep == 1], "symbol"))  %>% 
  binNet(.) 

glioma1x %>% 
  knitr::kable(.)

glioma2x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "U87 glioma" &
                  data.EGEOD18494$rep == 2], "symbol"))  %>% 
  binNet(.) 

glioma2x %>% 
  knitr::kable(.)

glioma3x <- 
exp.EGEOD18494.hif %>% 
  dplyr::select(c(data.EGEOD18494$codes[data.EGEOD18494$cell_line == "U87 glioma" &
                  data.EGEOD18494$rep == 3], "symbol"))  %>% 
  binNet(.) 

glioma3x %>% 
  knitr::kable(.)

# All glioma nets merged:

net <- reconstructNetwork(list(glioma1x, glioma2x, glioma3x), method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

# Individual nets of each replica:

net <- reconstructNetwork(glioma1x, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

net <- reconstructNetwork(glioma2x, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

net <- reconstructNetwork(glioma3x, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)

#' 
