#' ---
#' title: "BoolNet Inference  (GSE41491)"
#' output:
#'   pdf_document: default
#'   html_notebook: default
#'   github_document: 
#'     df_print: paged
#'     html_preview: FALSE
#'     keep_html: TRUE
#' ---
#' 

#' 
#' Hypoxia transcriptomic time-series data in three different cancer cell lines (GSE41491)
#' 
#' Tumour hypoxia exhibits a highly dynamic spatial and temporal distribution and is associated with increased malignancy and poor prognosis.
#' Assessment of time-dependent gene-expression changes in response to hypoxia may thus provide additional biological insights and help with patient prognosis.
#' 
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41491
#' 
#' Exponentially growing prostate (DU145), colon (HT29) and breast (MCF7) carcinoma cells were seeded on glass dishes in McCoy, DMEM or RPMI media, respectively with 10% FCS.
#' 
#' All analyses were performed in R (v2.11.1). Data processing was performed for each cell line separately using RMA (affy package, v1.26.1) for pre-processing, and updated Entrez GeneID annotation (hgu133plus2hsentrezgcdf package, v13.0.0).
#' 
## ----message=FALSE, warning=FALSE---------------------------------------------
packages_cran = c("igraph", "BoolNet", "BiocManager", "tidyverse", "fs", "ff")


# Install and load packages
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# For oligo and ArrayExpress First install:
#install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz',repos=NULL)

packages_bioconductor = c("Biobase", "GEOquery", "oligo", "ArrayExpress")

# Install and load packages
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(package.check, packages_bioconductor, packages_cran)


#' 
## -----------------------------------------------------------------------------

download_dir <- fs::path(".data_tmp")
if (!dir_exists(download_dir)) { dir_create(download_dir) }  

GSE41491 <-getGEO("GSE41491", destdir = download_dir, GSEMatrix = T)

expr.GSE41491 <- exprs(oligo::rma(GSE41491[[1]], normalize = TRUE))
data.GSE41491 <- pData(GSE41491[[1]])

data.GSE41491 <- data.frame(
                  codes = as.character(data.GSE41491$geo_accession),
                  cell_line = gsub("[(]","",gsub("[)]","",substr(as.character(data.GSE41491$source_name_ch1),18,25))),
                  time = substr(as.character(data.GSE41491$characteristics_ch1.1),18,19),
                  condition = substr(as.character(data.GSE41491$characteristics_ch1),12,13))

# #Affymetrix Human Genome U133 Plus 2.0 Array [Brainarray Version 15.0.0, HGU133Plus2_Hs_ENTREZG]
# #  load/install the package
# if(!require("hgu133plus2hsentrezg.db")){ install.packages('http://mbni.org/customcdf/13.0.0/entrezg.download/hgu133plus2hsentrezg.db_13.0.0.tar.gz',repos=NULL) }


library(readr)
GPL14877 <- read_delim(".data_tmp/GPL14877.soft", 
    "\t", escape_double = FALSE, trim_ws = TRUE, 
    skip = 393)

ENTREZID <- as.character(unique(GPL14877$SPOT_ID))
PROBES <- unique(rownames(expr.GSE41491))

library(hgu133plus2hsentrezg.db)
library(annotate)
SYMBOL <- getSYMBOL(PROBES, "hgu133plus2hsentrezg")
results <- as.data.frame(cbind(PROBES, SYMBOL))

# Convert the probes to Symbol names
anno.GSE41491 <- AnnotationDbi::select(hgu133plus2.db, 
       keys = ENTREZID, 
       columns=c("SYMBOL", "GENENAME"), 
       keytype="ENTREZID")

anno.GSE41491 <- merge(anno.GSE41491, results, by="SYMBOL", all.x = T) %>%
  filter(!is.na(SYMBOL)) %>%
  mutate_all(as.character) %>%
  rename_at(vars(c("SYMBOL", "ENTREZID", "GENENAME", "PROBES")), ~ c("symbol", "entrezid", "description", "probes"))

rm(download_dir, GSE41491, results, GPL14877, PROBES, SYMBOL, ENTREZID)


#' 
#' 
#' # Selecting the HIF Genes 
#'     
## -----------------------------------------------------------------------------
# Genes from Boolean Network:    
# HIF1a, HIF2a, p53, BNIP3, VEGF, cMyc, Oct4, cdc20, cycA, cycB, cycE, cycD, p27, Rb, E2F, cdh1, mdm2, BAD, BclX

# hif.symbols <- c("HIF1A", "HIF1", "PASD8", "MOP1", "EPAS1", "HIF2A", "HLF", "PASD2", "MOP2", "VEGFA", "VPF", "MVCD1", "VEGF-A", "TP53", "P53", "MYC", "C-Myc", "POU5F1", "OCT3", "OTF3", "CDC20", "P55CDC", "CDC20A", "CCNA1", "CCNA1", "CCNA2", "CCN1", "CCNA", "CCNB1", "CCNB", "CCNB1", "CCNB2", "CCNB2", "HsT17299", "CCND1", "PRAD1", "CCND2", "MPPH3", "CCNE1", "CCNE1", "CCNE", "PCCNE1", "CCNE2", "CCNE2", "CCNE", "PCCNE1", "CDKN1B", "KIP1", "P27KIP1", "CDKN4MEN4", "RB1", "PPP1R130", "Pp110", "E2F1", "RBAP-1", "BNIP3", "NIP3", "BCL2", "MCL1", "BCL2L3", "CDH1", "CD324", "UBE2C", "BAD", "BCL2L8", "VHL", "PVHL", "VHL1", "MDM2", "HDM2", "EP300", "P300", "KAT3B")

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
  

#' 
#' 
## -----------------------------------------------------------------------------
# Function to binarize according an consensus mean of probes, add the O2 state and rename columns 

binNet <- function(b){
  
  cols <- data.GSE41491$codes %in% names(b)
  
  binarizeTimeSeries(b[,-1], method="kmeans")$binarizedMeasurements  %>% 
  as.data.frame(.)  %>% 
  aggregate(., list(symbol = b$symbol), mean) %>% 
  mutate_at(vars(-symbol), funs(ifelse(. >= 0.5, 1, 0))) %>% 
  rbind(., c("O2", 1,,0,0,0,0,0,0,0)) %>% 
    rename_at(vars(data.GSE41491$codes[cols] ),
            ~paste0(substr(data.GSE41491$condition[cols],1,2),".",
                    data.GSE41491$time[cols],".",
                    substr(data.GSE41491$cell_line[cols],1,2))) %>% 
  column_to_rownames("symbol")
  
}

#' 
#' Exponentially growing prostate (DU145), colon (HT29) and breast (MCF7) carcinoma cells were seeded on glass dishes in McCoy, DMEM or RPMI media, respectively with 10% FCS.
#' 
## -----------------------------------------------------------------------------

breast1_MCF7 <- 
expr.GSE41491.hif %>% 
  dplyr::select(c("symbol", data.GSE41491$codes[data.GSE41491$cell_line == "MCF7"])) %>% arrange(symbol) %>% 
  rename_at(vars(data.GSE41491$codes[data.GSE41491$cell_line == "MCF7"]),
            ~paste0(data.GSE41491$condition[data.GSE41491$cell_line == "MCF7"],".",
                    data.GSE41491$time[data.GSE41491$cell_line == "MCF7"],".",
                    data.GSE41491$cell_line[data.GSE41491$cell_line == "MCF7"]))

knitr::kable(breast1_MCF7)

binarizeTimeSeries(breast1_MCF7[,-1], method="kmeans")$binarizedMeasurements  %>% 
  data.frame(.)  %>% 
  add_column(symbol = breast1_MCF7$symbol, .before=0) %>% 
  knitr::kable(.)
  
binarizeTimeSeries(breast1_MCF7[,-1], method="kmeans")$binarizedMeasurements  %>% 
  data.frame(.)  %>% 
  aggregate(., list(symbol = breast1_MCF7$symbol), mean) %>% 
  mutate_at(vars(-symbol), funs(ifelse(. >= 0.5, 1, 0))) %>% 
  rbind(., c("O2", 1,0,0,0,0,0,0,0)) %>% 
  knitr::kable(.)
  

#' 
#' # breast (MCF7) carcinoma cells 
#' 
## -----------------------------------------------------------------------------

breast_MCF7 <- 
expr.GSE41491.hif %>% 
  dplyr::select(c("symbol", data.GSE41491$codes[data.GSE41491$cell_line == "MCF7"]))  %>% 
  binNet(.) 
knitr::kable(breast_MCF7)


# Individual net

net <- reconstructNetwork(breast_MCF7, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)


#' 
#' # prostate (DU145) carcinoma cells 
#' 
## -----------------------------------------------------------------------------

prostate_DU145 <- 
expr.GSE41491.hif %>% 
  dplyr::select(c("symbol", data.GSE41491$codes[data.GSE41491$cell_line == "DU145"]))  %>% 
  binNet(.) 
knitr::kable(prostate_DU145)


net <- reconstructNetwork(prostate_DU145, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)


#' 
#' #  colon (HT29)  carcinoma cells 
#' 
## -----------------------------------------------------------------------------

colon_MCF7 <- 
expr.GSE41491.hif %>% 
  dplyr::select(c("symbol", data.GSE41491$codes[data.GSE41491$cell_line == "HT29"]))  %>% 
  binNet(.) 
knitr::kable(colon_MCF7)


# Individual net

net <- reconstructNetwork(colon_MCF7, method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)


#' 
#' 
#' # All three cells prostate (DU145), colon (HT29) and breast (MCF7) carcinoma 
#' 
## -----------------------------------------------------------------------------

net <- reconstructNetwork(list(breast_MCF7,prostate_DU145, colon_MCF7), method="bestfit",returnPBN=TRUE,readableFunctions=TRUE)
plotNetworkWiring(net)
print(net)


