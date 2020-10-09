#' ---
#' title: "Discretization of Hypoxia and Normoxia"
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

packages_cran = c("tidyverse", "BiTrinA")
  
#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from CRAN and loaded
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(package.check, packages_bioconductor, packages_cran)

#' 
## -----------------------------------------------------------------------------

load("..data/pipeline_EGEOD18494.RData")
load("..data/pipeline_GSE142867.RData")

#' 
#' 
#' # BiTrinA
#' 
#' 
#' Müssel C, Schmid F, Blätte TJ, Hopfensitz M, Lausser L, Kestler HA. BiTrinA--multiscale binarization and trinarization with quality analysis. Bioinformatics. 2016;32(3):465-468. doi:10.1093/bioinformatics/btv591. https://pubmed.ncbi.nlm.nih.gov/26468003/
#' 
#' https://cran.r-project.org/web/packages/BiTrinA/vignettes/Vignette.pdf
#' 
#' 
#' *Arguments*
#' 
#' - *method*: Chooses the BASC method to use (see details), i.e. either "A" or "B".
#' 
#' - *vect*: A real-valued vector of data to binarize.
#' 
#' - *tau*: This parameter adjusts the sensitivity and the specificity of the statistical testing procedure that rates the quality of the binarization. Defaults to 0.01.
#' 
#' - *Compute a series of step functions*: An initial step function is obtained by rearranging the original time series measurements in increasing order. Then, step functions with fewer discontinuities are calculated. BASC A calculates these step functions in such a way that each minimizes the Euclidean distance to the initial step function. BASC B obtains step functions from smoothened versions of the input function in a scale-space manner.
#' 
#' - *Find strongest discontinuity in each step function*: A strong discontinuity is a high jump size (derivative) in combination with a low approximation error.
#' 
#' - *Estimate location and variation of the strongest discontinuities*: Based on these estimates, data values can be excluded from further analyses.    
#'     
#'   
## -----------------------------------------------------------------------------

# Correcting an outlier on HIF1A:
exp.GSE142867.hif$GSM4246818 <- NULL
exp.GSE142867.hif.pivot <- exp.GSE142867.hif.pivot[exp.GSE142867.hif.pivot$codes !="GSM4246818", ]
pData.GSE142867 <- pData.GSE142867[pData.GSE142867$codes != "GSM4246818",]

result <- binarize.kMeans(exp.GSE142867.hif.pivot[exp.GSE142867.hif.pivot$Symbol =="HIF1A", "value"])
print(result@originalMeasurements)

print(result)
plot(result)

result <- binarize.BASC(exp.GSE142867.hif.pivot[exp.GSE142867.hif.pivot$Symbol =="HIF1A", "value"], method="A", tau=0.15)

print(result)
plot(result)

result <- binarize.BASC(exp.GSE142867.hif.pivot[exp.GSE142867.hif.pivot$Symbol =="HIF1A", "value"], method="B", tau=0.15)

print(result)
plot(result)


#' 
#' 
## -----------------------------------------------------------------------------
result <- binarize.kMeans(exp.GSE142867.hif.pivot[exp.GSE142867.hif.pivot$Symbol =="EPAS1", "value"])
print(result@originalMeasurements)
print(result)
plot(result)

result <- binarize.BASC(exp.GSE142867.hif.pivot[exp.GSE142867.hif.pivot$Symbol =="EPAS1", "value"], method="A", tau=0.15)

print(result)
plot(result)

result <- binarize.BASC(exp.GSE142867.hif.pivot[exp.GSE142867.hif.pivot$Symbol =="EPAS1", "value"], method="B", tau=0.15)

print(result)
plot(result)


#' 
#' # BinarizeMatrix of GSE142867
#' 
#' 
## -----------------------------------------------------------------------------

exp.GSE142867.hif.bin <- binarizeMatrix(exp.GSE142867.hif[,c(2:10)], 
               method = c("BASCA"), 
               adjustment = "none")

exp.GSE142867.hif.bin <- cbind(exp.GSE142867.hif.bin, exp.GSE142867.hif$Symbol)

pData.GSE142867$rep <- pData.GSE142867$title %>% 
                      as.character(.) %>%
                      strsplit( "rep " ) %>%
                      sapply( "[", 2 ) %>% str_replace_na(replacement="")


names(exp.GSE142867.hif.bin) <- c(paste0(substr(pData.GSE142867$condition,1,4),".", pData.GSE142867$time, ".", pData.GSE142867$rep), c("threshold", "p.value", "Symbol"))

# Selecting the probes with smaller p-value
exp.GSE142867.hif.bin <- data.frame(exp.GSE142867.hif.bin %>%
     group_by(Symbol) %>%
     slice(which.min(p.value)))

head(exp.GSE142867.hif.bin) %>% 
  knitr::kable(.)

#' 
## -----------------------------------------------------------------------------
exp.GSE142867.hif.mean <- exp.GSE142867.hif.bin %>%
 mutate(norm = rowMeans(dplyr::select(., starts_with("norm"))),
        hypo.d3 = rowMeans(dplyr::select(., starts_with("hypo.d3"))),
        hypo.d7 = rowMeans(dplyr::select(., starts_with("hypo.d7"))),
        hypo.d10 = rowMeans(dplyr::select(., starts_with("hypo.d10"))))  %>%
  dplyr::select(., -ends_with(c(".",".1",".2")))

exp.GSE142867.hif.mean %>% 
  knitr::kable(.)

#' 
#' 
## -----------------------------------------------------------------------------

exp.GSE142867.hif.pivot <- exp.GSE142867.hif.mean %>%
  group_by(Symbol) %>%
  pivot_longer(cols = starts_with(c("norm","hypo")), names_to = "codes", values_to = "value")

exp.GSE142867.hif.pivot$codes <- factor(exp.GSE142867.hif.pivot$codes,  levels =  c("norm", "hypo.d3" , "hypo.d7" , "hypo.d10"))

exp.GSE142867.hif.pivot$time <- as.numeric(exp.GSE142867.hif.pivot$codes)

head(exp.GSE142867.hif.pivot) %>% 
  knitr::kable(.)


#'    
#'    
## -----------------------------------------------------------------------------
ggplot(aes(x = factor(time), y = value, group = Symbol, color="red"),  
           data = exp.GSE142867.hif.pivot[exp.GSE142867.hif.pivot$Symbol %in% c("HIF1A", "EPAS1", "TP53", "CCND1", "MYC", "BAD"),]) +
  geom_point() + 
  geom_line() + 
  scale_x_discrete(breaks = c(1, 2, 3, 4), 
                 labels = c("Normoxia", "Hypoxia: Day 3" , "Hypoxia: Day 7" , "Hypoxia: Day10")) +
  xlab("Conditions") + ylab("Gene Expression") +
  theme(legend.position = "none", axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8)) +
  #geom_line(aes(linetype=Symbol, color=Symbol)) +
  facet_wrap(~ Symbol) 

ggplot(aes(x = factor(time), y = value, group = Symbol, color="red"),  
           data = exp.EGEOD18494.hif.mean.pivot[exp.EGEOD18494.hif.mean.pivot$Symbol %in% c("HIF1A", "EPAS1", "TP53", "CCND1", "MYC", "BAD"),]) +
  geom_point() + 
  geom_line() + 
  scale_x_discrete(breaks = c(1, 2, 3, 4), 
                 labels = c("Normoxia", "Hypoxia: 4h" , "Hypoxia: 8h" , "Hypoxia: 12h")) +
  xlab("Conditions") + ylab("Gene Expression") +
  theme(legend.position = "none", axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8)) +
  #geom_line(aes(linetype=Symbol, color=Symbol)) +
  facet_wrap(~ Symbol) 

#'    
## -----------------------------------------------------------------------------
ggplot(aes(x = factor(time), y = value, group = Symbol, color="red"),  
           data = exp.EGEOD18494.hif.mean.pivot[exp.EGEOD18494.hif.mean.pivot$Symbol %in% c("HIF1A", "EPAS1", "TP53", "CCND1", "MYC", "BAD"),]) +
  geom_point() + 
  geom_line() + 
  scale_x_discrete(breaks = c(1, 2, 3, 4), 
                 labels = c("Normoxia", "Hypoxia: 4h" , "Hypoxia: 8h" , "Hypoxia: 12h")) +
  xlab("Conditions") + ylab("Gene Expression") +
  theme(legend.position = "none", axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8)) +
  #geom_line(aes(linetype=Symbol, color=Symbol)) +
  facet_wrap(~ Symbol) 

#'    
#'    
## -----------------------------------------------------------------------------
result <- binarize.kMeans(exp.EGEOD18494.hif.pivot[exp.EGEOD18494.hif.pivot$Symbol =="HIF1A", "value"])
print(result@originalMeasurements)

print(result)
plot(result)

result <- binarize.BASC(exp.EGEOD18494.hif.pivot[exp.EGEOD18494.hif.pivot$Symbol =="HIF1A", "value"], method="A", tau=0.15)

print(result)
plot(result)

result <- binarize.BASC(exp.EGEOD18494.hif.pivot[exp.EGEOD18494.hif.pivot$Symbol =="HIF1A", "value"], method="B", tau=0.15)

print(result)
plot(result)


#' 
#' 
## -----------------------------------------------------------------------------
result <- binarize.kMeans(exp.EGEOD18494.hif.pivot[exp.EGEOD18494.hif.pivot$Symbol =="EPAS1", "value"])
print(result@originalMeasurements)
print(result)
plot(result)

result <- binarize.BASC(exp.EGEOD18494.hif.pivot[exp.EGEOD18494.hif.pivot$Symbol =="EPAS1", "value"], method="A", tau=0.15)

print(result)
plot(result)

result <- binarize.BASC(exp.EGEOD18494.hif.pivot[exp.EGEOD18494.hif.pivot$Symbol =="EPAS1", "value"], method="B", tau=0.15)

print(result)
plot(result)


#' 
#' # BinarizeMatrix of EGEOD18494
#' 
#' 
## -----------------------------------------------------------------------------
library(dplyr)

df <- exp.EGEOD18494.hif %>% dplyr::select(contains("GSM"))

exp.EGEOD18494.hif.bin <- binarizeMatrix(df, 
               method = c("BASCA"), 
               adjustment = "none")

exp.EGEOD18494.hif.bin <- cbind(exp.EGEOD18494.hif.bin, exp.EGEOD18494.hif$Symbol)

names <- exp.EGEOD18494.hif.bin %>% dplyr::select(contains("GSM")) %>% names()
pData.EGEOD18494 <- pData.EGEOD18494[order(pData.EGEOD18494$codes),]
pData.EGEOD18494$rep <- rep(1:3, n= length(pData.EGEOD18494$codes))
names(exp.EGEOD18494.hif.bin) <- c(paste0(substr(pData.EGEOD18494$condition,1,4),".", pData.EGEOD18494$time,".", substr(pData.EGEOD18494$cell_line,1,1) , pData.EGEOD18494$rep), c("threshold", "p.value", "Symbol"))


# Selecting the probes with smaller p-value
exp.EGEOD18494.hif.bin <- data.frame(exp.EGEOD18494.hif.bin %>%
     group_by(Symbol) %>%
     slice(which.min(p.value)))

head(exp.EGEOD18494.hif.bin)%>% 
  knitr::kable(.)


#' 
#' 
#' 
## -----------------------------------------------------------------------------
exp.EGEOD18494.hif.mean <- exp.EGEOD18494.hif.bin %>%
 mutate(norm = rowMeans(dplyr::select(., starts_with("norm.control"))),
        hypo.4h = rowMeans(dplyr::select(., starts_with("hypo.4"))),
        hypo.8h = rowMeans(dplyr::select(., starts_with("hypo.8"))),
        hypo.12h = rowMeans(dplyr::select(., starts_with("hypo.12"))))  %>%
  dplyr::select(., -ends_with(c("1","2","3")))

exp.EGEOD18494.hif.mean %>% 
  knitr::kable(.)

#' 
#' 
## -----------------------------------------------------------------------------

exp.EGEOD18494.hif.mean.pivot <- exp.EGEOD18494.hif.mean %>%
  group_by(Symbol) %>%
  pivot_longer(cols = starts_with(c("norm","hypo")), names_to = "codes", values_to = "value")

exp.EGEOD18494.hif.mean.pivot$codes <- factor(exp.EGEOD18494.hif.mean.pivot$codes,  levels =  c("norm", "hypo.4h" , "hypo.8h" , "hypo.12h"))

exp.EGEOD18494.hif.mean.pivot$time <- as.numeric(exp.EGEOD18494.hif.mean.pivot$codes)

head(exp.EGEOD18494.hif.mean.pivot) %>% 
  knitr::kable(.)


#'    
#'    
## -----------------------------------------------------------------------------
ggplot(aes(x = factor(time), y = value, group = Symbol, color="red"),  
           data = exp.EGEOD18494.hif.mean.pivot[exp.EGEOD18494.hif.mean.pivot$Symbol %in% c("HIF1A", "EPAS1", "TP53", "CCND1", "MYC", "BAD"),]) +
  geom_point() + 
  geom_line() + 
  scale_x_discrete(breaks = c(1, 2, 3, 4), 
                 labels = c("Normoxia", "Hypoxia: 4h" , "Hypoxia: 8h" , "Hypoxia: 12h")) +
  xlab("Conditions") + ylab("Gene Expression") +
  theme(legend.position = "none", axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8)) +
  #geom_line(aes(linetype=Symbol, color=Symbol)) +
  facet_wrap(~ Symbol) 

