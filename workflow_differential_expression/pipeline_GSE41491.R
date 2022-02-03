require("GEOquery")
GSEDATA <- getGEO("GSE41491", GSEMatrix=T, AnnotGPL=FALSE)
eset <- GSEDATA[[1]]
eset


# Get the annotation GPL id 
gpl <- getGEO('GPL14877', destdir=".")
geneProbes <- which(!is.na(Table(gpl)$ID))
probeids <- as.character(Table(gpl)$ID[geneProbes])
probes <- intersect(probeids, rownames(exprs(eset)))
expr.GSE41491 <- exprs(eset)[probes, ]
inds <- which(Table(gpl)$ID %in% probes)

library("Biobase")
library("AnnotationDbi")
library("org.Hs.eg.db")

symbol <- mapIds(org.Hs.eg.db, keys=as.character(Table(gpl)$SPOT_ID), 
                      column="SYMBOL", keytype="ENTREZID", multiVals="first")

ensgene <- mapIds(org.Hs.eg.db, keys=as.character(Table(gpl)$SPOT_ID), 
                  column="ENSEMBL", keytype="ENTREZID", multiVals="first")

description <- Table(gpl)$Description
anno.GSE41491 <- as.data.frame(cbind(probes,ensgene,symbol,description))
rownames(anno.GSE41491) <- NULL
#save(expr.GSE41491,anno.GSE41491, file = "data.GSE41491.RData")

library("tidyverse")
library("dplyr")
main_varsLabel <- c("geo_accession", "characteristics_ch1", "source_name_ch1", "characteristics_ch1.1")
data.GSE41491 <- pData(eset)[main_varsLabel] %>% 
  rename(codes = 'geo_accession',
         condition = 'characteristics_ch1',
         cell_line = 'source_name_ch1',
         time = 'characteristics_ch1.1') %>%
  droplevels() %>%
  mutate(condition = fct_recode(condition,
                                'normoxia'='treatment: none',
                                'hypoxia' = 'treatment: hypoxia'),
         time = fct_recode(time,
                           '0h' = 'time point (hr): 0',
                           '1h' = 'time point (hr): 1',
                           '2h' = 'time point (hr): 2',
                           '4h' = 'time point (hr): 4',
                           '8h' = 'time point (hr): 8',
                           '12h'= 'time point (hr): 12',
                           '16h'= 'time point (hr): 16',
                           '24h'= 'time point (hr): 24'),
         
         cell_line = fct_recode(cell_line,
                                'HT29' = 'cancer cell line (HT29)',
                                'DU145' = 'cancer cell line (DU145)',
                                'MCF7' = 'cancer cell line (MCF7)')) %>%
           
  mutate(condition = factor(condition, levels = c('normoxia','hypoxia')),
         time = factor(time, levels = c('0h','1h','2h','4h','8h','12h','16h','24h')),
         cell_line = factor(cell_line,levels = c('DU145','HT29','MCF7')))

row.names(data.GSE41491) <- NULL
# Arrange the expression data and the data description on same order:
expr.GSE41491 <- expr.GSE41491[,data.GSE41491$codes]

save(expr.GSE41491,anno.GSE41491,data.GSE41491, file = "data.GSE41491.RData")
