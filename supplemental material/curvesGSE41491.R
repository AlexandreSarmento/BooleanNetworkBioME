packages_cran = c("BiTrinA","igraph", "BoolNet", "BiocManager", "tidyverse", "fs")
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


setwd("~/R/dataAnalysis_UFRN/BooleanNetworkBioME/TBN_IBN")
load("~/R/dataAnalysis_UFRN/BooleanNetworkBioME/data/data.GSE41491.Rdata")

#minMaxNorm <- function(x) {
#  (x - min(x)) / (max(x) - min(x))
#}

#netnodes.symbols <- c("AKT1","ATR",
#"BAX","EGLN2","HIF1A","HIF1AN","MDM2",
#"PFKL",
# "PPP1R10","PTEN","TP53","VEGFA","INPP5D")

geneName <- "PIK3CA"
#plotTitle <- "HIF1a"
netnodes.symbols <- c(geneName)


netnodes.probes <- anno.GSE41491$probes[anno.GSE41491$symbol %in% netnodes.symbols]

# Select the probes and genes
# GSE41491
expr.GSE41491.nodes <- as.data.frame(expr.GSE41491) %>% 
  rownames_to_column('probes') %>% 
  filter(probes %in% netnodes.probes) %>% 
  merge(anno.GSE41491[anno.GSE41491$symbol %in% netnodes.symbols, c("probes","symbol")], by = "probes") %>% 
  group_by(symbol) %>%
  summarise_at(vars(-probes),funs(median(., na.rm=TRUE)))  %>%
  column_to_rownames(var = "symbol") %>%
  dplyr::select(c(data.GSE41491$codes[data.GSE41491$cell_line == "MCF7"]))


#expr.GSE41491.nodes.norm <- minMaxNorm(expr.GSE41491.nodes)
#expr.GSE41491.nodes.norm$symbol <- row.names(expr.GSE41491.nodes.norm)
row <- data.GSE41491$cell_line == "MCF7"
expr.GSE41491.nodes <- expr.GSE41491.nodes[, c(as.character(data.GSE41491$codes[row]))]
names(expr.GSE41491.nodes) <- c(paste0(substr(data.GSE41491$condition[row],1,4),".", data.GSE41491$time[row], "."))
expr.GSE41491.nodes$symbol <- row.names(expr.GSE41491.nodes)

expr.GSE41491.nodes.pivot <- expr.GSE41491.nodes %>%
  group_by(symbol) %>%
  pivot_longer(cols = starts_with(c("norm","hypo")), names_to = "codes", values_to = "value")

expr.GSE41491.nodes.pivot$codes <- factor(expr.GSE41491.nodes.pivot$codes,
                                          levels=c("norm.0h.","hypo.1h.", "hypo.2h.","hypo.4h.",
                                                   "hypo.8h.","hypo.12h.","hypo.16h.","hypo.24h."))

expr.GSE41491.nodes.pivot$time <- as.numeric(expr.GSE41491.nodes.pivot$codes)

p.MDA <- ggplot(aes(x = factor(time),y = value, group = symbol, color="red"),  
                data = expr.GSE41491.nodes.pivot[expr.GSE41491.nodes.pivot$symbol %in% netnodes.symbols,]) +
  geom_point() + 
  geom_line() + 
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7,8), 
                   labels = c("0h","1h","2h","4h","8h","12h","16h","24h")) +
  xlab("Conditions") + ylab("Gene Expression") +
  #ggtitle(plotTitle) +
  theme(legend.position = "none", axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8)) +
  #geom_line(aes(linetype=Symbol, color=Symbol)) +
  facet_wrap(~ symbol) 

p.MDA