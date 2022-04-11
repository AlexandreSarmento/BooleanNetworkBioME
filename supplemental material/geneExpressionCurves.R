packages_cran = c("pheatmap","ggplot2","BiTrinA","BiocManager", "tidyverse", "fs")
# Install and load packages
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
packages_bioconductor = c("Biobase","vsn", "hgu133plus2.db")
# Install and load packages
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("RNAAgeCalc",version = "1.4")

setwd("/storages/caico/home/alexandresarmento/R")
load("/storages/caico/home/alexandresarmento/R/DataBases/data.GSE41491.RData")

normalize = function(x){
  ans = (x - min(x)) / (max(x) - min(x))
  return(ans)
}

nodes.symbols <- c("AKT1","ATR","BAX","EGLN1","HIF1A","PIK3CA","GSK3B",
                   "HIF1AN","MDM2","PFKL","PPP1R10","PTEN",'EP300',
                   "TP53","VEGFA")

time_points <- c(0,1,2,4,8,12,16,24)

nodes.probes <- anno.GSE41491$probes[anno.GSE41491$symbol %in% nodes.symbols]

expr.GSE41491.nodes <- as.data.frame(expr.GSE41491) %>% 
  rownames_to_column('probes') %>% 
  filter(probes %in% nodes.probes) %>% 
  merge(anno.GSE41491[anno.GSE41491$symbol %in% nodes.symbols, c("probes","symbol")], by = "probes") %>% 
  group_by(symbol) %>%
  summarise_at(vars(-probes), funs(median(., na.rm=TRUE)))  %>%
  column_to_rownames(var = "symbol") %>%
  dplyr::select(c(data.GSE41491$codes[data.GSE41491$cell_line == "MCF7"]))

row <- data.GSE41491$cell_line == "MCF7"
expr.GSE41491.nodes <- expr.GSE41491.nodes[,c(as.character(data.GSE41491$codes[row]))]
names(expr.GSE41491.nodes) <- c(paste0(data.GSE41491$condition[row],".", data.GSE41491$time[row]))
write.csv(t(expr.GSE41491.nodes),"/storages/caico/home/alexandresarmento/R/MCF7_.csv")
expr.GSE41491.nodes <- normalize(expr.GSE41491.nodes)
expr.GSE41491.nodes$symbol <- row.names(expr.GSE41491.nodes)


expr.GSE41491.nodes.pivot <- expr.GSE41491.nodes %>%
  group_by(symbol) %>%
  pivot_longer(cols = starts_with(c("norm","hypo")), names_to = "codes", values_to = "value")


expr.GSE41491.nodes.pivot$codes <- factor(expr.GSE41491.nodes.pivot$codes,
                                          levels = c("normoxia.0h","hypoxia.1h","hypoxia.2h","hypoxia.4h",
                                                     "hypoxia.8h","hypoxia.12h","hypoxia.16h","hypoxia.24h"))


expr.GSE41491.nodes.pivot$time <- rep(time_points,times = 15)


p.MDA <- ggplot(aes(x = factor(time), y = value, group = symbol, color="red"),  
                data = expr.GSE41491.nodes.pivot[expr.GSE41491.nodes.pivot$symbol %in% nodes.symbols,]) +
  geom_point() + 
  geom_line() + 
  scale_x_discrete(breaks = time_points , 
                   labels = c("0","1","2","4","8","12","16","24")) +
  xlab("Time (hour)") + ylab("Gene Expression") +
  ggtitle("MCF7") +
  theme(legend.position = "none", axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8)) +
  #geom_line(aes(linetype=Symbol, color=Symbol)) +
  facet_wrap(~ symbol) 

p.MDA
