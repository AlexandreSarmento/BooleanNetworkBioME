packages_cran = c("BiocManager", "tidyverse", "fs")
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

setwd("/storages/caico/home/alexandresarmento/R/figures")
load("/storages/caico/home/alexandresarmento/R/DataBases/data.GSE41491.RData")

normalize = function(x){
  ans = (x - min(x)) / (max(x) - min(x))
  return(ans)
}

nodes.symbols <- c("CASP3","EP300","BBC3","BNIP3")

nodes.nickname <- c("CASP3","EP300","PUMA","BNIP3")
  
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
expr.GSE41491.nodes <- normalize(expr.GSE41491.nodes)
#expr.GSE41491.nodes$symbol <- row.names(expr.GSE41491.nodes)

experimentalTimePoints <- c(0,1,2,4,8,12,16,24) 
allTimePoints <- as.numeric(0:24) 

akt <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[1,], method = "natural")
AKT = akt(allTimePoints)

atr <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[2,], method = "natural")
ATR = atr(allTimePoints)

bax <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[3,], method = "natural")
BAX = bax(allTimePoints)

# puma <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[1,], method = "natural") 
# PUMA = puma(allTimePoints) 
# 
# bnip3 <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[2,], method = "natural") 
# BNIP3 = bnip3(allTimePoints) 
# 
# casp3 <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[3,], method = "natural") 
# CASP3 = casp3(allTimePoints) 

phd <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[7,], method = "natural")
PHD = phd(allTimePoints)

hif1a <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[8,], method = "natural")
HIF1A = hif1a(allTimePoints)

fih <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[9,], method = "natural")
FIH = fih(allTimePoints)

# p300 <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[4,], method = "natural") 
# EP300 = p300(allTimePoints) 

mdm2 <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[11,], method = "natural")
MDM2 = mdm2(allTimePoints)

pfkl <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[12,], method = "natural")
PFKL = pfkl(allTimePoints)

pnuts <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[13,], method = "natural")
PNUTS = pnuts(allTimePoints)

pten <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[14,], method = "natural")
PTEN = pten(allTimePoints)

p53 <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[15,], method = "natural")
TP53 = p53(allTimePoints)

vegf <- splinefun(x = experimentalTimePoints, y = expr.GSE41491.nodes[16,], method = "natural")
VEGF = vegf(allTimePoints)

# AKT,ATR,BAX,PUMA,BNIP3,CASP3,PHD,HIF1A,FIH,EP300,MDM2,PFKL,PNUTS,PTEN,TP53,VEGF
expr.GSE41491.nodes.interp <- rbind(PUMA,BNIP3,CASP3,EP300)
expr.GSE41491.nodes.interp <- data.frame(expr.GSE41491.nodes.interp)
colnames(expr.GSE41491.nodes.interp) <- c("norm.0h","hypo.1h","hypo.2h","hypo.3h","hypo.4h","hypo.5h",
                                          "hypo.6h","hypo.7h","hypo.8h","hypo.9h","hypo.10h","hypo.11h",
                                          "hypo.12h","hypo.13h","hypo.14h","hypo.15h","hypo.16h","hypo.17h",
                                          "hypo.18h","hypo.19h","hypo.20h","hypo.21h","hypo.22h","hypo.23h",
                                          "hypo.24h")
expr.GSE41491.nodes.interp$symbol <- row.names(expr.GSE41491.nodes.interp)


expr.GSE41491.nodes.pivot <- expr.GSE41491.nodes.interp %>%
  group_by(symbol) %>%
  pivot_longer(cols = starts_with(c("norm","hypo")), names_to = "codes", values_to = "value")


expr.GSE41491.nodes.pivot$codes <- factor(expr.GSE41491.nodes.pivot$codes,
                                          levels = c("norm.0h","hypo.1h","hypo.2h","hypo.3h","hypo.4h","hypo.5h",
                                                     "hypo.6h","hypo.7h","hypo.8h","hypo.9h","hypo.10h","hypo.11h",
                                                     "hypo.12h","hypo.13h","hypo.14h","hypo.15h","hypo.16h","hypo.17h",
                                                     "hypo.18h","hypo.19h","hypo.20h","hypo.21h","hypo.22h","hypo.23h",
                                                     "hypo.24h"))


expr.GSE41491.nodes.pivot$time <- as.numeric(expr.GSE41491.nodes.pivot$codes)


p.MDA <- ggplot(aes(x = factor(time), y = value, group = symbol, color="red"),  
                data = expr.GSE41491.nodes.pivot[expr.GSE41491.nodes.pivot$symbol %in% nodes.nickname,]) +
  geom_point() + 
  geom_line() + 
  scale_x_discrete(breaks = as.numeric(0:24), 
                   labels = c("0","1","2","3","4","5","6","7","8","9","10",
                              "11","12","13","14","15","16","17","18","19",
                              "20","21","22","23","24")) +
  xlab("Time (hour)") + ylab("Gene Expression") +
  ggtitle("MCF7") +
  theme(legend.position = "none", axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8)) +
  #geom_line(aes(linetype=Symbol, color=Symbol)) +
  facet_wrap(~ symbol) 

p.MDA
