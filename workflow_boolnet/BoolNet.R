#' ---
#' title: "BoolNet"
#' 
#' output: 
#'   pdf_document: default
#'   html_document: default
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
## -----------------------------------------------------------------------------
library(BoolNet)
setwd("/data4/terrematte/modelling/BoolNet")

#' 
#' 
#' # Boolean network from HIFaxis 10/09/2020
#' 
#' 
## -----------------------------------------------------------------------------
net <- loadNetwork("boolean_network_HIFaxis.bn")
net

#' 
#' 
## ----fig.height=6, fig.width=12-----------------------------------------------
attr.syn <- getAttractors(net, type = "synchronous")

# calculate number of different attractor lengths, 
# and plot attractors side by side in "table" mode
par(mfrow=c(1, length(table(sapply(attr.syn$attractors,
                          function(attr.syn)
                          {
                            length(attr.syn$involvedStates)
                          })))))

plotAttractors(attr.syn)

# plot attractors in "graph" mode
par(mfrow=c(1, length(attr.syn$attractors)))
plotAttractors(attr.syn, mode="graph")

# identify asynchronous attractors
attr.asyn <- getAttractors(net, type="asynchronous")

# plot attractors in "graph" mode
par(mfrow=c(1, length(attr.asyn$attractors)))
plotAttractors(attr.asyn, mode="graph")

#' 
## -----------------------------------------------------------------------------
plotStateGraph(attr.syn)

#' 
#' 
## -----------------------------------------------------------------------------
sim <- markovSimulation(net,
                        numIterations=1024,
                        returnTable=FALSE)
sim

#' 
## -----------------------------------------------------------------------------
plotNetworkWiring(net)

#' 
#' 
#' # Boolean network from 09/07/2020
#' 
#' 
## -----------------------------------------------------------------------------
net <- loadNetwork("boolean_network_2020_07_09.bn")
net

#' 
## ----fig.height=6, fig.width=12-----------------------------------------------
attr.syn <- getAttractors(net, type = "synchronous")

# calculate number of different attractor lengths, 
# and plot attractors side by side in "table" mode
par(mfrow=c(1, length(table(sapply(attr.syn$attractors,
                          function(attr.syn)
                          {
                            length(attr.syn$involvedStates)
                          })))))

plotAttractors(attr.syn)

# plot attractors in "graph" mode
par(mfrow=c(1, length(attr.syn$attractors)))
plotAttractors(attr.syn, mode="graph")

# # identify asynchronous attractors
# attr.asyn <- getAttractors(net, type="asynchronous")
# 
# # plot attractors in "graph" mode
# par(mfrow=c(1, length(attr.asyn$attractors)))
# plotAttractors(attr.asyn, mode="graph")

#' 
#' 
#' <!-- ```{r} -->
#' <!-- plotStateGraph(attr.syn) -->
#' <!-- ``` -->
#' 
#' 
## -----------------------------------------------------------------------------
sim <- markovSimulation(net,
                        numIterations=1024,
                        returnTable=FALSE)
sim

#' 
## -----------------------------------------------------------------------------
plotNetworkWiring(net)

#' 
#' 
#' # Boolean network from 10/06/2020
#' 
#' 
## -----------------------------------------------------------------------------
net <- loadNetwork("boolean_network_2020_06_10.bn")
net

#' 
## ----fig.height=6, fig.width=12-----------------------------------------------
attr.syn <- getAttractors(net, type = "synchronous")

# calculate number of different attractor lengths, 
# and plot attractors side by side in "table" mode
par(mfrow=c(1, length(table(sapply(attr.syn$attractors,
                          function(attr.syn)
                          {
                            length(attr.syn$involvedStates)
                          })))))

plotAttractors(attr.syn)

# # plot attractors in "graph" mode
# par(mfrow=c(1, length(attr.syn$attractors)))
# plotAttractors(attr.syn, mode="graph")

# # identify asynchronous attractors
# attr.asyn <- getAttractors(net, type="asynchronous")
# 
# # plot attractors in "graph" mode
# par(mfrow=c(1, length(attr.asyn$attractors)))
# plotAttractors(attr.asyn, mode="graph")

#' 
#' <!-- ```{r} -->
#' <!-- plotStateGraph(attr.syn) -->
#' <!-- ``` -->
#' 
#' 
#' <!-- ```{r} -->
#' <!-- sim <- markovSimulation(net, -->
#' <!--                         numIterations=1024, -->
#' <!--                         returnTable=FALSE) -->
#' <!-- sim -->
#' <!-- ``` -->
#' 
## -----------------------------------------------------------------------------
plotNetworkWiring(net)

#' 
#' 
#' 
