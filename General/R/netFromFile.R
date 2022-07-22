require(igraph)
require(WGCNA)
require(BioNet)
require(readxl)

setwd('G:/Projects/PhD Thesis/Results/Network Analysis/')
fn <- list.files(path = getwd(),pattern = '*.xlsx')

netData <- read_excel(paste0(fn))
netFile <- graph_from_data_frame(netData, directed = T)
E(netFile)$weight <- abs(netData$weightMat)
# netData <- removeSelfLoops(netData)
netFile <- simplify(netFile, remove.multiple = T, remove.loops = T)

saveNetwork(netFile,name = paste0(gsub('.xlsx','',fn)), file = paste0(gsub('.xlsx','',fn),'.sif'), type = 'sif')

Metabolite <- V(netFile)$name
nodeDegree <- igraph::degree(netFile, v=V(netFile),mode = 'out')
nodeBetweenness <- igraph::betweenness(netFile, v=V(netFile))
nodeCloseness <- igraph::closeness(netFile, vids = V(netFile),mode = 'out')
netTop <- data.frame(cbind(Metabolite, 
                           as.numeric(nodeDegree),
                           as.numeric(nodeBetweenness),
                           as.numeric(nodeCloseness)))
colnames(netTop) <- c('Metabolite','Degree','Betweenness','Closeness')
netTop <- netTop[order(as.numeric(netTop$Betweenness), decreasing = T),]
writexl::write_xlsx(netTop, 'NetworkTopology.xlsx')
