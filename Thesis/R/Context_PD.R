
library(igraph)
library(BioNet)
library(dplyr)

library(readr)
curatedGDAssoc <- read_delim('G:/Projects/PhD Thesis/Input Data/OMIC Data/Priori/Disease Priors/curated_gene_disease_associations.tsv', 
                                delim = '\t', escape_double = FALSE, 
                                trim_ws = TRUE)
curatedGDAssoc <- read_table2('G:/Projects/PhD Thesis/Input Data/OMIC Data/Priori/curated_gene_disease_associations.tsv')
BioGRID.Net <- loadNetwork.sif('G:/Projects/PhD Thesis/Input Data/OMIC Data/Priori/BioGRID/BioGRID.Network.Nodes.sif',
                               format = 'igraph')

curatedGDAssoc <- curatedGDAssoc[which((grepl('park',curatedGDAssoc$diseaseName,ignore.case = T)|
                                             grepl('park',curatedGDAssoc$diseaseType,ignore.case = T))==T),]$geneSymbol
curatedGDAssoc <- unique(curatedGDAssoc)

cluster_network <- cluster_walktrap(BioGRID.Net, steps = 4)
# Get the length of the modules

# We only consider those cluster with at least 10 genes in it
cluster_of_interest <- cluster_network[(sizes(cluster_network) >= 20)]

### Test if the cluster of interest include more disease genes than expected by chance

FisherTests <- data.frame(matrix(data = NA, nrow = length(cluster_of_interest), ncol = 2))
colnames(FisherTests) <- c('Cluster','FishersPValue')
FisherTests$Cluster <- names(cluster_of_interest)

for(i in 1:length(cluster_of_interest)){
  # N = size of the network, N_disease = Size of disease genes in the network, N_cluster = size of the cluster selected
  N <- length(V(BioGRID.Net)$name)
  N_disease <-  length(intersect(curatedGDAssoc, V(BioGRID.Net)$name))
  N_cluster <- length(cluster_of_interest[[i]])
  
  # Contingency table: cluster&disease, cluster&not_disease, not_cluster&disease, not_cluster&not_disease
  cluster_disease <- length(unique(intersect(curatedGDAssoc, cluster_of_interest[[i]])))
  cluster_not_disease <- N_cluster - cluster_disease
  not_cluster_disease <- N_disease - cluster_disease
  not_cluster_not_disease <- N - cluster_disease - cluster_not_disease - not_cluster_disease
  
  cont_table <- rbind(c(cluster_disease, cluster_not_disease),
                      c(not_cluster_disease, not_cluster_not_disease))
  
  # Fisher's exact test
  FisherTests[i,2] <- fisher.test(cont_table, alternative = "two.sided")
}
FisherTests <- FisherTests %>% filter(.$FishersPValue<0.05)
cluster_of_interest2 <- cluster_of_interest[names(cluster_of_interest)==FisherTests$Cluster]

library(purrr)
not_rel <- setdiff(V(BioGRID.Net)$name , cluster_of_interest2$`8`)
BioGRID.Net_1 <- delete.vertices(BioGRID.Net,not_rel)
saveNetwork(BioGRID.Net_1,name = 'BioGRID.Net.PD', file = 'BioGRID.Net.PD.sif', type = 'sif')

load('DrugModules.RData')
library(BioNet)
library(igraph)
library(tidyverse)
library(igraph)
library(BioNet)
library(dplyr)

library(readr)
# Load drug signature files
drug.signature <- readRDS('G:/Projects/TUBITAK/GSE164788_DGE_Sym.RDS')
# Load background network
back.Net <- loadNetwork.sif(sif.file = 'G:/Projects/PhD Thesis/Input Data/OMIC Data/Priori/BioGRID/BioGRID.Network.Nodes.sif')
df <- drug.signature$PVal %>% select(contains('0.3')) # Here we use 0.3uM dosage
cmpModules <- list()
for(i in 1:length(df)){
  print(colnames(df)[i])
  pVals <- df[,1] %>% as.data.frame()
  rownames(pVals) <- rownames(df) 
  colnames(pVals) <-  'pvals'
  pVals<- aggrPvals(as.matrix(pVals), order = 1, plot = FALSE) 
  
  subnet <- subNetwork(rownames(df), back.Net) 
  subnet <- rmSelfLoops(subnet) 
  
  fb <- fitBumModel(pVals, plot = F) 
  scores <- scoreNodes(network =  subnet, fb = fb, fdr = 0.005)
  module <- runFastHeinz(subnet, scores) 
  cmpModules[[paste0(colnames(df)[i])]] <- module

  saveNetwork(module, file = paste(colnames(df)[i],'_nodes', sep = ''), name = paste(colnames(df)[i],'_network.sif',sep = ''), type = 'sif')
}
save.image('G:/Projects/TUBITAK/DrugModules.RData')
