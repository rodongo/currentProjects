

setwd('G:/Projects/PhD Thesis/Input Data/OMIC Data/Priori/')
load('endoPhenosMining.RData')
# Mining literature for disease endophenotypes
# AD

# DisGeNet Data
library(readr)
curatedGDAssoc <- read_table2('G:/Thesis Data/OMIC Data/Priori/curated_gene_disease_associations.tsv')

ADGenes <- curatedGDAssoc[which((grepl('alz',curatedGDAssoc$diseaseName,ignore.case = T)|
                                   grepl('alz',curatedGDAssoc$diseaseType,ignore.case = T))==T),]

endos <- c('tau','amyloid','immun','inflammation',
           'oxidat','oxidative stress','free radicals','lysosom','vascular','endothelial',
           'mito','mitochond','proteostasis','protein degradation','synap')
endosNames <- c('Tauopathy','Amyloidopathy','Neuroinflammation','Neuroinflammation2',
           'Oxidative Stress','Oxidative Stress2','Oxidative Stress3','Lysosomal Dysfunction','Vascular Dysfunction','Vascular Dysfunction2',
           'Mitochondrial Dysfunction','Mitochondrial Dysfunction2','Proteastasis','Proteastasis2','Synaptic Dysfunction')

# 1: PubMed Mining
library(easyPubMed)
library(dplyr)
library(kableExtra)

# Get this from PUBMED
api_key <- 'adf37d133d67c8959656202115157db57407'

# my_query <- paste('[Alzheimer] AND [immun] AND [', ADGenes$geneSymbol[33],'] gene', sep = '')
# my_query <- get_pubmed_ids(my_query, api_key = api_key)
# 
# my_query <- 'Damiano Fantini[AU] AND "2018"[PDAT]'
# my_entrez_id <- get_pubmed_ids(my_query)
# my_abstracts_txt <- fetch_pubmed_data(my_entrez_id, format = "abstract")

ADGenes2 <- unique(ADGenes$geneSymbol)
endosList <- list()
for(g in 1:length(endos)){
  geneEvid <- matrix(data = NA, nrow = length(ADGenes2), ncol = 3)
  for(h in 1:length(ADGenes2)){
    print(ADGenes2[h])
    Sys.sleep(0.50)
    my_query <- paste('[Alzheimer] disease AND [',endos[g],']', 'AND[', ADGenes2[h],']gene', sep = ' ')
    my_query <- get_pubmed_ids(my_query, api_key = api_key)
    if(length(my_query$IdList)>0){
      geneEvid[h,1] <- paste(ADGenes2[h])
      geneEvid[h,2] <- paste(toString(my_query$IdList))
      geneEvid[h,3] <- length(my_query$IdList)
    }
    else{
      geneEvid[h,1] <- paste(ADGenes2[h])
      geneEvid[h,2] <- paste('No study found')
      geneEvid[h,3] <- 0
    }
  }
  geneEvid <- as.data.frame(geneEvid)
  colnames(geneEvid) <- c('AD Gene', 'PMID', 'Number of Publications')
  geneEvid <- geneEvid %>% dplyr::filter(geneEvid$PMID != 'No study found')
  endosList[[paste0(endos[g])]] <- geneEvid
  Sys.sleep(10)
}
names(endosList) <- endosNames

endosList %>% writexl::write_xlsx(path = 'AD_Endophen_PubMed.xlsx')

# 2: Gene Ontology
library(gprofiler2)
library(enrichR)


library(plyr)
dbs <- listEnrichrDbs()
Dbs <- c('HDSigDB_Human_2021',
         'GTEx_Aging_Signatures_2021',
         'PhenGenI_Association_2021',
         'MGI_Mammalian_Phenotype_Level_4_2021',
         'GO_Molecular_Function_2021',
         'GO_Cellular_Component_2021',
         'GO_Biological_Process_2021',
         'WikiPathway_2021_Human',
         'KEGG_2021_Human',
         'Allen_Brain_Atlas_10x_scRNA_2021',
         'PheWeb_2019')
enrichAnal <- enrichr(ADGenes2[33], Dbs)
enrichAnal <- lapply(enrichAnal,function(x){x[,c(1:3)]})
enrichAnal2 <- c()
for(i in 1:length(enrichAnal)){
  enrichAnal2 <- rbind(enrichAnal2,enrichAnal[[i]])
}
enrichAnal2 <- enrichAnal2[grepl('tau',enrichAnal2$Term, ignore.case = T),]


ADGenes2 <- unique(ADGenes$geneSymbol)
endosList2 <- list()
for(g in 2:length(endos)){
  geneEvid <- matrix(data = NA, nrow = length(ADGenes2), ncol = 2)
  for(h in 1:length(ADGenes2)){
    print(ADGenes2[h])
    Sys.sleep(0.50)
    enrichAnal <- enrichr(ADGenes2[h], Dbs)
    enrichAnal <- lapply(enrichAnal,function(x){x[,c(1:3)]})
    enrichAnal2 <- c()
    for(i in 1:length(enrichAnal)){
      enrichAnal2 <- rbind(enrichAnal2,enrichAnal[[i]])
    }
    if(nrow(enrichAnal2)!=0){
      enrichAnal2 <- enrichAnal2[grepl(paste0(endos[g]),enrichAnal2$Term, ignore.case = T),]
      nGO <- length(unique(enrichAnal2$Term))
      if(nrow(enrichAnal2)>0){
        geneEvid[h,1] <- paste(ADGenes2[h])
        geneEvid[h,2] <- paste(nGO)
      }
      else{
        geneEvid[h,1] <- paste(ADGenes2[h])
        geneEvid[h,2] <- paste('No study found')
      }
    }
    else{
      geneEvid[h,1] <- paste(ADGenes2[h])
      geneEvid[h,2] <- paste('No study found')
    }
    }
  geneEvid <- as.data.frame(geneEvid)
  colnames(geneEvid) <- c('AD Gene', '#GO Terms')
  geneEvid <- geneEvid %>% dplyr::filter(geneEvid$`#GO Terms` != 'No study found')
  endosList2[[paste0(endos[g])]] <- geneEvid
  Sys.sleep(10)
}
names(endosList2) <- endosNames

endosList2 %>% writexl::write_xlsx(path = 'AD_Endophen_ENRICHR.xlsx')

# g:Profiler

ADGenes2 <- unique(ADGenes$geneSymbol)
endosList3 <- list()
for(g in 2:length(endos)){
  geneEvid <- matrix(data = NA, nrow = length(ADGenes2), ncol = 2)
  for(h in 1:length(ADGenes2)){
    print(ADGenes2[h])
    Sys.sleep(0.50)
    gostres <- gost(query = ADGenes2[h],
                    # user_threshold = 0.05,
                    # correction_method = 'false_discovery_rate',
                    # custom_bg = as.character(modIDs$entrezgene_id),
                    organism = 'hsapiens')
    gostres2 <- gostres$result[,c(10,11)]
    
    if(!is.null(gostres2)){
      gostres2 <- gostres2[grepl(paste0(endos[g]),gostres2$term_name, ignore.case = T),]
      nGO <- length(unique(gostres2$term_name))
      if(nrow(gostres2)>0){
        geneEvid[h,1] <- paste(ADGenes2[h])
        geneEvid[h,2] <- paste(nGO)
      }
      else{
        geneEvid[h,1] <- paste(ADGenes2[h])
        geneEvid[h,2] <- paste('No study found')
      }
    }
    else{
      geneEvid[h,1] <- paste(ADGenes2[h])
      geneEvid[h,2] <- paste('No study found')
    }
  }
  geneEvid <- as.data.frame(geneEvid)
  colnames(geneEvid) <- c('AD Gene', '#GO Terms')
  geneEvid <- geneEvid %>% dplyr::filter(geneEvid$`#GO Terms` != 'No study found')
  endosList3[[paste0(endos[g])]] <- geneEvid
  Sys.sleep(10)
}
names(endosList3) <- endosNames

endosList3 %>% writexl::write_xlsx(path = 'AD_Endophen_gProfiler.xlsx')

# Mining GO
BiocManager::install('ontologyIndex')
BiocManager::install('ontologySimilarity')
library(ontologyIndex)
data(go)

library(ontologySimilarity)
data(gene_GO_terms)
names(gene_GO_terms)

endosList4 <- list()
for(k in 1:length(endos)){
  endoTerm <- endos[k]
  endosGO <- c()
  for(i in 1:length(ADGenes2)){
    if(ADGenes2[i] %in% names(gene_GO_terms)){
      # i <- 1
      GO_Terms <- gene_GO_terms[which(names(gene_GO_terms) == ADGenes2[i])]
      for(j in 1:length(GO_Terms[[1]])){
        # j <- 1
        GO_Term_Name <- go$name[go$id==GO_Terms[[1]][j]]
        if(grepl(paste0(endoTerm),GO_Term_Name, ignore.case = T)){
          endosGO <- c(endosGO,go$id[go$name == GO_Term_Name])
        }
      }
    }
  }
  endosList4[[paste0(endos[k])]] <- endosGO
}
# Identify associated genes
endosList5 <- list()
for(i in 1:length(endosList4)){
  # i <- 1
  endosI <- endosList4[[i]]
  assGenes <- c()
  for(j in 1:length(endosI[[1]])){
    # j <- 1
    goEndos <- endosI[[1]][j]
    for(k in 1:length(gene_GO_terms)){
      # k <- 1
      geneGoTerms <- gene_GO_terms[k][[1]]
      if(goEndos %in% geneGoTerms){
        assGenes <- c(assGenes, names(gene_GO_terms)[k])
      }
    }
  }
  endosList5[[paste0(names(endosList4)[i])]] <- assGenes
}
endosList42 <- lapply(endosList5,function(x){data.frame(x)})
endosList42 %>% writexl::write_xlsx(path = 'AD_Endophen_GO1.xlsx')
# Information from BioGRID
BiocManager::install('BioNet')
library(readr)
library(dplyr)
library(igraph)
library(BioNet)
BioGRID_Hsa <- read_delim('BioGRID-Hsa.txt', 
                          delim = '\t', escape_double = FALSE, 
                          trim_ws = TRUE)
phyInteract <- c('Affinity Capture-Luminescence','Affinity Capture-MS',
                 'Affinity Capture-RNA','Affinity Capture-Western',
                 'Biochemical Activity','Co-crystal Structure',
                 'Co-fractionation','Co-localization','Co-purification',
                 'Far Western','FRET','PCA','Protein-peptide','Protein-RNA',
                 'Proximity Label-MS','Reconstituted Complex','Two-hybrid')
BioGRID_Hsa2 <- BioGRID_Hsa %>% dplyr::filter(BioGRID_Hsa$EXPERIMENTAL_SYSTEM %in% phyInteract)
# BioGRID_Hsa2 <- BioGRID_Hsa2 %>% filter(BioGRID_Hsa2$ORGANISM_A_ID)
ppiNet <- cbind(as.character(BioGRID_Hsa2$OFFICIAL_SYMBOL_A), as.character(BioGRID_Hsa2$OFFICIAL_SYMBOL_B))

for (i in 1:nrow(ppiNet)){
  ppiNet[i, ] <- sort(ppiNet[i, ])
}
ppiNetSort <- ppiNet[!duplicated(ppiNet),]
ppiNetSort <- ftM2graphNEL(ppiNetSort, edgemode= 'undirected')

a <- ppiNetSort@edgeL
attributes(ppiNetSort@nodeData@data)<- a
saveNetwork(ppiNetSort, file = paste('BioGRID.Network.Nodes.txt', sep = ''), 
            name = paste('BioGRID.Network_network.sif',sep = ''), type = 'sif')

BioGRID.Net <- loadNetwork.sif('BioGRID.Network.Nodes.sif',
                                   format = 'igraph')
adNeigh <- c()
for(i in 1:length(ADGenes2)){
  neigh.nodes <- neighbors(BioGRID.Net, V(BioGRID.Net)$name==ADGenes2[i], mode='out')
  adNeigh <- c(adNeigh, as.character(as_ids(neigh.nodes)))
}
adNeigh <- unique(c(adNeigh,ADGenes2))

# length(intersect(ADGenes2,adNeigh))
endosList7 <- list()
k_vals <- c()
for(k in 1:length(endos)){
  endoTerm <- endos[k]
  endosGO <- c()
  for(i in 1:length(adNeigh)){
    if(adNeigh[i] %in% names(gene_GO_terms)){
      # i <- 1
      GO_Terms <- gene_GO_terms[which(names(gene_GO_terms) == adNeigh[i])]
      for(j in 1:length(GO_Terms[[1]])){
        # j <- 1
        GO_Term_Name <- go$name[go$id==GO_Terms[[1]][j]]
        if(grepl(paste0(endoTerm),GO_Term_Name, ignore.case = T)){
          endosGO <- c(endosGO,go$id[go$name == GO_Term_Name])
          k_vals <- c(k_vals,k)
        }
      }
    }
  }
  endosList7[[k]] <- endosGO
}
endosList8 <- endosList7 %>% keep( ~ !is.null(.) )
names(endosList8) <- endosNames[unique(k_vals)]
library(tidyverse)
# endosList8 <- endosList7 %>% keep( ~ !is.null(.) )
# Identify associated genes
endosList71 <- list()
for(i in 1:length(endosList8)){
  # i <- 1
  endosI <- endosList8[[i]]
  assGenes <- c()
  for(j in 1:length(endosI[[1]])){
    # j <- 1
    goEndos <- endosI[[1]][j]
    for(k in 1:length(gene_GO_terms)){
      # k <- 1
      geneGoTerms <- gene_GO_terms[k][[1]]
      if(goEndos %in% geneGoTerms){
        assGenes <- c(assGenes, names(gene_GO_terms)[k])
      }
    }
  }
  endosList71[[paste0(names(endosList8)[i])]] <- assGenes
}
endosList71 <- lapply(endosList71,function(x){data.frame(x)})
endosList72 <- list('Tauopathy' = data.frame(GeneID = endosList71$Tauopathy$x),
                    'Amyloidopathy' = data.frame(GeneID = endosList71$Amyloidopathy$x),
                    'Neuroinflammation' = data.frame(GeneID = endosList71$Neuroinflammation$x),
                    'Oxidative Stress' = data.frame(GeneID = endosList71$`Oxidative Stress`$x),
                    'Lysosomal Dysfunction' = data.frame(GeneID = endosList71$`Lysosomal Dysfunction`$x),
                    'Vascular Dysfunction' = data.frame(GeneID = unique(c(endosList71$`Vascular Dysfunction2`$x,
                                                                          endosList71$`Vascular Dysfunction`$x))),
                    'Mitochondrial Dysfunction' = data.frame(GeneID = endosList71$`Mitochondrial Dysfunction`),
                    'Synaptic Dysfunction' = data.frame(GeneID = endosList71$`Synaptic Dysfunction`$x))
endosList72 %>% writexl::write_xlsx(path = 'AD_Endophen_GO2.xlsx')
saveRDS(endosList72,'.RDS')
#########################################################################################################################################
# PubMed Again!

adNeigh2 <- toupper(adNeigh)
# adNeigh2[nchar(adNeigh2)<2] <- NULL

library(stringr)
# adNeigh2 <- str_extract_all(adNeigh2, '\\w{2,}')[[1]]
adNeigh2 <- adNeigh2[grepl('\\w{2,}',adNeigh2)==T]

ADGenes2 <- unique(ADGenes$geneSymbol)
endosList2 <- list()
for(g in 1:length(endos)){
  geneEvid <- data.frame(matrix(data = NA, nrow = length(adNeigh2), ncol = 3))
  colnames(geneEvid) <- c('AD.Gene', 'PMID', 'Frequency')
  for(h in 1:length(adNeigh2)){
    print(adNeigh2[h])
    # Sys.sleep(0.50)
    my_query <- paste('[Alzheimer] disease AND [',endos[g],']', 'AND[', adNeigh2[h],']gene', sep = ' ')
    my_query <- get_pubmed_ids(my_query, api_key = api_key)
    if(length(my_query$IdList)>0){
      geneEvid[h,1] <- paste(adNeigh2[h])
      geneEvid[h,2] <- paste(toString(my_query$IdList))
      geneEvid[h,3] <- length(my_query$IdList)
    }
    else{
      geneEvid[h,1] <- paste(adNeigh2[h])
      geneEvid[h,2] <- paste('No study found')
      geneEvid[h,3] <- 0
    }
  }
  geneEvid <- geneEvid %>% dplyr::filter(geneEvid$PMID != 'No study found')
  endosList2[[paste0(endos[g])]] <- geneEvid
  # Sys.sleep(10)
}
endosList2.1 <- endosList2[c(12:26)]
names(endosList2.1) <- endosNames

endosList22 <- list('Tauopathy' = endosList2.1$Tauopathy,
                    'Amyloidopathy' = endosList2.1$Amyloidopathy,
                    'Neuroinflammation' = rbind(endosList2.1$Neuroinflammation,endosList2.1$Neuroinflammation2),
                    'Oxidative Stress' = rbind(endosList2.1$`Oxidative Stress`,endosList2.1$`Oxidative Stress2`,endosList2.1$`Oxidative Stress3`),
                    'Lysosomal Dysfunction' = endosList2.1$`Lysosomal Dysfunction`,
                    'Vascular Dysfunction' = rbind(endosList2.1$`Vascular Dysfunction`,endosList2.1$`Vascular Dysfunction2`),
                    'Mitochondrial Dysfunction' = rbind(endosList2.1$`Mitochondrial Dysfunction`,endosList2.1$`Mitochondrial Dysfunction2`),
                    'Proteastasis' = rbind(endosList2.1$Proteastasis,endosList2.1$Proteastasis2),
                    'Synaptic Dysfunction' = endosList2.1$`Synaptic Dysfunction`)

endosList22.1 <- lapply(endosList22, function(x){x = x %>% filter(.$Frequency>4)})
endosList22.1 <- lapply(endosList22.1, function(x){x = x[!duplicated(x$AD.Gene),]})
endosList22.1 <- endosList22.1[-c(5,7,9)]
# endosList22.2 <- list()
# for(i in 1:length(endosList22.1)){
#   df <- endosList22.1[[i]]
#   colnames(df)[1] <- 'GeneID'
#   df[,'nChar'] <- nchar(df$GeneID)
#   endosList22.2[[paste0(names(endosList22.1)[i])]] <- df[df$nChar!=1,] %>% select(-nChar)
# }
# 
# endosList22.1 <- lapply(endosList22.2, function(x){x = x %>% mutate(Evidence = as.numeric(x$`Number of Publications`)*100/max(as.numeric(x$`Number of Publications`)),
#                                                                   Gene = toupper(x$GeneID)) %>% select(-GeneID)})
# # endosList22.1 <- lapply(endosList22.1, function(x){x = x %>% filter(.$Evidence > 30)})
# endosList22.1 <- lapply(endosList22.1, function(x){x = x[!duplicated(x$Gene),] %>% select(Gene, Evidence)})
# 


endosList22 %>% writexl::write_xlsx(path = 'AD_Endophen_PubMed2.xlsx')
endosList22.1 %>% writexl::write_xlsx(path = 'G:/Projects/PhD Thesis/Results/Endophenotypes/AD_Endophenotype2Gene_2.xlsx')

endosList22.2 <- lapply(endosList22.1, function(x){x = data.frame(GeneID = x$AD.Gene)})
saveRDS(endosList22.1, 'G:/Projects/PhD Thesis/Results/Endophenotypes/AD_Endophenotype2Gene_2.RDS')
modIDs <- readRDS('G:/ParkinsonsDatasets/Combined Datasets/modIDs.RDS')
# endosList2$Tauopathy$`AD Gene`
Intersect <- function(inList){
  x <- crossprod(table(stack(inList)))
  x[lower.tri(x)] <- NA
  diag(x) <- NA
  out <- na.omit(data.frame(as.table(x)))
  out[order(out$ind),]
}
endosList22.3 <- lapply(endosList22.1, function(x){x = x$AD.Gene})
endoSim <- Intersect(endosList22.3)
colnames(endoSim) <- c('Endophenotype.1','Endophenotype.2','NumberInCommon')

writexl::write_xlsx(endoSim, path = 'G:/Projects/PhD Thesis/Results/Endophenotypes/EndophenotypeGeneSimilarity.xlsx')
# endosList3 <- lapply(endosList22, function(x){length(intersect(x$`AD Gene`,modIDs$hgnc_symbol))})

save.image('endoPhenosMining.RData')
