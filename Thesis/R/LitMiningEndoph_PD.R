
setwd('G:/Projects/PhD Thesis/Input Data/OMIC Data/Priori/')
load('endoPhenosMining_PD.RData')
# Mining literature for disease endophenotypes
# AD

# DisGeNet Data
library(readr)
curatedGDAssoc <- read_table2('G:/Thesis Data/OMIC Data/Priori/curated_gene_disease_associations.tsv')

PDGenes <- curatedGDAssoc[which((grepl('park',curatedGDAssoc$diseaseName,ignore.case = T)|
                                   grepl('park',curatedGDAssoc$diseaseType,ignore.case = T))==T),]

endos <- c('synuc','immun','inflammation',
           'oxidat','lysosom',
           'mito','proteostasis','synap')
endosNames <- c('Synucleinopathy','Neuroinflammation','Neuroinflammation2',
                'Oxidative Stress','Lysosomal Dysfunction',
                'Mitochondrial Dysfunction','Proteastasis','Synaptic Dysfunction')

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

PDGenes2 <- unique(PDGenes$geneSymbol)
endosList <- list()
for(g in 1:length(endos)){
  geneEvid <- matrix(data = NA, nrow = length(PDGenes2), ncol = 3)
  for(h in 1:length(PDGenes2)){
    print(PDGenes2[h])
    Sys.sleep(0.50)
    my_query <- paste('[Alzheimer] disease AND [',endos[g],']', 'AND[', PDGenes2[h],']gene', sep = ' ')
    my_query <- get_pubmed_ids(my_query, api_key = api_key)
    if(length(my_query$IdList)>0){
      geneEvid[h,1] <- paste(PDGenes2[h])
      geneEvid[h,2] <- paste(toString(my_query$IdList))
      geneEvid[h,3] <- length(my_query$IdList)
    }
    else{
      geneEvid[h,1] <- paste(PDGenes2[h])
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
enrichAnal <- enrichr(PDGenes2[33], Dbs)
enrichAnal <- lapply(enrichAnal,function(x){x[,c(1:3)]})
enrichAnal2 <- c()
for(i in 1:length(enrichAnal)){
  enrichAnal2 <- rbind(enrichAnal2,enrichAnal[[i]])
}
enrichAnal2 <- enrichAnal2[grepl('tau',enrichAnal2$Term, ignore.case = T),]


PDGenes2 <- unique(PDGenes$geneSymbol)
endosList2 <- list()
for(g in 2:length(endos)){
  geneEvid <- matrix(data = NA, nrow = length(PDGenes2), ncol = 2)
  for(h in 1:length(PDGenes2)){
    print(PDGenes2[h])
    Sys.sleep(0.50)
    enrichAnal <- enrichr(PDGenes2[h], Dbs)
    enrichAnal <- lapply(enrichAnal,function(x){x[,c(1:3)]})
    enrichAnal2 <- c()
    for(i in 1:length(enrichAnal)){
      enrichAnal2 <- rbind(enrichAnal2,enrichAnal[[i]])
    }
    if(nrow(enrichAnal2)!=0){
      enrichAnal2 <- enrichAnal2[grepl(paste0(endos[g]),enrichAnal2$Term, ignore.case = T),]
      nGO <- length(unique(enrichAnal2$Term))
      if(nrow(enrichAnal2)>0){
        geneEvid[h,1] <- paste(PDGenes2[h])
        geneEvid[h,2] <- paste(nGO)
      }
      else{
        geneEvid[h,1] <- paste(PDGenes2[h])
        geneEvid[h,2] <- paste('No study found')
      }
    }
    else{
      geneEvid[h,1] <- paste(PDGenes2[h])
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

PDGenes2 <- unique(PDGenes$geneSymbol)
endosList3 <- list()
for(g in 2:length(endos)){
  geneEvid <- matrix(data = NA, nrow = length(PDGenes2), ncol = 2)
  for(h in 1:length(PDGenes2)){
    print(PDGenes2[h])
    Sys.sleep(0.50)
    gostres <- gost(query = PDGenes2[h],
                    # user_threshold = 0.05,
                    # correction_method = 'false_discovery_rate',
                    # custom_bg = as.character(modIDs$entrezgene_id),
                    organism = 'hsapiens')
    gostres2 <- gostres$result[,c(10,11)]
    
    if(!is.null(gostres2)){
      gostres2 <- gostres2[grepl(paste0(endos[g]),gostres2$term_name, ignore.case = T),]
      nGO <- length(unique(gostres2$term_name))
      if(nrow(gostres2)>0){
        geneEvid[h,1] <- paste(PDGenes2[h])
        geneEvid[h,2] <- paste(nGO)
      }
      else{
        geneEvid[h,1] <- paste(PDGenes2[h])
        geneEvid[h,2] <- paste('No study found')
      }
    }
    else{
      geneEvid[h,1] <- paste(PDGenes2[h])
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
  for(i in 1:length(PDGenes2)){
    if(PDGenes2[i] %in% names(gene_GO_terms)){
      # i <- 1
      GO_Terms <- gene_GO_terms[which(names(gene_GO_terms) == PDGenes2[i])]
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
for(i in 1:length(PDGenes2)){
  neigh.nodes <- neighbors(BioGRID.Net, V(BioGRID.Net)$name==PDGenes2[i], mode='out')
  adNeigh <- c(adNeigh, as.character(as_ids(neigh.nodes)))
}
adNeigh <- unique(adNeigh)

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
# PubMed Again
PDGenes2 <- unique(PDGenes$geneSymbol)
endosList2 <- list()
for(g in 1:length(endos)){
  geneEvid <- matrix(data = NA, nrow = length(adNeigh), ncol = 3)
  for(h in 1:length(adNeigh)){
    print(adNeigh[h])
    Sys.sleep(0.50)
    my_query <- paste('[Parkinson] disease AND [',endos[g],']', 'AND[', adNeigh[h],']gene', sep = ' ')
    my_query <- get_pubmed_ids(my_query, api_key = api_key)
    if(length(my_query$IdList)>0){
      geneEvid[h,1] <- paste(adNeigh[h])
      geneEvid[h,2] <- paste(toString(my_query$IdList))
      geneEvid[h,3] <- length(my_query$IdList)
    }
    else{
      geneEvid[h,1] <- paste(adNeigh[h])
      geneEvid[h,2] <- paste('No study found')
      geneEvid[h,3] <- 0
    }
  }
  geneEvid <- as.data.frame(geneEvid)
  colnames(geneEvid) <- c('AD Gene', 'PMID', 'Number of Publications')
  geneEvid <- geneEvid %>% dplyr::filter(geneEvid$PMID != 'No study found')
  endosList2[[paste0(endos[g])]] <- geneEvid
  Sys.sleep(10)
}
names(endosList2) <- endosNames

endosList22 <- list('Tauopathy' = endosList2$Tauopathy,
                    'Amyloidopathy' = endosList2$Amyloidopathy,
                    'Neuroinflammation' = rbind(endosList2$Neuroinflammation,endosList2$Neuroinflammation2),
                    'Oxidative Stress' = endosList2$`Oxidative Stress`,
                    'Lysosomal Dysfunction' = endosList2$`Lysosomal Dysfunction`,
                    'Vascular Dysfunction' = rbind(endosList2$`Vascular Dysfunction`,endosList2$`Vascular Dysfunction2`),
                    'Mitochondrial Dysfunction' = endosList2$`Mitochondrial Dysfunction`,
                    'Proteastasis' = endosList2$Proteastasis,
                    'Synaptic Dysfunction' = endosList2$`Synaptic Dysfunction`)

endosList22 %>% writexl::write_xlsx(path = 'AD_Endophen_PubMed2.xlsx')
modIDs <- readRDS('G:/ParkinsonsDatasets/Combined Datasets/modIDs.RDS')
# endosList2$Tauopathy$`AD Gene`
endosList3 <- lapply(endosList22, function(x){length(intersect(x$`AD Gene`,modIDs$hgnc_symbol))})
saveRDS(endosList3,'')
save.image('endoPhenosMining_PD.RData')
