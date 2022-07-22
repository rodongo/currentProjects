BiocManager::install('EBSeq')
BiocManager::install('edgeR')
BiocManager::install('GEOquery')
BiocManager::install('biomaRt')
BiocManager::install('limma')
BiocManager::install('factoextra')
BiocManager::install('gprofiler2')
BiocManager::install('WGCNA')
BiocManager::install('igraph')
BiocManager::install('GO.db')
BiocManager::install('impute')
BiocManager::install('preprocessCore')


library(edgeR)
library(GEOquery)
library(readxl)
library(biomaRt)
library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(WGCNA)
library(igraph)
library(tictoc)
library(PCAtools)
library(factoextra)
library(reshape2)
library(gprofiler2)
library(EBSeq)
library(limma)
setwd('G:/Projects/PhD Thesis/AD/FreshMicro Study/')
# load('ThesisEpigenetics.RData')
load('ThesisEpigenetics-3.RData')

source('F:/Functions/Thesis/R/geTMM.R')
source('F:/Functions/Thesis/R/pcaPlot.R')
source('F:/Functions/Thesis/R/plotLoadings3.R')


# AD Microglia Data
# Metadata
FreshMicro_individual_metadata <- read_csv('FreshMicro_individual_metadata.csv')
AD_Samples <- FreshMicro_individual_metadata %>% filter(.$diagnosis=='Alzheimer Disease')
CON_Samples <- FreshMicro_individual_metadata %>% filter(.$diagnosis=='control')

remSamples <- FreshMicro_individual_metadata[as.character(FreshMicro_individual_metadata$individualID)==c(FM_082,FM_099),] 


# %>% 
#   as.data.frame() %>% 
#   dplyr::filter(.$individualID ==c(FM_082,FM_099))
# ATAC-Seq data
# ATACseq_count_matrix$FM_002__ATACseq
ATACseq_count_matrix <- read_csv('ATACseq_count_matrix.csv')
ATACseq_PeakInfo <- read_csv('ATACseq_PeakInfo.csv')
colnames(ATACseq_count_matrix) <- gsub('__ATACseq','',colnames(ATACseq_count_matrix))
colnames(ATACseq_count_matrix)[1] <- 'PeakID'
FreshMicro_assay_ATACseq_metadata <- read_csv('FreshMicro_assay_ATACseq_metadata.csv')


# RNA-Seq data
FreshMicro_assay_RNAseq_metadata <- read_csv('FreshMicro_assay_RNAseq_metadata.csv')
RNAseq_count_matrix <- read_csv('RNAseq_count_matrix.csv')

colnames(RNAseq_count_matrix) <- gsub('__RNAseq','',colnames(RNAseq_count_matrix))
colnames(RNAseq_count_matrix)[1] <- 'GeneID'
RNAseq_GeneInfo <- read_csv('RNAseq_GeneInfo.csv')

# RNA-seq + ATAC-seq samples
commonSamples <- intersect(colnames(ATACseq_count_matrix),colnames(RNAseq_count_matrix))
combinedSampsAD <- intersect(commonSamples,AD_Samples$individualID)
combinedSampsCO <- intersect(commonSamples,CON_Samples$individualID)

# Data Normalization and Annotation: RNA-Seq- TMM and gene length,
# ATAC-seq - TMM!
# RNA-Seq
AD_Expr <- RNAseq_count_matrix %>% remove_rownames() %>%
  column_to_rownames(var = 'GeneID') %>%
  select(combinedSampsAD,combinedSampsCO) %>%
  as.data.frame()
geneLength <- readRDS(file = 'G:/Thesis Data/OMIC Data/Priori/geneLengthDF.RDS')
geneLength2 <- geneLength[match(rownames(AD_Expr),geneLength$ENSEMBLID),]
geneLength2 <- geneLength2[complete.cases(geneLength2),]

AD_Expr <- AD_Expr[rownames(AD_Expr) %in% geneLength2$ENSEMBLID,]
samGroups <- as.factor(c(rep('AD',length(combinedSampsAD)),rep('Control',length(combinedSampsCO))))
Expr.Norm <- geTMM(as.matrix(AD_Expr),sampleGrps = samGroups,geneLength2$GENELENGTH,'COUNT')

AD_Expr.Norm <- Expr.Norm %>% select(combinedSampsAD)

# PCA analysis: done independent of ATAC seq peaks availability
# key <- AD_Samples$diagnosis[match(comSampsAD,AD_Samples$individualID)]
# print(pcaPlot(AD_Expr.Norm,key,'AD - Microglia: Before Outlier Removal'))
# 

AD_Expr.Norm.PCA <- Expr.Norm %>% select(combinedSampsAD,combinedSampsCO)
Keys <- FreshMicro_individual_metadata$diagnosis[match(colnames(AD_Expr.Norm.PCA),
                                                       FreshMicro_individual_metadata$individualID)] 

png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_RNA-Seq-normcount-b4.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(AD_Expr.Norm.PCA,Keys,'AD - Microglia: Before Outlier Removal'))
dev.off()

rmOutSamps <- AD_Expr.Norm.PCA %>% select(-c(FM_082,FM_099))

Keys <- FreshMicro_individual_metadata$diagnosis[match(colnames(rmOutSamps),
                                                       FreshMicro_individual_metadata$individualID)]
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_RNA-Seq-normcount-after.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(rmOutSamps,Keys,'AD - Microglia: After Outlier Removal'))
dev.off()

# ATAC-seq


geTMM.Peak <- function(
  countsMat,#Raw counts data matrix
  sampleGrps,#A factor object with sample-group information
  peakIDs
){
  require(edgeR)
  countsMat <- as.matrix(countsMat)
  myCPM <- cpm(countsMat, log = F)
  thresh <- myCPM > 0.1
  keep <- rowSums(thresh) >= min(table(sampleGrps))
  # Subset the rows of countdata to keep the more highly expressed genes
  countsMat <- countsMat[keep,]
  peakIDs <- peakIDs[keep]
  # if not rpk, convert and normalize with gene length
  rpk <- countsMat
  # Create a DGEList object using edgeR
  DGE <- DGEList(counts=rpk, group=sampleGrps)
  # normalize for library size by cacluating scaling factor using TMM (default method)
  DGE <- calcNormFactors(DGE, method = 'TMM')
  # count per million read (normalized count)
  normCounts <- cpm(DGE)
  return(data.frame(Peaks = peakIDs, normCounts))
}
library(dplyr)
library(tidyverse)
ATACseq_count <- ATACseq_count_matrix

CO_SamplesUpdated <- intersect(combinedSampsCO,colnames(rmOutSamps))
AD_SamplesUpdated <- intersect(combinedSampsAD,colnames(rmOutSamps))
# sampGroups <- c(intersect(comSampsAD,colnames(rmOutSamps)),intersect(combinedSampsCO,colnames(rmOutSamps)))
sampGroups <- c(rep('CONTROL',length(CO_SamplesUpdated)),
                     rep('AD',length(AD_SamplesUpdated)))
ATACseq_count.Peak <- ATACseq_count %>% select(c(CO_SamplesUpdated,
                                                 AD_SamplesUpdated))
ATACseq_count.Peak.1 <- geTMM.Peak(
  ATACseq_count.Peak,#Raw counts data matrix
  sampGroups,#A factor object with sample-group information
  peakIDs = ATACseq_count$PeakID
)
# Peak Annotation 
ATACseq_PeakInf <- ATACseq_PeakInfo %>% 
  select(PeakID,Associated.Gene.Name) %>%
  mutate(ENSEMBLID = geneLength$ENSEMBLID[match(.$Associated.Gene.Name,geneLength$SYMBOL)])
ATACseq_count.Peak.1 <- ATACseq_count.Peak.1 %>%
  mutate(ENSEMBL = ATACseq_PeakInf$ENSEMBLID[match(.$Peaks,ATACseq_PeakInf$PeakID)])
ATACseq_count.Peak.1 <- ATACseq_count.Peak.1[!is.na(ATACseq_count.Peak.1$ENSEMBL),]

AD_Peak.Norm.PCA <- ATACseq_count.Peak.1 %>% select(-c(ENSEMBL,Peaks))
Keys <- FreshMicro_individual_metadata$diagnosis[match(colnames(AD_Peak.Norm.PCA),
                                                       FreshMicro_individual_metadata$individualID)] 

png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_ATAC-Seq-normpeaks.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(AD_Peak.Norm.PCA,Keys,'AD - Microglia'))
dev.off()

AD_Peak.Norm.PCA <- ATACseq_count.Peak.1 %>% select(-Peaks)

# print(pcaPlot(rmOutSamps,key$diagnosis,'AD - Microglia: After Outlier Removal'))

# Select most relevant genes
# AD_SamplesFinal <- intersect(colnames(rmOutSamps),comSampsAD)
GeneIDConverter <- readRDS('G:/Projects/Common/GeneCode_GeneInfo_v40_2.RDS') %>%
  select(ID3, ENSEMBL_ID)
GeneIDConverter <- GeneIDConverter[!duplicated(GeneIDConverter$ENSEMBL_ID),]
expr.Data <- rmOutSamps %>% 
  select(CO_SamplesUpdated,AD_SamplesUpdated) #%>%
#   mutate(GeneID = GeneIDConverter$ID3[match(rownames(rmOutSamps),GeneIDConverter$ENSEMBL_ID)])
# rownames(expr.Data) <- GeneIDConverter$ID3[match(rownames(expr.Data),GeneIDConverter$ENSEMBL_ID)]
# Focus on highly variable genes


IDs <- gconvert(query = rownames(expr.Data), organism = "hsapiens",
         target='HGNC', mthreshold = Inf, filter_na = TRUE)
# rownames(expr.Data.2.1) <- ENSEMBL.df.1$SYMBOL[match(rownames(expr.Data.2.1),ENSEMBL.df.1$ENSEMBL)]

# com.genes <- HGNC_ENSEMBL[HGNC_ENSEMBL$ENSEMBL %in% rownames(expr.Data),]
# 
# com.genes <- setdiff(rownames(expr.Data),HGNC_ENSEMBL$ENSEMBL)


CommonGenes <- IDs[,c(2,4)]
colnames(CommonGenes) <- c('ENSEMBL','SYMBOL')
CommonGenes <- rbind(CommonGenes,HGNC_ENSEMBL)
CommonGenes <- rbind(CommonGenes,ENSEMBL.df.1)
CommonGenes <- CommonGenes[!duplicated(CommonGenes$SYMBOL),]
saveRDS(CommonGenes, 'G:/Projects/Common/CommonGenes.RDS')
# IDs.1 <- gconvert(query = com.genes, organism = "hsapiens",
#                 target='HGNC', mthreshold = Inf, filter_na = TRUE)


expr.Data <- expr.Data %>% filter(rownames(expr.Data) %in% CommonGenes$ENSEMBL)
rownames(expr.Data) <- CommonGenes$SYMBOL[match(rownames(expr.Data),CommonGenes$ENSEMBL)]

dfdata.S <- scale(t(expr.Data), scale = T, center = T)
p <- pca(dfdata.S, transposed = T, scale = F, center = F)
explained_var <- p$variance/sum(p$variance)*100
sum(explained_var[1:19])

relGenes <- plotLoadings3(
  p, # Provide a PCA object input
  components = 1:19, # The components you wish to include
  rangeRetain = 0.15 # The percentage  to be considered
)

write.table(relGenes, 'G:/Projects/PhD Thesis/Results/WGCNA/WGCNA-Genes.txt', row.names = F)
expr.Data.2 <- expr.Data %>% filter(rownames(expr.Data) %in% relGenes)

# Patient Sub-groups: RNA-seq

df <- scale(t(expr.Data.2 %>% select(AD_SamplesUpdated)),center = T, scale = T) %>%
  as.data.frame() %>% t()
df <- hcut(t(df), k = 2,hc_method = 'ward.D', hc_metric = 'spearman',
                    hc_func = 'hclust')

png(filename = 'G:/Projects/PhD Thesis/Results/Figures/AD_FreshMicro_SubGroups.png', 
    width = 14.5, height = 7.5, units = 'in', res = 600)
fviz_dend(df, rect = TRUE, cex = 0.5, xlab = 'Samples',
          ylab = 'Distance')+
  ggtitle(paste0('Spearman Correlation on AD Samples: ', 
                 nrow(expr.Data.2),' genes'))
dev.off()

# Patient Sub-groups: ATAC-seq
AD_Peak.Norm_max <- aggregate(.~ENSEMBL, AD_Peak.Norm.PCA, max) %>%
  remove_missing() %>%
  column_to_rownames(var = 'ENSEMBL') %>%
  select(colnames(expr.Data))
ATAC_dfdata.S <- scale(t(AD_Peak.Norm_max), scale = T, center = T)
p <- pca(ATAC_dfdata.S, transposed = T, scale = F, center = F)
explained_var <- p$variance/sum(p$variance)*100
sum(explained_var[1:9])

relGenes <- plotLoadings3(
  p, # Provide a PCA object input
  components = 1:9, # The components you wish to include
  rangeRetain = 0.05 # The percentage  to be considered
)
ATAC_dfdata.S <- AD_Peak.Norm_max %>% filter(rownames(AD_Peak.Norm_max) %in% relGenes)
df <- scale(t(ATAC_dfdata.S),center = T, scale = T) %>%
  as.data.frame() %>% t()
df <- hcut(t(df), k=6,hc_method = 'ward.D', hc_metric = 'spearman',
           hc_func = 'hclust')

png(filename = 'G:/Projects/PhD Thesis/Results/Figures/AD_FreshMicro_SubGroups_ATAC.png', 
    width = 14.5, height = 7.5, units = 'in', res = 600)
fviz_dend(df, rect = TRUE, cex = 0.5, xlab = 'Samples',
          ylab = 'Distance')+
  ggtitle(paste0('Spearman Correlation on AD Samples: ', 
                 nrow(ATAC_dfdata.S),' genes'))
dev.off()
library(tibble)
# Save patient clusters
AD_Clusters <- df$cluster %>% 
  as.data.frame() %>%
  rownames_to_column(var = 'PatientID')
colnames(AD_Clusters)[2] <- 'ClusterLabel'
clustersAD <- unique(AD_Clusters$ClusterLabel)
AD_ClustersDF <- list()
for(i in 1:length(clustersAD)){
  AD_ClustersDF[[paste0(clustersAD[i])]] <- toString(AD_Clusters[which(AD_Clusters$ClusterLabel==clustersAD[i]),]$PatientID)
}

comSampsCON2 <- intersect(comSampsCON,colnames(rmOutSamps))
for(i in 1:length(clustersAD)){
  df <- list()
  AD_PatID <- AD_Clusters[which(AD_Clusters$ClusterLabel==clustersAD[i]),]$PatientID
  df[['AD']] <- rmOutSamps %>% select(AD_PatID) %>% rownames_to_column(var = 'GeneID')
  df[['Control']] <- rmOutSamps %>% select(CO_SamplesUpdated) %>% rownames_to_column(var = 'GeneID')
  df %>% writexl::write_xlsx(paste0('G:/Projects/PhD Thesis/Results/Tables/AD Clusters/Expression/AD_Cluster_RNAExpr_',
                                    clustersAD[i],'.xlsx'))
}

AD_Peak.Norm.PCA.1 <- AD_Peak.Norm.PCA %>% filter(.$ENSEMBL %in% HumanGEM.genes$NAME)
AD_Peak.Norm.PCA.2 <- aggregate(.~ENSEMBL,AD_Peak.Norm.PCA.1,max)
for(i in 1:length(clustersAD)){
  df <- list()
  AD_PatID <- AD_Clusters[which(AD_Clusters$ClusterLabel==clustersAD[i]),]$PatientID
  df[['AD']] <- AD_Peak.Norm.PCA.2 %>% select(c(ENSEMBL,AD_PatID)) #%>% rownames_to_column(var = 'GeneID')
  df[['Control']] <- AD_Peak.Norm.PCA.2 %>% select(c(ENSEMBL,CO_SamplesUpdated))
  df %>% writexl::write_xlsx(paste0('G:/Projects/PhD Thesis/Results/Tables/AD Clusters/Peak/AD_Cluster_Peak_',
                                    clustersAD[i],'.xlsx'))
}

# Cluster Specific Differential Expression Analysis: PPDE
DGE_PPDE <- function(AD, CON){
  df <- cbind(CON,AD)
  design <- factor(c(rep('control',ncol(CON)),rep('disease',ncol(AD))))
  df <- data.matrix(df)
  Sizes <- MedianNorm(df)
  EBSeq.method <- EBTest(Data = df,
                         Conditions= design,
                         sizeFactors=Sizes, maxround=50)
  # LIMMA
  df <- cbind(CON,AD)
  design <- c(rep('control',ncol(CON)),rep('disease',ncol(AD)))
  design <- model.matrix(~0 + factor(design))#Adjust this based on the experimental setup
  colnames(design) <- c('control', 'disease')#
  contrast <- makeContrasts(disease - control, levels = design)
  fit <- lmFit(log2(df+1), design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  degDF <- topTable(fit2, n = 'Inf', sort.by = 'none')
  
  df <- cbind(degDF$logFC, degDF$adj.P.Val, EBSeq.method$PPDE) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'GeneID')
  colnames(df)[2:4] <- c('FC','AdjPVal','PPDE')
  return(df)
}
clustersAD.DGE <- list()


# df <- cbind(control,disease)
# design <- factor(c(rep('control',ncol(control)),rep('disease',ncol(disease))))
# df <- data.matrix(df)
# Sizes <- MedianNorm(df)
# EBSeq.method <- EBTest(Data = df,
#                        Conditions= design,
#                        sizeFactors=Sizes, maxround=50)
# EBSeq.method.1 <- data.frame(EBSeq.method$PPDE)
# clustersAD.Expre <- list()
for(i in 1:length(clustersAD)){
  # df <- list()
  AD_PatID <- AD_Clusters[which(AD_Clusters$ClusterLabel==clustersAD[i]),]$PatientID
  disease <- rmOutSamps %>% select(AD_PatID)
  control <- rmOutSamps %>% select(CO_SamplesUpdated)
  df <- DGE_PPDE(disease,control)
  clustersAD.DGE[[paste0(clustersAD[i])]] <- df
  # clustersAD.Expre[[paste0(clustersAD[i])]] <- list(Controls = control, AD = disease)
  # df <- list(Controls = control, AD = disease)
  df %>% writexl::write_xlsx(paste0('G:/Projects/PhD Thesis/Results/Tables/AD Clusters/DGE_mRNA/AD_Cluster_RNA_',
                                    clustersAD[i],'.xlsx'))
}

library(gprofiler2)
cluster.RNA <- list()
for(i in 1:length(clustersAD.DGE)){
  # endos.df[i,1] <- names(endoPhen.1)[i]
  cluster.RNA[[paste0(names(clustersAD.DGE)[i])]] <- gost(query = clustersAD.DGE[[i]]$GeneID[clustersAD.DGE[[i]]$PPDE>.9],
                                                     organism = "hsapiens",
                                                     correction_method = 'bonferroni')$result[,-c(1,7,8,12:14)] %>%
    as_tibble() %>% filter(.$source =='WP')
  # endos.df[i,2] <- nrow(endoPhen.1[[i]])
}
# degs <- clustersAD.DGE$`1`[which(clustersAD.DGE$`1`$AdjPVal<0.05 &
#                   clustersAD.DGE$`1`$PPDE>(1-0.05)),]
# degs.PPDE <- clustersAD.DGE$`1`[clustersAD.DGE$`1`$PPDE<0.05,]
# degs.ADJ <- clustersAD.DGE$`1`[clustersAD.DGE$`1`$AdjPVal<0.05,]

AD_Peak.Norm.PCA.3 <- aggregate(.~ENSEMBL,AD_Peak.Norm.PCA.1,median)
clustersAD.DPE <- list()
for(i in 1:length(clustersAD)){
  # df <- list()
  AD_PatID <- AD_Clusters[which(AD_Clusters$ClusterLabel==clustersAD[i]),]$PatientID
  disease <- AD_Peak.Norm.PCA.2 %>% select(c(ENSEMBL,AD_PatID)) %>%
    remove_rownames() %>%
    column_to_rownames(var = 'ENSEMBL')
  control <- AD_Peak.Norm.PCA.2 %>% select(c(ENSEMBL,CO_SamplesUpdated)) %>%
    remove_rownames() %>%
    column_to_rownames(var = 'ENSEMBL')
  df <- DGE_PPDE(disease,control)
  clustersAD.DPE[[paste0(clustersAD[i])]] <- df
  writexl::write_xlsx(df, paste0('G:/Projects/PhD Thesis/Results/Tables/AD Clusters/DGE_Peaks/AD_Cluster_Peak_',
                                    clustersAD[i],'.xlsx'))
}

#######################################################################################################
# Co-expression Analysis: WGCNA
set.seed(6699)

# Perform initial clustering analysis to ensure that no more outliers exist


expr.Data.2.1 <- expr.Data.2



sampleTree <- hclust(dist(t(expr.Data.2)), method = 'average')

png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/WGCNA-initial-clustering.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(plot(sampleTree))
dev.off()


expr.Data.3 <- expr.Data.2 %>% select(-FM_098)

expr.Data.4 <- t(expr.Data.2)
sampleTree <- hclust(dist(t(expr.Data.3)), method = 'average')
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/WGCNA-after-clustering.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(plot(sampleTree))
dev.off()

# Sample meta-information

meta_info <- FreshMicro_individual_metadata %>% filter(.$individualID %in% as.character(colnames(expr.Data.3)))
meta_info$ageDeath <- as.numeric(gsub('\\+','',meta_info$ageDeath))
meta_info$Braak[c(nrow(meta_info)-1,nrow(meta_info))] <- 0
# First, we empirically estimate the power based on network scale free topology
powers <- c(1:20)
# Call the network topology analysis function
tic()
sft <- pickSoftThreshold(t(expr.Data.2), 
                         powerVector = powers, 
                         verbose = 5)
toc()

# cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
df <- data.frame(Power = sft$fitIndices[,1],
                 R2 = -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                 MeanK = sft$fitIndices[,5])
# Scale-free topology plot
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/WGCNA-scale-free-topology.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
ggplot(data = df)+
  aes(x=Power, y = R2)+
  geom_point() +
  geom_text(aes(label = row.names(df)),hjust=0.6, vjust=0, size = 5)+
  geom_hline(yintercept=.95, color = 'red')+
  labs(title = paste0('Scale-free Topology Plot'))+
  theme_bw()
dev.off()
# We choose power = 7
# Mean connectivity plot
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/WGCNA-mean-connectivity.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
ggplot(data = df)+
  aes(x=Power, y = MeanK)+
  geom_point() +
  geom_text(aes(label = row.names(df)),hjust=0.6, vjust=0, size = 5)+
  labs(title = paste0('Mean Connectivity'))+
  theme_bw()
dev.off()
# Construct an adjacency matrix
tic()
adjaMat <- adjacency(t(expr.Data.3), power = 7)
toc()
# Turn adjacency into topological overlap matrix
tic()
TOM <- TOMsimilarity(adjaMat)
dissTOM <- 1-TOM
toc()
# Clustering the TOM matrix
geneTree <- hclust(as.dist(dissTOM), method = 'average')
# Plot the resulting clustering tree (dendrogram)
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/WGCNA-TOM-clusters.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
plot(geneTree, xlab="", sub="", main = 'Gene clustering on TOM-based dissimilarity',
     labels = FALSE, hang = 0.04);
dev.off()
# We like large modules, so we set the minimum module size relatively high:
# minModuleSize <- 20
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 50);
# table(dynamicMods)

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/WGCNA-dendrogram.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
# Merge modules with similar patterns
# Calculate eigengenes
MEList <- moduleEigengenes(t(expr.Data.3), colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs);
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres <- 0.4 # or 60% similarity!
# Plot the cut line into the dendrogram
# abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge <- mergeCloseModules(t(expr.Data.3), 
                           dynamicColors, 
                           cutHeight = MEDissThres, 
                           verbose = 3)
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/WGCNA-dendrogram-merged.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# 
# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

HGNC_ENSEMBL <- read.table('G:/Projects/Common/Homo_sapiens.GRCh38.106.gtf', sep = '\t')

# ENSEMBL.df <- data.frame(matrix(data = NA, ncol = 2, nrow = nrow(HGNC_ENSEMBL)))
# colnames(ENSEMBL.df) <- c('ENSEMBL','SYMBOL')
# for(i in 1:nrow(HGNC_ENSEMBL)){
#   # i <- 2344919
#   id.info <- HGNC_ENSEMBL$V9[i]
#   
#   id.info <- str_split(id.info, pattern = ';')
#   id.info <- as.character(id.info[[1]])
#   geneName <- id.info[grepl('^ gene_name', id.info)==T]
#   if(length(geneName)>0){
#     ENSEMBL.df[i,1] <- str_split(id.info[grepl('^gene_id', id.info)==T], pattern = ' ')[[1]][2]
#     ENSEMBL.df[i,2] <- str_split(geneName, pattern = ' ')[[1]][3]
#   }
# }
ENSEMBL.df.1 <- ENSEMBL.df[!duplicated(ENSEMBL.df$ENSEMBL),]
ENSEMBL.df.1 <- ENSEMBL.df.1[complete.cases(ENSEMBL.df.1),]

writexl::write_xlsx(ENSEMBL.df.1, 'G:/Projects/Common/Ensembl2GeneSymbol.xlsx')
HGNC_ENSEMBL <- read_excel("HGNC_ENSEMBL.xlsx", 
                           sheet = "Sheet2")[,c(2,3,10)] %>%
  mutate(ENSEMBL = HGNC_ENSEMBL$`Ensembl ID(supplied by Ensembl)`,
         SYMBOL = HGNC_ENSEMBL$`Approved symbol`) %>%
  select(ENSEMBL,SYMBOL)

modules <- unique(moduleColors)
# Select module probes
probes <- rownames(adjaMat)
cytModules <- list()
for(i in 1:length(modules)){
  inModule<- is.finite(match(moduleColors, modules[i]))
  modProbes <- probes[inModule]
  # modGenes <- HGNC_ENSEMBL$SYMBOL[match(modProbes, HGNC_ENSEMBL$ENSEMBL)]
  # Select the corresponding Topological Overlap
  modTOM <- TOM[inModule, inModule]
  # dimnames(modTOM) <- list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cytModules[[paste0(modules[i])]] <- exportNetworkToCytoscape(modTOM,
                                      edgeFile = paste("G:/Projects/PhD Thesis/Results/WGCNA/Modules/CytoscapeInput-edges-", 
                                                       paste(modules[i], collapse="-"), ".txt", sep=""),
                                      nodeFile = paste("G:/Projects/PhD Thesis/Results/WGCNA/Modules/CytoscapeInput-nodes-", 
                                                       paste(modules[i], collapse="-"), ".txt", sep=""),
                                      weighted = TRUE,
                                      threshold = 0.05,
                                      nodeNames = modProbes,
                                      # altNodeNames = modGenes,
                                      nodeAttr = moduleColors[inModule])
}


write.table(cytModules$turquoise$edgeData,'G:/Projects/PhD Thesis/Results/WGCNA/turquoiseModule.txt', row.names = F)

# Save modules
modCols <- unique(moduleColors)
moduleGenes <- list()
df.Modules.WGCNA <- data.frame(matrix(data = NA, ncol = 2, nrow = length(cytModules)))
colnames(df.Modules.WGCNA) <- c('Module','Number of Genes')
for(i in 1:length(cytModules)){
  moduleGenes[[paste0(names(cytModules)[i])]] <- cytModules[[i]]$nodeData$nodeName
  df.Modules.WGCNA[i,] <- c(names(cytModules)[i], length(cytModules[[i]]$nodeData$nodeName))
}
writexl::write_xlsx(df.Modules.WGCNA, 'G:/Projects/PhD Thesis/Results/WGCNA/ModuleSizes.xlsx')
# ################################################################################################################
# Relating module genes to endophenotypes

# endoPhen <- readRDS('G:/Thesis Data/OMIC Data/Priori/endoBioGRID.RDS')
endoPhen.1 <- readRDS('G:/Projects/PhD Thesis/Results/Endophenotypes/AD_Endophenotype2Gene_2.RDS')
endoPhen.1 <- lapply(endoPhen.1, function(x){x = x %>% 
  filter(.$Frequency >5) %>% 
  select(AD.Gene)})

# Display the number of genes associated with each endophenotype
endos.df <- data.frame(matrix(data = NA, 
                              nrow = length(endoPhen.1), 
                              ncol = 2))
colnames(endos.df) <- c('Endophenotype','NumberOfAssociatedGenes')
for(i in 1:length(endoPhen.1)){
  endos.df[i,1] <- names(endoPhen.1)[i]
  endos.df[i,2] <- nrow(endoPhen.1[[i]])
}
as_tibble(endos.df)

# Perform enrichment analysis for each endophenotype gene list: confirmatory analysis
library(gprofiler2)
endos.list <- list()
for(i in 1:length(endoPhen.1)){
  # endos.df[i,1] <- names(endoPhen.1)[i]
  endos.list[[paste0(names(endoPhen.1)[i])]] <- gost(query = endoPhen.1[[i]]$AD.Gene,
                  organism = "hsapiens",
                  correction_method = 'bonferroni')$result[,-c(1,7,8,12:14)] %>%
    as_tibble()
  # endos.df[i,2] <- nrow(endoPhen.1[[i]])
}
# as_tibble(endos.df)
# endoPhen <- endoPhen.1[-c(4,5,7,9)]
endoPhen2 <- endoPhen.1
# 
moduleEnrich <- data.frame(matrix(data = NA, 
                                  ncol = length(moduleGenes), 
                                  nrow = length(endoPhen2)))
colnames(moduleEnrich) <- names(moduleGenes)
rownames(moduleEnrich) <- names(endoPhen2)
moduleEnrich.Fisher.Adj <- moduleEnrich
moduleEnrich.HyperG.Adj <- moduleEnrich
moduleEnrich.FisherEx <- moduleEnrich
moduleEnrich.HyperGeo <- moduleEnrich
universe.geneID <- rownames(expr.Data)
# universe.geneID <- universe.geneID[gsub(' ', 'NA',universe.geneID)]
# universe.geneID <- universe.geneID[!is.na(universe.geneID)]
# universe.geneID <- universe.geneID[universe.geneID != ""]

endoPhenUni <- unique(unlist(endoPhen2))
for(i in 1:length(moduleGenes)){
  # modGenesI <- unique(geneLength2$SYMBOL[match(moduleGenes[[i]],geneLength2$ENSEMBLID)]) 
  modGenesI <- moduleGenes[[i]] 
  hyperTest <- c()
  FisherTests <- c()
  for(j in 1:length(endoPhen2)){
    endoPhenJ <- endoPhen2[[j]]$AD.Gene
    # Hypergeometric Test
    q <- length(intersect(modGenesI,endoPhenJ))
    m <- length(endoPhenJ)
    n <- length(endoPhenUni)
    k <- length(modGenesI)
    hyperTest <- c(hyperTest, phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE))
    
    # Fishers Exact Test
    # nrow(expr.Data)
    N <- length(universe.geneID)
    N_disease <-  length(intersect(endoPhenJ, universe.geneID))
    N_cluster <- length(modGenesI)
    # endoPhen2[[j]]$
    # Contingency table: cluster&disease, cluster&not_disease, not_cluster&disease, not_cluster&not_disease
    cluster_disease <- length(unique(intersect(endoPhenJ, modGenesI)))
    cluster_not_disease <- N_cluster - cluster_disease
    not_cluster_disease <- N_disease - cluster_disease
    not_cluster_not_disease <- N - cluster_disease - cluster_not_disease - not_cluster_disease
    
    cont_table <- rbind(c(cluster_disease, cluster_not_disease),
                        c(not_cluster_disease, not_cluster_not_disease))
    
    # Fisher's exact test
    FisherTests <- c(FisherTests,fisher.test(cont_table, alternative = "two.sided")$p.value)
    # FisherTests <- fisher.test(cont_table, alternative = "two.sided")
    # FisherTests$p.value
  }
  moduleEnrich.HyperGeo[,i] <- hyperTest
  moduleEnrich.FisherEx[,i] <- FisherTests
  moduleEnrich.HyperG.Adj[,i] <- p.adjust(hyperTest,method="BH") # Benjamini-Hochberg Corrected
  moduleEnrich.Fisher.Adj[,i] <- p.adjust(FisherTests,method="BH") # Benjamini-Hochberg Corrected
}

# moduleEnrich.HyperGeo[moduleEnrich.HyperGeo>0.05] <- 1
# moduleEnrich.FisherEx[moduleEnrich.FisherEx>0.05] <- 1
# moduleEnrich.HyperG.Adj[moduleEnrich.HyperG.Adj>0.005] <- 1
# moduleEnrich.Fisher.Adj[moduleEnrich.Fisher.Adj>0.05] <- 1
library(tibble)
WGCNA.Res <- list()
WGCNA.Res[['HyperGeo.WG']] <- moduleEnrich.HyperGeo %>% rownames_to_column(var = 'Endophenotype')
WGCNA.Res[['FisherEx.WG']] <- moduleEnrich.FisherEx %>% rownames_to_column(var = 'Endophenotype')
WGCNA.Res[['HyperG.Adj.WG']] <- moduleEnrich.HyperG.Adj %>% rownames_to_column(var = 'Endophenotype')
WGCNA.Res[['Fisher.Adj.WG']] <- moduleEnrich.Fisher.Adj %>% rownames_to_column(var = 'Endophenotype')

WGCNA.Res %>% writexl::write_xlsx('G:/Projects/PhD Thesis/Results/WGCNA/WGCNA-Endophenotypes.xlsx')


cytModules.1 <- cytModules[c(10,12)]
saveRDS(cytModules.1,'G:/Projects/PhD Thesis/Results/WGCNA/WGCNA-Modules.RDS')
writexl::write_xlsx(moduleEnrich %>% rownames_to_column(var = 'Endophenotype'), 
                    'G:/Projects/PhD Thesis/Results/Tables/module2endophenotype-HPV_0005.xlsx')
writexl::write_xlsx(moduleEnrich.Adj %>% rownames_to_column(var = 'Endophenotype'), 
                    'G:/Projects/PhD Thesis/Results/Tables/module2endophenotype-HBH_0005.xlsx')

writexl::write_xlsx(moduleEnrich.Adj %>% rownames_to_column(var = 'Endophenotype'), 
                    'G:/Projects/PhD Thesis/Results/Tables/module2endophenotype-HBH_0005.xlsx')
writexl::write_xlsx(moduleEnrich.Fisher.Adj %>% rownames_to_column(var = 'Endophenotype'), 
                    'G:/Projects/PhD Thesis/Results/Tables/module2endophenotype-HBH_0005.xlsx')

# Co-expression Analysis 2: MEGENA
# library(devtools)
# install_github('songw01/MEGENA')
install.packages('MEGENA')
library(MEGENA)
expr.Data.4 <- expr.Data.3 #%>% mutate(GeneSymbol = geneLength$SYMBOL[match(rownames(expr.Data.3),geneLength$ENSEMBLID)])
# expr.Data.4 <- expr.Data.4[!which(expr.Data.4$GeneSymbol==' '),]

# expr.Data.4 <- expr.Data.4[!(is.na(expr.Data.4$GeneSymbol) | expr.Data.4$GeneSymbol=="") & !duplicated(expr.Data.4$GeneSymbol), ] %>%
#   remove_rownames() %>%
#   column_to_rownames(var = 'GeneSymbol') 

# %>%
#   t()
# meg.1 <- MEGENA::calculate.correlation(expr.Data.4)

# Calculate correlations with 50 permutations
set.seed(669588)
meg.1 <- calculate.correlation(expr.Data.3, doPerm = 500, 
                               method = 'spearman', FDR.cutoff = 0.05)
meg.2 <- calculate.PFN(meg.1)
g <- graph.data.frame(meg.2,directed = FALSE)
MEGENA.output <- do.MEGENA(g = g,remove.unsig = FALSE,
                           doPar = T,n.perm = 500)
output.summary <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = 0.05, hub.pvalue = 0.05,
                                       min.size = 50, max.size = 500,
                                       annot.table = NULL,id.col = NULL,symbol.col = NULL,
                                       output.sig = TRUE)
# MEGENA::plot_module(output.summary, PFN = meg.2)
# Perform Module - Endophenotype Enrichment Analysis

endoPhen.1 <- readRDS('G:/Projects/PhD Thesis/Results/Endophenotypes/AD_Endophenotype2Gene_2.RDS')
endoPhen.1 <- lapply(endoPhen.1, function(x){x = x %>% 
  filter(.$Evidence >quantile(.$Evidence,.5)) %>% 
  select(Gene)})

# endoPhen <- endoPhen.1[-c(4,5,7,9)]
endoPhen2 <- endoPhen.1
# for(i in 1:length(endoPhen)){
#   df <- endoPhen[[i]]
#   colnames(df)[1] <- 'GeneID'
#   df$nChar <- nchar(df$GeneID)
#   df <- df[df$nChar!=1,]
#   # df$GeneID2 <- df$GeneID1
#   # for(j in 1:length(df$GeneID)){
#   #   df$GeneID2[j] <- ifelse(nchar(df$GeneID1[j])!=1,df$GeneID1[j],'')
#   # }
#   # df <- df$GeneID2[-df$GeneID2=='']
#   
#   
#   endoPhen2[[paste0(names(endoPhen)[i])]] <- data.frame(GeneID = df$GeneID2[df$GeneID2!=''])
# }

MEGENA.Modules.Res <- output.summary$modules

df.Modules.MEGENA <- data.frame(matrix(data = NA, ncol = 2, nrow = length(MEGENA.Modules.Res)))
colnames(df.Modules.MEGENA) <- c('Module','Number of Genes')
for(i in 1:length(MEGENA.Modules.Res)){
  # MEGENA.Modules.Res$c1_6
  # moduleGenes[[paste0(names(cytModules)[i])]] <- cytModules[[i]]$nodeData$nodeName
  df.Modules.MEGENA[i,] <- c(names(MEGENA.Modules.Res)[i], length(MEGENA.Modules.Res[[i]]))
}
writexl::write_xlsx(df.Modules.MEGENA, 'G:/Projects/PhD Thesis/Results/MEGENA/ModuleSizes.xlsx')

moduleEnrich.ME <- data.frame(matrix(data = NA, 
                                  ncol = length(MEGENA.Modules.Res), 
                                  nrow = length(endoPhen2)))
colnames(moduleEnrich.ME) <- names(MEGENA.Modules.Res)
rownames(moduleEnrich.ME) <- names(endoPhen2)
moduleEnrich.Fisher.Adj.ME <- moduleEnrich.ME
moduleEnrich.HyperG.Adj.ME <- moduleEnrich.ME
moduleEnrich.FisherEx.ME <- moduleEnrich.ME
moduleEnrich.HyperGeo.ME <- moduleEnrich.ME

universe.geneID.1 <- universe.geneID[universe.geneID != ""]

for(i in 1:length(MEGENA.Modules.Res)){
  modGenesI <- MEGENA.Modules.Res[[i]] 
  hyperTest <- c()
  FisherTests <- c()
  for(j in 1:length(endoPhen2)){
    endoPhenJ <- endoPhen2[[j]]$AD.Gene
    q <- length(intersect(modGenesI,endoPhenJ))
    m <- length(endoPhenJ)
    n <- length(endoPhenUni)-m
    k <- length(modGenesI)
    hyperTest <- c(hyperTest, phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE))
    # df <- data.frame("cured" = c(60, 30), "noncured" = c(10, 25), row.names = c("treated", "nontreated"))
    # Fishers Exact Test
    # nrow(expr.Data)
    N <- length(universe.geneID.1)
    N_disease <-  length(intersect(endoPhenJ, universe.geneID.1))
    N_cluster <- length(modGenesI)
    # endoPhen2[[j]]$
    # Contingency table: cluster&disease, cluster&not_disease, not_cluster&disease, not_cluster&not_disease
    cluster_disease <- length(unique(intersect(endoPhenJ, modGenesI)))
    cluster_not_disease <- N_cluster - cluster_disease
    not_cluster_disease <- N_disease - cluster_disease
    not_cluster_not_disease <- N - cluster_disease - cluster_not_disease - not_cluster_disease
    
    cont_table <- rbind(c(cluster_disease, cluster_not_disease),
                        c(not_cluster_disease, not_cluster_not_disease))
    
    # Fisher's exact test
    FisherTests <- c(FisherTests,fisher.test(cont_table, alternative = "two.sided")$p.value)
    # FisherTests <- fisher.test(cont_table, alternative = "two.sided")
    # FisherTests$p.value
  }
  moduleEnrich.HyperGeo.ME[,i] <- hyperTest
  moduleEnrich.FisherEx.ME[,i] <- FisherTests
  moduleEnrich.HyperG.Adj.ME[,i] <- p.adjust(hyperTest,method="BH") # Benjamini-Hochberg Corrected
  moduleEnrich.Fisher.Adj.ME[,i] <- p.adjust(FisherTests,method="BH") # Benjamini-Hochberg Corrected
}

# moduleEnrich.HyperGeo.ME.1 <- apply(moduleEnrich.HyperGeo.ME, 1, function(x){x[numeric(x)>0.05] = 1})
# moduleEnrich.HyperGeo.ME <- moduleEnrich.HyperGeo.ME[,-apply(moduleEnrich.HyperGeo.ME,2,sum)==nrow(moduleEnrich.HyperGeo.ME)]
# moduleEnrich1[moduleEnrich1>0.005] <- 1
library(tibble)
MEGENA.Res <- list()
MEGENA.Res[['HyperGeo.ME']] <- moduleEnrich.HyperGeo.ME %>% rownames_to_column(var = 'Endophenotype')
MEGENA.Res[['FisherEx.ME']] <- moduleEnrich.FisherEx.ME %>% rownames_to_column(var = 'Endophenotype')
MEGENA.Res[['HyperG.Adj.ME']] <- moduleEnrich.HyperG.Adj.ME %>% rownames_to_column(var = 'Endophenotype')
MEGENA.Res[['Fisher.Adj.ME']] <- moduleEnrich.Fisher.Adj.ME %>% rownames_to_column(var = 'Endophenotype')

MEGENA.Res %>% writexl::write_xlsx('G:/Projects/PhD Thesis/Results/MEGENA/MEGENA-Endophenotypes.xlsx')

MEGENA.SigModules <- moduleEnrich.Fisher.Adj.ME %>% select(c1_46,c1_50,c1_196)


cytModules.2 <- MEGENA.Modules.Res[c(28,30,49)]
saveRDS(cytModules.2,'G:/Projects/PhD Thesis/Results/MEGENA/MEGENA-Modules.RDS')

# 
# Save significant modules
stat.Sig <- moduleEnrich.Fisher.Adj.ME %>% rownames_to_column(var = 'Endophenotype')
sig.Modules <- list()
for(i in 1:nrow(stat.Sig)){
  endos <- stat.Sig[i,]
  sig.Modules.1 <- list()
  for(j in 1:ncol(endos[,-1])){
    mod <- colnames(endos)[j]
    if(endos[j] < 0.01){
      sig.Modules.1[[paste0(mod)]] <- MEGENA.Modules.Res[names(MEGENA.Modules.Res)==mod]
    }
  }
  sig.Modules[[paste0(stat.Sig[i,1])]] <- sig.Modules.1
}
library(rlist)
sig.Modules <- list.clean(sig.Modules, fun = is.null, recursive = FALSE)
saveRDS(sig.Modules, 'G:/Projects/PhD Thesis/Results/Tables/sig.Modules_MEGENA.RDS')

writexl::write_xlsx(moduleEnrich1 %>% rownames_to_column(var = 'Endophenotype'), 
                    'G:/Projects/PhD Thesis/Results/Tables/module2endophenotype-HPV_0005-MEGENA.xlsx')
writexl::write_xlsx(moduleEnrich2 %>% rownames_to_column(var = 'Endophenotype'), 
                    'G:/Projects/PhD Thesis/Results/Tables/module2endophenotype-HBH_0005-MEGENA.xlsx')
endoPhen %>% writexl::write_xlsx('G:/Projects/PhD Thesis/Results/Tables/EndophenotypeGene.xlsx')
# # no coloring
# sbobj <- draw_sunburst_wt_fill(module.df = output.summary$module.table,
#                               feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
# sbobj
# # get some coloring (with log transform option)
# mdf <- output.summary$module.table
# mdf$heat.pvalue <- runif(nrow(mdf),0,0.1)
# sbobj <- draw_sunburst_wt_fill(module.df = mdf,feat.col = "heat.pvalue",log.transform = TRUE,
#                               fill.type = "continuous",
#                               fill.scale = scale_fill_gradient2(low = "white",mid = "white",high = "red",
#                                                                 midpoint = -log10(0.05),na.value = "white"),
#                               id.col = "module.id",parent.col = "module.parent")
# sbobj
# # get discrete coloring done
# mdf$category <- factor(sample(x = c("A","B"),size = nrow(mdf),replace = TRUE))
# sbobj <- draw_sunburst_wt_fill(module.df = mdf,feat.col = "category",
#                               fill.type = "discrete",
#                               fill.scale = scale_fill_manual(values = c("A" = "red","B" = "blue")),
#                               id.col = "module.id",parent.col = "module.parent")
# sbobj
# 
# pnet.obj <- plot_module(output = output.summary,PFN = g,subset.module = NULL,
#                         layout = "kamada.kawai",label.hubs.only = FALSE,
#                         gene.set = list("hub.set" = c("CD3E","CD2")),color.code = c("red"),
#                         output.plot = FALSE,out.dir = "modulePlot",col.names = c("grey","grey","grey"),
#                         hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)
# 
# print(pnet.obj[[1]])

# Bayesian RF-Gene Network Modelling using Priors (Genome Accessibility)
BiocManager::install('bnlearn')
BiocManager::install('R.matlab')

library(R.matlab)
library(bnlearn)
# Prepare Metabolic Gene Expression
library(readxl)
HumGEM <- read_excel('G:/Thesis Data/Human-GEM-main/model/Human-GEM.xlsx', 
                        sheet = 'GENES')$NAME
# humGEM <- readMat('G:/Thesis Data/Human-GEM-main/model/Human-GEM.mat')

expr.Data.5 <- Expr.Norm %>% select(comSampsAD,comSampsCON) %>% 
  filter(rownames(AD_Expr.Norm) %in% HumGEM)

# TF/RF
library(readr)
TF_Target_information <- read_delim('G:/Thesis Data/OMIC Data/Priori/TF/TF-Target-information.txt', 
                                    delim = '\t', escape_double = FALSE, 
                                    trim_ws = TRUE)
TF_Target_information <- TF_Target_information %>% 
  select(-tissue) %>%
  distinct() %>%
  mutate(ENSEMBLID = geneLength$ENSEMBLID[match(.$target,geneLength$SYMBOL)])
colnames(TF_Target_information) <- c('RF','Gene.Symbol', 'ENSEMBLID' )
TF_Target_information$Type <- 'TF'
# microRNAs
miRTarBase_MTI <- read_excel('G:/Thesis Data/OMIC Data/Priori/microRNA/miRTarBase_MTI.xlsx') %>%
  filter(.$`Species (miRNA)` == 'Homo sapiens') %>%
  mutate(ENSEMBLID = geneLength$ENSEMBLID[match(miRTarBase_MTI$`Target Gene`,geneLength$SYMBOL)])

miRNet_mir_gene_hsa_mirtarbase <- read_csv('G:/Thesis Data/OMIC Data/Priori/microRNA/miRNet-mir-gene-hsa-mirtarbase.csv')

hsa_MTI <- read_excel("G:/Projects/Common/TF Data/hsa_MTI.xlsx") %>%
  filter(.$`Support Type` == 'Functional MTI')
# miRNet_mir_gene_hsa_mirecords <- read_csv('G:/Thesis Data/OMIC Data/Priori/microRNA/miRNet-mir-gene-hsa-mirecords.csv')
# consolidate
miRTarBase_MTI.1 <- miRTarBase_MTI[,c(2,4,10)]
colnames(miRTarBase_MTI.1) <- c('miRNA','Gene.Symbol','ENSEMBLID')
miRNet_mir_gene_hsa_mirtarbase.1 <- miRNet_mir_gene_hsa_mirtarbase[,c(2,4,5)]
colnames(miRNet_mir_gene_hsa_mirtarbase.1) <- c('miRNA','Gene.Symbol','ENSEMBLID')

miRNAs <- rbind(miRTarBase_MTI.1,miRNet_mir_gene_hsa_mirtarbase.1)
# Remove duplicates
miRNAs <- miRNAs %>% dplyr::distinct()
colnames(miRNAs)[1] <- 'RF'
miRNAs$Type <- 'miRNA'

RF.data <- rbind(TF_Target_information,miRNAs) %>%
  as.data.frame() %>%
  distinct()
RF.data.1 <- RF.data[!is.na(RF.data$ENSEMBLID),]
tfs <- grepl('TF',RF.data.1$Type,ignore.case = T)
mirs <- grepl('miRNA',RF.data.1$Type,ignore.case = T)
RF.data.1$RF.ENS[tfs==T] <- geneLength$ENSEMBLID[match(RF.data.1$RF[tfs==T],geneLength$SYMBOL)]
RF.data.1$RF.ENS[mirs==T] <- RF.data.1$RF[mirs==T]


# Harmonize Regulatory Data
miRData <- hsa_MTI[,c(2,4)]
colnames()
# Prior Network

# ATACseq_count.1 <- ATACseq_count %>%
#   filter(ATACseq_count$TF_ENS %in% RF.data.1$RF.ENS)
# ATACseq_count.2 <- ATACseq_count.1[complete.cases(ATACseq_count.1$TF_ENS),]

# ####################################################################
# Filter for strong TF activity using Peak quantiles

# TF Effect Type

library(readr)
trrust_rawdata_human <- read_delim('G:/Projects/Common/TF Data/trrust_rawdata.human.tsv', 
                                   delim = '\t', escape_double = FALSE, 
                                   col_names = FALSE, trim_ws = TRUE) %>% 
  select(-X4)
#  %>%
# filter(.$X3 != 'Unknown')
colnames(trrust_rawdata_human) <- c('RF','Target','Type')

# trrust_rawdata_human$TF <- geneLength$ENSEMBLID[match(trrust_rawdata_human$TF,geneLength$SYMBOL)]
# trrust_rawdata_human$Target <- geneLength$ENSEMBLID[match(trrust_rawdata_human$Target,geneLength$SYMBOL)]

trrust_rawdata_human.1 <- trrust_rawdata_human[!(is.na(trrust_rawdata_human$TF)|is.na(trrust_rawdata_human$Target)),]

trrust_rawdata_human.1 <- trrust_rawdata_human %>%
  mutate(Interaction.Sign = case_when(.$Type=='Activation' ~ '1',
                              .$Type== 'Repression'~ '-1',
                              TRUE ~ .$Type)) %>%
  select(c(RF,Target,Interaction.Sign))

# Harmonize Regulatory Data
miRData <- hsa_MTI[,c(2,4)]
colnames(miRData) <- c('RF','Target')
miRData$Interaction.Sign <- '-1'

RFData <- rbind(trrust_rawdata_human.1,miRData)
writexl::write_xlsx(RFData, 'G:/Projects/Common/TF Data/RFData.xlsx')
# 
ATACseq_peakQ <- ATACseq_count.2 %>%
  select(-c(Peaks, TF_ENS)) %>%
  as.matrix() %>%
  quantile(probs = seq(0,1,.25))
ATACseq_peakQuant <- ATACseq_count.2 %>%
  select(-c(Peaks, TF_ENS)) %>%
  as.matrix()
for(i in 1:nrow(ATACseq_peakQuant)){
  for(j in 1:ncol(ATACseq_peakQuant)){
    ATACseq_peakQuant[i,j] <- if_else(ATACseq_peakQuant[i,j]<ATACseq_peakQ[[4]],
                                      0,
                                      round(ATACseq_peakQuant[i,j]/ATACseq_peakQ[[4]], digits = 0))
  }
}
ATACseq_peakQuant.1 <- data.frame(TF = ATACseq_count.2$TF_ENS,
                                  ATACseq_peakQuant)
# colnames(ATACseq_peakQuant.1) <- gsub('__ATACseq','',colnames(ATACseq_peakQuant.1))
ATACseq_peakQuant.2 <- ATACseq_peakQuant.1 %>% 
  select(c('TF',intersect(comSampsAD,AD_CON),
           intersect(comSampsCON,AD_CON)))
# Build an individual Prior Matrix
indPriorMat <- ATACseq_peakQuant.1[,c(1:2)] %>%
  filter(.$FM_002 != 0)
length(unique(indPriorMat$TF))
tfs <- unique(indPriorMat$TF)
tfs <- tfs[!is.na(tfs)] %>% as.data.frame()
write_tsv(tfs, 'Inferelator/Input/tf_names.tsv')
# Use Peaks with maximum effect
trrust_rawdata_human.1$Int.Sign <- as.numeric(trrust_rawdata_human.1$Int.Sign)
indPriorMat <- aggregate(.~indPriorMat$TF,indPriorMat,max)[,c(2,3)]
tf.targets <- RF.data.1 %>% filter(.$RF.ENS %in% indPriorMat$TF)
tf.targets <- data.frame(RF = tf.targets$RF.ENS, Target = tf.targets$ENSEMBLID)
tfs <- unique(tf.targets$RF)
targs <- unique(tf.targets$Target)
tf.df <- data.frame(matrix(data = NA, 
                           ncol = length(tfs), 
                           nrow = length(targs)))
colnames(tf.df) <- tfs
rownames(tf.df) <- targs
for(i in 1:ncol(tf.df)){#
  # i <- 1
  tf.i <- colnames(tf.df)[i]
  tf.i.w <- as.numeric(indPriorMat$FM_002[indPriorMat$TF==tf.i])
  for(j in 1:nrow(tf.df)){#
    # j <- 3
    tar.j <- rownames(tf.df)[j]
    # tar.j <- 'ENSG00000033011'
    # Extract the TF targets from the TF-targets db
    rf.dat <- tf.targets %>% filter(.$RF %in% tar.i)
    # Check if the target is in this list
    tf.df.1 <- rf.dat[(rf.dat$Target == tar.j),]
    
    if(nrow(tf.df.1)==0){
      tf.df.2 <- 0
    }
    else{
      # If there is a prior info then capture that with a magnitude equal to FC
      # and direction from repression/activation status
      tf.prior <- trrust_rawdata_human.1 %>% filter(.$TF==tf.df.1$RF)
      tf.prior <- ifelse(dim(tf.prior)[1] !=0,
                         tf.prior$Int.Sign,
                          1)# Take the sign if available otherwise assume a positive
      tf.df.2 <- tf.i.w*tf.prior
      # tf.df.2 <- if_else(tf.df.1$TF %in% trrust_rawdata_human.1$TF,
      #                    trrust_rawdata_human.1$Int.Sign[trrust_rawdata_human.1$TF==tf.df.1$TF]*indPriorMat[indPriorMat$TF==tf.df.1$TF,2],
      #                    indPriorMat[indPriorMat$TF==tf.df.1$TF,2])
    }
    tf.df[j,i] <- tf.df.2
  }
}
BiocManager::install('tictoc')
library(foreach)
library(doParallel)
library(tictoc)

registerDoParallel(cores = detectCores() - 3)
tf.i <- colnames(tf.df)
tr.j <- rownames(tf.df)
tf.df.2 <- foreach(i = tf.i, .combine = 'cbind') %:%
  foreach(j = tr.j, .combine = 'c') %dopar% {
    rf.dat <- tf.targets %>% filter(.$RF %in% tar.i)
    ifelse(j %in% rf.dat$Target, 1, 0)
  }
stopImplicitCluster()


library(doMC)

foreachVersion <- function(dataFrame, aVector, iterations, cores) {
  registerDoMC(cores) # unusual, but reasonable with doMC
  rows <- nrow(dataFrame)
  cols <- length(aVector)
  results <-
    foreach(i=1:iterations, .combine='rbind') %dopar% {
      # The value of the inner foreach loop is returned as
      # the value of the body of the outer foreach loop
      foreach(aElem=aVector, .combine='c') %do% {
        roadMapRow <- double(length=rows)
        roadMapRow[sample(rows,aElem)] <- 1
        sum(roadMapRow)
      }     
    }
  results <- as.data.frame(results)
  names(results) <- aVector
  results
}

# Build Prior Matrices
priorMats <- list()
for(k in 1:ncol(ATACseq_peakQuant.2)){
  indPriorMat <- ATACseq_peakQuant.2[,c(1,k)] 
  patID <- colnames(indPriorMat)[2]
  indPriorMat <- indPriorMat %>%
    filter(patID == 1)
  length(unique(indPriorMat$TF))
  tf.targets <- RF.data.1 %>% filter(.$RF.ENS %in% indPriorMat$TF)
  tf.targets <- data.frame(RF = tf.targets$RF.ENS, Target = tf.targets$ENSEMBLID)
  tfs <- unique(tf.targets$RF)
  targs <- unique(tf.targets$Target)
  tf.df <- data.frame(matrix(data = NA, 
                             ncol = length(tfs), 
                             nrow = length(targs)))
  colnames(tf.df) <- tfs
  rownames(tf.df) <- targs
  for(i in 1:nrow(tf.df)){
    tar.i <- rownames(tf.df)[i]
    for(j in 1:ncol(tf.df)){
      rf.j <- colnames(tf.df)[j]
      rf.dat <- tf.targets %>% filter(.$RF %in% rf.j)
      tf.df[i,j] <- ifelse(tar.i %in% rf.dat$Targets, 1, 0)
    }
  }
  priorMats[[paste0(colnames(indPriorMat)[2])]] <- tf.df
}

# ####################################################################
metGenPeaks <- ATACseq_PeakInf %>% 
  filter(.$ENSEMBLID %in% HumGEM)

metGenPeaks.1 <- ATACseq_count.2 %>%
  filter(ATACseq_count.2$Peaks %in% metGenPeaks$PeakID)
metGenPeaks.2 <- metGenPeaks.1 %>% 
  mutate(GeneID = metGenPeaks$ENSEMBLID[match(metGenPeaks.1$Peaks,metGenPeaks$PeakID)]) %>%
  select(-Peaks)
# colnames(metGenPeaks.2) <- gsub('__ATACseq','',colnames(metGenPeaks.2))
# metGenRF <- metGenPeaks %>% filter()

AD_CON <- intersect(c(AD_Samples$individualID,CON_Samples$individualID),colnames(metGenPeaks.2))
metGenPeaks.4 <- metGenPeaks.2 %>% select(-GeneID) %>%
  as.data.frame() %>%
  select(AD_CON)

Keys <- FreshMicro_individual_metadata$diagnosis[match(colnames(metGenPeaks.4),
                                                       FreshMicro_individual_metadata$individualID)] 

png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_ATAC-Seq-Peaks.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(metGenPeaks.4,Keys,'ATAC-Seq Peaks'))
dev.off()

# FM_138!
# Export Data to Excel

# We use regions associated with metabolic genes with the max peak
metGenPeaks.5 <- aggregate(.~GeneID,metGenPeaks.2,max)

ATAC_Data <- list()
ATAC_Data[['AD_Peaks']] <- data.frame(GeneID =metGenPeaks.5$GeneID, metGenPeaks.5 %>% select(intersect(comSampsAD,AD_CON)))
ATAC_Data[['CO_Peaks']] <- data.frame(GeneID =metGenPeaks.5$GeneID, metGenPeaks.5 %>% select(intersect(comSampsCON,AD_CON)))

AD_Expr.Norm <- Expr.Norm %>% select(comSampsAD,comSampsCON)
Keys <- FreshMicro_individual_metadata$diagnosis[match(colnames(AD_Expr.Norm),FreshMicro_individual_metadata$individualID)] 

png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_RNA-Seq-normcount.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(AD_Expr.Norm,Keys,'RNA-Seq Controls'))
dev.off()


mRNA_Data <- list()
mRNA_Data[['AD_Expr']] <- data.frame(GeneID = rownames(AD_Expr.Norm), AD_Expr.Norm %>% select(intersect(comSampsAD,AD_CON)))
mRNA_Data[['CO_Expr']] <- data.frame(GeneID = rownames(AD_Expr.Norm), AD_Expr.Norm %>% select(intersect(comSampsCON,AD_CON)))

expression.df <- mRNA_Data$AD_Expr %>% select(-GeneID)
colnames(expression.df)[1] <- ''
library(R.utils)
library(data.table)
fwrite(expression.df,'Inferelator/Input/expression.tsv',sep = '\t')
write.table(expression.df,'Inferelator/Input/expression.tsv',sep = '\t',
            row.names = T)
gzip('Inferelator/Input/expression.tsv','Inferelator/Input/expression.tsv.gz')

AD_Expr.Norm.1 <- AD_Expr.Norm %>% 
  select(c(intersect(comSampsAD,AD_CON),intersect(comSampsCON,AD_CON)))
meta_data <- data.frame(isTs = as.character(rep('FALSE',ncol(expression.df))),
                        is1stLast = as.character(rep('e',ncol(expression.df))),
                        PrevCol = as.character(rep('NA',ncol(expression.df))),
                        del.t = as.character(rep('NA',ncol(expression.df))),
                        condName = as.character(colnames(expression.df)))
fwrite(meta_data,'Inferelator/Input/meta_data.tsv',sep = '\t')
write.table(meta_data,'Inferelator/Input/meta_data.tsv',sep = '\t',
            row.names = F)


# Gold Standard: RF-gene

gold_standard <- read_delim("Inferelator/Input/gold_standard.tsv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
# write_tsv(expression.df, 'Inferelator/Input/expression.tsv')

ATAC_Data %>% writexl::write_xlsx('G:/Thesis Data/OMIC Data/Expression/AD/Processed/FreshMicro_ATAC_Data.xlsx')
mRNA_Data %>% writexl::write_xlsx('G:/Thesis Data/OMIC Data/Expression/AD/Processed/FreshMicro_mRNA_Data.xlsx')

# Which metabolic genes have RF activity in the chromosome accessibility data?
RF.Accessibility <- RF.data.1 %>% filter(.$ENSEMBLID %in% metGenPeaks$ENSEMBLID)
# Establish a gold standard matrix
rfs <- unique(RF.Accessibility$RF)
tars <- unique(RF.Accessibility$ENSEMBLID)
gs.Matrix <- data.frame(matrix(data = NA,
                               ncol = length(rfs),
                               nrow = length(tars)))

# Personalized Differential Expression: a) mRNA and b) ATAC-seq Peaks
library(readxl)
HumanGEM.genes <- read_excel("G:/Projects/NIH/Input/Animal Models/GEM Simulations/Mouse  - Alzheimer/Human-GEM.xlsx", 
                        sheet = "GENES")

metGenPeaks.6 <- metGenPeaks.5 %>% 
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# Number of peaks associated with metabolic genes in the HumanGEM
length(intersect(rownames(metGenPeaks.6),HumanGEM.genes$NAME))

ATAC_seq_cons <- data.frame(GeneID = rownames(metGenPeaks.6), metGenPeaks.6 %>% 
                              select(intersect(comSampsCON,AD_CON)) %>% 
                              as.matrix())
ATAC_seq_AD <- data.frame(GeneID = rownames(metGenPeaks.6),
                             metGenPeaks.6 %>% 
                               select(intersect(comSampsAD,AD_CON)) %>% 
                               as.matrix())
mRNA_seq_cons <- data.frame(GeneID = rownames(AD_Expr.Norm), 
                            AD_Expr.Norm %>% 
                              select(intersect(comSampsCON,AD_CON)) %>% 
                              as.matrix())
mRNA_seq_AD <- data.frame(GeneID = rownames(AD_Expr.Norm), 
                          AD_Expr.Norm %>% 
                            select(intersect(comSampsAD,AD_CON)) %>% 
                            as.matrix())

library(limma)
# library(preprocessCore)
persDEG <- function(cons, dis){
  require(dplyr)
  cons.1 <- data.matrix(cons[,-1])
  dis.1 <- data.matrix(dis[,-1])
  p.values <- data.frame(matrix(data = NA,ncol = ncol(dis.1),nrow = nrow(dis.1)))
  f.changes <- data.frame(matrix(data = NA,ncol = ncol(dis.1),nrow = nrow(dis.1)))
  for(i in 1:ncol(dis.1)){
    expDesign <- c(rep(1, ncol(cons)-1), rep(2, 1))
    exprData <- log2(abs(as.matrix(cbind(cons.1,dis.1[,i])))+1)
    design <- model.matrix(~0 + factor(expDesign))#Adjust this based on the exprimental setup
    colnames(design) <- c('control', 'disease')#
    contrast <- makeContrasts(disease - control, levels = design)
    fit <- lmFit(exprData, design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    degDF <- topTable(fit2, n = 'Inf', sort.by = 'none')[, c(1,4)]
    fc <- degDF[,1]
    for(j in 1:length(fc)){
      fc[j] <- ifelse(fc[j]>0, 2^fc[j],-1*2^abs(fc[j]))
    }
    f.changes[,i] <- fc
    p.values[,i] <- degDF[,2]
  }
  f.changes <- data.frame(GeneID = cons$GeneID, f.changes)
  colnames(f.changes)[-1] <- colnames(dis)[-1]
  p.values <- data.frame(GeneID = cons$GeneID, p.values)
  colnames(p.values)[-1] <- colnames(dis)[-1]
  return(list(p.vals = p.values, 
              f.changes = f.changes))
}


deg.mRNA <- persDEG(mRNA_seq_cons, mRNA_seq_AD)
deg.ATAC <- persDEG(ATAC_seq_cons, ATAC_seq_AD)

deg.mRNA %>% writexl::write_xlsx('G:/Projects/PhD Thesis/Input Data/OMIC Data/DEG Data/deg.mRNA.xlsx')
deg.ATAC %>% writexl::write_xlsx('G:/Projects/PhD Thesis/Input Data/OMIC Data/DEG Data/deg.ATAC.xlsx')

# Posterior Probability
BiocManager::install('BADER')
library(BADER)
set.seed(5512)
mRNA.Data <- cbind(mRNA_seq_cons[,-1],mRNA_seq_AD[,-1])
design <- factor(c(rep('control',ncol(mRNA_seq_cons)-1),rep('AD',ncol(mRNA_seq_AD)-1)))
results <- BADER(mRNA.Data,design,burn=10,reps=20)
results.1 <- cbind(results$logMeanA,results$logMeanB,results$logFoldChange,results$logDispA,
                 results$logDispB,results$diffProb)
rownames(results.1) <- rownames(mRNA.Data)

# ##################################
# EBSeq

# mRNA

# data2 <- AD_Expr %>% 
#   select(c(intersect(comSampsCON,AD_CON),intersect(comSampsAD,AD_CON)))
# data2 <- data.matrix(data2)

design <- factor(c(rep('control',ncol(mRNA_seq_cons)-1),rep('AD',ncol(mRNA_seq_AD)-1)))
mRNA.Data.1 <- data.matrix(mRNA.Data)
Sizes <- MedianNorm(mRNA.Data.1)
EBSeq.method <- EBTest(Data = mRNA.Data.1,
               Conditions= design,
               sizeFactors=Sizes, maxround=5)


DenNHist(EBSeq.method, GeneLevel = F)

data(GeneMat)
# LIMMA
design <- model.matrix(~0 + factor(design))#Adjust this based on the exprimental setup
colnames(design) <- c('control', 'disease')#
contrast <- makeContrasts(disease - control, levels = design)
fit <- lmFit(log2(mRNA.Data.1+1), design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
degDF <- topTable(fit2, n = 'Inf', sort.by = 'none')

degDF$FC <- NULL
for(i in 1:nrow(degDF)){
  degDF$FC[i] <- ifelse(degDF$logFC[i]>0,2^degDF$logFC[i],-1*2^abs(degDF$logFC[i]))
}
mRNA.DGE <- cbind(degDF$FC,EBSeq.method$PPDE) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'GeneID')
colnames(mRNA.DGE)[2:3] <- c('FC','PPDE')
writexl::write_xlsx(mRNA.DGE, 'G:/Projects/PhD Thesis/Input Data/OMIC Data/DEG Data/mRNA.DGE.xlsx')

# ATAC-Seq
# Estimate posterior probabilities
ATAC.Data <- cbind(ATAC_seq_cons[,-1],ATAC_seq_AD[,-1])
design <- factor(c(rep('control',ncol(ATAC_seq_cons)-1),rep('AD',ncol(ATAC_seq_AD)-1)))
ATAC.Data <- data.matrix(ATAC.Data)
Sizes <- MedianNorm(ATAC.Data)
EBSeq.method <- EBTest(Data = ATAC.Data,
                       Conditions= design,
                       sizeFactors=Sizes, maxround=5)

# limma differential gene expression analysis
design <- model.matrix(~0 + factor(design))#Adjust this based on the exprimental setup
colnames(design) <- c('control', 'disease')#
contrast <- makeContrasts(disease - control, levels = design)
fit <- lmFit(log2(ATAC.Data+1), design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
degDF <- topTable(fit2, n = 'Inf', sort.by = 'none')

ATAC.DGE <- cbind(degDF$logFC,EBSeq.method$PPDE) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'GeneID')
colnames(ATAC.DGE)[2:3] <- c('FC','PPDE')
writexl::write_xlsx(ATAC.DGE, 'G:/Projects/PhD Thesis/Input Data/OMIC Data/DEG Data/ATAC.DGE.xlsx')

# AD Genes
library(readr)
AD_associated_genes <- read_delim('G:/Thesis Data/OMIC Data/Priori/Open Targets/MONDO_0004975-associated-diseases.tsv', 
                                                delim = '\t', escape_double = FALSE, 
                                                trim_ws = TRUE)
brain_nutrition_uptake <- read_delim('G:/Thesis Data/GPMM-master/data/brain_nutrition_uptake.txt', 
                                     delim = '\t', escape_double = FALSE, 
                                     trim_ws = TRUE)

save.image('ThesisEpigenetics-3.RData')
