setwd('G:/Projects/NIH/Mouse/Input/Animal Models/Tx Data/')
load('Mouse.RData')

# BiocManager::install('edgeR')
# BiocManager::install('AnnotationDbi')
# BiocManager::install('biomaRt')
# BiocManager::install('org.Mm.eg.db')
# BiocManager::install('VarfromPDB')
# BiocManager::install('RISmed')
# BiocManager::install('stringi')
# BiocManager::install('readr')
# BiocManager::install('GEOquery')
# BiocManager::install('readxl')
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(VarfromPDB)
library(RISmed)
library(stringi)
library(readr)
library(tidyverse)
library(edgeR)
library(GEOquery)
library(readxl)
library(EBSeq)
library(biomaRt)

dgeData <- list()
geneLength <- readRDS(file = 'G:/Alzheimer Disease/Tx Data/geneLengthDF.RDS')
# geneLength2 <- geneLength
# geneLength2 <- geneLength2 %>%mutate(ENTREZID = mapIds(org.Mm.eg.db,keys = geneLength2$SYMBOL,column = 'ENTREZID',keytype = 'SYMBOL'))%>%
#   as.data.frame()
# A function to normalize RNA-seq data for within-sample and between-sample
# differences.
geTMM <- function(
  countsMat,#Raw counts data matrix
  sampleGrps,#A factor object with sample-group information
  geneLength,#Length of each gene in the raw count matrix
  preNorm
){
  require(edgeR)
  countsMat <- as.matrix(countsMat)
  myCPM <- cpm(countsMat, log = F)
  thresh <- myCPM > 0.1
  keep <- rowSums(thresh) >= min(table(sampleGrps))
  # Subset the rows of countdata to keep the more highly expressed genes
  countsMat <- countsMat[keep,]
  geneLength <- geneLength[keep]
  
  # if not rpk, convert and normalize with gene length
  if(preNorm=='FPKM'){
    rpk <- countsMat
  }
  else{
    rpk <- countsMat/(geneLength/10^3)
  }
  # Create a DGEList object using edgeR
  DGE <- DGEList(counts=rpk, group=sampleGrps)
  # normalize for library size by cacluating scaling factor using TMM (default method)
  DGE <- calcNormFactors(DGE, method = 'TMM')
  # count per million read (normalized count)
  normCounts <- cpm(DGE)
  return(as.data.frame(normCounts))
}
processRawSeq <- function(GSE_ID, myData,filepaths,sepPat){
  dir.create(paste0(GSE_ID,'/data'))
  untar(paste0(myData), exdir = paste0(GSE_ID,'/data'))
  cels <- list.files(paste0(GSE_ID,'/data/'), pattern = 'txt')
  sapply(paste(paste0(GSE_ID,'/data/'), cels, sep = '/'), gunzip)
  
  files <- list.files(paste0(GSE_ID,'/data/'), pattern='*.txt')
  labs <- sapply(strsplit(files, split=paste0(sepPat), fixed=TRUE), function(x) (x[1]))
  cov <- list()
  for (i in files) {
    filepath <- file.path(paste0(filepaths,i,sep=''))
    cov[[i]] <- read.table(filepath, sep = "\t", header=F, stringsAsFactors=FALSE)
    colnames(cov[[i]]) <- c('ENSEMBL', i)
  }
  df <-Reduce(function(x,y) merge(x = x, y = y, by = 'ENSEMBL'), cov)
  
  colnames(df)[-1] <- labs
  return(df)
}
library(limma)
# library(preprocessCore)
persDEG <- function(cons, dis){
  require(dplyr)
  expDesign <- c(rep(1, ncol(cons)-1), rep(2, ncol(dis)-1))
  # cons <- data %>% dplyr::select(cons)
  # dis <- data %>% dplyr::select(dis)
  
  exprData <- log2(abs(as.matrix(cbind(cons[,-1],dis[,-1]))))
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
  return(data.frame(GeneID = cons$GeneID, PPDE = degDF[,2], FC = fc))
}

DGE_PPDE <- function(AD, CON,logTrans=FALSE){
  # df <- log2(cbind(CON,AD)+1)
  # design <- factor(c(rep('control',ncol(CON)),rep('disease',ncol(AD))))
  # df <- data.matrix(df)
  # Sizes <- MedianNorm(df)
  # EBSeq.method <- EBTest(Data = df,
  #                        Conditions= design,
  #                        sizeFactors=Sizes, maxround=5)
  # LIMMA
  df <- cbind(CON,AD)
  design <- c(rep('control',ncol(CON)),rep('disease',ncol(AD)))
  design <- model.matrix(~0 + factor(design))#Adjust this based on the exprimental setup
  colnames(design) <- c('control', 'disease')#
  contrast <- makeContrasts(disease - control, levels = design)
  fit <- ifelse(logTrans==TRUE,
         lmFit(df, design),
         lmFit(log2(df+1), design))
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  degDF <- topTable(fit2, n = 'Inf', sort.by = 'none')
  
  df <- cbind(degDF$logFC, degDF$adj.P.Val, EBSeq.method$PPDE) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'GeneID')
  colnames(df)[2:4] <- c('FC','AdjPVal','PPDE')
  return(df)
}
# Gene Expression Data
# GSE163857
gset <- getGEO('GSE163857', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE163857/')

gset.pheno <- gset[[2]]@phenoData@data
# gset.expr <- exprs(gset[[1]])
# gset.feat <- fData(gset[[1]])

GSE163857 <- read_csv('GSE163857/GSE163857_Mouse_Microglia_counts.csv')
colnames(GSE163857)[-1] <- sapply(colnames(GSE163857)[-1],function(x){gsub('_merged.fastq.gz','',x)})
GSE163857 <- GSE163857 %>% remove_rownames %>% column_to_rownames(var='ENSID')

colnames(GSE163857) <- gset.pheno$geo_accession[match(gset.pheno$description.1, colnames(GSE163857))]


e3Con <- gset.pheno$geo_accession[grepl('^(E3)',gset.pheno$`genotype:ch1`) & grepl('^(Targeted)',gset.pheno$`background:ch1`)]
e4Con <- gset.pheno$geo_accession[grepl('^(E4)',gset.pheno$`genotype:ch1`) & grepl('^(Targeted)',gset.pheno$`background:ch1`)]

e3FAD <- gset.pheno$geo_accession[grepl('^(E3)',gset.pheno$`genotype:ch1`) & grepl('^(FAD)',gset.pheno$`background:ch1`)]
e4FAD <- gset.pheno$geo_accession[grepl('^(E4)',gset.pheno$`genotype:ch1`) & grepl('^(FAD)',gset.pheno$`background:ch1`)]

rownames(GSE163857) <- sapply(rownames(GSE163857), function(x){strsplit(x,'.',fixed = T)[[1]][1]})

geneLengths <- geneLength[match(rownames(GSE163857), geneLength$ENSEMBLID),]
avaiLengths <- complete.cases(geneLengths)
geneLengths <- geneLengths[avaiLengths,]
GSE163857 <- GSE163857[avaiLengths,]

# Normalization
e3e4Con <- GSE163857[,c(e3Con,e4Con)]
# e3e4Con1 <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(e3e4Con),geneLengths$ENSEMBLID)], GSE163857)
samGroups <- as.factor(c(rep('control',length(e3Con)),rep('treatment',length(e4Con))))
e3e4Con <- geTMM(as.matrix(e3e4Con),samGroups,as.numeric(geneLengths$GENELENGTH),'COUNT')
e3e4Con <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(e3e4Con),geneLengths$ENSEMBLID)], e3e4Con)


# e3e4Con <- e3e4Con[!is_empty(e3e4Con$GeneID),]


# data_frame_mod <- e3e4Con[!which(e3e4Con$GeneID==' '), ] 

GSE163857_MIC_E3E4_WT <- list()
GSE163857_MIC_E3E4_WT$E3 <- data.frame(GeneID = e3e4Con$GeneID, e3e4Con %>%dplyr::select(e3Con))
GSE163857_MIC_E3E4_WT$E4 <- data.frame(GeneID = e3e4Con$GeneID, e3e4Con %>%dplyr::select(e4Con))

GSE163857_MIC_E3E4_WT %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE163857_MIC_E3E4_WT.xlsx')

df <- aggregate(.~GeneID, e3e4Con, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,e3Con)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e4Con)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE163857_MIC_E3E4_WT')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE163857_MIC_E3E4_WT.xlsx')
# #####################
e3e4FAD <- GSE163857[,c(e3FAD,e4FAD)]
e3e4Con1 <- e3e4FAD
samGroups <- as.factor(c(rep('control',length(e3FAD)),rep('treatment',length(e4FAD))))
e3e4FAD <- geTMM(as.matrix(e3e4FAD),samGroups,as.numeric(geneLengths$GENELENGTH),'COUNTS')
e3e4FAD <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(e3e4FAD), geneLengths$ENSEMBLID)], e3e4FAD)

GSE163857_MIC_E3E4_FAD <- list()
GSE163857_MIC_E3E4_FAD$E3FAD <- data.frame(GeneID = e3e4FAD$GeneID, e3e4FAD %>%dplyr::select(e3FAD))
GSE163857_MIC_E3E4_FAD$E4FAD <- data.frame(GeneID = e3e4FAD$GeneID, e3e4FAD %>%dplyr::select(e4FAD))

GSE163857_MIC_E3E4_FAD %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE163857_MIC_E3E4_FAD.xlsx')


df <- aggregate(.~GeneID, e3e4FAD, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,e3FAD)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e4FAD)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE163857_MIC_E3E4_FAD')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE163857_MIC_E3E4_FAD.xlsx')
# ############################################################

# GSE158152
# ################################################################################################################

gset <- getGEO('GSE158152', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE158152/')
gset.pheno <- gset[[1]]@phenoData@data
GSE158152 <- read_excel("G:/Projects/NIH/Input/Animal Models/Tx Data/GSE158152/GSE158152_dst150_processed.xlsx", 
                                         sheet = "normalized_cpm")#[,-c(2:5)]
# geneLengths <- geneLength[match(GSE158152$feature_id, geneLength$ENSEMBLID),]
GSE158152 <- GSE158152 %>% remove_rownames %>% 
  column_to_rownames(var='feature_id')%>% select(-c(name, meta, source, symbol))
colnames(GSE158152) <- gset.pheno$geo_accession[match(gset.pheno$title, colnames(GSE158152))]

geneLengths <- geneLength[match(rownames(GSE158152), geneLength$ENSEMBLID),]
avaiLengths <- complete.cases(geneLengths)
geneLengths <- geneLengths[avaiLengths,]
GSE158152 <- GSE158152[avaiLengths,]

# groups
WT <- gset.pheno$geo_accession[gset.pheno$description == 'Microglia from WT animal']
Hom <- gset.pheno$geo_accession[gset.pheno$description == 'Microglia from APP-SAA homozygous animal']
Het <- gset.pheno$geo_accession[gset.pheno$description == 'Microglia from APP-SAA heterozygous animal']

# WT vs Hom and WT vs Het
GSE158152_HomExpr <- GSE158152 %>% select(c(WT,Hom))
samGroups <- as.factor(c(rep('control',length(WT)),rep('treatment',length(Hom))))
GSE158152_HomExpr <- geTMM(GSE158152_HomExpr,samGroups,as.numeric(geneLengths$GENELENGTH),'FPKM')
GSE158152_HomExpr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE158152_HomExpr), geneLengths$ENSEMBLID)], GSE158152_HomExpr)

GSE158152_MIC_FAD_APPHom <- list()
GSE158152_MIC_FAD_APPHom$WT <- data.frame(GeneID = GSE158152_HomExpr$GeneID, GSE158152_HomExpr %>%dplyr::select(WT))
GSE158152_MIC_FAD_APPHom$Hom <- data.frame(GeneID = GSE158152_HomExpr$GeneID, GSE158152_HomExpr %>%dplyr::select(Hom))

GSE158152_MIC_FAD_APPHom %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE158152_MIC_FAD_APPHom.xlsx')

# cons <- GSE158152_HomExpr %>% dplyr::select(c(GeneID,WT))
# dis <- GSE158152_HomExpr %>% dplyr::select(c(GeneID,Hom))
df <- aggregate(.~GeneID, GSE158152_HomExpr, max)
cons <- df %>% dplyr::select(c(GeneID,WT)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,Hom)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE158152_MIC_FAD_APPHom')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE158152_MIC_FAD_APPHom.xlsx')

GSE158152_HetExpr <- GSE158152 %>% select(c(WT,Het))
samGroups <- as.factor(c(rep('control',length(WT)),rep('treatment',length(Het))))
GSE158152_HetExpr <- geTMM(GSE158152_HetExpr,samGroups,as.numeric(geneLengths$GENELENGTH),'FPKM')
GSE158152_HetExpr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE158152_HetExpr), geneLengths$ENSEMBLID)], GSE158152_HetExpr)

GSE158152_MIC_FAD_APPHet <- list()
GSE158152_MIC_FAD_APPHet$WT <- data.frame(GeneID = GSE158152_HetExpr$GeneID, GSE158152_HetExpr %>%dplyr::select(WT))
GSE158152_MIC_FAD_APPHet$Het <- data.frame(GeneID = GSE158152_HetExpr$GeneID, GSE158152_HetExpr %>%dplyr::select(Het))

GSE158152_MIC_FAD_APPHet %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE158152_MIC_FAD_APPHet.xlsx')

# cons <- GSE158152_HetExpr %>% dplyr::select(c(GeneID,WT))
# dis <- GSE158152_HetExpr %>% dplyr::select(c(GeneID,Het))
df <- aggregate(.~GeneID, GSE158152_HetExpr, max)
cons <- df %>% dplyr::select(c(GeneID,WT)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,Het)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE158152_MIC_FAD_APPHet')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE158152_MIC_FAD_APPHet.xlsx')

# #######################################################################################################################################
# gset <- getGEO('GSE158153', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE158153/')
# gset <- gset[[1]]@phenoData@data
# GSE158153 <- read_excel("G:/Alzheimer Disease/Tx Data/GSE158153/GSE158153_dst168_processed.xlsx", 
#                                          sheet = "normalized_cpm")[,-c(2:5)]
# 
# GSE158153 <- GSE158153 %>% remove_rownames %>% column_to_rownames(var='feature_id')
# colnames(GSE158153) <- gset$geo_accession[match(gset$title, colnames(GSE158153))]
# 
# exprData <- edgeRNorm(GSE158153)
# rownames(exprData$counts) <- mapIds(org.Mm.eg.db,keys = rownames(exprData$counts),column = 'SYMBOL',keytype = 'ENSEMBL')
# 
# exprData <- exprData$counts %>% as.data.frame() %>%
#   rownames_to_column('GeneID')
# nas <- grepl('^(NA)', exprData$GeneID)
# exprData2 <- exprData[!nas=='TRUE',]
# # groups
# WT <- gset$geo_accession[gset$description == 'Methoxy-negative microglia from WT animal']
# HomNeg <- gset$geo_accession[gset$description == 'Methoxy-negative microglia from APP-SAA homozygous animal']
# HomPos <- gset$geo_accession[gset$description == 'Methoxy-positive microglia from APP-SAA homozygous animal']
# 
# GSE158153_HomNeg <- list()
# GSE158153_HomNeg$WT <- data.frame(GeneID = exprData2$GeneID, exprData2 %>%dplyr::select(WT))
# GSE158153_HomNeg$HomNeg <- data.frame(GeneID = exprData2$GeneID, exprData2 %>%dplyr::select(HomNeg))
# 
# GSE158153_HomNeg %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/GSE158153_HomNeg.xlsx')
# 
# GSE158153_HomPos <- list()
# GSE158153_HomPos$WT <- data.frame(GeneID = exprData2$GeneID, exprData2 %>%dplyr::select(WT))
# GSE158153_HomPos$HomPos <- data.frame(GeneID = exprData2$GeneID, exprData2 %>%dplyr::select(HomPos))
# 
# GSE158153_HomPos %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/GSE158153_HomPos.xlsx')

# ###########################################################################################################################
# GSE168137
# WT vs 4 time points
gset <- getGEO('GSE168137', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE168137/')
gset.pheno <- gset[[1]]@phenoData@data
library(readr)
GSE168137 <- read_table('GSE168137/GSE168137_countList.txt')

GSE168137 <- GSE168137 %>% remove_rownames %>% column_to_rownames(var='gene_id')
colnames(GSE168137) <- gsub(';','_',colnames(GSE168137))
colnames(GSE168137) <- gsub('Male','male',colnames(GSE168137),ignore.case = F)
colnames(GSE168137) <- gsub('Female','female',colnames(GSE168137),ignore.case = F)
colnames(GSE168137) <- gsub('5xFAD_BL6','5xFAD',colnames(GSE168137))

colnames(GSE168137) <- gset.pheno$geo_accession[match(colnames(GSE168137),gset.pheno$title)]
colnames(GSE168137)[c(66,106)] <- c('GSM5129656','GSM5129696')

rownames(GSE168137) <- sapply(rownames(GSE168137), function(x) {strsplit(x,'.',fixed = T)[[1]][1]})
geneLengths <- geneLength[match(rownames(GSE168137), geneLength$ENSEMBLID),]
avaiLengths <- complete.cases(geneLengths)
geneLengths <- geneLengths[avaiLengths,]
GSE168137 <- GSE168137[avaiLengths,]
# groups
# 4Months
WT_CO_4 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'cortex' & gset.pheno$`time (month):ch1`=='4' & grepl('^(BL6)',gset.pheno$title)==1]
WT_HI_4 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'hippocampus'& gset.pheno$`time (month):ch1`=='4' & grepl('^(BL6)',gset.pheno$title)==1]
e5xFAD_CO_4 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'cortex'& gset.pheno$`time (month):ch1`=='4' & grepl('^(5xFAD)',gset.pheno$title)==1]
e5xFAD_HI_4 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'hippocampus'& gset.pheno$`time (month):ch1`=='4' & grepl('^(5xFAD)',gset.pheno$title)==1]
# 8Months
WT_CO_8 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'cortex' & gset.pheno$`time (month):ch1`=='8' & grepl('^(BL6)',gset.pheno$title)==1]
WT_HI_8 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'hippocampus'& gset.pheno$`time (month):ch1`=='8' & grepl('^(BL6)',gset.pheno$title)==1]
e5xFAD_CO_8 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'cortex'& gset.pheno$`time (month):ch1`=='8' & grepl('^(5xFAD)',gset.pheno$title)==1]
e5xFAD_HI_8 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'hippocampus'& gset.pheno$`time (month):ch1`=='8' & grepl('^(5xFAD)',gset.pheno$title)==1]
# 12Months
WT_CO_12 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'cortex' & gset.pheno$`time (month):ch1`=='12' & grepl('^(BL6)',gset.pheno$title)==1]
WT_HI_12 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'hippocampus'& gset.pheno$`time (month):ch1`=='12' & grepl('^(BL6)',gset.pheno$title)==1]
e5xFAD_CO_12 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'cortex'& gset.pheno$`time (month):ch1`=='12' & grepl('^(5xFAD)',gset.pheno$title)==1]
e5xFAD_HI_12 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'hippocampus'& gset.pheno$`time (month):ch1`=='12' & grepl('^(5xFAD)',gset.pheno$title)==1]
# 18Months
WT_CO_18 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'cortex' & gset.pheno$`time (month):ch1`=='18' & grepl('^(BL6)',gset.pheno$title)==1]
WT_HI_18 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'hippocampus'& gset.pheno$`time (month):ch1`=='18' & grepl('^(BL6)',gset.pheno$title)==1]
e5xFAD_CO_18 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'cortex'& gset.pheno$`time (month):ch1`=='18' & grepl('^(5xFAD)',gset.pheno$title)==1]
e5xFAD_HI_18 <- gset.pheno$geo_accession[gset.pheno$`tissue:ch1` == 'hippocampus'& gset.pheno$`time (month):ch1`=='18' & grepl('^(5xFAD)',gset.pheno$title)==1]

GSE168137_CO_4_Expr <- GSE168137 %>% select(c(WT_CO_4,e5xFAD_CO_4))
samGroups <- as.factor(c(rep('control',length(WT_CO_4)),rep('treatment',length(e5xFAD_CO_4))))
GSE168137_CO_4_Expr <- geTMM(GSE168137_CO_4_Expr,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNTS')
GSE168137_CO_4_Expr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE168137_CO_4_Expr), geneLengths$ENSEMBLID)], GSE168137_CO_4_Expr)

GSE168137_CO_FAD_4M <- list()
GSE168137_CO_FAD_4M$WT_CO_4 <- data.frame(GeneID = GSE168137_CO_4_Expr$GeneID, GSE168137_CO_4_Expr %>%dplyr::select(WT_CO_4))
GSE168137_CO_FAD_4M$e5xFAD_CO_4 <- data.frame(GeneID = GSE168137_CO_4_Expr$GeneID, GSE168137_CO_4_Expr %>%dplyr::select(e5xFAD_CO_4))

GSE168137_CO_FAD_4M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE168137_CO_FAD_4M.xlsx')

# cons <- GSE168137_CO_4_Expr %>% dplyr::select(c(GeneID,WT_CO_4))
# dis <- GSE168137_CO_4_Expr %>% dplyr::select(c(GeneID,e5xFAD_CO_4))

df <- aggregate(.~GeneID, GSE168137_CO_4_Expr, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT_CO_4)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e5xFAD_CO_4)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE168137_CO_FAD_4M')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE168137_CO_FAD_4M.xlsx')


GSE168137_HI_4_Expr <- GSE168137 %>% select(c(WT_HI_4,e5xFAD_HI_4))
samGroups <- as.factor(c(rep('control',length(WT_HI_4)),rep('treatment',length(e5xFAD_HI_4))))
GSE168137_HI_4_Expr <- geTMM(GSE168137_HI_4_Expr,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNTS')
GSE168137_HI_4_Expr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE168137_HI_4_Expr), geneLengths$ENSEMBLID)], GSE168137_HI_4_Expr)

GSE168137_HI_FAD_4M <- list()
GSE168137_HI_FAD_4M$WT_HI_4 <- data.frame(GeneID = GSE168137_HI_4_Expr$GeneID, GSE168137_HI_4_Expr %>%dplyr::select(WT_HI_4))
GSE168137_HI_FAD_4M$e5xFAD_HI_4 <- data.frame(GeneID = GSE168137_HI_4_Expr$GeneID, GSE168137_HI_4_Expr %>%dplyr::select(e5xFAD_HI_4))

GSE168137_HI_FAD_4M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE168137_HI_FAD_4M.xlsx')

# cons <- GSE168137_HI_4_Expr %>% dplyr::select(c(GeneID,WT_HI_4))
# dis <- GSE168137_HI_4_Expr %>% dplyr::select(c(GeneID,e5xFAD_HI_4))

df <- aggregate(.~GeneID, GSE168137_HI_4_Expr, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT_HI_4)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e5xFAD_HI_4)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')

# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE168137_HI_FAD_4M')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE168137_HI_FAD_4M.xlsx')


GSE168137_CO_8_Expr <- GSE168137 %>% select(c(WT_CO_8,e5xFAD_CO_8))
samGroups <- as.factor(c(rep('control',length(WT_CO_8)),rep('treatment',length(e5xFAD_CO_8))))
GSE168137_CO_8_Expr <- geTMM(GSE168137_CO_8_Expr,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNTS')
GSE168137_CO_8_Expr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE168137_CO_8_Expr), geneLengths$ENSEMBLID)], GSE168137_CO_8_Expr)

GSE168137_CO_FAD_8M <- list()
GSE168137_CO_FAD_8M$WT_CO_8 <- data.frame(GeneID = GSE168137_CO_8_Expr$GeneID, GSE168137_CO_8_Expr %>%dplyr::select(WT_CO_8))
GSE168137_CO_FAD_8M$e5xFAD_CO_8 <- data.frame(GeneID = GSE168137_CO_8_Expr$GeneID, GSE168137_CO_8_Expr %>%dplyr::select(e5xFAD_CO_8))

GSE168137_CO_FAD_8M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE168137_CO_FAD_8M.xlsx')

# cons <- GSE168137_CO_8_Expr %>% dplyr::select(c(GeneID,WT_CO_8))
# dis <- GSE168137_CO_8_Expr %>% dplyr::select(c(GeneID,e5xFAD_CO_8))

df <- aggregate(.~GeneID, GSE168137_CO_8_Expr, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT_CO_8)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e5xFAD_CO_8)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE168137_CO_FAD_8M')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE168137_CO_FAD_8M.xlsx')


GSE168137_HI_8_Expr <- GSE168137 %>% select(c(WT_HI_8,e5xFAD_HI_8))
samGroups <- as.factor(c(rep('control',length(WT_HI_8)),rep('treatment',length(e5xFAD_HI_8))))
GSE168137_HI_8_Expr <- geTMM(GSE168137_HI_8_Expr,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNTS')
GSE168137_HI_8_Expr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE168137_HI_8_Expr), geneLengths$ENSEMBLID)], GSE168137_HI_8_Expr)

GSE168137_HI_FAD_8M <- list()
GSE168137_HI_FAD_8M$WT_HI_8 <- data.frame(GeneID = GSE168137_HI_8_Expr$GeneID, GSE168137_HI_8_Expr %>%dplyr::select(WT_HI_8))
GSE168137_HI_FAD_8M$e5xFAD_HI_8 <- data.frame(GeneID = GSE168137_HI_8_Expr$GeneID, GSE168137_HI_8_Expr %>%dplyr::select(e5xFAD_HI_8))

GSE168137_HI_FAD_8M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE168137_HI_FAD_8M.xlsx')

# cons <- GSE168137_HI_8_Expr %>% dplyr::select(c(GeneID,WT_HI_8))
# dis <- GSE168137_HI_8_Expr %>% dplyr::select(c(GeneID,e5xFAD_HI_8))

df <- aggregate(.~GeneID, GSE168137_HI_8_Expr, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT_HI_8)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e5xFAD_HI_8)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')

# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE168137_HI_FAD_8M')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE168137_HI_FAD_8M.xlsx')


GSE168137_CO_12_Expr <- GSE168137 %>% select(c(WT_CO_12,e5xFAD_CO_12))
samGroups <- as.factor(c(rep('control',length(WT_CO_12)),rep('treatment',length(e5xFAD_CO_12))))
GSE168137_CO_12_Expr <- geTMM(GSE168137_CO_12_Expr,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNTS')
GSE168137_CO_12_Expr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE168137_CO_12_Expr), geneLengths$ENSEMBLID)], GSE168137_CO_12_Expr)

GSE168137_CO_FAD_12M <- list()
GSE168137_CO_FAD_12M$WT_CO_12 <- data.frame(GeneID = GSE168137_CO_12_Expr$GeneID, GSE168137_CO_12_Expr %>%dplyr::select(WT_CO_12))
GSE168137_CO_FAD_12M$e5xFAD_CO_12 <- data.frame(GeneID = GSE168137_CO_12_Expr$GeneID, GSE168137_CO_12_Expr %>%dplyr::select(e5xFAD_CO_12))

GSE168137_CO_FAD_12M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE168137_CO_FAD_12M.xlsx')

# cons <- GSE168137_CO_12_Expr %>% dplyr::select(c(GeneID,WT_CO_12))
# dis <- GSE168137_CO_12_Expr %>% dplyr::select(c(GeneID,e5xFAD_CO_12))


df <- aggregate(.~GeneID, GSE168137_CO_12_Expr, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT_CO_12)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e5xFAD_CO_12)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')

# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE168137_CO_FAD_12M')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE168137_CO_FAD_12M.xlsx')

GSE168137_HI_12_Expr <- GSE168137 %>% select(c(WT_HI_12,e5xFAD_HI_12))
samGroups <- as.factor(c(rep('control',length(WT_HI_12)),rep('treatment',length(e5xFAD_HI_12))))
GSE168137_HI_12_Expr <- geTMM(GSE168137_HI_12_Expr,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNTS')
GSE168137_HI_12_Expr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE168137_HI_12_Expr), geneLengths$ENSEMBLID)], GSE168137_HI_12_Expr)

GSE168137_HI_FAD_12M <- list()
GSE168137_HI_FAD_12M$WT_HI_12 <- data.frame(GeneID = GSE168137_HI_12_Expr$GeneID, GSE168137_HI_12_Expr %>%dplyr::select(WT_HI_12))
GSE168137_HI_FAD_12M$e5xFAD_HI_12 <- data.frame(GeneID = GSE168137_HI_12_Expr$GeneID, GSE168137_HI_12_Expr %>%dplyr::select(e5xFAD_HI_12))

GSE168137_HI_FAD_12M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE168137_HI_FAD_12M.xlsx')

# cons <- GSE168137_HI_12_Expr %>% dplyr::select(c(GeneID,WT_HI_12))
# dis <- GSE168137_HI_12_Expr %>% dplyr::select(c(GeneID,e5xFAD_HI_12))

df <- aggregate(.~GeneID, GSE168137_HI_12_Expr, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT_HI_12)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e5xFAD_HI_12)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE168137_HI_FAD_12M')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE168137_HI_FAD_12M.xlsx')


GSE168137_CO_18_Expr <- GSE168137 %>% select(c(WT_CO_18,e5xFAD_CO_18))
samGroups <- as.factor(c(rep('control',length(WT_CO_18)),rep('treatment',length(e5xFAD_CO_18))))
GSE168137_CO_18_Expr <- geTMM(GSE168137_CO_18_Expr,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNTS')
GSE168137_CO_18_Expr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE168137_CO_18_Expr), geneLengths$ENSEMBLID)], GSE168137_CO_18_Expr)

GSE168137_CO_FAD_18M <- list()
GSE168137_CO_FAD_18M$WT_CO_18 <- data.frame(GeneID = GSE168137_CO_18_Expr$GeneID, GSE168137_CO_18_Expr %>%dplyr::select(WT_CO_18))
GSE168137_CO_FAD_18M$e5xFAD_CO_18 <- data.frame(GeneID = GSE168137_CO_18_Expr$GeneID, GSE168137_CO_18_Expr %>%dplyr::select(e5xFAD_CO_18))

GSE168137_CO_FAD_18M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE168137_CO_FAD_18M.xlsx')

# cons <- GSE168137_CO_18_Expr %>% dplyr::select(c(GeneID,WT_CO_18))
# dis <- GSE168137_CO_18_Expr %>% dplyr::select(c(GeneID,e5xFAD_CO_18))

df <- aggregate(.~GeneID, GSE168137_CO_18_Expr, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT_CO_18)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e5xFAD_CO_18)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE168137_CO_FAD_18M')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE168137_CO_FAD_18M.xlsx')

GSE168137_HI_18_Expr <- GSE168137 %>% select(c(WT_HI_18,e5xFAD_HI_18))
samGroups <- as.factor(c(rep('control',length(WT_HI_18)),rep('treatment',length(e5xFAD_HI_18))))
GSE168137_HI_18_Expr <- geTMM(GSE168137_HI_18_Expr,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNTS')
GSE168137_HI_18_Expr <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE168137_HI_18_Expr), geneLengths$ENSEMBLID)], GSE168137_HI_18_Expr)

GSE168137_HI_FAD_18M <- list()
GSE168137_HI_FAD_18M$WT_HI_18 <- data.frame(GeneID = GSE168137_HI_18_Expr$GeneID, GSE168137_HI_18_Expr %>%dplyr::select(WT_HI_18))
GSE168137_HI_FAD_18M$e5xFAD_HI_18 <- data.frame(GeneID = GSE168137_HI_18_Expr$GeneID, GSE168137_HI_18_Expr %>%dplyr::select(e5xFAD_HI_18))

GSE168137_HI_FAD_18M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE168137_HI_FAD_18M.xlsx')

# cons <- GSE168137_HI_18_Expr %>% dplyr::select(c(GeneID,WT_HI_18))
# dis <- GSE168137_HI_18_Expr %>% dplyr::select(c(GeneID,e5xFAD_HI_18))

df <- aggregate(.~GeneID, GSE168137_HI_18_Expr, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT_HI_18)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,e5xFAD_HI_18)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE168137_HI_FAD_18M')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE168137_HI_FAD_18M.xlsx')


# ##########################################################################
# GSE117868: DROP!

gset <- getGEO('GSE117868', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE117868/')
gset.pheno <- gset[[1]]@phenoData@data
gset.pheno$title <- gsub('RNA-Seq ','',gset.pheno$title)
library(readr)
GSE117868 <- read_table('GSE117868/GSE117868_FPKM_table_CKP25.txt') %>% remove_rownames() %>%
  column_to_rownames(var = 'Genes')

WT_CK <- gset.pheno$geo_accession[grepl('^(WT_CK_)',gset.pheno$title)==1]
WT_CKP25 <- gset.pheno$geo_accession[grepl('^(WT_CKP25)',gset.pheno$title)==1]
PU.1cKO_CKP25 <- gset.pheno$geo_accession[grepl('^(PU.1cKO_CKP25)',gset.pheno$title)==1]

geneLengths <- geneLength[match(rownames(GSE117868), geneLength$SYMBOL),]
avaiLengths <- complete.cases(geneLengths)
geneLengths <- geneLengths[avaiLengths,]
GSE117868 <- GSE117868[avaiLengths,]

colnames(GSE117868) <- c(PU.1cKO_CKP25,WT_CK,WT_CKP25)

GSE117868_Expr_CK <- GSE117868 %>% select(c(WT_CK,WT_CKP25))
samGroups <- as.factor(c(rep('control',length(WT_CK)),rep('treatment',length(WT_CKP25))))
GSE117868_Expr_CK <- geTMM(GSE117868_Expr_CK,samGroups,as.numeric(geneLengths$GENELENGTH),'FPKM')
GSE117868_Expr_CK <- data.frame(GeneID = rownames(GSE117868_Expr_CK), GSE117868_Expr_CK)

GSE117868_MIC_CK <- list()
GSE117868_MIC_CK$WT_CK <- data.frame(GeneID = GSE117868_Expr_CK$GeneID, GSE117868_Expr_CK %>%dplyr::select(WT_CK))
GSE117868_MIC_CK$WT_CKP25 <- data.frame(GeneID = GSE117868_Expr_CK$GeneID, GSE117868_Expr_CK %>%dplyr::select(WT_CKP25))

GSE117868_MIC_CK %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE117868_MIC_CK.xlsx')

# dgeData <- persDEG(WT_CK,WT_CKP25, GSE117868_Expr_CK)
dgeData <- DGE_PPDE(WT_CKP25, WT_CK)
writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE117868_MIC_CK.xlsx')

GSE117868_Expr_KO <- GSE117868 %>% select(c(WT_CK,PU.1cKO_CKP25))
samGroups <- as.factor(c(rep('control',length(WT_CK)),rep('treatment',length(PU.1cKO_CKP25))))
GSE117868_Expr_KO <- geTMM(GSE117868_Expr_KO,samGroups,as.numeric(geneLengths$GENELENGTH),'FPKM')
GSE117868_Expr_KO <- data.frame(GeneID = rownames(GSE117868_Expr_KO), GSE117868_Expr_KO)

GSE117868_MIC_KO <- list()
GSE117868_MIC_KO$WT_CK <- data.frame(GeneID = GSE117868_Expr_KO$GeneID, GSE117868_Expr_KO %>%dplyr::select(WT_CK))
GSE117868_MIC_KO$PU.1cKO_CKP25 <- data.frame(GeneID = GSE117868_Expr_KO$GeneID, GSE117868_Expr_KO %>%dplyr::select(PU.1cKO_CKP25))

GSE117868_MIC_KO %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE117868_MIC_KO.xlsx')

# dgeData <- persDEG(WT_CK,PU.1cKO_CKP25, GSE117868_Expr_KO)
dgeData <- DGE_PPDE(PU.1cKO_CKP25, WT_CK)
writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE117868_MIC_KO.xlsx')

# ##############################################################################
# # GSE165306: scRNA Seq
# gset <- getGEO('GSE165306', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE165306/')
# gset <- gset[[1]]@phenoData@data

# ##############################################################################
# GSE165306: Not well controlled!
gset <- getGEO('GSE165306', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE165306/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- exprs(gset[[1]])
gset.feat <- fData(gset[[1]])
gset.expr <- read_csv('GSE165306/GSE165306_Bulk_counts.csv')
library(readr)
GSE135999 <- read_delim("GSE135999/GSE135999_series_matrix.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      comment = "!", trim_ws = TRUE)%>% remove_rownames() %>%
  column_to_rownames(var = 'ID_REF')

# ens <- useMart('ENSEMBL_MART_ENSEMBL')
# # ensembl <- useEnsembl(biomart = 'ensembl', 
# #                       dataset = 'hsapiens_gene_ensembl', 
# #                       mirror = 'useast')
# ensData <- useDataset('mmusculus_gene_ensembl', ens)
# 
# ipro <- getBM(attributes=c('entrezgene_id','hgnc_symbol'), 
#               filters='entrezgene_id',
#               values=GSE135999$ID_REF, 
#               mart=ensData)
# ##############################################################################
# GSE181903: 
gset <- getGEO('GSE181903', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE181903/')
gset.pheno <- gset[[1]]@phenoData@data
gset.pheno$title <- str_sub(gset.pheno$title,1,nchar(gset.pheno$title)-2)
gset.expr <- exprs(gset[[1]])
gset.feat <- fData(gset[[1]])

GSE181903 <- read_csv('GSE181903/GSE181903_readcount_genename.csv') %>% remove_rownames() %>%
  column_to_rownames(var = 'gene_id') %>%
  dplyr::select(-c(gene_name,gene_chr,gene_start,gene_end,gene_strand,gene_length,gene_biotype,gene_description,tf_family))
colnames(GSE181903) <- gset.pheno$geo_accession[match(colnames(GSE181903),gset.pheno$description.1)]

WT <- gset.pheno$geo_accession[gset.pheno$title=='WT']
R47H <- gset.pheno$geo_accession[gset.pheno$title=='R47H']
R47H_Tau <- gset.pheno$geo_accession[gset.pheno$title=='R47H_Tau']

geneLengths <- geneLength[match(rownames(GSE181903), geneLength$ENSEMBLID),]
avaiLengths <- complete.cases(geneLengths)
geneLengths <- geneLengths[avaiLengths,]
GSE181903 <- GSE181903[avaiLengths,]

GSE181903_R47H <- GSE181903 %>%dplyr::select(c(WT,R47H))
samGroups <- as.factor(c(rep('control',length(WT)),rep('treatment',length(R47H))))
GSE181903_R47H <- geTMM(GSE181903_R47H,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNT')
GSE181903_R47H <- data.frame(GeneID =geneLengths$SYMBOL[match(rownames(GSE181903_R47H),geneLengths$ENSEMBLID)] , GSE181903_R47H)

GSE181903_MI_R47H <- list()
GSE181903_MI_R47H$WT <- data.frame(GeneID = GSE181903_R47H$GeneID, GSE181903_R47H %>%dplyr::select(WT))
GSE181903_MI_R47H$R47H <- data.frame(GeneID = GSE181903_R47H$GeneID, GSE181903_R47H %>%dplyr::select(R47H))

GSE181903_MI_R47H %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE181903_MI_R47H.xlsx')

# cons <- GSE181903_R47H %>% dplyr::select(c(GeneID,WT))
# dis <- GSE181903_R47H %>% dplyr::select(c(GeneID,R47H))

df <- aggregate(.~GeneID, GSE181903_R47H, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,R47H)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE181903_MI_R47H')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE181903_MI_R47H.xlsx')

GSE181903_R47H_Tau <- GSE181903 %>%dplyr::select(c(WT,R47H_Tau))
samGroups <- as.factor(c(rep('control',length(WT)),rep('treatment',length(R47H_Tau))))
GSE181903_R47H_Tau <- geTMM(GSE181903_R47H_Tau,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNT')
GSE181903_R47H_Tau <- data.frame(GeneID =geneLengths$SYMBOL[match(rownames(GSE181903_R47H_Tau),geneLengths$ENSEMBLID)], GSE181903_R47H_Tau)

GSE181903_MI_R47H_Tau <- list()
GSE181903_MI_R47H_Tau$WT <- data.frame(GeneID = GSE181903_R47H_Tau$GeneID, GSE181903_R47H_Tau %>%dplyr::select(WT))
GSE181903_MI_R47H_Tau$R47H_Tau <- data.frame(GeneID = GSE181903_R47H_Tau$GeneID, GSE181903_R47H_Tau %>%dplyr::select(R47H_Tau))

GSE181903_MI_R47H_Tau %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE181903_MI_R47H_Tau.xlsx')

# cons <- GSE181903_R47H_Tau %>% dplyr::select(c(GeneID,WT))
# dis <- GSE181903_R47H_Tau %>% dplyr::select(c(GeneID,R47H_Tau))

df <- aggregate(.~GeneID, GSE181903_R47H_Tau, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,R47H_Tau)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE181903_MI_R47H_Tau')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE181903_MI_R47H_Tau.xlsx')

# GSE165111: Agilent microarray!
library(plyr)
# dir.create('GSE165111')
gset <- getGEO('GSE165111', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE165111/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- read_excel('GSE165111/GSE165111.xlsx') #%>% remove_rownames() %>%
  #column_to_rownames(var = 'GeneID')
# B Cell
# Hippo
WT_HP <- gset.pheno$geo_accession[grepl('^(HIPP_WT)',gset.pheno$source_name_ch1)==1]
AD_HP <- gset.pheno$geo_accession[grepl('^(HIPP_3XTgAD)',gset.pheno$source_name_ch1)==1]

GSE165111_HP_B <- list()
GSE165111_HP_B$WT_CK <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(WT_HP))
GSE165111_HP_B$PU.1cKO_CKP25 <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(AD_HP))

GSE165111_HP_B %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE165111_HP_B.xlsx')

# cons <- gset.expr %>% dplyr::select(c(GeneID,WT_HP))
# dis <- gset.expr %>% dplyr::select(c(GeneID,AD_HP))

df <- aggregate(.~GeneID, gset.expr, max)
cons <- df %>% dplyr::select(c(GeneID,WT_HP)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,AD_HP)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE165111_HP_B')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE165111_HP_B.xlsx')
# Brain
WT_BR <- gset.pheno$geo_accession[grepl('^(HIPP_WT)',gset.pheno$source_name_ch1)==1]
AD_BR <- gset.pheno$geo_accession[grepl('^(BRAIN_3XTgAD)',gset.pheno$source_name_ch1)==1]

GSE165111_BR_B <- list()
GSE165111_BR_B$WT_BR <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(WT_BR))
GSE165111_BR_B$AD_BR <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(AD_BR))

GSE165111_BR_B %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE165111_BR_B.xlsx')

# cons <- gset.expr %>% dplyr::select(c(GeneID,WT_BR))
# dis <- gset.expr %>% dplyr::select(c(GeneID,AD_BR))

df <- aggregate(.~GeneID, gset.expr, max)
cons <- df %>% dplyr::select(c(GeneID,WT_BR)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,AD_BR)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE165111_BR_B')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE165111_BR_B.xlsx')

# B Cell KO
# Hippo
WT_HP_KO <- gset.pheno$geo_accession[grepl('^(BRAIN_WT)',gset.pheno$source_name_ch1)==1]
AD_HP_KO <- gset.pheno$geo_accession[grepl('^(HIPP_JHT_3XTgAD)',gset.pheno$source_name_ch1)==1]

GSE165111_HP_BKO <- list()
GSE165111_HP_BKO$WT_HP_KO <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(WT_HP_KO))
GSE165111_HP_BKO$AD_HP_KO <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(AD_HP_KO))

GSE165111_HP_BKO %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE165111_HP_BKO.xlsx')

# cons <- gset.expr %>% dplyr::select(c(GeneID,WT_HP_KO))
# dis <- gset.expr %>% dplyr::select(c(GeneID,AD_HP_KO))

df <- aggregate(.~GeneID, gset.expr, max)
cons <- df %>% dplyr::select(c(GeneID,WT_HP_KO)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,AD_HP_KO)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE165111_HP_BKO')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE165111_HP_BKO.xlsx')
# Brain
WT_BR_KO <- gset.pheno$geo_accession[grepl('^(BRAIN_WT)',gset.pheno$source_name_ch1)==1]
AD_BR_KO <- gset.pheno$geo_accession[grepl('^(BRAIN_JHT_3XTgAD)',gset.pheno$source_name_ch1)==1]

GSE165111_BR_BKO <- list()
GSE165111_BR_BKO$WT_BR_KO <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(WT_BR_KO))
GSE165111_BR_BKO$AD_BR_KO <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(AD_BR_KO))

GSE165111_BR_BKO %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE165111_BR_BKO.xlsx')

# cons <- gset.expr %>% dplyr::select(c(GeneID,WT_BR_KO))
# dis <- gset.expr %>% dplyr::select(c(GeneID,AD_BR_KO))

df <- aggregate(.~GeneID, gset.expr, max)
cons <- df %>% dplyr::select(c(GeneID,WT_BR_KO)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,AD_BR_KO)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE165111_BR_BKO')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE165111_BR_BKO.xlsx')
# GSE135999: Agilent Microarray!
dir.create('GSE135999')
gset <- getGEO('GSE135999', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE135999/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- read_excel('GSE135999/GSE135999.xlsx') #%>% remove_rownames() %>%
  #column_to_rownames(var = 'GeneID')
# Cortex
WT_CX <- gset.pheno$geo_accession[grepl('^(Cortex_WT_CTR)',gset.pheno$source_name_ch1)==1]
AD_CX <- gset.pheno$geo_accession[grepl('^(Cortex_AD_REPLICATE)',gset.pheno$source_name_ch1)==1]
GSE135999_CX_APPS1 <- list()
GSE135999_CX_APPS1$WT_CX <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(WT_CX))
GSE135999_CX_APPS1$AD_CX <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(AD_CX))

GSE135999_CX_APPS1 %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE135999_CX_APPS1.xlsx')

# cons <- gset.expr %>% dplyr::select(c(GeneID,WT_CX))
# dis <- gset.expr %>% dplyr::select(c(GeneID,AD_CX))

df <- aggregate(.~GeneID, gset.expr, max)
cons <- df %>% dplyr::select(c(GeneID,WT_CX)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,AD_CX)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE135999_CX_APPS1')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE135999_CX_APPS1.xlsx')

# Hippocampus
WT_HP <- gset.pheno$geo_accession[grepl('^(Hippocampus_WT_CTR)',gset.pheno$source_name_ch1)==1]
AD_HP <- gset.pheno$geo_accession[grepl('^(Hippocampus_AD_REPLICATE)',gset.pheno$source_name_ch1)==1]

GSE135999_HP_APPS1 <- list()
GSE135999_HP_APPS1$WT_HP <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(WT_HP))
GSE135999_HP_APPS1$AD_HP <- data.frame(GeneID = row.names(gset.expr), gset.expr %>%dplyr::select(AD_HP))

GSE135999_HP_APPS1 %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE135999_HP_APPS1.xlsx')

# cons <- gset.expr %>% dplyr::select(c(GeneID,WT_HP))
# dis <- gset.expr %>% dplyr::select(c(GeneID,AD_HP))
# dgeData <- persDEG(cons,dis)
df <- aggregate(.~GeneID, gset.expr, max)
cons <- df %>% dplyr::select(c(GeneID,WT_HP)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,AD_HP)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dgeData[[paste0('GSE135999_HP_APPS1')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE135999_HP_APPS1.xlsx')

# GSE179078: Available as counts
dir.create('GSE179078')
gset <- getGEO('GSE179078', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE179078/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- read_csv('GSE179078/GSE179078_cw2001-WT_F5KO_DG_STAR_Gene_Counts.csv') %>% 
  remove_rownames() %>% column_to_rownames(var = 'Gene_ID')
gset.expr <- gset.expr[,c(1,3,5,7,9,2,4,6,8,10)]
WT <- gset.pheno$geo_accession[grepl('^(Wild)',gset.pheno$title)==1]
Fndc5_KO <- gset.pheno$geo_accession[grepl('^(Fndc5)',gset.pheno$title)==1]
colnames(gset.expr) <- c(WT,Fndc5_KO)

geneLengths <- geneLength[match(rownames(gset.expr),toupper(geneLength$SYMBOL)),]
avaiLengths <- complete.cases(geneLengths)
geneLengths <- geneLengths[avaiLengths,]
GSE179078 <- gset.expr[avaiLengths,] #%>% remove_rownames() %>% column_to_rownames(var = 'Geneid')

GSE179078_DG <- GSE179078 %>% select(c(WT,Fndc5_KO))
samGroups <- as.factor(c(rep('control',length(WT)),rep('treatment',length(Fndc5_KO))))
GSE179078_DG <- geTMM(GSE179078_DG,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNT')
GSE179078_DG <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE179078_DG),toupper(geneLengths$SYMBOL))], GSE179078_DG)

GSE179078_DG_Fndc5_KO <- list()
GSE179078_DG_Fndc5_KO$WT <- data.frame(GeneID = GSE179078_DG$GeneID, GSE179078_DG %>%dplyr::select(WT))
GSE179078_DG_Fndc5_KO$Fndc5_KO <- data.frame(GeneID = GSE179078_DG$GeneID, GSE179078_DG %>%dplyr::select(Fndc5_KO))

GSE179078_DG_Fndc5_KO %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE179078_DG_Fndc5_KO.xlsx')

# cons <- GSE179078_DG %>% dplyr::select(c(GeneID,WT))
# dis <- GSE179078_DG %>% dplyr::select(c(GeneID,Fndc5_KO))

df <- aggregate(.~GeneID, GSE179078_DG, max)
cons <- df %>% dplyr::select(c(GeneID,WT)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,Fndc5_KO)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE179078_DG_Fndc5_KO')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE179078_DG_Fndc5_KO.xlsx')

# gset.feat <- fData(gset[[1]])
# GSE174314
dir.create('GSE174314')
gset <- getGEO('GSE174314', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE174314/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- read_csv('GSE174314/GSE174314_Raw_counts_matrix.csv') %>% select(-c(GeneName, NewName))
gset.pheno$title <- gsub(' ','_',gset.pheno$title)
colnames(gset.expr)[-1] <- gset.pheno$geo_accession[match(colnames(gset.expr)[-1],gset.pheno$title)]

WT <- gset.pheno$geo_accession[gset.pheno$`genotype:ch1`=='WT']
APP <- gset.pheno$geo_accession[gset.pheno$`genotype:ch1`=='APP/PS1']

geneLengths <- geneLength[match(gset.expr$Geneid,geneLength$ENSEMBLID),]
avaiLengths <- complete.cases(geneLengths)
geneLengths <- geneLengths[avaiLengths,]
GSE174314 <- gset.expr[avaiLengths,] %>% remove_rownames() %>% column_to_rownames(var = 'Geneid')

GSE174314_FCX <- GSE174314 %>% select(c(WT,APP))
samGroups <- as.factor(c(rep('control',length(WT)),rep('treatment',length(APP))))
GSE174314_FCX <- geTMM(GSE174314_FCX,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNT')
GSE174314_FCX <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE174314_FCX),geneLengths$ENSEMBLID)], GSE174314_FCX)

GSE174314_FCX_APPS1 <- list()
GSE174314_FCX_APPS1$WT <- data.frame(GeneID = GSE174314_FCX$GeneID, GSE174314_FCX %>%dplyr::select(WT))
GSE174314_FCX_APPS1$APP <- data.frame(GeneID = GSE174314_FCX$GeneID, GSE174314_FCX %>%dplyr::select(APP))

GSE174314_FCX_APPS1 %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE174314_FCX_APPS1.xlsx')

# cons <- GSE174314_FCX %>% dplyr::select(c(GeneID,WT))
# dis <- GSE174314_FCX %>% dplyr::select(c(GeneID,APP))
# 

df <- aggregate(.~GeneID, GSE174314_FCX, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,APP)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE174314_FCX_APPS1')]] <- DGE_PPDE(dis, cons)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE174314_FCX_APPS1.xlsx')

# # GSE176032: Low sample size
# dir.create('GSE176032')
# gset <- getGEO('GSE176032', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE176032/')
# gset.pheno <- gset[[1]]@phenoData@data
# gset.expr <- exprs(gset[[1]])
# gset.feat <- fData(gset[[1]])

# GSE163289: available as CPM, Not relevant!
# dir.create('GSE163289')
# gset <- getGEO('GSE163289', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE163289/')
# gset.pheno <- gset[[1]]@phenoData@data
# gset.expr <- read.delim('GSE163289/GSE163289_allSamples_cpm.tsv')
# 
# WT <- gset.pheno$geo_accession[gset.pheno$source_name_ch1=='5xFAD_Microglia']
# KO <- gset.pheno$geo_accession[gset.pheno$source_name_ch1=='Il3-/-5xFAD_Microglia']
# 
# colnames(gset.expr)[-1] <- c(KO, WT)
# 
# gset.expr <- gset.expr %>% as.data.frame() %>%
#   remove_rownames() %>% column_to_rownames(var = 'X')
# 
# geneLengths <- geneLength[match(row.names(gset.expr),geneLength$SYMBOL),]
# avaiLengths <- complete.cases(geneLengths)
# geneLengths <- geneLengths[avaiLengths,]
# GSE163289 <- gset.expr[avaiLengths,] %>% remove_rownames() %>% column_to_rownames(var = 'Geneid')
# 
# GSE174314_FCX <- GSE163289 %>% select(c(WT,APP))
# samGroups <- as.factor(c(rep('control',length(WT)),rep('treatment',length(APP))))
# GSE174314_FCX <- geTMM(GSE174314_FCX,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNT')
# GSE174314_FCX <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE174314_FCX),geneLengths$ENSEMBLID)], GSE174314_FCX)


# gset.feat <- fData(gset[[1]])ata
# GSE158962: scRNA Seq
dir.create('GSE158962')
gset <- getGEO('GSE158962', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE158962/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- read_excel('GSE158962/GSE158962_NFF-2020_READS-update-v2.xlsx')
# gset.feat <- fData(gset[[1]])
# GSE163725: Not appropriate
dir.create('GSE163725')
gset <- getGEO('GSE163725', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE163725/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- exprs(gset[[1]])
gset.feat <- fData(gset[[1]])
# # GSE169083: Groups to compare?
# dir.create('GSE169083')
gset <- getGEO('GSE169083', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE169083/')
gset.pheno <- gset[[1]]@phenoData@data
# 
# filepaths <- 'G:/Alzheimer Disease/Tx Data/GSE169083/data/'
# myData <- 'G:/Alzheimer Disease/Tx Data/GSE169083/GSE169083_RAW.tar'
# sepPat <- '_'
# gseID <- 'GSE169083'
# 
# 
# GSE169083 <- processRawSeq(gseID,myData,filepaths,sepPat)
# 
# geneLengths <- geneLength[match(GSE169083$ENSEMBL,geneLength$ENSEMBLID),]
# avaiLengths <- complete.cases(geneLengths)
# geneLengths <- geneLengths[avaiLengths,]
# GSE169083 <- GSE169083[avaiLengths,] %>% remove_rownames() %>% column_to_rownames(var = 'ENSEMBL')
# 
# GSE169083_FCX <- GSE169083 %>% select(c(WT,APP))
# samGroups <- as.factor(c(rep('control',length(WT)),rep('treatment',length(APP))))
# GSE174314_FCX <- geTMM(GSE174314_FCX,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNT')
# GSE174314_FCX <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE174314_FCX),geneLengths$ENSEMBLID)], GSE174314_FCX)

# GSE129296: Problems with annotation
dir.create('GSE129296')
gset <- getGEO('GSE129296', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE129296/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- exprs(gset[[1]]) %>% as.data.frame()%>% rownames_to_column(var = 'ID')

# BiocManager::install("pd.clariom.s.mouse")
# library(pd.clariom.s.mouse)
# data(pmSequence)
# data:pmSequence
gset.feat <- read_delim('GSE129296/GPL23038-70510.txt', 
                        delim = "\t", escape_double = FALSE, 
                        comment = "#", trim_ws = TRUE, skip = 1)
gset.feat2 <- data.frame(PROBE = gset.feat$probeset_id,
                         MRNA = sapply(strsplit(gset.feat$mrna_assignment, 
                                                               split='//', fixed=TRUE), 
                                                      function(x) (x[3])))
gset.feat2$Symbol <- sapply(strsplit(gset.feat2$MRNA, 
                                     split='(', fixed=TRUE), 
                            function(x) (x[2]))
library(data.table)
gset.feat2$Symbol <- NA
for(i in 1:nrow(gset.feat2)){
  t2 <- tstrsplit(gset.feat2$MRNA[i], "),")[[1]]
  t2 <- tstrsplit(t2, " ")
  if(grepl('\\(',t2[length(t2)])==FALSE){
    if(grepl('\\)',t2[length(t2)])==TRUE){
      t2 <- tstrsplit(gset.feat2$MRNA[i], "),")
      if(length(t2)>1){
        t2 <- t2[[2]]
        t2 <- tstrsplit(t2, " ")
      }
    }
  }
  t2 <- t2[length(t2)]
  gset.feat2$Symbol[i] <- gsub('\\(','',t2)
}
library(plyr)
source('F:/Functions/Thesis/R/pcaPlot.R')
gset.expr2 <- merge(gset.feat2,gset.expr,by.x='PROBE',by.y='ID') %>% 
  dplyr::select(-c(MRNA,PROBE)) %>%
  as.data.frame() %>%
  filter(.$Symbol %in% geneLength$SYMBOL) %>%
  as.data.frame()
gset.expr2 <- aggregate(.~Symbol, gset.expr2, max) %>%
  as.data.frame()
WT_3M <- gset.pheno$geo_accession[grepl('^(Microglia WT 3m)',gset.pheno$title)==1]
WT_12M <- gset.pheno$geo_accession[grepl('^(Microglia WT 12m)',gset.pheno$title)==1]
# APP
APP_3M <- gset.pheno$geo_accession[grepl('^(Microglia APP 3m)',gset.pheno$title)==1]

GSE129296_MIC_APPS1_3M <- list()
GSE129296_MIC_APPS1_3M$WT_3M <- data.frame(GeneID = gset.expr2$Symbol, gset.expr2 %>%dplyr::select(WT_3M))
GSE129296_MIC_APPS1_3M$APP_3M <- data.frame(GeneID = gset.expr2$Symbol, gset.expr2 %>%dplyr::select(APP_3M))

GSE129296_MIC_APPS1_3M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE129296_MIC_APPS1_3M.xlsx')

# cons <- gset.expr2 %>% dplyr::select(c(Symbol,WT_3M))
# colnames(cons)[1] <- 'GeneID'
# dis <- gset.expr2 %>% dplyr::select(c(Symbol,APP_3M))

df <- aggregate(.~Symbol, gset.expr2, max)[-1,]
cons <- gset.expr2 %>% dplyr::select(c(Symbol,WT_3M)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Symbol')%>% as.matrix()
dis <- gset.expr2 %>% dplyr::select(c(Symbol,APP_3M)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Symbol')%>% as.matrix()
dis <- 2^dis
cons <- 2^cons

Keys <- c(rep('control',ncol(cons)),rep('AD',ncol(dis)))
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_RNA-Seq-normcount-b4.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(cbind(cons,dis),Keys,'AD - Microglia: Before Outlier Removal'))
dev.off()

# rownames(dis) <- gset.expr2$Symbol
# rownames(cons) <- gset.expr2$Symbol
dgeData[[paste0('GSE129296_MIC_APPS1_3M')]] <- DGE_PPDE(2^dis, 2^cons, logTrans=FALSE)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE129296_MIC_APPS1_3M.xlsx')


APP_12M <- gset.pheno$geo_accession[grepl('^(Microglia APP 12m)',gset.pheno$title)==1]

GSE129296_MIC_APPS1_12M <- list()
GSE129296_MIC_APPS1_12M$WT_12M <- data.frame(GeneID = gset.expr2$Symbol, gset.expr2 %>%dplyr::select(WT_12M))
GSE129296_MIC_APPS1_12M$APP_12M <- data.frame(GeneID = gset.expr2$Symbol, gset.expr2 %>%dplyr::select(APP_12M))

GSE129296_MIC_APPS1_12M %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE129296_MIC_APPS1_12M.xlsx')

cons <- gset.expr2 %>% dplyr::select(c(Symbol,WT_12M)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Symbol')
dis <- gset.expr2 %>% dplyr::select(c(Symbol,APP_12M)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Symbol')
dis <- 2^dis
cons <- 2^cons
Keys <- c(rep('control',ncol(cons)),rep('AD',ncol(dis)))
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_RNA-Seq-normcount-b4.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(cbind(cons,dis),Keys,'AD - Microglia: Before Outlier Removal'))
dev.off()
dgeData[[paste0('GSE129296_MIC_APPS1_12M')]] <- DGE_PPDE(dis, cons, logTrans=TRUE)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE129296_MIC_APPS1_12M.xlsx')

# TAU
TAU_3M <- gset.pheno$geo_accession[grepl('^(Microglia TAU 3m)',gset.pheno$title)==1]

GSE129296_MIC_TAU_3M <- list()
GSE129296_MIC_TAU_3M$WT_3M <- data.frame(GeneID = gset.expr2$Symbol, gset.expr2 %>%dplyr::select(WT_3M))
GSE129296_MIC_TAU_3M$TAU_3M <- data.frame(GeneID = gset.expr2$Symbol, gset.expr2 %>%dplyr::select(TAU_3M))

GSE129296_MIC_TAU_3M %>% writexl::write_xlsx(path = 'G:/Projects/NIH/Input/Animal Models/Tx Data/Processed/Expression/GSE129296_MIC_TAU_3M.xlsx')

cons <- gset.expr2 %>% dplyr::select(c(Symbol,WT_3M)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Symbol')
dis <- gset.expr2 %>% dplyr::select(c(Symbol,TAU_3M)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Symbol')
dis <- 2^dis
cons <- 2^cons
Keys <- c(rep('control',ncol(cons)),rep('AD',ncol(dis)))
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_RNA-Seq-normcount-b4.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(cbind(cons,dis),Keys,'AD - Microglia: Before Outlier Removal'))
dev.off()
dgeData[[paste0('GSE129296_MIC_TAU_3M')]] <- DGE_PPDE(dis, cons, logTrans=TRUE)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE129296_MIC_TAU_3M.xlsx')

TAU_12M <- gset.pheno$geo_accession[grepl('^(Microglia TAU12m)',gset.pheno$title)==1]

GSE129296_MIC_TAU_12M <- list()
GSE129296_MIC_TAU_12M$WT_12M <- data.frame(GeneID = gset.expr2$Symbol, gset.expr2 %>%dplyr::select(WT_12M))
GSE129296_MIC_TAU_12M$TAU_12M <- data.frame(GeneID = gset.expr2$Symbol, gset.expr2 %>%dplyr::select(TAU_12M))

GSE129296_MIC_TAU_12M %>% writexl::write_xlsx(path = 'G:/Projects/NIH/Input/Animal Models/Tx Data/Processed/Expression/GSE129296_MIC_TAU_12M.xlsx')

cons <- gset.expr2 %>% dplyr::select(c(Symbol,WT_12M)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Symbol')
dis <- gset.expr2 %>% dplyr::select(c(Symbol,TAU_12M)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Symbol')
dis <- 2^dis
cons <- 2^cons
Keys <- c(rep('control',ncol(cons)),rep('AD',ncol(dis)))
png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_RNA-Seq-normcount-b4.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(cbind(cons,dis),Keys,'AD - Microglia: Before Outlier Removal'))
dev.off()
dgeData[[paste0('GSE129296_MIC_TAU_12M')]] <- DGE_PPDE(dis, cons, logTrans=FALSE)
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE129296_MIC_TAU_12M.xlsx')

# GSE149661
dir.create('GSE149661')
gset <- getGEO('GSE149661', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE149661/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- read_table2('GSE149661/GSE149661_HTSeq_counts_all_samples.txt')
gset.pheno$title <- gsub('_Brain','',gset.pheno$title)
colnames(gset.expr)[-1] <- gsub('_counts','',colnames(gset.expr)[-1])
colnames(gset.expr)[-1] <- gset.pheno$geo_accession[match(colnames(gset.expr)[-1],gset.pheno$title)]

gset.pheno$title <- paste(sapply(strsplit(gset.pheno$title, 
                                          split='_', fixed=TRUE), 
                                 function(x) (x[2])))
WT <- gset.pheno$geo_accession[gset.pheno$title=='WTPBS']
APP <- gset.pheno$geo_accession[gset.pheno$title=='APPPS1PBS']

geneLengths <- geneLength[match(gset.expr$Ensemble_ID,geneLength$ENSEMBLID),]
avaiLengths <- complete.cases(geneLengths)
geneLengths <- geneLengths[avaiLengths,]
GSE149661 <- gset.expr[avaiLengths,] %>% remove_rownames() %>% column_to_rownames(var = 'Ensemble_ID')

GSE149661_BR <- GSE149661 %>% select(c(WT,APP))
samGroups <- as.factor(c(rep('control',length(WT)),rep('treatment',length(APP))))
GSE149661_BR <- geTMM(GSE149661_BR,samGroups,as.numeric(geneLengths$GENELENGTH),'COUNT')
GSE149661_BR <- data.frame(GeneID = geneLengths$SYMBOL[match(rownames(GSE149661_BR),geneLengths$ENSEMBLID)], GSE149661_BR)

GSE149661_BR_APPS1 <- list()
GSE149661_BR_APPS1$WT <- data.frame(GeneID = GSE149661_BR$GeneID, GSE149661_BR %>%dplyr::select(WT))
GSE149661_BR_APPS1$APP <- data.frame(GeneID = GSE149661_BR$GeneID, GSE149661_BR %>%dplyr::select(APP))

GSE149661_BR_APPS1 %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE149661_BR_APPS1.xlsx')

# cons <- GSE149661_BR %>% dplyr::select(c(GeneID,WT))
# dis <- GSE149661_BR %>% dplyr::select(c(GeneID,APP))

df <- aggregate(.~GeneID, GSE149661_BR, max)[-1,]
cons <- df %>% dplyr::select(c(GeneID,WT)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
dis <- df %>% dplyr::select(c(GeneID,APP)) %>%
  remove_rownames() %>%
  column_to_rownames(var = 'GeneID')
# dgeData <- persDEG(cons,dis)
dgeData[[paste0('GSE149661_BR_APPS1')]] <- DGE_PPDE(dis, cons)

dgeData %>% writexl::write_xlsx('G:/Projects/NIH/Results/DGE/MouseDGE.xlsx')
# writexl::write_xlsx(dgeData,path = 'G:/Projects/NIH/Results/DGE/GSE149661_BR_APPS1.xlsx')

# # GSE134706: Human
# dir.create('GSE134706')
# gset <- getGEO('GSE134706', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE134706/')
# gset.pheno <- gset[[1]]@phenoData@data
# gset.expr <- read_table2('GSE134706/GSE134706_mouse_mm10_raw_.txt')

# GSE101112: drop``
dir.create('GSE101112')
gset <- getGEO('GSE101112', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE101112/')
gset.pheno <- gset[[1]]@phenoData@data
gset.expr <- read_delim('GSE101112/GSE101112_SAMP8_FPKM.txt', 
                        delim = '\t', escape_double = FALSE, 
                        trim_ws = TRUE)
gset.expr2 <- gset.expr %>% mutate(GeneSymbol = sapply(strsplit(gset.expr$`Annotation/Divergence`, 
                                                                split='|', fixed=TRUE), 
                                                       function(x) (x[1])) )
gset.expr3 <- data.frame(GeneSymbol=gset.expr2$GeneSymbol, gset.expr2[,-c(1:8,32)] %>% select(-contains('J147')))

Control <- colnames(gset.expr3[2:6])
AD <- colnames(gset.expr3[7:18])

geneLengths <- geneLength[match(gset.expr3$GeneSymbol,geneLength$SYMBOL),]
avaiLengths <- complete.cases(geneLengths)
geneLengths <- geneLengths[avaiLengths,]
GSE101112 <- gset.expr3[avaiLengths,] %>% remove_rownames() %>% column_to_rownames(var = 'GeneSymbol')

GSE101112_HP <- GSE101112 %>% select(c(Control,AD))
samGroups <- as.factor(c(rep('control',length(Control)),rep('treatment',length(AD))))
GSE101112_HP <- geTMM(GSE101112_HP,samGroups,as.numeric(geneLengths$GENELENGTH),'FPKM')
GSE101112_HP <- data.frame(GeneID = rownames(GSE101112_HP), GSE101112_HP)

GSE101112_HP_SAMP8 <- list()
GSE101112_HP_SAMP8$YOUNG <- data.frame(GeneID = GSE101112_HP$GeneID, GSE101112_HP %>%dplyr::select(Control))
GSE101112_HP_SAMP8$OLD <- data.frame(GeneID = GSE101112_HP$GeneID, GSE101112_HP %>%dplyr::select(AD))

GSE101112_HP_SAMP8 %>% writexl::write_xlsx(path = 'G:/Alzheimer Disease/Tx Data/Processed/Expression/GSE101112_HP_SAMP8.xlsx')

save.image('Mouse.RData')
