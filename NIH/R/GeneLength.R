load('G:/Alzheimer Disease/Tx Data/GeneLengths.RData')


data <- read.table('G:/Alzheimer Disease/MouseAlz/Gene Annotation/gencode.vM28.annotation.gtf', header = FALSE, sep = '\t')
# data1 <- subset(data, data$V2=='ENSEMBL')
data$GeneLength <- abs(data$V5-data$V4)

data$GeneID <- sapply(data$V9, function(x){strsplit(x,';',fixed = T)[[1]][1]})
data$TxID <- sapply(data$V9, function(x){strsplit(x,';',fixed = T)[[1]][2]})
data$TxID <- gsub('transcript_id ','',data$TxID)
data$GeneID <- gsub('gene_id ','',data$GeneID)

library(dplyr)
library(tidyverse)

data2 <- data %>% select(c(V2, V3, GeneLength,GeneID,TxID)) %>%
  filter(.$V3 %in% c('gene','exon','transcript'))

IDs <- unique(data2$GeneID)
geneLengthDF <- data.frame()
for(i in 1:length(IDs)){
  dataIDs <- data2 %>% filter(.$GeneID==IDs[i]) %>% filter(.$V3 %in% 'exon')
  if(length(unique(dataIDs$V2))>1){
    dataIDs <- dataIDs %>% filter(.$V2 %in% 'HAVANA')
    geneLength <- sum(dataIDs$GeneLength)
  }
  else{
    geneLength <- mean(dataIDs$GeneLength)
  }
  txid <- data2[which(data2$GeneID==IDs[3]),] %>% filter(.$V3 %in% 'exon')
  txid <- txid %>% filter(str_detect(.$TxID, "^ ENS"))
  if(length(unique(txid$TxID))>1){
    geneLengthDF1 <- data.frame(matrix(data = NA, nrow = length(unique(txid$TxID)), ncol = 3))
    for(j in 1:length(unique(txid$TxID))){
      geneLengthDF1[j,] <- c(IDs[i],unique(txid$TxID)[j],geneLength)
    }
  }
  else{
    geneLengthDF1 <- data.frame(IDs[i], txid$TxID, geneLength)
  }
  geneLengthDF <- rbind(geneLengthDF,geneLengthDF1)
}
colnames(geneLengthDF) <- c('ENSEMBLID','TranscriptID','GeneLength')

library(biomaRt)

ensembl <- useMart('ensembl', dataset='mmusculus_gene_ensembl')
geneInfo <- getBM(attributes=c('ensembl_gene_id_version','ensembl_gene_id',
                               'external_gene_name'),
                  filters = 'ensembl_gene_id_version',
                  values = geneLengthDF2$ENSEMBLID,
                  mart = ensembl)
colnames(geneInfo) <- c('ENSEMBLID','ENSEMBLID_VER','SYMBOL')
library(plyr)
geneLengthDF_Final <- join(geneLengthDF, geneInfo,by = 'ENSEMBLID')
colnames(geneLengthDF_Final) <- c('ENSEMBLID_VER','TranscriptID','GENELENGTH','ENSEMBLID','SYMBOL')
saveRDS(geneLengthDF_Final,file = 'G:/Alzheimer Disease/Tx Data/geneLengthDF.RDS')


library(org.Mm.eg.db)
geneLengthDF2 <- geneLengthDF
geneLengthDF2$GeneName <- mapIds(org.Mm.eg.db,keys = geneLengthDF2$TranscriptID,column = 'SYMBOL',keytype = 'ENSEMBLTRANS')


TxID <- data$TxID[grepl('^gene_type',as.character(data$TxID))]
data <- subset(data, data$V3=='exon')
# data1 <- data[!duplicated(data$GeneID),]

data1 <- data[,c(11,10)]
data2 <- aggregate(.~GeneID, data1,sum)

mouseGeneLength <- data.frame(GeneID = data$GeneID, GeneName = data$GeneName, GeneLength = data$GeneLength)
saveRDS(mouseGeneLength, file = 'G:/Alzheimer Disease/MouseAlz/Gene Annotation/mouseGeneLength.RDS')
data <- read.table('G:/Alzheimer Disease/Hsa/Gene Annotation/gencode.v19.annotation.gtf', header = FALSE, sep = '\t')
data$GeneLength <- abs(data$V5-data$V4)
data$GeneID <- sapply(data$V9, function(x){strsplit(x,';',fixed = T)[[1]][1]})
data$GeneName <- sapply(data$V9, function(x){strsplit(x,';',fixed = T)[[1]][3]})
data$GeneName <- gsub('gene_name ','',data$GeneName)
data$GeneID <- gsub('gene_id ','',data$GeneID)
data <- subset(data, data$V3=='exon')

data1 <- data[,c(11,10)]
data2 <- aggregate(.~GeneID, data1,sum)



library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
library(GenomicFeatures)
library(vroom)
# library(xlsx)
hg19.ens <- makeTxDbFromUCSC(genome="hg19", tablename="ensGene")
exonic <- exonsBy(hg19.ens, by="gene")
red.exonic <- reduce(exonic)
exon.lengths <- sum(width(red.exonic))
exon.lengths.df <- data.frame(exon.lengths)

save.image('G:/Alzheimer Disease/Tx Data/GeneLengths.RData')
