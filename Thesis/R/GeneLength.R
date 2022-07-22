

data <- read.table('G:/Thesis Data/OMIC Data/Priori/gencode.v39.annotation.gtf', header = FALSE, sep = '\t')
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

ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
geneInfo <- getBM(attributes=c('ensembl_gene_id_version','ensembl_gene_id',
                               'external_gene_name'),
                  filters = 'ensembl_gene_id_version',
                  values = geneLengthDF$ENSEMBLID,
                  mart = ensembl)
colnames(geneInfo) <- c('ENSEMBLID_VER','ENSEMBLID','SYMBOL')
library(plyr)
geneLengthDF2 <- geneLengthDF
colnames(geneLengthDF2)[1] <- c('ENSEMBLID_VER')
geneLengthDF_Final <- join(geneLengthDF2, geneInfo,by = 'ENSEMBLID_VER')
colnames(geneLengthDF_Final) <- c('ENSEMBLID_VER','TranscriptID','GENELENGTH','ENSEMBLID','SYMBOL')
saveRDS(geneLengthDF_Final,file = 'G:/Thesis Data/OMIC Data/Priori/geneLengthDF.RDS')
save.image('G:/Thesis Data/OMIC Data/Expression/AD/HSaGeneLeng.RData')
