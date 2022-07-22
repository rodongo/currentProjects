
library(GEOquery)
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
setwd('G:/Projects/NIH/Rat/Input/Transcriptome/')
RatTranscriptome <- readRDS('RatTranscriptome.RDS')
source('F:/Functions/Thesis/R/pcaPlot.R')


# RatTranscriptome$GSE161233
GSE161233 <- RatTranscriptome$GSE161233

GSE161233.PCA <- GSE161233$Expr %>% 
  select(GSE161233$Controls,GSE161233$AD) %>% 
  as.matrix()
Keys <- c(rep('control',length(GSE161233$Controls)),rep('AD',length(GSE161233$AD)))

png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_RNA-Seq-normcount-b4.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(GSE161233.PCA,Keys,'AD - Microglia: Before Outlier Removal'))
dev.off()


# RatTranscriptome$GSE129051
GSE129051 <- RatTranscriptome$GSE129051

GSE129051.PCA <- GSE129051$Expr %>% 
  select(c(GSE129051$Controls,GSE129051$AD)) %>%
  as.data.frame()
Keys <- c(rep('control',length(GSE129051$Controls)),rep('AD',length(GSE129051$AD)))


pcaANAL <- prcomp(t(2^GSE129051.PCA), scale. = T)
pcaPlotDF <- data.frame(pcaANAL$x[,1:2]) 

ggplot(data = pcaPlotDF, aes(x = PC1, y = PC2, label = row.names(pcaPlotDF),
                             color = Key))+
  geom_point()+
  geom_text(aes(label = row.names(pcaPlotDF)),hjust=0.6, vjust=0, size = 5)+
  labs(title = paste0(title))+
  theme_bw()

png(filename = paste('G:/Projects/PhD Thesis/Results/Figures/PCA_RNA-Seq-normcount-b4.png'),
    width = 14.5, height = 10, units = 'in', res = 600)
print(pcaPlot(GSE161233.PCA,Keys,'AD - Microglia: Before Outlier Removal'))
dev.off()
# RatTranscriptome$GSE85952
GSE85952 <- RatTranscriptome$GSE85952
# RatTranscriptome$GSE14499

GSE14499 <- RatTranscriptome$GSE14499
