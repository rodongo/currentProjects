
library(GEOquery)
library(dplyr)
library(tidyverse)
library(tidyr)
library(readr)
setwd('G:/Projects/NIH/Rat/Input/Transcriptome/')

RatTranscriptome <- list()
# GSE161233
GSE161233 <- getGEO('GSE161233', destdir = 'GSE161233/', getGPL = TRUE)
meta.data <- pData(GSE161233[[1]])
GSE161233.expr <- read_csv("GSE161233/GSE161233_norm_counts_hotflash_hippo_old.csv") 

colnames(GSE161233.expr)[-1] <- meta.data$geo_accession[match(colnames(GSE161233.expr)[-1],meta.data$title)]
colnames(GSE161233.expr)[1] <- 'GeneID'
controls <- meta.data$geo_accession[meta.data$`treatment:ch1`=='SHAM']
AD <- meta.data$geo_accession[meta.data$`treatment:ch1`=='OVX']

GSE161233.expr <- GSE161233.expr %>% select(c(GeneID,controls, AD))
GSE161233 <- list(Expr = GSE161233.expr,
                  Controls = controls,
                  AD = AD)
RatTranscriptome[['GSE161233']] <- GSE161233

# GSE161142
GSE161142 <- getGEO('GSE161142', destdir = 'GSE161142/', getGPL = TRUE)
meta.data <- pData(GSE161142[[1]])
GSE161142.expr <- read_csv("GSE161142/GSE161142_norm_counts_P3_hippo_old.csv")
colnames(GSE161142.expr)[-1] <- meta.data$geo_accession[match(colnames(GSE161142.expr)[-1],meta.data$title)]
colnames(GSE161142.expr)[1] <- 'GeneID'

# GSE149540
GSE149540 <- getGEO('GSE149540', destdir = 'GSE149540/', getGPL = TRUE)
meta.data <- pData(GSE149540[[1]])
GSE161142.expr <- read_delim("GSE149540/GSE149540_fpkm.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
colnames(GSE161142.expr)[-1] <- meta.data$geo_accession[match(colnames(GSE161142.expr)[-1],meta.data$title)]
colnames(GSE161142.expr)[1] <- 'GeneID'

# GSE101112
GSE101112 <- getGEO('GSE101112', destdir = 'GSE101112/', getGPL = TRUE)
meta.data <- pData(GSE101112[[2]])
GSE101112.expr <- read_delim("GSE149540/GSE149540_fpkm.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
colnames(GSE161142.expr)[-1] <- meta.data$geo_accession[match(colnames(GSE161142.expr)[-1],meta.data$title)]
colnames(GSE161142.expr)[1] <- 'GeneID'

# GSE129051
GSE129051 <- getGEO('GSE129051', destdir = 'GSE129051/', getGPL = TRUE)
meta.data <- pData(GSE129051[[1]])
feat.data <- fData(GSE129051[[1]])
GSE129051.expr <- exprs(GSE129051[[1]])
GSE129051.expr <- GSE129051.expr %>% as.data.frame() %>%
  mutate(GeneID = feat.data$GENE_SYMBOL[match(rownames(GSE129051.expr),feat.data$ID)])
GSE129051.expr <- aggregate(.~GeneID, GSE129051.expr, max)
controls <- meta.data$geo_accession[grepl('^Sham',meta.data$title)==T]
AD <- meta.data$geo_accession[grepl('^AD',meta.data$title)==T]
GSE129051 <- list(Expr = GSE129051.expr,
                  Controls = controls,
                  AD = AD)
RatTranscriptome[['GSE129051']] <- GSE129051

# GSE129015
GSE129015 <- getGEO('GSE129015', destdir = 'GSE129015/', getGPL = TRUE)
meta.data <- pData(GSE129015[[1]])
feat.data <- fData(GSE129015[[1]])
GSE129015.expr <- exprs(GSE129015[[1]])

# GSE85952
GSE85952 <- getGEO('GSE85952', destdir = 'GSE85952/', getGPL = TRUE)
meta.data <- pData(GSE85952[[1]])
feat.data <- fData(GSE85952[[1]])
GSE85952.expr <- exprs(GSE85952[[1]])
GSE85952.expr <- GSE85952.expr %>% as.data.frame() %>%
  mutate(GeneID = feat.data$Symbol[match(rownames(GSE85952.expr),feat.data$ID)])
GSE85952.expr <- aggregate(.~GeneID, GSE85952.expr, max)
controls <- meta.data$geo_accession[grepl('^cont',meta.data$title)==T]
AD.T40 <- meta.data$geo_accession[grepl('^T40',meta.data$title)==T]
AD.T42 <- meta.data$geo_accession[grepl('^T42',meta.data$title)==T]
AD.XL40 <- meta.data$geo_accession[grepl('^XL40',meta.data$title)==T]
AD.F42 <- meta.data$geo_accession[grepl('^F42',meta.data$title)==T]

GSE85952 <- list(Expr = GSE85952.expr,
                  Controls = controls,
                 AD.T40 = AD.T40, AD.T42 = AD.T42, 
                 AD.XL40 = AD.XL40, AD.F42 = AD.F42)
RatTranscriptome[['GSE85952']] <- GSE85952


# GSE62686
GSE62686 <- getGEO('GSE62686', destdir = 'GSE62686/', getGPL = TRUE)
meta.data <- pData(GSE62686[[1]])
feat.data <- fData(GSE62686[[1]])
GSE62686.expr <- exprs(GSE62686[[1]])

# GSE37618
GSE37618 <- getGEO('GSE37618', destdir = 'GSE37618/', getGPL = TRUE)
meta.data <- pData(GSE37618[[1]])
feat.data <- fData(GSE37618[[1]])
GSE37618.expr <- exprs(GSE37618[[1]])

# GSE27620
GSE27620 <- getGEO('GSE27620', destdir = 'GSE27620/', getGPL = TRUE)
meta.data <- pData(GSE27620[[1]])
feat.data <- fData(GSE27620[[1]])
GSE27620.expr <- exprs(GSE27620[[1]])
GSE27620.expr <- GSE27620.expr %>% as.data.frame() %>%
  mutate(GeneID = feat.data$Symbol[match(rownames(GSE27620.expr),feat.data$ID)])
GSE27620.expr <- aggregate(.~GeneID, GSE85952.expr, max)
controls <- meta.data$geo_accession[grepl('^cont',meta.data$title)==T]
AD.T40 <- meta.data$geo_accession[grepl('^T40',meta.data$title)==T]
AD.T42 <- meta.data$geo_accession[grepl('^T42',meta.data$title)==T]
AD.XL40 <- meta.data$geo_accession[grepl('^XL40',meta.data$title)==T]
AD.F42 <- meta.data$geo_accession[grepl('^F42',meta.data$title)==T]

# GSE14499
GSE14499 <- getGEO('GSE14499', destdir = 'GSE14499/', getGPL = TRUE)
meta.data <- pData(GSE14499[[1]])
feat.data <- fData(GSE14499[[1]])
GSE14499.expr <- exprs(GSE14499[[1]])
GSE14499.expr <- GSE14499.expr %>% as.data.frame() %>%
  mutate(GeneID = feat.data$`Gene Symbol`[match(rownames(GSE14499.expr),feat.data$ID)])
GSE14499.expr <- aggregate(.~GeneID, GSE14499.expr, max)


controls.EC <- meta.data$geo_accession[grepl('^Entorhinal cortex, non-transgenic animals, sham lesion, biological',meta.data$title)==T]
AD.EC <- meta.data$geo_accession[grepl('^Entorhinal cortex, APP transgenic animals, GFP treated, biological',meta.data$title)==T]

controls.HP <- meta.data$geo_accession[grepl('Hippocampus, non-transgenic animals, sham lesion, biological',meta.data$title)==T]
AD.HP <- meta.data$geo_accession[grepl('^Hippocampus, APP transgenic animals, GFP treated, biological',meta.data$title)==T]


GSE14499 <- list(Expr = GSE14499.expr,
                 controls.EC = controls.EC,controls.HP = controls.HP,
                 AD.EC = AD.EC, AD.HP = AD.HP)
RatTranscriptome[['GSE14499']] <- GSE14499

# GSE9990
GSE9990 <- getGEO('GSE9990', destdir = 'GSE9990/', getGPL = TRUE)
meta.data <- pData(GSE9990[[1]])
feat.data <- fData(GSE9990[[1]])
GSE9990.expr <- exprs(GSE9990[[1]])

saveRDS(RatTranscriptome, file = 'RatTranscriptome.RDS')
