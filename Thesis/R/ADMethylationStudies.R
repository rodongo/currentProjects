
BiocManager::install('methylKit')
BiocManager::install('genomation')
BiocManager::install('GenomicRanges')
# BiocManager::install('methylKit')
library(readxl)
library(GEOquery)
library(tidyverse)

library(methylKit)
# Annotation package
library(genomation)
library(GenomicRanges)

Human_GEM <- read_excel("G:/Projects/Common/Human-GEM.xlsx") %>%
  filter(.$SUBSYSTEM == 'Exchange/demand reactions')

setwd('G:/Projects/PhD Thesis/AD/Methylation Studies/')

# GSE66351
GSE66351 <- getGEO('GSE66351', destdir = 'GSE66351/')

list.files(path = '',pattern = '*.txt')
file.list <- list(
  "/sw/courses/epigenomics/DNAmethylation/biseq_data/P6_1.bismark.cov.gz",
  "/sw/courses/epigenomics/DNAmethylation/biseq_data/P6_4.bismark.cov.gz",
  "/sw/courses/epigenomics/DNAmethylation/biseq_data/P8_3.bismark.cov.gz",
  "/sw/courses/epigenomics/DNAmethylation/biseq_data/P8_6.bismark.cov.gz")

# read the listed files into a methylRawList object making sure the other
# parameters are filled in correctly.
myobj <- methRead(file.list,
                  sample.id=list("Luminal_1","Luminal_2","Basal_1","Basal_2"),
                  pipeline = "bismarkCoverage",
                  assembly="mm10",
                  treatment=c(1,1,0,0),
                  mincov = 10
)

# check number of samples
myobj

# What type of data is stored here?
head(myobj[[1]])

# GSE109627
GSE109627 <- getGEO('GSE109627', destdir = 'GSE109627/')

# GSE134379
GSE134379 <- getGEO('GSE134379', destdir = 'GSE134379/')

# GSE80970
GSE80970 <- getGEO('GSE80970', destdir = 'GSE80970/')

# GSE105109
GSE105109 <- getGEO('GSE105109', destdir = 'GSE105109/')

# GSE59685
GSE59685 <- getGEO('GSE59685', destdir = 'GSE59685/')

