setwd('G:/Thesis Data/OMIC Data/Expression/')
H5Fopen('myThesisData.h5')
BiocManager::install('rhdf5')
library(rhdf5)

library(edgeR)
library(GEOquery)
library(readxl)

library(biomaRt)

h5createFile('myThesisData.h5')
h5createGroup('myThesisData.h5','GSE145153')
h5createGroup('myThesisData.h5','GSE110720')
h5createGroup('myThesisData.h5','GSE129041')
h5createGroup('myThesisData.h5','GSE174367')
h5createGroup('myThesisData.h5','GSE134379')
# h5createGroup('myThesisData.h5','GSE145153')
# h5createGroup('myThesisData.h5','GSE145153')
h5ls('myThesisData.h5')
# PD
# GSE145153

dir.create('GSE145153')
gset <- getGEO('GSE145153', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE145153/')
gset <- getGEO('GSE135036', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE145153/')

gset.pheno <- gset[[1]]@phenoData@data
gset.feat <- gset[[1]]@featureData@data

gset.expr <- gset[[1]]@assayData$exprs

datFiles <- list(gset,gset.pheno,gset.feat,gset.expr)
h5write(gset.expr, 'myThesisData.h5','GSE145153/gset.expr')
# GSE110720
dir.create('GSE110720')
gset <- getGEO('GSE110720', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE110720/')
gset <- getGEO('GSE135036', GSEMatrix =TRUE, getGPL=FALSE,destdir = 'GSE145153/')

gset.pheno <- gset[[1]]@phenoData@data
gset.feat <- gset[[1]]@featureData@data

gset.expr <- gset[[1]]@assayData$exprs


# AD
BiocManager::install('rtracklayer')
library(rtracklayer)
# GSE129041
setwd('G:/Thesis Data/OMIC Data/Expression/AD/GSE129041/GSE129041_RAW/')
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
gr_narrowPeak <- import(file_narrowPeak, format = "BED",
                        extraCols = extraCols_narrowPeak)


# GSE153657
library(readr)
GSE153657 <- read_delim("G:/Thesis Data/OMIC Data/Expression/AD/GSE153658/GSE153657_gene.results.merged.tpm.txt/GSE153657_gene.results.merged.tpm.txt", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE)


# GSE174367

GSE174367 <- H5Fopen('G:/Thesis Data/OMIC Data/Expression/AD/GSE174367/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5')
GSE174367 <- h5read('G:/Thesis Data/OMIC Data/Expression/AD/GSE174367/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5')



library(readr)
GSE174367_snATAC_seq_cell_meta <- read_csv("G:/Thesis Data/OMIC Data/Expression/AD/GSE174367/GSE174367_snATAC-seq_cell_meta.csv/GSE174367_snATAC-seq_cell_meta.csv")
h5write(GSE174367_snATAC_seq_cell_meta, 'myThesisData.h5','GSE174367/GSE174367_snATAC_seq_cell_meta')



# GSE134379
GSE134379_AD <- read_delim("AD/GSE134379/GSE134379_processedSamples_cblAD.txt/GSE134379_processedSamples_cblAD.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)
h5write(GSE134379_AD, 'myThesisData.h5','GSE134379/GSE134379_AD')

h5closeAll()

# Saving multiple objects to an HDF5 file (h5save)
h5save(A, B, D, file="newfile2.h5")
h5dump("newfile2.h5")

h5delete(file = "myhdf5file.h5", name = "df")
