BiocManager::install('org.Mm.eg.db')
library(org.Mm.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(VarfromPDB)
library(RISmed)
library(stringi)

library(readr)
GSE163857 <- read_csv('G:/Alzheimer Disease/Tx Data/GSE163857/GSE163857_Mouse_Microglia_counts.csv')

ENSEMBLEID<- mapIds(org.Mm.eg.db,keys = GSE163857$ENSID,column = 'SYMBOL',keytype = 'ENSEMBL')
mm <- columns(org.Mm.eg.db)

x <- org.Mm.egENSEMBLTRANS
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the Ensembl gene IDs for the first five proteins
  xx[1:5]
  # Get the first one
  xx[[1]]
}
#For the reverse map ENSEMBLTRANS2EG:
# Convert to a list
xx <- as.list(org.Mm.egENSEMBLTRANS2EG)
if(length(xx) > 0){
  # Gets the entrez gene IDs for the first five Ensembl IDs
  xx[1:5]
  # Get the first one
  xx[[1]]
}



x <- org.Mm.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the Ensembl gene IDs for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
#For the reverse map ENSEMBL2EG:
# Convert to a list
xx <- as.list(org.Mm.egENSEMBL2EG)
if(length(xx) > 0){
  # Gets the entrez gene IDs for the first five Ensembl IDs
  xx[1:5]
  # Get the first one
  xx[[1]]
}